#include <libpq-fe.h>

#include "format.hpp"
#include "middle.hpp"
#include "options.hpp"
#include "osmtypes.hpp"
#include "output-gazetteer.hpp"
#include "pgsql.hpp"
#include "util.hpp"
#include "wkb.hpp"

#include <cstring>
#include <iostream>
#include <memory>

void output_gazetteer_t::delete_unused_classes(char osm_type, osmid_t osm_id)
{
    if (!m_options.append) {
        return;
    }

    if (m_style.has_data()) {
        std::string cls = m_style.class_list();
        m_copy.delete_object(osm_type, osm_id, cls);
    } else {
        /* unconditional delete of all places */
        m_copy.delete_object(osm_type, osm_id);
    }
}

void output_gazetteer_t::start()
{
    int srid = m_options.projection->target_srs();

    /* (Re)create the table unless we are appending */
    if (!m_options.append) {
        pg_conn_t conn{m_options.database_options.conninfo()};

        /* Drop any existing table */
        conn.exec("DROP TABLE IF EXISTS place CASCADE");

        /* Create the new table */

        std::string sql =
            "CREATE TABLE place ("
            "  osm_id " POSTGRES_OSMID_TYPE " NOT NULL,"
            "  osm_type char(1) NOT NULL,"
            "  class TEXT NOT NULL,"
            "  type TEXT NOT NULL,"
            "  name HSTORE,"
            "  admin_level SMALLINT,"
            "  address HSTORE,"
            "  extratags HSTORE," +
            "  geometry Geometry(Geometry,{}) NOT NULL"_format(srid) + ")";
        if (m_options.tblsmain_data) {
            sql += " TABLESPACE " + m_options.tblsmain_data.get();
        }

        conn.exec(sql);

        std::string index_sql =
            "CREATE INDEX place_id_idx ON place USING BTREE (osm_type, osm_id)";
        if (m_options.tblsmain_index) {
            index_sql += " TABLESPACE " + m_options.tblsmain_index.get();
        }
        conn.exec(index_sql);
    }
}

void output_gazetteer_t::commit() { m_copy.sync(); }

void output_gazetteer_t::process_node(osmium::Node const &node)
{
    m_copy.prepare();
    m_style.process_tags(node);
    delete_unused_classes('N', node.id());

    /* Are we interested in this item? */
    if (m_style.has_data()) {
        auto wkb = m_builder.get_wkb_node(node.location());
        m_style.copy_out(node, wkb, m_copy);
    }
}

void output_gazetteer_t::process_way(osmium::Way *way)
{
    m_copy.prepare();
    m_style.process_tags(*way);
    delete_unused_classes('W', way->id());

    /* Are we interested in this item? */
    if (m_style.has_data()) {
        /* Fetch the node details */
        m_mid->nodes_get_list(&(way->nodes()));

        /* Get the geometry of the object */
        geom::osmium_builder_t::wkb_t geom;
        if (way->is_closed()) {
            geom = m_builder.get_wkb_polygon(*way);
        }
        if (geom.empty()) {
            auto wkbs = m_builder.get_wkb_line(way->nodes(), 0.0);
            if (wkbs.empty()) {
                delete_unused_full('W', way->id());
                return;
            }

            geom = wkbs[0];
        }

        m_style.copy_out(*way, geom, m_copy);
    }
}

void output_gazetteer_t::process_relation(osmium::Relation const &rel)
{
    m_copy.prepare();

    auto const &tags = rel.tags();
    char const *type = tags["type"];
    if (!type) {
        delete_unused_full('R', rel.id());
        return;
    }

    bool is_waterway = strcmp(type, "waterway") == 0;

    if (strcmp(type, "associatedStreet") == 0 ||
        !(strcmp(type, "boundary") == 0 || strcmp(type, "multipolygon") == 0 ||
          is_waterway)) {
        delete_unused_full('R', rel.id());
        return;
    }

    m_style.process_tags(rel);
    delete_unused_classes('R', rel.id());

    /* Are we interested in this item? */
    if (!m_style.has_data()) {
        return;
    }

    /* get the boundary path (ways) */
    osmium_buffer.clear();
    auto num_ways = m_mid->rel_way_members_get(rel, nullptr, osmium_buffer);

    if (num_ways == 0) {
        delete_unused_full('R', rel.id());
        return;
    }

    for (auto &w : osmium_buffer.select<osmium::Way>()) {
        m_mid->nodes_get_list(&(w.nodes()));
    }

    auto geoms = is_waterway
                     ? m_builder.get_wkb_multiline(osmium_buffer, 0.0)
                     : m_builder.get_wkb_multipolygon(rel, osmium_buffer);

    if (geoms.empty()) {
        delete_unused_full('R', rel.id());
    } else {
        m_style.copy_out(rel, geoms[0], m_copy);
    }
}
