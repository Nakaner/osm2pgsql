#ifndef OSM2PGSQL_OUTPUT_GAZETTEER_HPP
#define OSM2PGSQL_OUTPUT_GAZETTEER_HPP

#include <memory>
#include <string>

#include <osmium/memory/buffer.hpp>

#include "db-copy-mgr.hpp"
#include "gazetteer-style.hpp"
#include "osmium-builder.hpp"
#include "osmtypes.hpp"
#include "output.hpp"
#include "util.hpp"

class output_gazetteer_t : public output_t
{
    output_gazetteer_t(output_gazetteer_t const *other,
                       std::shared_ptr<middle_query_t> const &cloned_mid,
                       std::shared_ptr<db_copy_thread_t> const &copy_thread)
    : output_t(cloned_mid, other->m_options), m_copy(copy_thread),
      m_builder(other->m_options.projection, true),
      osmium_buffer(PLACE_BUFFER_SIZE, osmium::memory::Buffer::auto_grow::yes)
    {}

public:
    output_gazetteer_t(std::shared_ptr<middle_query_t> const &mid,
                       options_t const &options,
                       std::shared_ptr<db_copy_thread_t> const &copy_thread)
    : output_t(mid, options), m_copy(copy_thread),
      m_builder(options.projection, true),
      osmium_buffer(PLACE_BUFFER_SIZE, osmium::memory::Buffer::auto_grow::yes)
    {
        m_style.load_style(options.style);
    }

    std::shared_ptr<output_t>
    clone(std::shared_ptr<middle_query_t> const &mid,
          std::shared_ptr<db_copy_thread_t> const &copy_thread) const override
    {
        return std::shared_ptr<output_t>(
            new output_gazetteer_t(this, mid, copy_thread));
    }

    void start() override;
    void stop(osmium::thread::Pool *) override {}
    void commit() override;

    void enqueue_ways(pending_queue_t &, osmid_t, size_t, size_t &) override {}
    void pending_way(osmid_t, int) override {}

    void enqueue_relations(pending_queue_t &, osmid_t, size_t,
                           size_t &) override
    {}
    void pending_relation(osmid_t, int) override {}

    void node_add(osmium::Node const &node) override { process_node(node); }

    void way_add(osmium::Way *way) override { process_way(way); }

    void relation_add(osmium::Relation const &rel) override
    {
        process_relation(rel);
    }

    void node_modify(osmium::Node const &node) override { process_node(node); }

    void way_modify(osmium::Way *way) override { process_way(way); }

    void relation_modify(osmium::Relation const &rel) override
    {
        process_relation(rel);
    }

    void node_delete(osmid_t id) override { m_copy.delete_object('N', id); }

    void way_delete(osmid_t id) override { m_copy.delete_object('W', id); }

    void relation_delete(osmid_t id) override { m_copy.delete_object('R', id); }

private:
    enum
    {
        PLACE_BUFFER_SIZE = 4096
    };

    /// Delete all places that are not covered by the current style results.
    void delete_unused_classes(char osm_type, osmid_t osm_id);
    void process_node(osmium::Node const &node);
    void process_way(osmium::Way *way);
    void process_relation(osmium::Relation const &rel);

    void delete_unused_full(char osm_type, osmid_t osm_id)
    {
        if (m_options.append) {
            m_copy.delete_object(osm_type, osm_id);
        }
    }

    gazetteer_copy_mgr_t m_copy;
    gazetteer_style_t m_style;

    geom::osmium_builder_t m_builder;
    osmium::memory::Buffer osmium_buffer;
};

#endif // OSM2PGSQL_OUTPUT_GAZETTEER_HPP
