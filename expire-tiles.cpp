/*
 * Dirty tile list generation
 *
 * Steve Hill <steve@nexusuk.org>
 *
 * Please refer to the OpenPisteMap expire_tiles.py script for a demonstration
 * of how to make use of the output:
 * http://subversion.nexusuk.org/projects/openpistemap/trunk/scripts/expire_tiles.py
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <string>

#include <osmium/geom/tile.hpp>

#include "expire-tiles.hpp"
#include "options.hpp"
#include "reprojection.hpp"
#include "table.hpp"
#include "wkb.hpp"

#define EARTH_CIRCUMFERENCE		40075016.68
#define HALF_EARTH_CIRCUMFERENCE	(EARTH_CIRCUMFERENCE / 2)

tile_output_t::tile_output_t(const char *filename)
: outfile(fopen(filename, "a"))
{
    if (outfile == nullptr) {
        fprintf(stderr, "Failed to open expired tiles file (%s).  Tile expiry "
                        "list will not be written!\n",
                strerror(errno));
    }
}

tile_output_t::~tile_output_t()
{
    if (outfile) {
        fclose(outfile);
    }
}

void tile_output_t::output_dirty_tile(uint32_t x, uint32_t y, uint32_t zoom)
{
    if (outfile) {
        fprintf(outfile, "%i/%i/%i\n", zoom, x, y);
        ++outcount;
        if (outcount % 1000 == 0) {
            fprintf(stderr, "\rWriting dirty tile list (%iK)", outcount / 1000);
        }
    }
}

void expire_tiles::output_and_destroy(const char *filename, uint32_t minzoom)
{
    tile_output_t output_writer(filename);
    output_and_destroy<tile_output_t>(output_writer, minzoom);
}

expire_tiles::expire_tiles(uint32_t max, double bbox,
                           const std::shared_ptr<reprojection> &proj)
: max_bbox(bbox), maxzoom(max), projection(proj)
{
    if (maxzoom > 0) {
        map_width = 1U << maxzoom;
        tile_width = EARTH_CIRCUMFERENCE / map_width;
        last_tile_x = map_width + 1;
        last_tile_y = map_width + 1;
    }
}

uint64_t expire_tiles::xy_to_quadkey(uint32_t x, uint32_t y, uint32_t zoom)
{
    uint64_t quadkey = 0;
    // the two highest bits are the bits of zoom level 1, the third and fourth bit are level 2, …
    for (uint32_t z = 0; z < zoom; z++) {
        quadkey |= ((x & (1ULL << z)) << z);
        quadkey |= ((y & (1ULL << z)) << (z + 1));
    }
    return quadkey;
}

xy_coord_t expire_tiles::quadkey_to_xy(uint64_t quadkey_coord, uint32_t zoom)
{
    xy_coord_t result;
    for (uint32_t z = zoom; z > 0; --z) {
        /* The quadkey contains Y and X bits interleaved in following order: YXYX...
         * We have to pick out the bit representing the y/x bit of the current zoom
         * level and then shift it back to the right on its position in a y-/x-only
         * coordinate.*/
        result.y = result.y + static_cast<uint32_t>(
                                  (quadkey_coord & (1ULL << (2 * z - 1))) >> z);
        result.x = result.x +
                   static_cast<uint32_t>(
                       (quadkey_coord & (1ULL << (2 * (z - 1)))) >> (z - 1));
    }
    return result;
}

bool expire_tiles::valid_tile_coord(uint32_t coord)
{
    // The check for coord < 0 is not necessary because coord is unsigned.
    return coord <= map_width;
}

double expire_tiles::normalise_tile_coord(double coord)
{
    // The check for coord < 0 is not necessary because coord is unsigned.
    if (coord > map_width) {
        return (2 << maxzoom) - 1;
    } else if (coord < 0) {
        return 0;
    }
    return coord;
}

void expire_tiles::expire_tile(uint32_t x, uint32_t y)
{
    if (!valid_tile_coord(x) || !valid_tile_coord(y)) {
        return;
    }
    // Only try to insert to tile into the set if the last inserted tile
    // is different from this tile.
    if (last_tile_x != x || last_tile_y != y) {
        m_dirty_tiles.insert(xy_to_quadkey(x, y, maxzoom));
        last_tile_x = x;
        last_tile_y = y;
    }
}

void expire_tiles::expire_line_segment(double x1, double y1, double x2,
                                       double y2)
{
    assert(x1 <= x2);
    assert(x2 - x1 <= map_width / 2);
    if (x1 == x2 && y1 == y2) {
        // The line is degenerated and only a point.
        return;
    }
    // The following if block ensures that x2-x1 does not cause an
    // underflow which could cause a division by zero.
    if ((x2 - x1 < 1) && ((static_cast<int>(x2) == static_cast<int>(x1)) ||
                          (x2 - x1 < 0.00000001))) {
        /**
         * Case 1: The linestring is parallel to a meridian or does not
         * cross a tile border. Therefore we can treat it as a vertical
         * linestring.
         *
         * Case 2: This linestring is almost parallel (very small error).
         * We just treat it as a parallel of a meridian. The resulting
         * error is negligible.
         */
        if (y2 < y1) {
            // swap coordinates
            double temp = y2;
            y2 = y1;
            y1 = temp;
        }
        double x_buffer_west =
            normalise_tile_coord(x1 - TILE_EXPIRY_LEEWAY);
        expire_vertical_line(x_buffer_west, y1, y2);
        // Expire parallels of this line with a distance of
        // TILE_EXPIRY_LEEWAY. If it is not necessary because the parallels
        // run through the same tiles, we don't call expire_vertical_line()
        // again.
        double x_buffer_east =
            normalise_tile_coord(x1 + TILE_EXPIRY_LEEWAY);
        if (static_cast<uint32_t>(x_buffer_west) !=
            static_cast<uint32_t>(x_buffer_east)) {
            expire_vertical_line(x_buffer_east, y1, y2);
        }
        return;
    }
    // build the buffer
    // length of the vector from (x1,y1) to (x2,y2)
    double segment_length = sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    double x_norm = (x2 - x1) / segment_length;
    double y_norm = (y2 - y1) / segment_length;
    double x_buffer = 0.1 * x_norm;
    double y_buffer = 0.1 * y_norm;
    // normal vector to the right: (-y,x), to the left: (y,-x)
    // left parallel:
    expire_line(x1 - x_buffer - y_buffer, y1 - y_buffer + x_buffer,
                x2 + x_buffer - y_buffer, y2 + y_buffer + x_buffer);
    // right parallel:
    expire_line(x1 - x_buffer + y_buffer, y1 - y_buffer - x_buffer,
                x2 + x_buffer + y_buffer, y2 + y_buffer - x_buffer);
}

void expire_tiles::expire_line(double x1, double y1, double x2, double y2)
{
    assert(x1 < x2);
    // y(x) = m * x + c with incline as m and y_intercept as c
    double incline = (y2 - y1) / (x2 - x1);
    double y_intercept = y2 - incline * x2;

    // Check if it is a horizontal line and lies fully outside the domain of
    // definition.
    if (incline == 0 && (y_intercept < 0 || y_intercept > map_width)) {
        return;
    }

    // If x2 is smaller or equal than 0, x1 will be also smaller and the whole
    // line outside the domain of the Web Mercator projection.
    if (x2 <= 0) {
        return;
    }
    // Check if x1 is smaller than 0. Set x to 0 and y1 to the correct value
    // if it is so.
    if (x1 < 0) {
        x1 = 0;
        y1 = y_intercept;
    }
    // It is not a problem if any of x1, x2, y1 and y2 is larger than map_width
    // because expire_tile() will just return without doing anything if we want
    // to expire a tile with x > map_width. It is only a problem if one of them
    // is smaller than 0 because expire_tile() accepts unsigned integers.
    if (y1 < 0) {
        y1 = 0;
        // 0 = incline * x + y_intercept
        x1 = -y_intercept / incline;
    }
    if (y2 < 0) {
        y2 = 0;
        x2 = -y_intercept / incline;
    }

    // mark start tile as expired
    expire_tile(static_cast<uint32_t>(x1), static_cast<uint32_t>(y1));
    // expire all tiles the line enters by crossing their western edge
    for (uint32_t x = static_cast<uint32_t>(x1 + 1);
         x <= static_cast<uint32_t>(x2); ++x) {
        double y = incline * x + y_intercept;
        expire_tile(x, static_cast<uint32_t>(y));
    }
    // the same for all tiles which are entered by crossing their northern edge
    double min_y = std::min(y1, y2);
    double max_y = std::max(y1, y2);
    for (uint32_t y = static_cast<uint32_t>(min_y + 1);
         y <= static_cast<uint32_t>(max_y); ++y) {
        double x = (y - y_intercept) / incline;
        if (y2 > y1) {
            // line going to the south on its way from (x1,y1) to (x2,y2)
            expire_tile(static_cast<uint32_t>(x), y);
        } else {
            expire_tile(static_cast<uint32_t>(x), y - 1);
        }
    }
}

void expire_tiles::expire_vertical_line(double x, double y1, double y2)
{
    assert(y1 < y2); // line in correct order and not collapsed
    // mark the tile of the southern end and its buffer as expired
    from_bbox(x, y1, x, y1);
    // mark all tiles above it as expired until we reach the northern end of
    // the line
    for (uint32_t y = static_cast<uint32_t>(y1 + 1);
         y < static_cast<uint32_t>(y2); y++) {
        expire_tile(static_cast<uint32_t>(x), y);
    }
    // mark the tile at the northern end and its buffer as expired
    from_bbox(x, y2, x, y2);
}

void expire_tiles::from_line_lon_lat(double lon_a, double lat_a, double lon_b,
                                     double lat_b)
{
    double tile_x_a;
    double tile_y_a;
    double tile_x_b;
    double tile_y_b;
    projection->coords_to_tile(&tile_x_a, &tile_y_a, lon_a, lat_a, map_width);
    projection->coords_to_tile(&tile_x_b, &tile_y_b, lon_b, lat_b, map_width);
    // swap ends of this segment if necessary because we go from left to right
    if (tile_x_a > tile_x_b) {
        std::swap(tile_x_a, tile_x_b);
        std::swap(tile_y_a, tile_y_b);
    }
    if (tile_x_b - tile_x_a > map_width / 2) {
        // line crosses 180th meridian → split the line at its intersection
        // with this meridian
        double y_split; // y-coordinate of intersection point
        if (tile_x_b == map_width && tile_x_a == 0) {
            // The line is part of the 180th meridian. We have to treat this in
            // a special way, otherwise there will
            // be a division by 0 in the following code.
            expire_line_segment(0, tile_y_a, 0, tile_y_b);
            return;
        }
        // This line runs from western to eastern hemisphere over the 180th meridian
        // use intercept theorem to get the intersection point of the line and
        // the 180th meridian
        // x-distance between left point and 180th meridian
        double x_distance = map_width + tile_x_a - tile_x_b;
        // apply intercept theorem: (y2-y1)/(y_split-y1) = (x2-x1)/(x_split-x1)
        y_split = tile_y_a + (tile_y_b - tile_y_a) * (tile_x_a / x_distance);
        expire_line_segment(0, y_split, tile_x_a, tile_y_a);
        expire_line_segment(tile_x_b, tile_y_b, map_width, y_split);
    } else {
        expire_line_segment(tile_x_a, tile_y_a, tile_x_b, tile_y_b);
    }
}

void expire_tiles::from_point(double lon, double lat)
{
    double tile_x;
    double tile_y;
    projection->coords_to_tile(&tile_x, &tile_y, lon, lat, map_width);
    from_bbox(tile_x, tile_y, tile_x, tile_y);
}

void expire_tiles::from_bbox_lon_lat(double min_x, double min_y, double max_x,
                                     double max_y)
{
    double x_min;
    double y_min;
    double x_max;
    double y_max;
    projection->coords_to_tile(&x_min, &y_max, min_x, min_y, map_width);
    projection->coords_to_tile(&x_max, &y_min, max_x, max_y, map_width);
    from_bbox(x_min, y_min, x_max, y_max);
}

void expire_tiles::from_bbox(double min_x, double min_y, double max_x,
                             double max_y)
{
    min_x -= TILE_EXPIRY_LEEWAY;
    min_y -= TILE_EXPIRY_LEEWAY;
    max_x += TILE_EXPIRY_LEEWAY;
    max_y += TILE_EXPIRY_LEEWAY;
    from_bbox_without_buffer(
        static_cast<uint32_t>(min_x), static_cast<uint32_t>(min_y),
        static_cast<uint32_t>(max_x), static_cast<uint32_t>(max_y));
}

void expire_tiles::from_bbox_without_buffer(uint32_t min_x, uint32_t min_y,
                                            uint32_t max_x, uint32_t max_y)
{
    for (uint32_t iterator_x = min_x; iterator_x <= max_x; ++iterator_x) {
        for (uint32_t iterator_y = min_y; iterator_y <= max_y; ++iterator_y) {
            expire_tile(iterator_x, iterator_y);
        }
    }
}


void expire_tiles::from_wkb(const char *wkb, osmid_t osm_id)
{
    if (maxzoom == 0) {
        return;
    }

    auto parse = ewkb::parser_t(wkb);

    int header = parse.read_header();
    switch (header) {
    case ewkb::wkb_point:
        from_wkb_point(&parse);
        break;
    case ewkb::wkb_line:
        from_wkb_line(&parse);
        break;
    case ewkb::wkb_polygon:
        from_wkb_polygon(&parse, osm_id);
        break;
    case ewkb::wkb_multi_line: {
        auto num = parse.read_length();
        for (unsigned i = 0; i < num; ++i) {
            parse.read_header();
            from_wkb_line(&parse);
        }
        break;
    }
    case ewkb::wkb_multi_polygon: {
        auto num = parse.read_length();
        for (unsigned i = 0; i < num; ++i) {
            parse.read_header();
            from_wkb_polygon(&parse, osm_id);
        }
        break;
    }
    default:
        fprintf(stderr, "OSM id %" PRIdOSMID
                        ": Unknown geometry type %d, cannot expire.\n",
                osm_id, header);
    }
}

void expire_tiles::from_wkb_point(ewkb::parser_t *wkb)
{
    auto c = wkb->read_point();
    from_point(c.x, c.y);
}

void expire_tiles::from_wkb_line(ewkb::parser_t *wkb)
{
    auto sz = wkb->read_length();

    if (sz == 0) {
        return;
    }

    if (sz == 1) {
        from_wkb_point(wkb);
    } else {
        auto prev = wkb->read_point();
        for (size_t i = 1; i < sz; ++i) {
            auto cur = wkb->read_point();
            from_line_lon_lat(prev.x, prev.y, cur.x, cur.y);
            prev = cur;
        }
    }
}

void expire_tiles::from_wkb_polygon(ewkb::parser_t *wkb, osmid_t osm_id)
{
    auto num_rings = wkb->read_length();
    assert(num_rings > 0);

    auto start = wkb->save_pos();

    auto num_pt = wkb->read_length();
    auto initpt = wkb->read_point();

    osmium::geom::Coordinates min{initpt}, max{initpt};

    // get bounding box of the polygon
    for (size_t i = 1; i < num_pt; ++i) {
        auto c = wkb->read_point();
        if (c.x < min.x)
            min.x = c.x;
        if (c.y < min.y)
            min.y = c.y;
        if (c.x > max.x)
            max.x = c.x;
        if (c.y > max.y)
            max.y = c.y;
    }
    wkb->rewind(start);
    /* Bounding boxes wider than half of the circumfence of the earth are
       treated as evil polygons because
       (1) they currently do not exist in OSM due to (2),
       (2) most software does not handle them correctly,
       (3) it is not unsafe if they are not expired.
       We had to split them at the antimeridian if we want to handle them
       properly. */
    if (max.x - min.x > max_bbox || max.y - min.y > max_bbox) {
        // expire all rings as if they were only lines
        for (unsigned ring = 0; ring < num_rings; ++ring) {
            wkb->rewind(start);
            from_wkb_line(wkb);
            return;
        }
    }
    // reproject coordinates of bounding box
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    // min and max are swapped when calling projection->coords_to_tile()
    // because
    projection->coords_to_tile(&min_x, &min_y, min.x, max.y, map_width);
    projection->coords_to_tile(&max_x, &max_y, max.x, min.y, map_width);

    // If the polygon does not cross the border between two tile columns in
    // maxzoom, it can be simply expired by expiring its bounding box.
    if (static_cast<uint32_t>(min_x) == static_cast<uint32_t>(max_x)) {
        from_bbox(min_x, min_y, max_x, max_y);
    }

    wkb->rewind(start);
    // expire interior of outer ring and a few tiles more
    // arguments swapped because of different direction of x and y axis
    intersecting_tiles_t tiles(min_x, max_x, map_width, TILE_EXPIRY_LEEWAY);
    for (unsigned ring = 0; ring < num_rings; ++ring) {
        auto ring_size = wkb->read_length();
        if (ring_size <= 1 && ring == 0) {
            // outer ring degenerated, ignore the whole polygon
            return;
        } else if (ring_size <= 3) {
            /* degenerated inner rings don't reduce the number of expired tiles
             * We don't have to care for them. */
            continue;
        } else {
            auto prev = wkb->read_point();
            for (size_t i = 1; i < ring_size; ++i) {
                auto cur = wkb->read_point();
                // reproject the coordinates
                double tile_x_a;
                double tile_y_a;
                double tile_x_b;
                double tile_y_b;
                projection->coords_to_tile(&tile_x_a, &tile_y_a, prev.x, prev.y,
                                           map_width);
                projection->coords_to_tile(&tile_x_b, &tile_y_b, cur.x, cur.y,
                                           map_width);
                // ring == 0 is an outer ring, all other rings are inner rings
                tiles.evaluate_segment(tile_x_a, tile_y_a, tile_x_b, tile_y_b,
                                       (ring == 0));
                prev = cur;
            }
        }
    }
    // mark tiles as expired
    tiles.sort_bounds();
    do {
        while (tiles.column_has_intervals()) {
            std::unique_ptr<std::pair<uint32_t, uint32_t>> interval =
                tiles.get_next_pair();
            //TODO handling if interval is nullptr
            if (!interval) {
                continue;
            }
            //TODO handling last column
            from_bbox_without_buffer(tiles.get_current_x(), interval->first,
                                     tiles.get_current_x(), interval->second);
        }
    } while (tiles.move_to_next_column());
}

/*
 * Expire tiles based on an osm element.
 * What type of element (node, line, polygon) osm_id refers to depends on
 * sql_conn. Each type of table has its own sql_conn and the prepared statement
 * get_wkb refers to the appropriate table.
 *
 * The function returns -1 if expiry is not enabled. Otherwise it returns the number
 * of elements that refer to the osm_id.

 */
int expire_tiles::from_db(table_t* table, osmid_t osm_id) {
    //bail if we dont care about expiry
    if (maxzoom == 0)
        return -1;

    //grab the geom for this id
    auto wkbs = table->get_wkb_reader(osm_id);

    //dirty the stuff
    const char* wkb = nullptr;
    while ((wkb = wkbs.get_next())) {
        auto binwkb = ewkb::parser_t::wkb_from_hex(wkb);
        from_wkb(binwkb.c_str(), osm_id);
    }

    //return how many rows were affected
    return wkbs.get_count();
}

void expire_tiles::merge_and_destroy(expire_tiles &other)
{
    if (map_width != other.map_width) {
        throw std::runtime_error(
            (boost::format("Unable to merge tile expiry sets when "
                           "map_width does not match: %1% != %2%.") %
             map_width % other.map_width)
                .str());
    }

    if (tile_width != other.tile_width) {
        throw std::runtime_error(
            (boost::format("Unable to merge tile expiry sets when "
                           "tile_width does not match: %1% != %2%.") %
             tile_width % other.tile_width)
                .str());
    }

    if (m_dirty_tiles.size() == 0) {
        m_dirty_tiles = std::move(other.m_dirty_tiles);
    } else {
        m_dirty_tiles.insert(other.m_dirty_tiles.begin(),
                             other.m_dirty_tiles.end());
    }

    other.m_dirty_tiles.clear();
}
