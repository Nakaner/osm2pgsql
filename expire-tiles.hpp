#ifndef EXPIRE_TILES_H
#define EXPIRE_TILES_H

#include <memory>
#include <unordered_set>

#include "intersecting_tiles.hpp"
#include "osmtypes.hpp"

class reprojection;
class table_t;
class tile;
namespace ewkb {
class parser_t;
}

/**
 * \brief Simple struct for the x and y index of a tile ID.
 */
struct xy_coord_t
{
    uint32_t x;
    uint32_t y;
    xy_coord_t() : x(0), y(0) {}
};

/**
 * Implementation of the output of the tile expiry list to a file.
 */
class tile_output_t
{
    FILE *outfile;
    uint32_t outcount = 0;

public:
    tile_output_t(const char *filename);

    ~tile_output_t();

    /**
     * Output dirty tile.
     *
     * \param x x index
     * \param y y index
     * \param zoom zoom level of the tile
     */
    void output_dirty_tile(uint32_t x, uint32_t y, uint32_t zoom);
};

struct expire_tiles
{
    expire_tiles(uint32_t maxzoom, double maxbbox,
                 const std::shared_ptr<reprojection> &projection);

    /**
     * Expire the tile including a small buffer around it where a point is located.
     *
     * \param lon longitude of the point
     * \param lat latitude of the point
     */
    void from_point(double lon, double lat);

    /**
     * Expire tiles intersecting this bounding box.
     *
     * This method is similar to from_bbox(double, double, double, double) but
     * accepts coordinates in the coordinate system of your database,
     * transforms them to tile IDs and then calls
     * from_bbox(double, double, double, double). Please note that the
     * direction of the y axis is different between most projections and tile
     * IDs.
     *
     * \param min_lon x coordinate of lower left corner
     * \param min_lat y coordinate of lower left corner
     * \param max_lon x coordinate of upper right corner
     * \param max_lat y coordinate of upper right corner
     */
    void from_bbox_lon_lat(double min_x, double min_y, double max_x, double max_y);

    /**
     * Expire the tiles intersecting this bounding box. A buffer will be added to
     * the bounding box.
     *
     * \param min_lon x coordinate of upper left corner in tile IDs
     * \param min_lat y coordinate of upper left corner in tile IDs
     * \param max_lon x coordinate of lower right corner in tile IDs
     * \param max_lat y coordinate of lower right corner in tile IDs
     */
    void from_bbox(double min_x, double min_y, double max_x, double max_y);

    /**
     * Expire the tiles intersecting this bounding box. No buffer will be added.
     */
    void from_bbox_without_buffer(uint32_t min_x, uint32_t min_y,
                                  uint32_t max_x, uint32_t max_y);

    /**
     * Expire a line segment including a buffer.
     *
     * Input coordinates are the coordinates in the projection of the database.
     *
     * This method checks if the line segment crosses the 180th meridian and splits
     * it if necessary.
     *
     * \param lon_a longitude of start point
     * \param lat_a latitude of start point
     * \param lon_b longitude of end point
     * \param lat_b latitude of end point
     */
    void from_line_lon_lat(double lon_a, double lat_a, double lon_b,
                           double lat_b);

    /**
     * Expire all tiles a line segment intersects with including a small buffer.
     *
     * Coordinates (x and y) are in tile IDs (but double).
     * The start point must have a smaller or equal x index than the end point.
     *
     * The difference between x2 and x1 must be smaller than half of the
     * circumfence of the earth.
     *
     * \param x1 x index of the west end of the segment
     * \param y1 y index of the west end of the segment
     * \param x2 x index of the east end of the segment
     * \param y2 y index of the east end of the segment
     */
    void expire_line_segment(double x1, double y1, double x2, double y2);

    /**
     * Expire all tiles a line from (x1,y1) to (x2,y2) intersects. A buffer is
     * not included.
     *
     * Coordinates (x and y) are in tile IDs (but double).
     * The start point must have a smaller than the x index than the end point.
     *
     * \param x1 x index of the west end of the segment
     * \param y1 y index of the west end of the segment
     * \param x2 x index of the east end of the segment
     * \param y2 y index of the east end of the segment
     */
    void expire_line(double x1, double y1, double x2, double y2);

    /**
     * Expire a line segment which runs straight from south to north or runs
     * nearly in that direction. A buffer is not included.
     *
     * Coordinates (x and y) are in tile IDs (but double).
     * The start point must have a smaller y index than the end point.
     *
     * If a input coordinate is smaller than 0, 0 will be used instead.
     * If a input coordinate is larger than map_width, map_width will be used
     * instead.
     *
     * \param x x index of start and end point
     * \param y1 y index of the start point
     * \param y2 y index of the end point
     */
    void expire_vertical_line(double x, double y1, double y2);

    void from_wkb(const char* wkb, osmid_t osm_id);
    int from_db(table_t* table, osmid_t osm_id);

    /**
     * Write the list of expired tiles to a file.
     *
     * You will probably use tile_output_t as template argument for production code
     * and another class which does not write to a file for unit tests.
     *
     * \param filename name of the file
     * \param minzoom minimum zoom level
     */
    void output_and_destroy(const char *filename, uint32_t minzoom);

    /**
     * Output expired tiles on all requested zoom levels.
     *
     * \tparam TILE_WRITER class which implements the method
     * output_dirty_tile(uint32_t x, uint32_t y, uint32_t zoom) which usually writes the tile ID to a file
     * (production code) or does something else (usually unit tests)
     *
     * \param minzoom minimum zoom level
     */
    template <class TILE_WRITER>
    void output_and_destroy(TILE_WRITER &output_writer, uint32_t minzoom)
    {
        assert(minzoom <= maxzoom);
        // build a sorted vector of all expired tiles
        std::vector<uint64_t> tiles_maxzoom(m_dirty_tiles.begin(),
                                            m_dirty_tiles.end());
        std::sort(tiles_maxzoom.begin(), tiles_maxzoom.end());
        /* Loop over all requested zoom levels (from maximum down to the minimum zoom level).
         * Tile IDs of the tiles enclosing this tile at lower zoom levels are calculated using
         * bit shifts.
         *
         * last_quadkey is initialized with a value which is not expected to exist
         * (larger than largest possible quadkey). */
        uint64_t last_quadkey = 1ULL << (2 * maxzoom);
        for (std::vector<uint64_t>::const_iterator it = tiles_maxzoom.cbegin();
             it != tiles_maxzoom.cend(); ++it) {
            for (uint32_t dz = 0; dz <= maxzoom - minzoom; dz++) {
                // scale down to the current zoom level
                uint64_t qt_current = *it >> (dz * 2);
                /* If dz > 0, there are propably multiple elements whose quadkey
                 * is equal because they are all sub-tiles of the same tile at the current
                 * zoom level. We skip all of them after we have written the first sibling.
                 */
                if (qt_current == last_quadkey >> (dz * 2)) {
                    continue;
                }
                xy_coord_t xy = quadkey_to_xy(qt_current, maxzoom - dz);
                output_writer.output_dirty_tile(xy.x, xy.y, maxzoom - dz);
            }
            last_quadkey = *it;
        }
    }

    /**
    * merge the list of expired tiles in the other object into this
    * object, destroying the list in the other object.
    */
    void merge_and_destroy(expire_tiles &other);

    /**
     * Helper method to convert a tile ID (x and y) into a quadkey
     * using bitshifts.
     *
     * Quadkeys are interleaved this way: YXYX…
     *
     * \param x x index
     * \param y y index
     * \param zoom zoom level
     * \returns quadtree ID as integer
     */
    static uint64_t xy_to_quadkey(uint32_t x, uint32_t y, uint32_t zoom);

    /**
     * Convert a quadkey into a tile ID (x and y) using bitshifts.
     *
     * Quadkeys coordinates are interleaved this way: YXYX…
     *
     * \param quadkey the quadkey to be converted
     * \param zoom zoom level
     */
    static xy_coord_t quadkey_to_xy(uint64_t quadkey, uint32_t zoom);

private:
    /**
     *  How many tiles worth of space to leave either side of a changed feature
     **/
    static constexpr double TILE_EXPIRY_LEEWAY = 0.1;

    /**
     * Check if a coordinate (x or y) of a tile at maxzoom zoom level is valid.
     *
     * This method checks if the coordinate is within the bounds for tile IDs
     * at this zoom level and sets it to the bounds if it is too large.
     *
     * Because coord is unsigned, no check for coord < 0 is done.
     *
     * \coord coordinate to be checked
     *
     * \returns normalised coordinates
     */
    bool valid_tile_coord(uint32_t coord);

    /**
     * Normalise the coordinate (x or y) of a tile at maxzoom zoom level.
     *
     * This method checks if the coordinate is within the bounds for tile IDs
     * at this zoom level and sets it to the bounds if it is too large.
     *
     * \returns normalised coordinates
     */
    double normalise_tile_coord(double coord);

    /**
     * Expire a single tile.
     *
     * \param x x index of the tile to be expired.
     * \param y y index of the tile to be expired.
     */
    void expire_tile(uint32_t x, uint32_t y);

    void from_wkb_point(ewkb::parser_t *wkb);
    void from_wkb_line(ewkb::parser_t *wkb);
    void from_wkb_polygon(ewkb::parser_t *wkb, osmid_t osm_id);

    /**
     * Evaluate a segment of the outer ring of a polygon.
     *
     * Mark all tiles as expired which are inside the bounding box of this polygon.
     * If there are already tiles with the same x index which are expired,
     * all tiles between them and the tiles of the bounding box of this segment will be
     * expired, too.
     *
     * All parameters should be the coordinates of the output SRS. This method
     * will transform them to tile IDs.
     */
    void evaluate_segment(intersecting_tiles_t &tiles, double x1, double y1,
                          double x2, double y2);

    double tile_width;
    double max_bbox;
    uint32_t map_width;
    uint32_t maxzoom;
    std::shared_ptr<reprojection> projection;

    /**
     * x coordinate of the tile which has been added as last tile to the unordered set
     */
    uint32_t last_tile_x;

    /**
     * y coordinate of the tile which has been added as last tile to the unordered set
     */
    uint32_t last_tile_y;

    /**
     * manages which tiles have been marked as empty
     *
     * This set stores the IDs of the tiles at the maximum zoom level. We don't
     * store the IDs of the expired tiles of lower zoom levels. They are calculated
     * on the fly at the end.
     *
     * Tile IDs are converted into so-called quadkeys as used by Bing Maps.
     * https://msdn.microsoft.com/en-us/library/bb259689.aspx
     * A quadkey is generated by interleaving the x and y index in following order:
     * YXYX...
     *
     * Example:
     * x = 3 = 0b011, y = 5 = 0b101
     * results in the quadkey 0b100111.
     *
     * Bing Maps itself uses the quadkeys as a base-4 number converted to a string.
     * We interpret this IDs as simple 64-bit integers due to performance reasons.
     */
    std::unordered_set<uint64_t> m_dirty_tiles;
};

#endif
