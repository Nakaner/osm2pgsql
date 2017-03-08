/*
 * intersecting_tiles.hpp
 *
 *  Created on: 07.03.2017
 *      Author: michael
 */

#ifndef INTERSECTING_TILES_HPP_
#define INTERSECTING_TILES_HPP_

#include <cstdint>
#include <memory>
#include <vector>

using crossings_t = std::vector<uint32_t>;
using stripe_t = std::vector<crossings_t>;

/**
 * This class stores the information which tiles are inside a polygon.
 *
 * The heart of this class are two vectors. Each element of these vectors
 * represents one column in the tile grid (i.e. stripes from south to north).
 * The vectors themselves contain vectors with doubles.
 *
 * We check for each segment whether the interior of the polygon is above the
 * segment or below. If the interior is above, we add the maximum y value of
 * the bounding box of this segment to the maximums vector. If the interior of
 * the polygon is below the this segment, we add the minimum y value of the
 * bounding box of this segment to the minimums vector. That's the way it works
 * in short. For all edge cases see the comments in the source file.
 *
 * This class is also able to handle inner rings.
 *
 * After reading and handling all segments of all rings, the vectors have to be
 * sorted by calling sort_bounds(). After that you can iterate over the columns
 * and get intervals (a pair of minimum and maximum tile ID) which are inside
 * the polygon. Use move_to_next_column(), column_has_intervals() and
 * get_next_pair() for that purpose.
 */
class intersecting_tiles_t
{
    uint32_t offset_x;
    uint32_t max_tile_id;
    const double TILE_EXPIRY_LEEWAY;

    static constexpr double PI = 3.14159265358979323846;

    /**
     * vector containing the minimum y values of all columns
     */
    stripe_t min_bounds;
    /**
     * vector containing the maximum y values of all columns
     */
    stripe_t max_bounds;

    /**
     * x index of last column whose interval was queried using get_next_pair()
     */
    uint32_t current_x = 0;

    /**
     * index of the last minimum value returned via get_next_pair()
     */
    uint32_t next_idx_min = 0;

    /**
     * index of the last maximum value returned via get_next_pair()
     */
    uint32_t next_idx_max = 0;

    /**
     * Add a minimum for a tile column.
     *
     * TILE_EXPIRY_LEEWAY will be deduced from the min argument.
     *
     * \param x x coordinate of the stripe
     * \param min minimum y coordinate
     */
    void add_minimum(uint32_t x, double min);

    /**
     * Add a maximum for a tile column.
     *
     * TILE_EXPIRY_LEEWAY will be added to the max argument.
     *
     * \param x x coordinate of the stripe
     * \param max maximum y coordinate
     */
    void add_maximum(uint32_t x_index, double max);

    /**
     * Examinate add a maxium or a minimum to the bounds vectors.
     *
     * This method examinates on its own if it is a maximum or a minimum.
     * The direction of the segment (left to right or right to left) does not matter.
     *
     * \param x x index of the tile
     * \param y1 y index of one end of the segment
     * \param y2 y index of the other end of the segment
     */
    void add_minimum_or_maximum(uint32_t x, double y1, double y2, bool interior_above);

    /**
     * Get next minimum bound in the column pointed by current_x.
     *
     * \returns next minimum bound
     *
     * \throws std::runtime_error if no minimum is found.
     */
    uint32_t get_next_minimum();

    /**
     * Get next maximum bound in the column pointed by current_x.
     *
     * \returns next maximum bound
     *
     * \throws std::runtime_error if no maximum is found.
     */
    uint32_t get_next_maximum();

    /**
     * Helper method to check whether a returned tile ID is valid.
     *
     * \param value tile ID to be checked
     */
    bool is_valid_bounds_entry(uint32_t value);

public:
    intersecting_tiles_t(double x_min, double x_max, uint32_t map_width,
                         double leeway);

    intersecting_tiles_t() = delete;

    /**
     * Sort the bounds vector.
     */
    void sort_bounds();

    /**
     * Get the index in the x_stripes vector corresponding to this x coordinate (tile ID).
     */
    uint32_t get_x_index(uint32_t x);

    /**
     * Evaluate a segment of the outer ring.
     *
     * All coordinates are in tile units.
     *
     * \param x1 x coordinate of the start point
     * \param y1 y coordinate of the start point
     * \param x2 x coordinate of the end point
     * \param y2 y coordinate of the end point
     * \param outer_ring true if the segment belongs to an inner ring
     */
    void evaluate_segment(double x1, double y1, double x2, double y2,
                          bool outer_ring);

    /**
     * Get next pair of minimum and maximum.
     *
     * \returns a unique_ptr with minimum and maximum. An empty unique_ptr is
     * returned if either the no minimum or no maximum is found.
     */
    std::unique_ptr<std::pair<uint32_t, uint32_t>> get_next_pair();

    /**
     * Check if there are intervals left in the current column which have not been
     * queried yet.
     */
    bool column_has_intervals();

    /**
     * Move to the next tile column for output.
     *
     * This increases current_x.
     *
     * \returns False if there is no "next column" (end of bounds vector reached),
     * true otherwise.
     */
    bool move_to_next_column();

    /**
     * Get current column of the output.
     */
    uint32_t get_current_x();

    /**
     * Get the side of the ring segment where the interior of the polygon is.
     * Use this for outer rings only.
     *
     * If the angel is exactly ± PI/2, false will be returned.
     *
     * All coordinates are in tile units.
     *
     * \param x1 x coordinate of the start point
     * \param y1 y coordinate of the start point
     * \param x2 x coordinate of the end point
     * \param y2 y coordinate of the end point
     *
     * \returns true if the interior is above, false otherwise.
     */
    static bool interior_side_above(double x1, double y1, double x2, double y2);
};

#endif /* INTERSECTING_TILES_HPP_ */
