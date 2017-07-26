#include "intersecting_tiles.hpp"
#include <algorithm>
#include <assert.h>

constexpr double PI = 3.14159265358979323846;

intersecting_tiles_t::intersecting_tiles_t(double x_min, double x_max,
                                           uint32_t map_width, double leeway)
: offset_x(static_cast<uint32_t>(x_min - leeway)), max_tile_id(map_width),
  TILE_EXPIRY_LEEWAY(leeway),
  min_bounds(static_cast<uint32_t>(x_max + leeway) -
             static_cast<uint32_t>(x_min - leeway) + 1),
  max_bounds(static_cast<uint32_t>(x_max + leeway) -
             static_cast<uint32_t>(x_min - leeway) + 1)
{
}

uint32_t intersecting_tiles_t::get_x_index(uint32_t x) { return x - offset_x; }

void intersecting_tiles_t::add_minimum(uint32_t x, double min)
{
    min_bounds.at(get_x_index(x)).emplace_back(static_cast<uint32_t>(min -
                  TILE_EXPIRY_LEEWAY));
}

void intersecting_tiles_t::add_maximum(uint32_t x, double max)
{
    max_bounds.at(get_x_index(x)).emplace_back(static_cast<uint32_t>(max +
                  TILE_EXPIRY_LEEWAY));
}

void intersecting_tiles_t::evaluate_segment(double x1, double y1, double x2,
                                            double y2, bool outer_ring)
{
    // Segments which do not cross the border of two columns get a special
    // treatment. They would only introduce fake minimum or fake maximum
    // entries which lack a corresponding maximum/minimum entry to be a
    // proper pair. That's why we add both a minimum and a maximum value for
    // them. Being a inner or outer segment and its direction does not matter
    // here. Most buildings will fall in this category.
    if (static_cast<uint32_t>(std::min(x1, x2) - TILE_EXPIRY_LEEWAY)
        == static_cast<uint32_t>(std::max(x1, x2) + TILE_EXPIRY_LEEWAY)) {
        add_minimum(static_cast<uint32_t>(x1), std::min(y1, y2));
        add_maximum(static_cast<uint32_t>(x1), std::max(y1, y2));
        return;
    }
    // get side where the interior of the polygon is
    bool interior_above = interior_side_above(x1, y1, x2, y2);

    // swap (x1,y1) and (x2,y2) if to ensure that the segment runs from west
    // to east. It does not matter any more what the original direction was
    // and whether it was a inner or outer ring. The correct direction matters
    // for the for loop a few lines below.
    if (x2 < x1) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }
    // iterate over x stripes this segment crosses
    const uint32_t start = static_cast<uint32_t>(x1 - TILE_EXPIRY_LEEWAY);
    const uint32_t end = static_cast<uint32_t>(x2 + TILE_EXPIRY_LEEWAY);
    for (uint32_t x = start; x <= end; ++x) {
        add_minimum_or_maximum(x, y1, y2, interior_above);
        // If the segment crosses the left and the right border of the column,
        // we have to add the minimum/maximum twice.
        if (x != start && x != end) {
            add_minimum_or_maximum(x, y1, y2, interior_above);
        }
    }
}

void intersecting_tiles_t::add_minimum_or_maximum(uint32_t x, double y1, double y2,
                                                  bool interior_above)
{
    if (interior_above) {
        // The second argument is the maximum of y1 and y2, not the minimum,
        // because the y axis points downwards.
        add_maximum(x, std::max(y1, y2));
    } else {
        add_minimum(x, std::min(y1, y2));
    }
}

void intersecting_tiles_t::sort_bounds()
{
    // sort bounds first
    for (stripe_t::iterator it = min_bounds.begin(); it != min_bounds.end();
            ++it) {
        std::sort(it->begin(), it->end());
    }
    for (stripe_t::iterator it = max_bounds.begin(); it != max_bounds.end();
            ++it) {
        std::sort(it->begin(), it->end());
    }

    // look for overlapping intervals and remove them
    for (size_t column = 0; column < std::min(min_bounds.size(),
         max_bounds.size()); ++column) {
        for (size_t entry = 1; entry < std::min(min_bounds.at(column).size(),
             max_bounds.at(column).size()); ++entry) {
            if (min_bounds.at(column).at(entry) <=
                max_bounds.at(column).at(entry - 1)) {
                // interval overlapping
                min_bounds.at(column).at(entry) =
                        std::min(min_bounds.at(column).at(entry),
                                min_bounds.at(column).at(entry - 1));
                max_bounds.at(column).at(entry) =
                        std::max(max_bounds.at(column).at(entry),
                                max_bounds.at(column).at(entry - 1));
                // invalidate previous element in both vectors
                min_bounds.at(column).at(entry - 1) =
                        std::numeric_limits<uint32_t>::max();
                max_bounds.at(column).at(entry - 1) =
                        std::numeric_limits<uint32_t>::max();
            }
        }
    }

    // sort again
    for (stripe_t::iterator it = min_bounds.begin(); it != min_bounds.end();
         ++it) {
        std::sort(it->begin(), it->end());
    }
    for (stripe_t::iterator it = max_bounds.begin(); it != max_bounds.end();
         ++it) {
        std::sort(it->begin(), it->end());
    }
}

uint32_t intersecting_tiles_t::get_next_minimum()
{
    uint32_t minimum = min_bounds.at(current_x).at(next_idx_min);
    ++next_idx_min;
    return minimum;
}

uint32_t intersecting_tiles_t::get_next_maximum()
{
    uint32_t maximum = max_bounds.at(current_x).at(next_idx_max);
    ++next_idx_max;
    return maximum;
}

bool intersecting_tiles_t::is_valid_bounds_entry(uint32_t value)
{
    return value < max_tile_id;
}

bool intersecting_tiles_t::move_to_next_column()
{
    current_x++;
    next_idx_min = 0;
    next_idx_max = 0;
    return current_x < min_bounds.size() && current_x < max_bounds.size();
}

bool intersecting_tiles_t::column_has_intervals()
{
    return next_idx_min < min_bounds.at(current_x).size()
            && next_idx_max < max_bounds.at(current_x).size();
}

std::unique_ptr<std::pair<uint32_t, uint32_t>>
intersecting_tiles_t::get_next_pair()
{
    if (!column_has_intervals()) {
        // is only thrown if column_has_intervals() has not been called before
        throw std::runtime_error("The end of this tile column has been reached already.");
    }
    uint32_t minimum = get_next_minimum();
    uint32_t maximum = get_next_maximum();
    if (is_valid_bounds_entry(maximum) && is_valid_bounds_entry(minimum)) {
        return std::unique_ptr<std::pair<uint32_t, uint32_t>>(
            new std::pair<uint32_t, uint32_t>(minimum, maximum));
    } else {
        return std::unique_ptr<std::pair<uint32_t, uint32_t>>();
    }
}

uint32_t intersecting_tiles_t::get_current_x()
{
    return current_x + offset_x;
}

/*static*/ bool intersecting_tiles_t::interior_side_above(double x1, double y1,
                                                          double x2, double y2)
{
    // y1 and y2 must be swapped because our y axis points downwards.
    return (std::atan2(y1 - y2, x2 - x1) < PI / 2 &&
            std::atan2(y1 - y2, x2 - x1) > -PI / 2);
}
