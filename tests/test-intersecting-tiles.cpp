#include "intersecting_tiles.hpp"
#include <boost/format.hpp>
#include <set>

/**
 * This represents a y interval. The first element is the x index, the second
 * is the minimum y and the third the maximum y value.
 */
using tile_interval_t = std::array<uint32_t, 3>;

namespace {

void run_test(const char *test_name, void (*testfunc)())
{
    try {
        fprintf(stderr, "%s\n", test_name);
        testfunc();
    } catch (const std::exception &e) {
        fprintf(stderr, "%s\n", e.what());
        fprintf(stderr, "FAIL\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "PASS\n");
}
#define RUN_TEST(x) run_test(#x, &(x))
#define ASSERT_EQ(a, b)                                                        \
    {                                                                          \
        if (!((a) == (b))) {                                                   \
            throw std::runtime_error(                                          \
                (boost::format("Expecting %1% == %2%, but %3% != %4%") % #a %  \
                 #b % (a) % (b))                                               \
                    .str());                                                   \
        }                                                                      \
    }

void add_intervals_to_result_set(std::set<tile_interval_t> &results,
                                 intersecting_tiles_t &tiles)
{
    tiles.sort_bounds();
    do {
        while (tiles.column_has_intervals()) {
            std::unique_ptr<std::pair<uint32_t, uint32_t>> interval =
                tiles.get_next_pair();
            if (!interval) {
                continue;
            }
            results.insert(tile_interval_t{tiles.get_current_x(),
                                           interval->first, interval->second});
        }
    } while (tiles.move_to_next_column());
}

void test_intersecting_tiles()
{
    // rectangle, lower left: 2.4, 1.6; upper right: 2.6, 1.4
    intersecting_tiles_t tiles(2.4, 2.6, 4, 0.1);
    tiles.evaluate_segment(2.4, 1.6, 2.6, 1.6, true);
    tiles.evaluate_segment(2.6, 1.6, 2.6, 1.4, true);
    tiles.evaluate_segment(2.6, 1.4, 2.4, 1.4, true);
    tiles.evaluate_segment(2.4, 1.4, 2.4, 1.6, true);
    std::set<tile_interval_t> minmaxs;

    add_intervals_to_result_set(minmaxs, tiles);

    ASSERT_EQ(minmaxs.size(), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 1, 1}), 1);
}

void test_intersecting_tiles_two_columns()
{
    // rectangle, lower left: 2.4, 1.6; upper right: 2.6, 1.4
    intersecting_tiles_t tiles(2.4, 3.6, 4, 0.1);
    tiles.evaluate_segment(2.4, 1.6, 3.6, 1.6, true);
    tiles.evaluate_segment(3.6, 1.6, 3.6, 1.4, true);
    tiles.evaluate_segment(3.6, 1.4, 2.4, 1.4, true);
    tiles.evaluate_segment(2.4, 1.4, 2.4, 1.6, true);
    std::set<tile_interval_t> minmaxs;

    add_intervals_to_result_set(minmaxs, tiles);

    ASSERT_EQ(minmaxs.size(), 2);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 1, 1}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 1, 1}), 1);
}

void test_intersecting_tiles_more_columns()
{
    // rectangle, lower left: 2.4, 1.6; upper right: 2.6, 1.4
    intersecting_tiles_t tiles(2.5, 8.0, 16, 0.1);
    tiles.evaluate_segment(2.5, 4.8, 3.3, 6.0, true);
    tiles.evaluate_segment(3.3, 6.0, 6.8, 5.6, true);
    tiles.evaluate_segment(6.8, 5.6, 8.0, 2.6, true);
    tiles.evaluate_segment(8.0, 2.6, 6.6, 1.7, true);
    tiles.evaluate_segment(6.6, 1.7, 6.8, 3.5, true);
    tiles.evaluate_segment(6.8, 3.5, 3.8, 5.2, true);
    tiles.evaluate_segment(3.8, 5.2, 3.4, 1.8, true);
    tiles.evaluate_segment(3.4, 1.8, 2.5, 4.8, true);
    std::set<tile_interval_t> minmaxs;

    add_intervals_to_result_set(minmaxs, tiles);

    ASSERT_EQ(minmaxs.size(), 7);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 1, 6}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 1, 6}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{4, 3, 6}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{5, 3, 6}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{6, 1, 6}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{7, 1, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{8, 1, 5}), 1);
}

void test_intersecting_tiles_u_shape()
{
    // rectangle, lower left: 2.4, 1.6; upper right: 2.6, 1.4
    intersecting_tiles_t tiles(1.3, 5.7, 8, 0.1);
    tiles.evaluate_segment(1.3, 3.7, 2.5, 5.6, true);
    tiles.evaluate_segment(2.5, 5.6, 5.5, 4.5, true);
    tiles.evaluate_segment(5.5, 4.5, 5.3, 4.2, true);
    tiles.evaluate_segment(5.3, 4.2, 2.7, 4.7, true);
    tiles.evaluate_segment(2.7, 4.7, 2.2, 1.6, true);
    tiles.evaluate_segment(2.2, 1.6, 5.7, 0.9, true);
    tiles.evaluate_segment(5.7, 0.9, 5.6, 0.4, true);
    tiles.evaluate_segment(5.6, 0.4, 1.8, 1.4, true);
    tiles.evaluate_segment(1.8, 1.4, 1.3, 3.7, true);
    std::set<tile_interval_t> minmaxs;

    add_intervals_to_result_set(minmaxs, tiles);

    ASSERT_EQ(minmaxs.size(), 8);
    ASSERT_EQ(minmaxs.count(tile_interval_t{1, 0, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 0, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 0, 1}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 4, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{4, 0, 1}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{4, 4, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{5, 0, 1}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{5, 4, 5}), 1);
}

void test_intersecting_tiles_inner_ring()
{
    // rectangle, lower left: 2.4, 1.6; upper right: 2.6, 1.4
    intersecting_tiles_t tiles(0.6, 5.8, 8, 0.1);
    tiles.evaluate_segment(0.6, 0.3, 1.6, 5.2, true);
    tiles.evaluate_segment(1.6, 5.2, 5.5, 4.7, true);
    tiles.evaluate_segment(5.5, 4.7, 5.8, 0.2, true);
    tiles.evaluate_segment(5.8, 0.2, 0.6, 0.3, true);
    tiles.evaluate_segment(1.5, 0.7, 5.4, 0.7, false);
    tiles.evaluate_segment(5.4, 0.7, 5.3, 4.3, false);
    tiles.evaluate_segment(5.3, 4.3, 1.8, 4.2, false);
    tiles.evaluate_segment(1.8, 4.2, 1.5, 0.7, false);
    std::set<tile_interval_t> minmaxs;

    add_intervals_to_result_set(minmaxs, tiles);

    ASSERT_EQ(minmaxs.size(), 9);
    ASSERT_EQ(minmaxs.count(tile_interval_t{0, 0, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{1, 0, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 0, 0}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{2, 4, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 0, 0}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{3, 4, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{4, 0, 0}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{4, 4, 5}), 1);
    ASSERT_EQ(minmaxs.count(tile_interval_t{5, 0, 5}), 1);
}
} // anonymous namespace

int main(int argc, char *argv[])
{
    //try each test if any fail we will exit
    RUN_TEST(test_intersecting_tiles);
    RUN_TEST(test_intersecting_tiles_two_columns);
    RUN_TEST(test_intersecting_tiles_more_columns);
    RUN_TEST(test_intersecting_tiles_u_shape);
    RUN_TEST(test_intersecting_tiles_inner_ring);
    //passed
    return 0;
}
