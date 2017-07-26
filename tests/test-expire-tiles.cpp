#include "expire-tiles.hpp"
#include "options.hpp"
#include "wkb.hpp"

#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdexcept>
#include <boost/format.hpp>
#include <set>

#define EARTH_CIRCUMFERENCE (40075016.68)

namespace {

void run_test(const char* test_name, void (*testfunc)())
{
    try
    {
        fprintf(stderr, "%s\n", test_name);
        testfunc();
    }
    catch(const std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        fprintf(stderr, "FAIL\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "PASS\n");
}
#define RUN_TEST(x) run_test(#x, &(x))
#define ASSERT_EQ(a, b) { if (!((a) == (b))) { throw std::runtime_error((boost::format("Expecting %1% == %2%, but %3% != %4%") % #a % #b % (a) % (b)).str()); } }

struct xyz {
    uint32_t z;
    uint32_t x, y;
    xyz(uint32_t z_, uint32_t x_, uint32_t y_) : z(z_), x(x_), y(y_) {}
    bool operator==(const xyz &other) const
    {
        return ((z == other.z) && (x == other.x) && (y == other.y));
  }
  bool operator<(const xyz &other) const {
    return ((z < other.z) ||
            ((z == other.z) &&
             ((x < other.x) ||
              ((x == other.x) &&
               (y < other.y)))));
  }
  void to_bbox(double &x0, double &y0,
               double &x1, double &y1) const {
    const double datum = 0.5 * (1 << z);
    const double scale = EARTH_CIRCUMFERENCE / (1 << z);
    x0 = (x - datum) * scale;
    y0 = (datum - (y + 1)) * scale;
    x1 = ((x + 1) - datum) * scale;
    y1 = (datum - y) * scale;
  }

  /**
     * Convert the x and y coordinate of the upper left corner of a tile into
     * the coordinate of the center of the tile.
     */
  void to_centroid(double &x0, double &y0) const
  {
      x0 = x + 0.5;
      y0 = y + 0.5;
  }
};

static std::shared_ptr<reprojection> defproj(reprojection::create_projection(PROJ_SPHERE_MERC));
static std::shared_ptr<reprojection> latlonproj(reprojection::create_projection(PROJ_LATLONG));

std::ostream &operator<<(std::ostream &out, const xyz &tile) {
  out << tile.z << "/" << tile.x << "/" << tile.y;
  return out;
}

struct tile_output_set
{
    tile_output_set(uint32_t min) : min_zoom(min) {}

    ~tile_output_set() = default;

    void output_dirty_tile(uint32_t x, uint32_t y, uint32_t zoom)
    {
        m_tiles.insert(xyz(zoom, x, y));
    }

    bool contains(xyz other) { return m_tiles.count(other) == 1; }

    void print()
    {
        for (xyz element : m_tiles) {
            fprintf(stderr, "%u %u %u\n", element.z, element.x, element.y);
        }
    }

    std::set<xyz> m_tiles;
    uint32_t min_zoom;
};

void test_xy_to_quadkey_z3()
{
    uint64_t quadkey_expected = 0x27;
    uint64_t quadkey2 = expire_tiles::xy_to_quadkey(3, 5, 3);
    ASSERT_EQ(quadkey2, quadkey_expected);
    xy_coord_t xy = expire_tiles::quadkey_to_xy(quadkey_expected, 3);
    ASSERT_EQ(xy.x, 3);
    ASSERT_EQ(xy.y, 5);
}

void test_xy_to_quadkey_z16()
{
    uint64_t quadkey_expected = 0xffffffff;
    uint64_t quadkey2 = expire_tiles::xy_to_quadkey(65535, 65535, 16);
    ASSERT_EQ(quadkey2, quadkey_expected);
    xy_coord_t xy = expire_tiles::quadkey_to_xy(quadkey_expected, 16);
    ASSERT_EQ(xy.x, 65535);
    ASSERT_EQ(xy.y, 65535);
}

/**
 * This test prevents problems which occur if 32-bit integers are used
 * instead of 64-bit integers.
 */
void test_xy_to_quadkey_z18()
{
    uint64_t quadkey_expected = 0xfffffffff;
    uint64_t quadkey2 = expire_tiles::xy_to_quadkey(262143, 262143, 18);
    ASSERT_EQ(quadkey2, quadkey_expected);
    xy_coord_t xy = expire_tiles::quadkey_to_xy(quadkey_expected, 18);
    ASSERT_EQ(xy.x, 262143);
    ASSERT_EQ(xy.y, 262143);
    quadkey_expected = 0x3fffffff0;
    quadkey2 = expire_tiles::xy_to_quadkey(131068, 131068, 18);
    ASSERT_EQ(quadkey2, quadkey_expected);
    xy = expire_tiles::quadkey_to_xy(quadkey_expected, 18);
    ASSERT_EQ(xy.x, 131068);
    ASSERT_EQ(xy.y, 131068);
}

void test_expire_simple_z1() {
    uint32_t minzoom = 1;
    expire_tiles et(minzoom, 20000, defproj);
    tile_output_set set(minzoom);

    // as big a bbox as possible at the origin to dirty all four
    // quadrants of the world.
    et.from_bbox_lon_lat(-10000, -10000, 10000, 10000);
    et.output_and_destroy<tile_output_set>(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 4);
    std::set<xyz>::iterator itr = set.m_tiles.begin();
    ASSERT_EQ(*itr, xyz(1, 0, 0));
    ++itr;
    ASSERT_EQ(*itr, xyz(1, 0, 1));
    ++itr;
    ASSERT_EQ(*itr, xyz(1, 1, 0));
    ++itr;
    ASSERT_EQ(*itr, xyz(1, 1, 1));
    ++itr;
}

void test_expire_simple_z3() {
    uint32_t minzoom = 3;
    expire_tiles et(minzoom, 20000, defproj);
    tile_output_set set(minzoom);

    // as big a bbox as possible at the origin to dirty all four
    // quadrants of the world.
    et.from_bbox_lon_lat(-10000, -10000, 10000, 10000);
    et.output_and_destroy<tile_output_set>(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 4);
    std::set<xyz>::iterator itr = set.m_tiles.begin();
    ASSERT_EQ(*itr, xyz(3, 3, 3));
    ++itr;
    ASSERT_EQ(*itr, xyz(3, 3, 4));
    ++itr;
    ASSERT_EQ(*itr, xyz(3, 4, 3));
    ++itr;
    ASSERT_EQ(*itr, xyz(3, 4, 4));
    ++itr;
}

void test_expire_simple_z18() {
    uint32_t minzoom = 18;
    expire_tiles et(18, 20000, defproj);
    tile_output_set set(minzoom);

    // dirty a smaller bbox this time, as at z18 the scale is
    // pretty small.
    et.from_bbox_lon_lat(-1, -1, 1, 1);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 4);
    std::set<xyz>::iterator itr = set.m_tiles.begin();
    ASSERT_EQ(*itr, xyz(18, 131071, 131071));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131071, 131072));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131072, 131071));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131072, 131072));
    ++itr;
}

/**
 * Test tile expiry on two zoom levels.
 */
void test_expire_simple_z17_18()
{
    uint32_t minzoom = 17;
    expire_tiles et(18, 20000, defproj);
    tile_output_set set(minzoom);

    // dirty a smaller bbox this time, as at z18 the scale is
    // pretty small.
    et.from_bbox_lon_lat(-1, -1, 1, 1);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 8);
    std::set<xyz>::iterator itr = set.m_tiles.begin();
    ASSERT_EQ(*itr, xyz(17, 65535, 65535));
    ++itr;
    ASSERT_EQ(*itr, xyz(17, 65535, 65536));
    ++itr;
    ASSERT_EQ(*itr, xyz(17, 65536, 65535));
    ++itr;
    ASSERT_EQ(*itr, xyz(17, 65536, 65536));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131071, 131071));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131071, 131072));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131072, 131071));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131072, 131072));
    ++itr;
}

/**
 * Similar to test_expire_simple_z17_18 but now all z18 tiles are children
 * of the same z17 tile.
 */
void test_expire_simple_z17_18_one_superior_tile()
{
    uint32_t minzoom = 17;
    expire_tiles et(18, 20000, defproj);
    tile_output_set set(minzoom);

    et.from_bbox_lon_lat(-163, 140, -140, 164);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 5);
    std::set<xyz>::iterator itr = set.m_tiles.begin();
    ASSERT_EQ(*itr, xyz(17, 65535, 65535));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131070, 131070));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131070, 131071));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131071, 131070));
    ++itr;
    ASSERT_EQ(*itr, xyz(18, 131071, 131071));
    ++itr;
}

/**
 * Test expire_line() method on zoom level 12.
 *
 * This test might probably fail if one of the tests of expire_tile() fails
 * because expire_line() calls expire_tiles().
 */
void test_expire_line_z12()
{
    uint32_t minzoom = 12;
    expire_tiles et(12, 0.1, latlonproj);
    et.expire_line(2116.3, 1416.3, 2118.5, 1417.5);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 4);
    ASSERT_EQ(set.contains(xyz(12, 2116, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1417)), true);
    ASSERT_EQ(set.contains(xyz(12, 2118, 1417)), true);
}

/**
 * Test expire_line() method on zoom level 12.
 *
 * This test might probably fail if one of the tests of expire_tile() fails
 * because expire_line() calls expire_tiles().
 */
void test_expire_line_z12_long()
{
    uint32_t minzoom = 12;
    expire_tiles et(12, 0.1, latlonproj);
    et.expire_line(2116.3, 1416.3, 2119.3, 1419.6);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 7);
    ASSERT_EQ(set.contains(xyz(12, 2116, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2116, 1417)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1417)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1418)), true);
    ASSERT_EQ(set.contains(xyz(12, 2118, 1418)), true);
    ASSERT_EQ(set.contains(xyz(12, 2118, 1419)), true);
    ASSERT_EQ(set.contains(xyz(12, 2119, 1419)), true);
}

/**
 * Test expire_line() method on zoom level 12 with a horizontal line.
 *
 * This test might probably fail if one of the tests of expire_tile() fails
 * because expire_line() calls expire_tiles().
 */
void test_expire_line_z12_horizontal()
{
    uint32_t minzoom = 12;
    expire_tiles et(12, 0.1, latlonproj);
    et.expire_line(2116.3, 1416.3, 2119.3, 1416.3);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 4);
    ASSERT_EQ(set.contains(xyz(12, 2116, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2118, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2119, 1416)), true);
}

/**
 * Test expire_line() method on zoom level 19.
 *
 * This test might probably fail if one of the tests of expire_tile() fails
 * because expire_line() calls expire_tiles().
 */
void test_expire_line_z19()
{
    uint32_t minzoom = 19;
    expire_tiles et(19, 0.1, latlonproj);
    et.expire_line(274374.3, 180067.5, 274376.5, 180066.3);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 4);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180067)), true);
    ASSERT_EQ(set.contains(xyz(19, 274375, 180067)), true);
    ASSERT_EQ(set.contains(xyz(19, 274375, 180066)), true);
    ASSERT_EQ(set.contains(xyz(19, 274376, 180066)), true);
}

/**
 * Test expire_line() method on zoom level 12 with a horizontal line.
 *
 * This test might probably fail if one of the tests of expire_tile() fails
 * because expire_line() calls expire_tiles().
 */
void test_expire_vertical_line_z19()
{
    uint32_t minzoom = 19;
    expire_tiles et(19, 0.1, latlonproj);
    et.expire_vertical_line(274374.3, 180063.3, 180067.5);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 5);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180067)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180066)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180065)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180064)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180063)), true);
}

void expire_line_segment()
{
    uint32_t minzoom = 12;
    expire_tiles et(12, 0.1, latlonproj);
    et.expire_line_segment(2116.3, 1416.3, 2118.5, 1417.5);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 4);
    ASSERT_EQ(set.contains(xyz(12, 2116, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1416)), true);
    ASSERT_EQ(set.contains(xyz(12, 2117, 1417)), true);
    ASSERT_EQ(set.contains(xyz(12, 2118, 1417)), true);
}

void expire_line_segment_vertical()
{
    uint32_t minzoom = 19;
    expire_tiles et(19, 0.1, latlonproj);
    et.expire_line_segment(274374.3, 180063.3, 274374.3, 180067.5);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 5);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180067)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180066)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180065)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180064)), true);
    ASSERT_EQ(set.contains(xyz(19, 274374, 180063)), true);
}

/**
 * Test from_line_lon_lat() with a line segment crossing the antimeridian.
 */
void test_from_line_lon_lat_crossing()
{
    uint32_t minzoom = 8;
    expire_tiles et(8, 0.1, latlonproj);
    et.from_line_lon_lat(179.1332, -16.4748, -179.1969, -17.7244);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 3);
    ASSERT_EQ(set.contains(xyz(8, 0, 140)), true);
    ASSERT_EQ(set.contains(xyz(8, 255, 139)), true);
    ASSERT_EQ(set.contains(xyz(8, 255, 140)), true);
}

/**
 * Test from_line_lon_lat() with a line segment whose ends have to be swapped.
 */
void test_from_line_lon_lat_wrong_order()
{
    uint32_t minzoom = 6;
    expire_tiles et(6, 0.1, latlonproj);
    et.from_line_lon_lat(86.3316, 34.9294, 78.1798, 28.6021);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);
    ASSERT_EQ(set.m_tiles.size(), 4);
    ASSERT_EQ(set.contains(xyz(6, 45, 26)), true);
    ASSERT_EQ(set.contains(xyz(6, 46, 26)), true);
    ASSERT_EQ(set.contains(xyz(6, 46, 25)), true);
    ASSERT_EQ(set.contains(xyz(6, 47, 25)), true);
}

/**
 * Test expire_from_wkb_polygon() method with a polygon with six corners, three
 * are rectangular.
 */
void test_expire_from_wkb_polygon_no_inner_z16()
{
    uint32_t minzoom = 16;
    // OSM way #8048087
    std::string wkb =
        "0103000020110F0000010000001100000057A95F38907B2C41F0E8C3BAC4F757416FFB"
        "127BCD7C2C41C5C89D9690F7574169A4E922A97D2C41D60AE4206CF757417DAF8A811B"
        "7E2C4194FD742F5AF757411E040521637E2C41740DE18952F757410090DFBA157F2C41"
        "63A2768D46F7574190CCAAF1937F2C411EB84E9D43F75741C0568F1B03802C411CFC02"
        "4842F757413D245E2F44802C41CE32C1AE41F75741393BFAB143802C419D7AC7A944F7"
        "574143597ADA5C802C415E993D7847F7574106C6AFB484802C418B9A952B48F757418A"
        "8E4FA624812C410A5C6C22D2F75741B7A1C90DBC802C41A614557AE0F75741FCD0140D"
        "7A7F2C418BD7FF0213F85741E979BB80EF7E2C41B56C66EF24F8574157A95F38907B2C"
        "41F0E8C3BAC4F75741";
    std::string binwkb = ewkb::parser_t::wkb_from_hex(wkb);
    expire_tiles et(16, 20000, defproj);
    et.from_wkb(binwkb.c_str(), 1);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 6);
    ASSERT_EQ(set.contains(xyz(16, 34294, 22492)), true);
    ASSERT_EQ(set.contains(xyz(16, 34294, 22493)), true);
    ASSERT_EQ(set.contains(xyz(16, 34294, 22494)), true);
    ASSERT_EQ(set.contains(xyz(16, 34295, 22492)), true);
    ASSERT_EQ(set.contains(xyz(16, 34295, 22493)), true);
    ASSERT_EQ(set.contains(xyz(16, 34295, 22494)), true);
}

/**
 * Test expire_from_wkb_polygon() method with a polygon with six corners, three
 * are rectangular.
 */
void test_expire_from_wkb_polygon_no_inner_z12()
{
    uint32_t minzoom = 12;
    // OSM way #8048087
    std::string wkb =
        "0103000020110F0000010000001100000057A95F38907B2C41F0E8C3BAC4F757416FF"
        "B127BCD7C2C41C5C89D9690F7574169A4E922A97D2C41D60AE4206CF757417DAF8A81"
        "1B7E2C4194FD742F5AF757411E040521637E2C41740DE18952F757410090DFBA157F2"
        "C4163A2768D46F7574190CCAAF1937F2C411EB84E9D43F75741C0568F1B03802C411C"
        "FC024842F757413D245E2F44802C41CE32C1AE41F75741393BFAB143802C419D7AC7A"
        "944F7574143597ADA5C802C415E993D7847F7574106C6AFB484802C418B9A952B48F7"
        "57418A8E4FA624812C410A5C6C22D2F75741B7A1C90DBC802C41A614557AE0F75741F"
        "CD0140D7A7F2C418BD7FF0213F85741E979BB80EF7E2C41B56C66EF24F8574157A95F"
        "38907B2C41F0E8C3BAC4F75741";
    std::string binwkb = ewkb::parser_t::wkb_from_hex(wkb);
    expire_tiles et(12, 20000, defproj);
    et.from_wkb(binwkb.c_str(), 1);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 1);
    ASSERT_EQ(set.contains(xyz(12, 2143, 1405)), true);
}

/**
 * Test expire_from_wkb_polygon() method with a polygon with one outer ring
 * like in test_expire_from_wkb_polygon_no_inner() but also one inner ring
 * with four corners.
 *
 * The EWKB string was created by running following PostGIS query:
 * ```sql
 * SELECT ST_GeomFromText('POLYGON ((13.4989 52.3512, 13.5727 52.3512, 13.5727
 * 52.3836, 13.4989 52.3836, 13.49666 52.36135, 13.48731 52.35558, 13.4989
 * 52.3512), (13.5053 52.3563, 13.5053 52.3811, 13.5679 52.3811, 13.5679
 * 52.3563, 13.5053 52.3563))', 4326);
 * ```
 */
void test_expire_from_wkb_polygon_with_inner()
{
    uint32_t minzoom = 14;
    std::string wkb =
        "0103000020E61000000200000007000000E5F21FD26FFF2A40772D211FF42C4A40A1D"
        "634EF38252B40772D211FF42C4A40A1D634EF38252B40705F07CE19314A40E5F21FD2"
        "6FFF2A40705F07CE19314A40757632384AFE2A406E3480B7402E4A400C0742B280F92"
        "A404C8E3BA5832D4A40E5F21FD26FFF2A40772D211FF42C4A40050000006D567DAEB6"
        "022B404BC8073D9B2D4A406D567DAEB6022B40B84082E2C7304A40FBCBEEC9C3222B4"
        "0B84082E2C7304A40FBCBEEC9C3222B404BC8073D9B2D4A406D567DAEB6022B404BC8"
        "073D9B2D4A40";
    std::string binwkb = ewkb::parser_t::wkb_from_hex(wkb);
    expire_tiles et(14, 20000, latlonproj);
    et.from_wkb(binwkb.c_str(), 1);
    tile_output_set set(minzoom);
    et.output_and_destroy(set, minzoom);

    ASSERT_EQ(set.m_tiles.size(), 12);
    ASSERT_EQ(set.contains(xyz(14, 8805, 5384)), true);
    ASSERT_EQ(set.contains(xyz(14, 8805, 5385)), true);
    ASSERT_EQ(set.contains(xyz(14, 8806, 5383)), true);
    ASSERT_EQ(set.contains(xyz(14, 8806, 5384)), true);
    ASSERT_EQ(set.contains(xyz(14, 8806, 5385)), true);
    ASSERT_EQ(set.contains(xyz(14, 8807, 5383)), true);
    ASSERT_EQ(set.contains(xyz(14, 8807, 5385)), true);
    ASSERT_EQ(set.contains(xyz(14, 8808, 5383)), true);
    ASSERT_EQ(set.contains(xyz(14, 8808, 5385)), true);
    ASSERT_EQ(set.contains(xyz(14, 8809, 5383)), true);
    ASSERT_EQ(set.contains(xyz(14, 8809, 5384)), true);
    ASSERT_EQ(set.contains(xyz(14, 8809, 5385)), true);
}

std::set<xyz> generate_random(uint32_t zoom, size_t count)
{
    size_t num = 0;
    std::set<xyz> set;
    const uint32_t coord_mask = (1U << zoom) - 1;

    while (num < count) {
        xyz item(zoom, static_cast<uint32_t>(rand()) & coord_mask,
                 static_cast<uint32_t>(rand()) & coord_mask);
        if (set.count(item) == 0) {
            set.insert(item);
            ++num;
        }
    }
    return set;
}

void assert_tilesets_equal(const std::set<xyz> &a,
                           const std::set<xyz> &b) {
  ASSERT_EQ(a.size(), b.size());
  std::set<xyz>::const_iterator a_itr = a.begin();
  std::set<xyz>::const_iterator b_itr = b.begin();
  while ((a_itr != a.end()) &&
         (b_itr != b.end())) {
    ASSERT_EQ(*a_itr, *b_itr);
    ++a_itr;
    ++b_itr;
  }
}

void expire_centroids(const std::set<xyz> &check_set,
                      expire_tiles &et) {
  for (std::set<xyz>::const_iterator itr = check_set.begin();
       itr != check_set.end(); ++itr) {
    double x0 = 0.0, y0 = 0.0;
    itr->to_centroid(x0, y0);
    et.from_bbox(x0, y0, x0, y0);
  }
}

// tests that expiring a set of tile centroids means that
// those tiles get expired.
void test_expire_set() {
    uint32_t zoom = 18;
    for (int i = 0; i < 100; ++i) {
        expire_tiles et(zoom, 20000, defproj);
        tile_output_set set(zoom);

        std::set<xyz> check_set = generate_random(zoom, 100);
        expire_centroids(check_set, et);

        et.output_and_destroy(set, zoom);

        assert_tilesets_equal(set.m_tiles, check_set);
  }
}

// this tests that, after expiring a random set of tiles
// in one expire_tiles object and a different set in
// another, when they are merged together they are the
// same as if the union of the sets of tiles had been
// expired.
void test_expire_merge() {
    uint32_t zoom = 18;

    for (int i = 0; i < 100; ++i) {
        expire_tiles et(zoom, 20000, defproj);
        expire_tiles et1(zoom, 20000, defproj);
        expire_tiles et2(zoom, 20000, defproj);
        tile_output_set set(zoom);

        std::set<xyz> check_set1 = generate_random(zoom, 100);
        expire_centroids(check_set1, et1);

        std::set<xyz> check_set2 = generate_random(zoom, 100);
        expire_centroids(check_set2, et2);

        et.merge_and_destroy(et1);
        et.merge_and_destroy(et2);

        std::set<xyz> check_set;
        std::set_union(check_set1.begin(), check_set1.end(), check_set2.begin(),
                       check_set2.end(),
                       std::inserter(check_set, check_set.end()));

        et.output_and_destroy(set, zoom);

        assert_tilesets_equal(set.m_tiles, check_set);
    }
}

// tests that merging two identical sets results in
// the same set. this guarantees that we check some
// pathways of the merging which possibly could be
// skipped by the random tile set in the previous
// test.
void test_expire_merge_same() {
    uint32_t zoom = 18;

    for (int i = 0; i < 100; ++i) {
        expire_tiles et(zoom, 20000, defproj);
        expire_tiles et1(zoom, 20000, defproj);
        expire_tiles et2(zoom, 20000, defproj);
        tile_output_set set(zoom);

        std::set<xyz> check_set = generate_random(zoom, 100);
        expire_centroids(check_set, et1);
        expire_centroids(check_set, et2);

        et.merge_and_destroy(et1);
        et.merge_and_destroy(et2);

        et.output_and_destroy(set, zoom);

        assert_tilesets_equal(set.m_tiles, check_set);
  }
}

// makes sure that we're testing the case where some
// tiles are in both.
void test_expire_merge_overlap() {
    uint32_t zoom = 18;

    for (int i = 0; i < 100; ++i) {
        expire_tiles et(zoom, 20000, defproj);
        expire_tiles et1(zoom, 20000, defproj);
        expire_tiles et2(zoom, 20000, defproj);
        tile_output_set set(zoom);

        std::set<xyz> check_set1 = generate_random(zoom, 100);
        expire_centroids(check_set1, et1);

        std::set<xyz> check_set2 = generate_random(zoom, 100);
        expire_centroids(check_set2, et2);

        std::set<xyz> check_set3 = generate_random(zoom, 100);
        expire_centroids(check_set3, et1);
        expire_centroids(check_set3, et2);

        et.merge_and_destroy(et1);
        et.merge_and_destroy(et2);

        std::set<xyz> check_set;
        std::set_union(check_set1.begin(), check_set1.end(), check_set2.begin(),
                       check_set2.end(),
                       std::inserter(check_set, check_set.end()));
        std::set_union(check_set1.begin(), check_set1.end(), check_set3.begin(),
                       check_set3.end(),
                       std::inserter(check_set, check_set.end()));

        et.output_and_destroy(set, zoom);

        assert_tilesets_equal(set.m_tiles, check_set);
  }
}

// checks that the set union still works when we expire
// large contiguous areas of tiles (i.e: ensure that we
// handle the "complete" flag correctly).
void test_expire_merge_complete() {
    uint32_t zoom = 18;

    for (int i = 0; i < 100; ++i) {
        expire_tiles et(zoom, 20000, defproj);
        expire_tiles et0(zoom, 20000, defproj);
        expire_tiles et1(zoom, 20000, defproj);
        expire_tiles et2(zoom, 20000, defproj);
        tile_output_set set(zoom);
        tile_output_set set0(zoom);

        // et1&2 are two halves of et0's box
        et0.from_bbox_lon_lat(-10000, -10000, 10000, 10000);
        et1.from_bbox_lon_lat(-10000, -10000, 0, 10000);
        et2.from_bbox_lon_lat(0, -10000, 10000, 10000);

        et.merge_and_destroy(et1);
        et.merge_and_destroy(et2);

        et.output_and_destroy(set, zoom);
        et0.output_and_destroy(set0, zoom);

        assert_tilesets_equal(set.m_tiles, set0.m_tiles);
  }
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    srand(0);

    //try each test if any fail we will exit
    RUN_TEST(test_xy_to_quadkey_z3);
    RUN_TEST(test_xy_to_quadkey_z16);
    RUN_TEST(test_xy_to_quadkey_z18);
    RUN_TEST(test_expire_simple_z1);
    RUN_TEST(test_expire_simple_z3);
    RUN_TEST(test_expire_simple_z18);
    RUN_TEST(test_expire_simple_z17_18);
    RUN_TEST(test_expire_simple_z17_18_one_superior_tile);
    RUN_TEST(test_expire_line_z12);
    RUN_TEST(test_expire_line_z12_long);
    RUN_TEST(test_expire_line_z12_horizontal);
    RUN_TEST(test_expire_line_z19);
    RUN_TEST(test_expire_vertical_line_z19);
    RUN_TEST(test_expire_from_wkb_polygon_no_inner_z16);
    RUN_TEST(test_expire_from_wkb_polygon_no_inner_z12);
    RUN_TEST(test_expire_from_wkb_polygon_with_inner);
    RUN_TEST(expire_line_segment);
    RUN_TEST(expire_line_segment_vertical);
    RUN_TEST(test_from_line_lon_lat_crossing);
    RUN_TEST(test_from_line_lon_lat_wrong_order);
    RUN_TEST(test_expire_set);
    RUN_TEST(test_expire_merge);
    RUN_TEST(test_expire_merge_same);
    RUN_TEST(test_expire_merge_overlap);
    RUN_TEST(test_expire_merge_complete);

    //passed
    return 0;
}
