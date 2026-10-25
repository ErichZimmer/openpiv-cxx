#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <initializer_list>

#include "core/vector_field.h"

using namespace openpiv::core;

TEST_CASE("vector_field_test - initialization")
{
    const grid_data fields[] = {grid_data{4, 3}, grid_data{size{4, 3}}};

    for (const auto& field : fields)
    {
        for (const auto* candidate : {
                 &field.u1, &field.v1, &field.s2n1, &field.peak1,
                 &field.u2, &field.v2, &field.s2n2, &field.peak2,
                 &field.u3, &field.v3, &field.s2n3, &field.peak3})
        {
            REQUIRE(candidate->size() == size{4, 3});
            for (const auto value : *candidate)
                CHECK(std::isnan(value.v));
        }

        REQUIRE(field.peak_ind.size() == size{4, 3});
        REQUIRE(field.flag.size() == size{4, 3});
        for (size_t i = 0; i < field.u.pixel_count(); ++i)
        {
            CHECK(field.peak_ind[i].v == 0);
            CHECK(field.flag[i].v == static_cast<uint8_t>(FLAG::VALID));
            CHECK(field.u[i].v == 0.0);
            CHECK(field.v[i].v == 0.0);
            CHECK(field.s2n[i].v == 0.0);
            CHECK(field.peak[i].v == 0.0);
            CHECK(field.p2p[i].v == 0.0);
        }
    }
}

TEST_CASE("vector_field_test - independent storage and copy")
{
    grid_data field{4, 3};
    double value = 1.0;
    for (const auto member : {
             &grid_data::u1, &grid_data::v1, &grid_data::s2n1, &grid_data::peak1,
             &grid_data::u2, &grid_data::v2, &grid_data::s2n2, &grid_data::peak2,
             &grid_data::u3, &grid_data::v3, &grid_data::s2n3, &grid_data::peak3})
    {
        (field.*member)[0] = value++;
    }
    field.peak_ind[0] = 2;
    field.flag[0] = static_cast<uint8_t>(FLAG::SUBSTITUTED);

    grid_data copy{field};
    field.u[0] = 100.0;
    field.v[0] = 200.0;
    field.s2n[0] = 300.0;
    field.peak[0] = 400.0;

    value = 1.0;
    for (const auto member : {
             &grid_data::u1, &grid_data::v1, &grid_data::s2n1, &grid_data::peak1,
             &grid_data::u2, &grid_data::v2, &grid_data::s2n2, &grid_data::peak2,
             &grid_data::u3, &grid_data::v3, &grid_data::s2n3, &grid_data::peak3})
    {
        CHECK((field.*member)[0].v == value);
        CHECK((copy.*member)[0].v == value);
        CHECK(std::isnan((copy.*member)[1].v));
        (copy.*member)[0] = -value;
        CHECK((field.*member)[0].v == value);
        ++value;
    }

    CHECK(copy.u[0].v == 0.0);
    CHECK(copy.v[0].v == 0.0);
    CHECK(copy.s2n[0].v == 0.0);
    CHECK(copy.peak[0].v == 0.0);
    CHECK(copy.peak_ind[0].v == 2);
    CHECK(copy.flag[0].v == static_cast<uint8_t>(FLAG::SUBSTITUTED));
    copy.peak_ind[0] = 3;
    copy.flag[0] = static_cast<uint8_t>(FLAG::INVALID);
    CHECK(field.peak_ind[0].v == 2);
    CHECK(field.flag[0].v == static_cast<uint8_t>(FLAG::SUBSTITUTED));
}
