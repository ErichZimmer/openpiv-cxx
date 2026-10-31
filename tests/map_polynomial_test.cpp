#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "interp/map_polynomial.h"

using openpiv::interp::poly_weights;

namespace
{
    // Independent reference: the original Lagrange product over every node.
    double reference_weight(double r, int32_t k, int32_t i)
    {
        const double xi = i - (k - 0.5);
        double weight = 1.0;
        for (int32_t j = 0; j < 2 * k; ++j)
        {
            if (j != i)
            {
                const double xj = j - (k - 0.5);
                weight *= (r - xj) / (xi - xj);
            }
        }
        return weight;
    }
}

TEST_CASE("poly_weights - stencils match the original Lagrange basis")
{
    for (int32_t k = 1; k <= 5; ++k)
    {
        std::vector<double> w(2 * k);
        std::vector<double> reflected(2 * k);
        for (int32_t step = 0; step <= 4096; ++step)
        {
            const double r = -0.5 + static_cast<double>(step) / 4096.0;
            INFO("k=" << k << ", r=" << r);
            poly_weights(r, k, w);
            poly_weights(-r, k, reflected);

            for (int32_t i = 0; i < 2 * k; ++i)
            {
                REQUIRE(std::abs(w[i] - reference_weight(r, k, i)) < 2e-14);
                REQUIRE(std::abs(w[i] - reflected[2 * k - 1 - i]) < 2e-14);
            }
            REQUIRE(std::abs(std::accumulate(w.begin(), w.end(), 0.0) - 1.0) < 2e-14);
        }
    }
}

TEST_CASE("poly_weights - exact nodes select the corresponding sample")
{
    for (int32_t k = 1; k <= 5; ++k)
    {
        std::vector<double> w(2 * k);
        for (int32_t node = 0; node < 2 * k; ++node)
        {
            poly_weights(node - (k - 0.5), k, w);
            for (int32_t i = 0; i < 2 * k; ++i)
                REQUIRE(std::abs(w[i] - (i == node ? 1.0 : 0.0)) < 2e-14);
        }
    }
}

TEST_CASE("poly_weights - reproduce polynomials through the stencil degree")
{
    for (int32_t k = 1; k <= 5; ++k)
    {
        std::vector<double> w(2 * k);
        const double scale = k - 0.5;
        for (int32_t step = 0; step <= 128; ++step)
        {
            const double r = -0.5 + static_cast<double>(step) / 128.0;
            poly_weights(r, k, w);
            for (int32_t degree = 0; degree < 2 * k; ++degree)
            {
                INFO("k=" << k << ", r=" << r << ", degree=" << degree);
                double value = 0.0;
                for (int32_t i = 0; i < 2 * k; ++i)
                    value += w[i] * std::pow((i - scale) / scale, degree);
                REQUIRE(std::abs(value - std::pow(r / scale, degree)) < 2e-14);
            }
        }
    }
}

TEST_CASE("poly_weights - reuse buffers and preserve elements beyond the stencil")
{
    for (int32_t k = 1; k <= 5; ++k)
    {
        std::vector<double> w(2 * k + 3, 123.0);
        const auto* data = w.data();
        const auto capacity = w.capacity();
        for (double r : {-0.5, -0.2, 0.0, 0.3, 0.5})
        {
            poly_weights(r, k, w);
            REQUIRE(w.data() == data);
            REQUIRE(w.size() == static_cast<size_t>(2 * k + 3));
            REQUIRE(w.capacity() == capacity);
            for (int32_t i = 0; i < 2 * k; ++i)
                REQUIRE(std::abs(w[i] - reference_weight(r, k, i)) < 2e-14);
            for (size_t i = 2 * k; i < w.size(); ++i)
                REQUIRE(w[i] == 123.0);
        }
    }
}

TEST_CASE("poly_weights - reject invalid orders and undersized buffers before writing")
{
    for (int32_t k : {std::numeric_limits<int32_t>::min(), -1, 0, 6,
                      std::numeric_limits<int32_t>::max()})
    {
        std::vector<double> w(10, 123.0);
        const auto original = w;
        REQUIRE_THROWS_AS(poly_weights(0.0, k, w), std::runtime_error);
        REQUIRE(w == original);
    }

    for (int32_t k = 1; k <= 5; ++k)
    {
        for (int32_t length = 0; length < 2 * k; ++length)
        {
            std::vector<double> w(length, 123.0);
            w.reserve(2 * k); // Capacity alone does not provide writable elements.
            const auto original = w;
            REQUIRE_THROWS_AS(poly_weights(0.0, k, w), std::runtime_error);
            REQUIRE(w == original);
        }
    }
}
