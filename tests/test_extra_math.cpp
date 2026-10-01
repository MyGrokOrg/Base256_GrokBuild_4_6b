#include <base256/base256.hpp>

#include <cmath>
#include <iostream>

using base256::float256;

static int g_fails = 0;
static int g_pass  = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        if (cond) {                                                            \
            ++g_pass;                                                          \
        } else {                                                               \
            ++g_fails;                                                         \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  "        \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

int main() {
    CHECK(base256::asinh(float256::zero()).iszero());
    CHECK(base256::asinh(-float256::zero()).signbit());
    CHECK(base256::asinh(float256::inf()).isinf());
    {
        const float256 sh = base256::sinh(float256::one());
        CHECK(std::fabs((base256::asinh(sh) - float256::one()).to_double()) < 1e-12);
    }
    CHECK(base256::acosh(float256::one()).iszero());
    CHECK(base256::acosh(float256::half()).isnan());
    {
        const float256 ch = base256::cosh(float256(2));
        CHECK(std::fabs((base256::acosh(ch) - float256(2)).to_double()) < 1e-12);
    }
    CHECK(base256::atanh(float256::zero()).iszero());
    CHECK(base256::atanh(float256::one()).isinf());
    CHECK(!base256::atanh(float256::one()).signbit());
    CHECK(base256::atanh(-float256::one()).isinf() && base256::atanh(-float256::one()).signbit());
    {
        const float256 th = base256::tanh(float256::half());
        CHECK(std::fabs((base256::atanh(th) - float256::half()).to_double()) < 1e-12);
    }
    {
        const float256 huge = base256::ldexp(float256::one(), 200000);
        CHECK(base256::asinh(huge).isfinite());
        CHECK(base256::acosh(huge).isfinite());
    }

    CHECK(base256::fdim(float256(5), float256(2)) == float256(3));
    CHECK(base256::fdim(float256(2), float256(5)).iszero());
    CHECK(!base256::fdim(float256(2), float256(5)).signbit());
    CHECK(base256::fdim(float256::qnan(), float256(1)).isnan());

    CHECK(base256::logb(float256(8)) == float256(3));
    CHECK(base256::logb(float256::zero()).isinf() && base256::logb(float256::zero()).signbit());
    CHECK(base256::logb(float256::inf()).isinf() && !base256::logb(float256::inf()).signbit());
    CHECK(base256::nearbyint(float256(2.5)) == float256(2));
    CHECK(base256::nearbyint(float256(1.5)) == float256(2));

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
