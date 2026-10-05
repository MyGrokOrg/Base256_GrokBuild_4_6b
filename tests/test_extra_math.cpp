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

    // exp2 / log2 are exact on integers and powers of two.
    CHECK(base256::exp2(float256::zero()) == float256::one());
    CHECK(base256::exp2(float256(10)) == float256(1024));
    CHECK(base256::exp2(float256(-10)) == float256::one() / float256(1024));
    CHECK(base256::exp2(float256(0.5)) == float256::sqrt2() ||
          base256::abs(base256::exp2(float256::half()) - float256::sqrt2()) <
              float256::epsilon() * float256(8));
    CHECK(base256::log2(float256(1024)) == float256(10));
    CHECK(base256::log2(float256::one()) == float256::zero());
    CHECK(base256::log2(float256::half()) == -float256::one());
    CHECK(base256::log2(base256::ldexp(float256::one(), -80)) == float256(-80));
    CHECK(base256::log2(float256::denorm_min()) == float256(base256::ilogb(float256::denorm_min())));
    CHECK(base256::exp2(float256(262143)).isfinite());
    CHECK(base256::exp2(float256(262144)).isinf());
    CHECK(base256::log2(float256::zero()).isinf() && base256::log2(float256::zero()).signbit());
    CHECK(base256::log2(-float256::one()).isnan());

    // Tiny sinh must not cancel to zero. Huge finite tanh must not be NaN.
    {
        const float256 tiny = base256::ldexp(float256::one(), -200);
        const float256 sh = base256::sinh(tiny);
        CHECK(sh == tiny || base256::abs(sh - tiny) <= base256::ldexp(tiny, -20));
        CHECK(base256::cosh(tiny) == float256::one());
        const float256 huge = base256::ldexp(float256::one(), 20);
        const float256 th = base256::tanh(huge);
        CHECK(th == float256::one());
        const float256 nth = base256::tanh(-huge);
        CHECK(nth == -float256::one());
        CHECK(base256::tanh(float256::inf()) == float256::one());
        CHECK(base256::tanh(-float256::inf()) == -float256::one());
        const float256 half_th = base256::tanh(float256::half());
        CHECK(std::fabs(half_th.to_double() - std::tanh(0.5)) < 1e-15);
    }

    // erf / erfc
    CHECK(base256::erf(float256::zero()).iszero() && !base256::erf(float256::zero()).signbit());
    CHECK(base256::erf(-float256::zero()).iszero() && base256::erf(-float256::zero()).signbit());
    CHECK(base256::erf(float256::inf()) == float256::one());
    CHECK(base256::erf(-float256::inf()) == -float256::one());
    CHECK(base256::erfc(float256::inf()).iszero() && !base256::erfc(float256::inf()).signbit());
    CHECK(base256::erfc(-float256::inf()) == float256::two());
    CHECK(base256::erfc(float256::zero()) == float256::one());
    {
        const float256 z = float256::one();
        CHECK(base256::erf(z) + base256::erfc(z) == float256::one());
        // Negative side is 1 - (−erf) = 1+|erf|, which is not a Sterbenz
        // subtraction, so the sum is 1 within an ulp rather than bitwise.
        const float256 neg_id = base256::erf(-z) + base256::erfc(-z);
        CHECK(base256::abs(neg_id - float256::one()) <= float256::epsilon());
        CHECK(std::fabs(base256::erf(z).to_double() - std::erf(1.0)) < 1e-14);
        CHECK(std::fabs(base256::erfc(z).to_double() - std::erfc(1.0)) < 1e-14);
        const float256 three = float256(3);
        CHECK(base256::erf(three) + base256::erfc(three) == float256::one());
        CHECK(std::fabs(base256::erf(three).to_double() - std::erf(3.0)) < 1e-14);
        CHECK(std::fabs(base256::erfc(three).to_double() - std::erfc(3.0)) < 1e-15);
        CHECK(base256::erf(float256(20)) == float256::one());
        CHECK(base256::erfc(float256(20)) > float256::zero());
        CHECK(base256::erfc(float256(20)).isfinite());
        CHECK(base256::erf(-float256(8)) == -base256::erf(float256(8)));
    }

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
