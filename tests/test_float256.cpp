#include <base256/base256.hpp>

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using base256::float256;
using namespace base256::literals;

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

#define CHECK_EQ(a, b)                                                         \
    do {                                                                       \
        const auto _va = (a);                                                  \
        const auto _vb = (b);                                                  \
        if (_va == _vb) {                                                      \
            ++g_pass;                                                          \
        } else {                                                               \
            ++g_fails;                                                         \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  "        \
                      << #a << " == " << #b << "\n    lhs=" << _va             \
                      << "\n    rhs=" << _vb << "\n";                          \
        }                                                                      \
    } while (0)

int main() {
    using L = std::numeric_limits<float256>;

    // --- specials ----------------------------------------------------------
    CHECK(float256::zero().iszero());
    CHECK((-float256::zero()).iszero());
    CHECK(float256::zero().signbit() == false);
    CHECK((-float256::zero()).signbit() == true);
    CHECK(float256::inf().isinf());
    CHECK((-float256::inf()).isinf());
    CHECK(float256::qnan().isnan());
    CHECK(!float256::qnan().isfinite());
    CHECK(L::has_infinity);
    CHECK(L::digits == 237);
    CHECK(L::digits10 == 71);
    CHECK(L::is_iec559);

    // --- integer / double round-trip ---------------------------------------
    CHECK_EQ(float256(0).to_double(), 0.0);
    CHECK_EQ(float256(1).to_double(), 1.0);
    CHECK_EQ(float256(2).to_double(), 2.0);
    CHECK_EQ(float256(-3).to_double(), -3.0);
    CHECK_EQ(float256(42.0).to_double(), 42.0);
    CHECK_EQ(float256(-0.5).to_double(), -0.5);
    CHECK_EQ(float256(3.141592653589793).to_double(), 3.141592653589793);

    // --- constants vs double -----------------------------------------------
    CHECK(std::fabs(float256::pi().to_double() - 3.141592653589793) < 1e-15);
    CHECK(std::fabs(float256::e().to_double()  - 2.718281828459045) < 1e-15);
    CHECK(std::fabs(float256::ln2().to_double() - 0.6931471805599453) < 1e-15);

    // --- 1.0 encoding ------------------------------------------------------
    {
        auto b = float256::one().raw_bits();
        CHECK(b[0] == 0 && b[1] == 0 && b[2] == 0);
        CHECK(b[3] == 0x3ffff00000000000ULL);
    }
    {
        auto b = float256::two().raw_bits();
        CHECK(b[3] == 0x4000000000000000ULL);
    }

    // --- arithmetic: small integers ----------------------------------------
    CHECK_EQ((float256(1) + float256(1)).to_double(), 2.0);
    CHECK_EQ((float256(10) - float256(3)).to_double(), 7.0);
    CHECK_EQ((float256(6) * float256(7)).to_double(), 42.0);
    CHECK_EQ((float256(84) / float256(2)).to_double(), 42.0);
    CHECK_EQ((float256(1) / float256(2)).to_double(), 0.5);

    // --- signed zero -------------------------------------------------------
    CHECK((float256::zero() + float256::zero()).signbit() == false);
    CHECK(((-float256::zero()) + (-float256::zero())).signbit() == true);
    CHECK(float256::zero() == -float256::zero());

    // --- inf / nan propagation ---------------------------------------------
    CHECK((float256::inf() + float256(1)).isinf());
    CHECK((float256::inf() - float256::inf()).isnan());
    CHECK((float256::inf() * float256::zero()).isnan());
    CHECK((float256(1) / float256::zero()).isinf());
    CHECK((float256::zero() / float256::zero()).isnan());
    CHECK(!(float256::qnan() == float256::qnan()));
    CHECK(!(float256::qnan() < float256(1)));

    // --- precision beyond double -------------------------------------------
    // ulp(1.0) for double is 2^-52; for binary256 it is 2^-236.
    {
        const float256 one = float256::one();
        const float256 delta = base256::ldexp(float256::one(), -80);
        const float256 sum = one + delta;
        CHECK(sum != one);
        CHECK((sum - one) == delta);
        // double cannot represent 1 + 2^-80
        CHECK(1.0 + std::ldexp(1.0, -80) == 1.0);
    }

    // --- 1/3 * 3 recovers 1 ------------------------------------------------
    {
        const float256 third = float256::one() / float256(3);
        const float256 back  = third * float256(3);
        CHECK(back == float256::one());
    }

    // --- parse / format ----------------------------------------------------
    {
        auto p = float256::parse("1.0");
        CHECK(p.has_value());
        CHECK(*p == float256::one());
        auto p2 = float256::parse("2");
        CHECK(p2.has_value() && *p2 == float256::two());
        auto pi = float256::parse("3.14159265358979323846");
        CHECK(pi.has_value());
        CHECK(std::fabs(pi->to_double() - 3.141592653589793) < 1e-15);
        auto infs = float256::parse("-inf");
        CHECK(infs.has_value() && infs->isinf() && infs->signbit());
    }

    // 0.1 + 0.2 is closer to 0.3 than IEEE-754 binary64
    {
        auto a = float256::parse("0.1");
        auto b = float256::parse("0.2");
        auto c = float256::parse("0.3");
        CHECK(a && b && c);
        const float256 s = *a + *b;
        const float256 err256 = base256::abs(s - *c);
        const double err64 = std::fabs((0.1 + 0.2) - 0.3);
        CHECK(err256.to_double() < err64);
        CHECK(err256.to_double() < 1e-20);
    }

    // user-defined literal
    {
        auto x = 2.5_f256;
        CHECK(x.to_double() == 2.5);
        auto y = "0.125"_f256;
        CHECK(y.to_double() == 0.125);
    }

    // --- sqrt --------------------------------------------------------------
    {
        const float256 s = base256::sqrt(float256(4));
        CHECK(s == float256(2));
        const float256 s2 = base256::sqrt(float256(9));
        CHECK(s2 == float256(3));
        const float256 s0 = base256::sqrt(float256::zero());
        CHECK(s0.iszero());
        CHECK(base256::sqrt(-float256::one()).isnan());
        const float256 r2 = base256::sqrt(float256(2));
        const float256 back = r2 * r2;
        CHECK(base256::abs(back - float256(2)).to_double() < 1e-15);
        // high-precision: |sqrt(2)^2 - 2| should be at most a few ulps of 2
        const float256 ulp2 = base256::ldexp(float256::one(), 1 - 236);
        CHECK(base256::abs(back - float256(2)) <= ulp2 * float256(8));
    }

    // --- comparisons / nextafter -------------------------------------------
    CHECK(float256(1) < float256(2));
    CHECK(float256(-1) < float256(1));
    CHECK(float256::max() > float256::one());
    {
        const float256 n = base256::nextafter(float256::one(), float256(2));
        CHECK(n > float256::one());
        CHECK(base256::nextafter(n, float256::one()) == float256::one());
    }

    // --- ldexp / frexp -----------------------------------------------------
    {
        int e = 0;
        const float256 m = base256::frexp(float256(8), &e);
        CHECK(e == 4);
        CHECK(m.to_double() == 0.5);
        CHECK(base256::ldexp(m, e) == float256(8));
    }

    // --- floor / ceil / trunc ----------------------------------------------
    CHECK(base256::floor(float256(1.5)) == float256(1));
    CHECK(base256::ceil(float256(1.5)) == float256(2));
    CHECK(base256::trunc(float256(-1.5)) == float256(-1));
    CHECK(base256::floor(float256(-1.5)) == float256(-2));

    // --- exp / log ---------------------------------------------------------
    {
        const float256 e1 = base256::exp(float256::one());
        CHECK(std::fabs(e1.to_double() - 2.718281828459045) < 1e-12);
        const float256 back = base256::log(e1);
        CHECK(std::fabs(back.to_double() - 1.0) < 1e-12);
        const float256 l2 = base256::log(float256(2));
        CHECK(std::fabs((l2 - float256::ln2()).to_double()) < 1e-12);
    }

    // --- sin / cos ---------------------------------------------------------
    {
        CHECK(std::fabs(base256::sin(float256::zero()).to_double()) < 1e-30);
        CHECK(std::fabs(base256::cos(float256::zero()).to_double() - 1.0) < 1e-15);
        const float256 s = base256::sin(float256::pi() / float256(2));
        CHECK(std::fabs(s.to_double() - 1.0) < 1e-10);
        const float256 c = base256::cos(float256::pi());
        CHECK(std::fabs(c.to_double() + 1.0) < 1e-10);
    }

    // --- subnormals exist --------------------------------------------------
    {
        const float256 dmin = float256::denorm_min();
        CHECK(dmin.is_subnormal());
        CHECK(dmin > float256::zero());
        CHECK(base256::nextafter(float256::zero(), float256::one()) == dmin);
    }

    // --- max / min / epsilon -----------------------------------------------
    CHECK(float256::max().isnormal());
    CHECK(float256::min().isnormal());
    CHECK((float256::one() + float256::epsilon()) != float256::one());
    CHECK((float256::one() + float256::epsilon() * float256::half()) == float256::one() ||
          (float256::one() + float256::epsilon() * float256::half()) != float256::one());
    // 1 + eps/2 ties-to-even should stay 1 (lsb of 1.0 is 0)
    CHECK((float256::one() + float256::epsilon() * float256::half()) == float256::one());

    // --- to_string ---------------------------------------------------------
    {
        CHECK(float256(1).to_string(1).rfind("1", 0) == 0);
        const std::string s2 = float256(2).to_string(6);
        CHECK(s2.find("2.00000") != std::string::npos || s2.find("2e+") != std::string::npos);
        const std::string spi = float256::pi().to_string(20);
        CHECK(spi.find("3.141592653589793238") != std::string::npos);
        const std::string s01 = float256::parse("0.1")->to_string(20);
        CHECK(s01.find("1.0000000000000000000e-1") != std::string::npos);
    }

    // --- hex round-trip ----------------------------------------------------
    {
        const auto h = float256::pi().to_hex_string();
        auto p = float256::parse(h);
        CHECK(p.has_value());
        CHECK(*p == float256::pi());
    }

    std::cout << "passed=" << g_pass << " failed=" << g_fails << "\n";
    return g_fails ? 1 : 0;
}
