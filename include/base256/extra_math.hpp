#pragma once

/// Extra elementary functions layered on math.hpp.
#include "math.hpp"

namespace base256 {

[[nodiscard]] inline float256 nearbyint(float256 x) {
    // Default rounding mode is roundTiesToEven, same as rint.
    return rint(x);
}

[[nodiscard]] inline float256 fdim(float256 x, float256 y) {
    if (x.isnan() || y.isnan()) return float256::qnan();
    if (x > y) return x - y;
    return float256::zero(false);
}

// logb lives in float256.hpp (constexpr). A second definition here made
// <base256/base256.hpp> ill-formed.

namespace detail_spec {

// erf series on |x| <= 2. Odd, so a signed argument is fine.
[[nodiscard]] inline float256 erf_series(float256 x) {
    const float256 x2 = x * x;
    float256 term = x;
    float256 sum  = x;
    for (int n = 1; n <= 250; ++n) {
        term = -term * x2 * float256::from_int(static_cast<std::int64_t>(2 * n - 1))
               / (float256::from_int(static_cast<std::int64_t>(n)) *
                  float256::from_int(static_cast<std::int64_t>(2 * n + 1)));
        const float256 next = sum + term;
        if (next == sum) break;
        sum = next;
    }
    return sum * (float256::two() / sqrt(float256::pi()));
}

// erfc(x) for x > 2 via the continued fraction
//   exp(x^2) sqrt(pi) erfc(x) = 1 / (x + (1/2)/(x + (2/2)/(x + (3/2)/(x + ...))))
[[nodiscard]] inline float256 erfc_cf(float256 x) {
    const float256 x2 = x * x;
    if (x2.isinf() || x > float256::from_int(600)) return float256::zero();
    const float256 tiny = ldexp(float256::one(), -100);
    const float256 tol  = ldexp(float256::one(), -200);
    float256 f = x;
    float256 C = f;
    float256 D = float256::zero();
    for (int j = 1; j <= 500; ++j) {
        const float256 a = float256::from_int(static_cast<std::int64_t>(j)) * float256::half();
        D = x + a * D;
        if (D.iszero()) D = tiny;
        C = x + a / C;
        if (C.iszero()) C = tiny;
        D = float256::one() / D;
        const float256 delta = C * D;
        f = f * delta;
        if (abs(delta - float256::one()) <= tol) break;
    }
    const float256 num = exp(-x2);
    if (num.iszero()) return float256::zero();
    return num / (sqrt(float256::pi()) * f);
}

} // namespace detail_spec

[[nodiscard]] inline float256 erf(float256 x) {
    if (x.isnan()) return x;
    if (x.isinf()) return copysign(float256::one(), x);
    if (x.iszero()) return x;
    if (abs(x) <= float256::two()) return detail_spec::erf_series(x);
    const float256 y = float256::one() - detail_spec::erfc_cf(abs(x));
    return copysign(y, x);
}

[[nodiscard]] inline float256 erfc(float256 x) {
    if (x.isnan()) return x;
    if (x.isinf()) return x.signbit() ? float256::two() : float256::zero();
    if (abs(x) <= float256::two()) return float256::one() - detail_spec::erf_series(x);
    if (x.signbit()) return float256::two() - detail_spec::erfc_cf(abs(x));
    return detail_spec::erfc_cf(x);
}

[[nodiscard]] inline float256 asinh(float256 x) {
    if (x.isnan() || x.isinf() || x.iszero()) return x;
    const bool s = x.signbit();
    const float256 a = abs(x);
    float256 y;
    // asinh(x) = log(x + sqrt(x^2+1)); for large x, log(2x).
    const float256 thresh = ldexp(float256::one(), 64);
    if (a > thresh) {
        y = log(a) + float256::ln2();
    } else {
        y = log1p(a + a * a / (float256::one() + sqrt(a * a + float256::one())));
    }
    return s ? -y : y;
}

[[nodiscard]] inline float256 acosh(float256 x) {
    if (x.isnan()) return x;
    if (x < float256::one()) return float256::qnan();
    if (x == float256::one()) return float256::zero();
    if (x.isinf()) return x;
    const float256 thresh = ldexp(float256::one(), 64);
    if (x > thresh) return log(x) + float256::ln2();
    // acosh(x) = log(x + sqrt(x^2-1)) = log1p(x-1 + sqrt((x-1)(x+1)))
    const float256 xm1 = x - float256::one();
    return log1p(xm1 + sqrt(xm1 * (x + float256::one())));
}

[[nodiscard]] inline float256 atanh(float256 x) {
    if (x.isnan()) return x;
    if (x.iszero()) return x;
    const float256 a = abs(x);
    if (a > float256::one()) return float256::qnan();
    if (a == float256::one()) return float256::inf(x.signbit());
    // atanh(x) = 1/2 log((1+x)/(1-x)) = 1/2 log1p(2x/(1-x))
    const float256 y = log1p((a + a) / (float256::one() - a)) * float256::half();
    return x.signbit() ? -y : y;
}

} // namespace base256
