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

[[nodiscard]] inline float256 logb(float256 x) {
    if (x.isnan()) return x;
    if (x.iszero()) return float256::inf(true);
    if (x.isinf()) return float256::inf(false);
    return float256::from_int(static_cast<std::int64_t>(ilogb(x)));
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
