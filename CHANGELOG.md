# Changelog

## 4.6.3

- **`asinh` / `acosh` / `atanh`** with large-argument `log(2|x|)` reduction so they stay finite past `2^64`.
- **`fdim(x, y)`** returns `max(x-y, +0)`. **`logb`** is the floating `ilogb`. **`nearbyint`** matches `rint` (ties to even).

## 4.6.2

- **`pow(±0, y)`** keeps the sign of zero when `y` is an odd integer. `pow(−0, 3)` is `−0`; `pow(−0, −1)` is `−∞`.
- **`atan2`** on infinities returns `±π/4` or `±3π/4`. `atan2(y, +∞)` is `±0` and `atan2(y, −∞)` is `±π`.
- **`sin` / `cos` / `tan`** at exact multiples of this library's π no longer produce a spurious `−0`. `sin(π)` is `+0`, `cos(π/2)` is `+0`, `tan(π/2)` is `+∞`.
- **`rint`** rounds to nearest, ties to even. **`modf`** splits a value into a signed integer and a signed fraction. **`nextup`** / **`nextdown`** wrap `nextafter`. **`remainder`** uses ties-to-even.
- **`hypot(∞, NaN)`** returns `+∞`, matching IEEE 754.

## 4.6.1

- **`to_string`** no longer trusts a double `log10` estimate that lands on the wrong side of a power of ten (`0.01` used to print as `1e-1`).
- **`hypot`** scales into `[0.5, 1)` before squaring, so large finite inputs stay finite.
- **`fmin` / `fmax`** of two zeros follow IEEE 754.

## 4.6.0

Initial release: header-only IEEE 754 binary256 with correctly rounded `+ − × ÷`, FMA, `sqrt`, decimal and hex I/O, and elementary functions.
