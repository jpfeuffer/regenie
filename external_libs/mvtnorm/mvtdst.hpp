/*
 * C++ port of Alan Genz's MVTDST Fortran code for computing
 * non-central multivariate t and normal probabilities.
 *
 * Original algorithm: QRSVN described in
 *   "Comparison of Methods for the Computation of Multivariate
 *    t-Probabilities", Alan Genz and Frank Bretz,
 *    J. Comp. Graph. Stat. 11 (2002), pp. 950-971.
 *
 * Ported to C++ with:
 *   - All global/SAVE/COMMON state eliminated (thread-safe)
 *   - OpenMP parallelization of the QMC lattice rule loop
 *   - Eigen-based triangular mat-vec in the integrand
 *   - No Fortran dependencies
 *
 * License: Same as original (public domain scientific code).
 */

#ifndef MVTDST_HPP_
#define MVTDST_HPP_

// MSVC's <cmath> only defines M_PI/M_SQRT1_2/etc. (BSD/POSIX extensions, not
// standard C++) when this is defined first; must precede the very first
// inclusion of <cmath> in the translation unit, which for every caller of
// this header is this one, right below.
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <random>
#include <numeric>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

// Eigen: used for vectorized dot products in the integrand
#ifdef EIGEN_WORLD_VERSION
// Eigen already included by the including translation unit
#else
#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#define MVTDST_HAS_EIGEN 1
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#define MVTDST_HAS_EIGEN 1
#endif
#endif

#if defined(EIGEN_WORLD_VERSION) && !defined(MVTDST_HAS_EIGEN)
#define MVTDST_HAS_EIGEN 1
#endif

namespace mvtdst {

// ─── Constants ──────────────────────────────────────────────────────────────
static constexpr double PI      = 3.141592653589793238462643383;
static constexpr double TWOPI   = 6.283185307179586476925286767;
static constexpr double SQTWPI  = 2.506628274631000502415765285;
static constexpr double RTWO    = 1.414213562373095048801688724;
static constexpr int    NL      = 1000;
static constexpr int    PLIM    = 28;
static constexpr int    KLIM    = 100;
static constexpr int    FLIM    = 5000;
static constexpr int    MINSMP  = 8;

// ─── Context struct (replaces SAVE/COMMON blocks) ───────────────────────────
struct MvtContext {
    int nu = 0;
    int nd = 0;
    double snu = 0.0;
    std::vector<int> infi;
    std::vector<double> a, b, dl, cov, y;
};

// ─── Forward declarations ───────────────────────────────────────────────────

// Normal CDF Φ(z) — Schonfelder approximation accurate to 1e-15
inline double mvphi(double z);

// Normal quantile (inverse CDF) — AS241 algorithm
inline double mvphnv(double p);

// Student-t CDF
inline double mvstdt(int nu, double t);

// t-density
inline double mvtdns(int nu, double x);

// Bivariate normal P(X > sh, Y > sk)
double mvbvu(double sh, double sk, double r);

// Bivariate normal probability with limit flags
double mvbvn(const double* lower, const double* upper, const int* infin, double correl);

// Bivariate t probability P(X < dh, Y < dk)
double mvbvtl(int nu, double dh, double dk, double r);

// Bivariate t probability with limit flags
double mvbvt(int nu, const double* lower, const double* upper, const int* infin, double correl);

// Chi-squared quantile helper
double mvchnv(int n, double p);

// Main entry point
void mvtdst(int n, int nu, double* lower, double* upper, int* infin,
            double* correl, double* delta, int maxpts, double abseps,
            double releps, double* error, double* value, int* inform);

// ═══════════════════════════════════════════════════════════════════════════════
//                          IMPLEMENTATION
// ═══════════════════════════════════════════════════════════════════════════════

// ─── Normal CDF ─────────────────────────────────────────────────────────────
// Phi(z) = erfc(-z/sqrt(2)) / 2.  Apple/glibc erfc is ~15-20 cycles
// (hardware fexp on ARM64) vs ~70+ cycles for 25-step Clenshaw + exp.
inline double mvphi(double z) {
    return 0.5 * std::erfc(-z * M_SQRT1_2);
}

// ─── Normal quantile (inverse CDF) — AS241 ─────────────────────────────────
inline double mvphnv(double p) {
    constexpr double SPLIT1 = 0.425;
    constexpr double SPLIT2 = 5.0;
    constexpr double CONST1 = 0.180625;
    constexpr double CONST2 = 1.6;

    // Coefficients for P close to 0.5
    constexpr double A0 = 3.3871328727963666080e0;
    constexpr double A1 = 1.3314166789178437745e2;
    constexpr double A2 = 1.9715909503065514427e3;
    constexpr double A3 = 1.3731693765509461125e4;
    constexpr double A4 = 4.5921953931549871457e4;
    constexpr double A5 = 6.7265770927008700853e4;
    constexpr double A6 = 3.3430575583588128105e4;
    constexpr double A7 = 2.5090809287301226727e3;
    constexpr double B1 = 4.2313330701600911252e1;
    constexpr double B2 = 6.8718700749205790830e2;
    constexpr double B3 = 5.3941960214247511077e3;
    constexpr double B4 = 2.1213794301586595867e4;
    constexpr double B5 = 3.9307895800092710610e4;
    constexpr double B6 = 2.8729085735721942674e4;
    constexpr double B7 = 5.2264952788528545610e3;

    // Coefficients for P not close to 0, 0.5 or 1
    constexpr double C0 = 1.42343711074968357734e0;
    constexpr double C1 = 4.63033784615654529590e0;
    constexpr double C2 = 5.76949722146069140550e0;
    constexpr double C3 = 3.64784832476320460504e0;
    constexpr double C4 = 1.27045825245236838258e0;
    constexpr double C5 = 2.41780725177450611770e-1;
    constexpr double C6 = 2.27238449892691845833e-2;
    constexpr double C7 = 7.74545014427834407640e-4;
    constexpr double D1 = 2.05319162663775882187e0;
    constexpr double D2 = 1.67638483018380162940e0;
    constexpr double D3 = 6.89767334985100004550e-1;
    constexpr double D4 = 1.48103976427480074590e-1;
    constexpr double D5 = 1.51986665636164571966e-2;
    constexpr double D6 = 5.47593808499534494600e-4;
    constexpr double D7 = 1.05075007164441684324e-9;

    // Coefficients for P near 0 or 1
    constexpr double E0 = 6.65790464350110377720e0;
    constexpr double E1 = 5.46378491116411436990e0;
    constexpr double E2 = 1.78482653991729133580e0;
    constexpr double E3 = 2.96560571828504891230e-1;
    constexpr double E4 = 2.65321895265761230930e-2;
    constexpr double E5 = 1.24266094738807843860e-3;
    constexpr double E6 = 2.71155556874348757815e-5;
    constexpr double E7 = 2.01033439929228813265e-7;
    constexpr double F1 = 5.99832206555887937690e-1;
    constexpr double F2 = 1.36929880922735805310e-1;
    constexpr double F3 = 1.48753612908506148525e-2;
    constexpr double F4 = 7.86869131145613259100e-4;
    constexpr double F5 = 1.84631831751005468180e-5;
    constexpr double F6 = 1.42151175831644588870e-7;
    constexpr double F7 = 2.04426310338993978564e-15;

    double q = (2.0 * p - 1.0) / 2.0;
    double r, result;

    if (std::abs(q) <= SPLIT1) {
        r = CONST1 - q * q;
        result = q * (((((((A7*r + A6)*r + A5)*r + A4)*r + A3)*r + A2)*r + A1)*r + A0)
                    / (((((((B7*r + B6)*r + B5)*r + B4)*r + B3)*r + B2)*r + B1)*r + 1.0);
    } else {
        r = std::min(p, 1.0 - p);
        if (r > 0.0) {
            r = std::sqrt(-std::log(r));
            if (r <= SPLIT2) {
                r -= CONST2;
                result = (((((((C7*r + C6)*r + C5)*r + C4)*r + C3)*r + C2)*r + C1)*r + C0)
                       / (((((((D7*r + D6)*r + D5)*r + D4)*r + D3)*r + D2)*r + D1)*r + 1.0);
            } else {
                r -= SPLIT2;
                result = (((((((E7*r + E6)*r + E5)*r + E4)*r + E3)*r + E2)*r + E1)*r + E0)
                       / (((((((F7*r + F6)*r + F5)*r + F4)*r + F3)*r + F2)*r + F1)*r + 1.0);
            }
        } else {
            result = 9.0;
        }
        if (q < 0.0) result = -result;
    }
    return result;
}

// ─── Student-t CDF ──────────────────────────────────────────────────────────
inline double mvstdt(int nu, double t) {
    if (nu < 1) {
        return mvphi(t);
    } else if (nu == 1) {
        return (1.0 + 2.0 * std::atan(t) / PI) / 2.0;
    } else if (nu == 2) {
        return (1.0 + t / std::sqrt(2.0 + t * t)) / 2.0;
    } else {
        double tt = t * t;
        double csthe = static_cast<double>(nu) / (nu + tt);
        double polyn = 1.0;
        for (int j = nu - 2; j >= 2; j -= 2) {
            polyn = 1.0 + static_cast<double>(j - 1) * csthe * polyn / j;
        }
        double result;
        if (nu % 2 == 1) {
            double rn = static_cast<double>(nu);
            double ts = t / std::sqrt(rn);
            result = (1.0 + 2.0 * (std::atan(ts) + ts * csthe * polyn) / PI) / 2.0;
        } else {
            double snthe = t / std::sqrt(static_cast<double>(nu) + tt);
            result = (1.0 + snthe * polyn) / 2.0;
        }
        if (result < 0.0) result = 0.0;
        return result;
    }
}

// ─── t-density ──────────────────────────────────────────────────────────────
inline double mvtdns(int nu, double x) {
    if (nu > 0) {
        double prod = 1.0 / std::sqrt(static_cast<double>(nu));
        for (int i = nu - 2; i >= 1; i -= 2) {
            prod = prod * static_cast<double>(i + 1) / i;
        }
        if (nu % 2 == 0) {
            prod = prod / 2.0;
        } else {
            prod = prod / PI;
        }
        return prod / std::pow(1.0 + x * x / nu, (nu + 1) / 2.0);
    } else {
        if (std::abs(x) < 10.0) return std::exp(-x * x / 2.0) / SQTWPI;
        return 0.0;
    }
}

// ─── Bivariate normal P(X > sh, Y > sk) ────────────────────────────────────
inline double mvbvu(double sh, double sk, double r) {
    // Gauss-Legendre Points and Weights
    static const double W[10][3] = {
        {0.1713244923791705e+00, 0.4717533638651177e-01, 0.1761400713915212e-01},
        {0.3607615730481384e+00, 0.1069393259953183e+00, 0.4060142980038694e-01},
        {0.4679139345726904e+00, 0.1600783285433464e+00, 0.6267204833410906e-01},
        {0.0,                    0.2031674267230659e+00, 0.8327674157670475e-01},
        {0.0,                    0.2334925365383547e+00, 0.1019301198172404e+00},
        {0.0,                    0.2491470458134029e+00, 0.1181945319615184e+00},
        {0.0,                    0.0,                    0.1316886384491766e+00},
        {0.0,                    0.0,                    0.1420961093183821e+00},
        {0.0,                    0.0,                    0.1491729864726037e+00},
        {0.0,                    0.0,                    0.1527533871307259e+00}
    };
    static const double X[10][3] = {
        {-0.9324695142031522e+00, -0.9815606342467191e+00, -0.9931285991850949e+00},
        {-0.6612093864662647e+00, -0.9041172563704750e+00, -0.9639719272779138e+00},
        {-0.2386191860831970e+00, -0.7699026741943050e+00, -0.9122344282513259e+00},
        { 0.0,                    -0.5873179542866171e+00, -0.8391169718222188e+00},
        { 0.0,                    -0.3678314989981802e+00, -0.7463319064601508e+00},
        { 0.0,                    -0.1252334085114692e+00, -0.6360536807265150e+00},
        { 0.0,                     0.0,                    -0.5108670019508271e+00},
        { 0.0,                     0.0,                    -0.3737060887154196e+00},
        { 0.0,                     0.0,                    -0.2277858511416451e+00},
        { 0.0,                     0.0,                    -0.7652652113349733e-01}
    };

    int ng, lg;
    if (std::abs(r) < 0.3) {
        ng = 0; lg = 3;
    } else if (std::abs(r) < 0.75) {
        ng = 1; lg = 6;
    } else {
        ng = 2; lg = 10;
    }

    double h = sh, k = sk;
    double hk = h * k;
    double bvn = 0.0;

    if (std::abs(r) < 0.925) {
        double hs = (h * h + k * k) / 2.0;
        double asr = std::asin(r);
        for (int i = 0; i < lg; ++i) {
            double sn = std::sin(asr * (X[i][ng] + 1.0) / 2.0);
            bvn += W[i][ng] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
            sn = std::sin(asr * (-X[i][ng] + 1.0) / 2.0);
            bvn += W[i][ng] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
        }
        bvn = bvn * asr / (2.0 * TWOPI) + mvphi(-h) * mvphi(-k);
    } else {
        if (r < 0.0) {
            k = -k;
            hk = -hk;
        }
        if (std::abs(r) < 1.0) {
            double as = (1.0 - r) * (1.0 + r);
            double a = std::sqrt(as);
            double bs = (h - k) * (h - k);
            double c = (4.0 - hk) / 8.0;
            double d = (12.0 - hk) / 16.0;
            bvn = a * std::exp(-(bs / as + hk) / 2.0)
                    * (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);
            if (hk > -160.0) {
                double b = std::sqrt(bs);
                bvn = bvn - std::exp(-hk / 2.0) * std::sqrt(TWOPI) * mvphi(-b / a) * b
                      * (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
            }
            a = a / 2.0;
            for (int i = 0; i < lg; ++i) {
                double xs = (a * (X[i][ng] + 1.0)) * (a * (X[i][ng] + 1.0));
                double rs = std::sqrt(1.0 - xs);
                bvn += a * W[i][ng] *
                       (std::exp(-bs / (2.0 * xs) - hk / (1.0 + rs)) / rs
                      - std::exp(-(bs / xs + hk) / 2.0) * (1.0 + c * xs * (1.0 + d * xs)));
                xs = as * (-X[i][ng] + 1.0) * (-X[i][ng] + 1.0) / 4.0;
                rs = std::sqrt(1.0 - xs);
                bvn += a * W[i][ng] * std::exp(-(bs / xs + hk) / 2.0)
                       * (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs
                        - (1.0 + c * xs * (1.0 + d * xs)));
            }
            bvn = -bvn / TWOPI;
        }
        if (r > 0.0) {
            bvn += mvphi(-std::max(h, k));
        } else {
            bvn = -bvn + std::max(0.0, mvphi(-h) - mvphi(-k));
        }
    }
    return bvn;
}

// ─── Bivariate normal probability with limit flags ──────────────────────────
inline double mvbvn(const double* lower, const double* upper, const int* infin, double correl) {
    if (infin[0] == 2 && infin[1] == 2) {
        return mvbvu(lower[0], lower[1], correl)
             - mvbvu(upper[0], lower[1], correl)
             - mvbvu(lower[0], upper[1], correl)
             + mvbvu(upper[0], upper[1], correl);
    } else if (infin[0] == 2 && infin[1] == 1) {
        return mvbvu(lower[0], lower[1], correl) - mvbvu(upper[0], lower[1], correl);
    } else if (infin[0] == 1 && infin[1] == 2) {
        return mvbvu(lower[0], lower[1], correl) - mvbvu(lower[0], upper[1], correl);
    } else if (infin[0] == 2 && infin[1] == 0) {
        return mvbvu(-upper[0], -upper[1], correl) - mvbvu(-lower[0], -upper[1], correl);
    } else if (infin[0] == 0 && infin[1] == 2) {
        return mvbvu(-upper[0], -upper[1], correl) - mvbvu(-upper[0], -lower[1], correl);
    } else if (infin[0] == 1 && infin[1] == 0) {
        return mvbvu(lower[0], -upper[1], -correl);
    } else if (infin[0] == 0 && infin[1] == 1) {
        return mvbvu(-upper[0], lower[1], -correl);
    } else if (infin[0] == 1 && infin[1] == 1) {
        return mvbvu(lower[0], lower[1], correl);
    } else if (infin[0] == 0 && infin[1] == 0) {
        return mvbvu(-upper[0], -upper[1], correl);
    } else {
        return 1.0;
    }
}

// ─── Bivariate t probability P(X < dh, Y < dk) ─────────────────────────────
inline double mvbvtl(int nu, double dh, double dk, double r) {
    double snu = std::sqrt(static_cast<double>(nu));
    double ors = 1.0 - r * r;
    double hrk = dh - r * dk;
    double krh = dk - r * dh;

    double xnhk, xnkh;
    if (std::abs(hrk) + ors > 0.0) {
        xnhk = hrk * hrk / (hrk * hrk + ors * (nu + dk * dk));
        xnkh = krh * krh / (krh * krh + ors * (nu + dh * dh));
    } else {
        xnhk = 0.0;
        xnkh = 0.0;
    }

    // Fortran SIGN(one, x) = copysign(1.0, x)
    int hs = (dh - r * dk >= 0.0) ? 1 : -1;
    int ks = (dk - r * dh >= 0.0) ? 1 : -1;

    double bvt;
    if (nu % 2 == 0) {
        // Even degrees of freedom
        bvt = std::atan2(std::sqrt(ors), -r) / TWOPI;
        double gmph = dh / std::sqrt(16.0 * (nu + dh * dh));
        double gmpk = dk / std::sqrt(16.0 * (nu + dk * dk));
        double btnckh = 2.0 * std::atan2(std::sqrt(xnkh), std::sqrt(1.0 - xnkh)) / PI;
        double btpdkh = 2.0 * std::sqrt(xnkh * (1.0 - xnkh)) / PI;
        double btnchk = 2.0 * std::atan2(std::sqrt(xnhk), std::sqrt(1.0 - xnhk)) / PI;
        double btpdhk = 2.0 * std::sqrt(xnhk * (1.0 - xnhk)) / PI;
        for (int j = 1; j <= nu / 2; ++j) {
            bvt += gmph * (1.0 + ks * btnckh);
            bvt += gmpk * (1.0 + hs * btnchk);
            btnckh += btpdkh;
            btpdkh = 2.0 * j * btpdkh * (1.0 - xnkh) / (2.0 * j + 1.0);
            btnchk += btpdhk;
            btpdhk = 2.0 * j * btpdhk * (1.0 - xnhk) / (2.0 * j + 1.0);
            gmph = gmph * (2.0 * j - 1.0) / (2.0 * j * (1.0 + dh * dh / nu));
            gmpk = gmpk * (2.0 * j - 1.0) / (2.0 * j * (1.0 + dk * dk / nu));
        }
    } else {
        // Odd degrees of freedom
        double qhrk = std::sqrt(dh * dh + dk * dk - 2.0 * r * dh * dk + nu * ors);
        double hkrn = dh * dk + r * nu;
        double hkn  = dh * dk - nu;
        double hpk  = dh + dk;
        bvt = std::atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - nu * hpk * qhrk) / TWOPI;
        if (bvt < -1e-15) bvt += 1.0;
        double gmph = dh / (TWOPI * snu * (1.0 + dh * dh / nu));
        double gmpk = dk / (TWOPI * snu * (1.0 + dk * dk / nu));
        double btnckh = std::sqrt(xnkh);
        double btpdkh = btnckh;
        double btnchk = std::sqrt(xnhk);
        double btpdhk = btnchk;
        for (int j = 1; j <= (nu - 1) / 2; ++j) {
            bvt += gmph * (1.0 + ks * btnckh);
            bvt += gmpk * (1.0 + hs * btnchk);
            btpdkh = (2.0 * j - 1.0) * btpdkh * (1.0 - xnkh) / (2.0 * j);
            btnckh += btpdkh;
            btpdhk = (2.0 * j - 1.0) * btpdhk * (1.0 - xnhk) / (2.0 * j);
            btnchk += btpdhk;
            gmph = 2.0 * j * gmph / ((2.0 * j + 1.0) * (1.0 + dh * dh / nu));
            gmpk = 2.0 * j * gmpk / ((2.0 * j + 1.0) * (1.0 + dk * dk / nu));
        }
    }
    return bvt;
}

// ─── Bivariate t probability with limit flags ───────────────────────────────
inline double mvbvt(int nu, const double* lower, const double* upper, const int* infin, double correl) {
    if (nu < 1) {
        return mvbvn(lower, upper, infin, correl);
    }
    if (infin[0] == 2 && infin[1] == 2) {
        return mvbvtl(nu, upper[0], upper[1], correl)
             - mvbvtl(nu, upper[0], lower[1], correl)
             - mvbvtl(nu, lower[0], upper[1], correl)
             + mvbvtl(nu, lower[0], lower[1], correl);
    } else if (infin[0] == 2 && infin[1] == 1) {
        return mvbvtl(nu, -lower[0], -lower[1], correl)
             - mvbvtl(nu, -upper[0], -lower[1], correl);
    } else if (infin[0] == 1 && infin[1] == 2) {
        return mvbvtl(nu, -lower[0], -lower[1], correl)
             - mvbvtl(nu, -lower[0], -upper[1], correl);
    } else if (infin[0] == 2 && infin[1] == 0) {
        return mvbvtl(nu, upper[0], upper[1], correl)
             - mvbvtl(nu, lower[0], upper[1], correl);
    } else if (infin[0] == 0 && infin[1] == 2) {
        return mvbvtl(nu, upper[0], upper[1], correl)
             - mvbvtl(nu, upper[0], lower[1], correl);
    } else if (infin[0] == 1 && infin[1] == 0) {
        return mvbvtl(nu, -lower[0], upper[1], -correl);
    } else if (infin[0] == 0 && infin[1] == 1) {
        return mvbvtl(nu, upper[0], -lower[1], -correl);
    } else if (infin[0] == 1 && infin[1] == 1) {
        return mvbvtl(nu, -lower[0], -lower[1], correl);
    } else if (infin[0] == 0 && infin[1] == 0) {
        return mvbvtl(nu, upper[0], upper[1], correl);
    } else {
        return 1.0;
    }
}

// ─── Limits computation ─────────────────────────────────────────────────────
inline void mvlims(double a, double b, int infin, double& lower, double& upper) {
    lower = 0.0;
    upper = 1.0;
    if (infin >= 0) {
        if (infin != 0) lower = mvphi(a);
        if (infin != 1) upper = mvphi(b);
    }
    upper = std::max(upper, lower);
}

// ─── Chi-squared quantile ───────────────────────────────────────────────────
namespace detail {

inline double mvchnc(double lkn, int n, double p, double r) {
    constexpr double LRP = -0.22579135264472743235;
    constexpr double EPS = 1e-14;

    double rr = r * r;
    double chi;

    if (n < 2) {
        chi = 2.0 * mvphi(-r);
    } else if (n < 100) {
        double rn = 1.0;
        for (int i = n - 2; i >= 2; i -= 2) {
            rn = 1.0 + rr * rn / i;
        }
        double rr2 = rr / 2.0;
        if (n % 2 == 0) {
            chi = std::exp(std::log(rn) - rr2);
        } else {
            chi = std::exp(LRP + std::log(r * rn) - rr2) + 2.0 * mvphi(-r);
        }
    } else {
        double rr2 = rr / 2.0;
        double al = static_cast<double>(n) / 2.0;
        chi = std::exp(-rr2 + al * std::log(rr2) + lkn + std::log(2.0) * (n - 2) / 2.0);
        if (rr2 < al + 1.0) {
            // Incomplete Gamma series
            double dl = chi;
            for (int i = 1; i <= 1000; ++i) {
                dl = dl * rr2 / (al + i);
                chi += dl;
                if (std::abs(dl * rr2 / (al + i + 1.0 - rr2)) < EPS) break;
            }
            chi = 1.0 - chi / al;
        } else {
            // Incomplete Gamma continued fraction
            double bi = rr2 + 1.0 - al;
            double ci = 1.0 / EPS;
            double di = bi;
            chi = chi / bi;
            for (int i = 1; i <= 250; ++i) {
                double ai = i * (al - i);
                bi += 2.0;
                ci = bi + ai / ci;
                if (ci == 0.0) ci = EPS;
                di = bi + ai / di;
                if (di == 0.0) di = EPS;
                double dl = ci / di;
                chi *= dl;
                if (std::abs(dl - 1.0) < EPS) break;
            }
        }
    }

    double df = (p - chi) / std::exp(lkn + (n - 1) * std::log(r) - rr / 2.0);
    return r - df * (1.0 - df * (r - static_cast<double>(n - 1) / r) / 2.0);
}

} // namespace detail

inline double mvchnv(int n, double p) {
    constexpr double LRP = -0.22579135264472743235;

    if (n <= 1) {
        return -mvphnv(p / 2.0);
    }
    if (p >= 1.0) return 0.0;

    double r;
    if (n == 2) {
        r = std::sqrt(-2.0 * std::log(p));
    } else {
        double lkn = 0.0;
        for (int i = n - 2; i >= 2; i -= 2) {
            lkn -= std::log(static_cast<double>(i));
        }
        if (n % 2 == 1) lkn += LRP;

        if (n >= -5.0 * std::log(1.0 - p) / 4.0) {
            r = 2.0 / (9.0 * n);
            r = n * std::pow(-mvphnv(p) * std::sqrt(r) + 1.0 - r, 3);
            if (r > 2.0 * n + 6.0) {
                r = 2.0 * (lkn - std::log(p)) + (n - 2) * std::log(r);
            }
        } else {
            r = std::exp((std::log((1.0 - p) * n) - lkn) * 2.0 / n);
        }
        r = std::sqrt(r);
        double ro = r;
        r = detail::mvchnc(lkn, n, p, r);
        if (std::abs(r - ro) > 1e-6) {
            ro = r;
            r = detail::mvchnc(lkn, n, p, r);
            if (std::abs(r - ro) > 1e-6) r = detail::mvchnc(lkn, n, p, r);
        }
    }
    return r;
}

// ─── Swap helpers ───────────────────────────────────────────────────────────
inline void mvsswp(double& x, double& y) {
    double t = x; x = y; y = t;
}

inline void mvswap(int p, int q, double* a, double* b, double* d,
                   int* infin, int n, double* c) {
    // Swap rows and columns P and Q in situ (0-based p, q with p <= q)
    // Note: p, q are 1-based in this implementation to match Fortran indexing
    mvsswp(a[p-1], a[q-1]);
    mvsswp(b[p-1], b[q-1]);
    mvsswp(d[p-1], d[q-1]);
    int tmp = infin[p-1];
    infin[p-1] = infin[q-1];
    infin[q-1] = tmp;

    int jj = (p * (p - 1)) / 2;
    int ii = (q * (q - 1)) / 2;
    mvsswp(c[jj + p - 1], c[ii + q - 1]);
    for (int j = 1; j <= p - 1; ++j) {
        mvsswp(c[jj + j - 1], c[ii + j - 1]);
    }
    jj += p;
    for (int i = p + 1; i <= q - 1; ++i) {
        mvsswp(c[jj + p - 1], c[ii + i - 1]);
        jj += i;
    }
    ii += q;
    for (int i = q + 1; i <= n; ++i) {
        mvsswp(c[ii + p - 1], c[ii + q - 1]);
        ii += i;
    }
}

// ─── MVSORT: Sort integration limits and compute Cholesky factor ────────────
inline void mvsort(int n, double* lower, double* upper, double* delta,
                   double* correl, int* infin, double* y, bool pivot,
                   int& nd, double* a, double* b, double* dl, double* cov,
                   int* infi, int& inform) {
    constexpr double EPS = 1e-10;
    inform = 0;
    int ij = 0, ii = 0;
    nd = n;

    for (int i = 1; i <= n; ++i) {
        a[i-1] = 0.0;
        b[i-1] = 0.0;
        dl[i-1] = 0.0;
        infi[i-1] = infin[i-1];
        if (infi[i-1] < 0) {
            nd--;
        } else {
            if (infi[i-1] != 0) a[i-1] = lower[i-1];
            if (infi[i-1] != 1) b[i-1] = upper[i-1];
            dl[i-1] = delta[i-1];
        }
        for (int j = 1; j <= i - 1; ++j) {
            cov[ij] = correl[ii];
            ij++;
            ii++;
        }
        cov[ij] = 1.0;
        ij++;
    }

    // Move doubly infinite limits to innermost positions
    if (nd > 0) {
        for (int i = n; i >= nd + 1; --i) {
            if (infi[i-1] >= 0) {
                for (int j = 1; j <= i - 1; ++j) {
                    if (infi[j-1] < 0) {
                        mvswap(j, i, a, b, dl, infi, n, cov);
                        break;
                    }
                }
            }
        }

        // Sort remaining limits and determine Cholesky factor
        ii = 0;
        int jl = nd;
        for (int i = 1; i <= nd; ++i) {
            double demin = 1.0;
            int jmin = i;
            double cvdiag = 0.0;
            ij = ii;
            double epsi = EPS * i;
            if (!pivot) jl = i;

            for (int j = i; j <= jl; ++j) {
                if (cov[ij + j - 1] > epsi) {
                    double sumsq = std::sqrt(cov[ij + j - 1]);
                    double sum = dl[j-1];
                    for (int k = 1; k <= i - 1; ++k) {
                        sum += cov[ij + k - 1] * y[k-1];
                    }
                    double aj = (a[j-1] - sum) / sumsq;
                    double bj = (b[j-1] - sum) / sumsq;
                    double d_val, e_val;
                    mvlims(aj, bj, infi[j-1], d_val, e_val);
                    if (demin >= e_val - d_val) {
                        jmin = j;
                        double amin = aj;
                        double bmin = bj;
                        demin = e_val - d_val;
                        cvdiag = sumsq;
                        (void)amin; (void)bmin; // used below after swap
                    }
                }
                ij += j;
            }

            if (jmin > i) {
                mvswap(i, jmin, a, b, dl, infi, n, cov);
            }
            if (cov[ii + i - 1] < -epsi) {
                inform = 3;
            }
            cov[ii + i - 1] = cvdiag;

            // Recompute amin, bmin after potential swap
            double amin_f, bmin_f;
            {
                ij = ii;
                double sum = dl[i-1];
                for (int k = 1; k <= i - 1; ++k) {
                    sum += cov[ij + k - 1] * y[k-1];
                }
                if (cvdiag > 0.0) {
                    amin_f = (a[i-1] - sum) / cvdiag;
                    bmin_f = (b[i-1] - sum) / cvdiag;
                } else {
                    amin_f = 0.0;
                    bmin_f = 0.0;
                }
            }

            // Compute Ith column of Cholesky factor
            if (cvdiag > 0.0) {
                int il = ii + i;
                for (int l = i + 1; l <= nd; ++l) {
                    cov[il + i - 1] = cov[il + i - 1] / cvdiag;
                    ij = ii + i;
                    for (int j = i + 1; j <= l; ++j) {
                        cov[il + j - 1] -= cov[il + i - 1] * cov[ij + i - 1];
                        ij += j;
                    }
                    il += l;
                }

                // Expected Y
                if (demin > epsi) {
                    y[i-1] = 0.0;
                    if (infi[i-1] != 0) y[i-1] = mvtdns(0, amin_f);
                    if (infi[i-1] != 1) y[i-1] = y[i-1] - mvtdns(0, bmin_f);
                    y[i-1] /= demin;
                } else {
                    if (infi[i-1] == 0) y[i-1] = bmin_f;
                    if (infi[i-1] == 1) y[i-1] = amin_f;
                    if (infi[i-1] == 2) y[i-1] = (amin_f + bmin_f) / 2.0;
                }
                for (int j = 1; j <= i; ++j) {
                    ii++;
                    cov[ii-1] = cov[ii-1] / cvdiag;
                }
                a[i-1] /= cvdiag;
                b[i-1] /= cvdiag;
                dl[i-1] /= cvdiag;
            } else {
                // Zero diagonal case
                int il = ii + i;
                for (int l = i + 1; l <= nd; ++l) {
                    cov[il + i - 1] = 0.0;
                    il += l;
                }

                // Permute limits and rows if necessary
                for (int j = i - 1; j >= 1; --j) {
                    if (std::abs(cov[ii + j - 1]) > epsi) {
                        a[i-1] /= cov[ii + j - 1];
                        b[i-1] /= cov[ii + j - 1];
                        dl[i-1] /= cov[ii + j - 1];
                        if (cov[ii + j - 1] < 0.0) {
                            mvsswp(a[i-1], b[i-1]);
                            if (infi[i-1] != 2) infi[i-1] = 1 - infi[i-1];
                        }
                        for (int l = 1; l <= j; ++l) {
                            cov[ii + l - 1] /= cov[ii + j - 1];
                        }
                        for (int l = j + 1; l <= i - 1; ++l) {
                            if (cov[(l-1)*l/2 + j] > 0.0) {
                                ij = ii;
                                for (int k = i - 1; k >= l; --k) {
                                    for (int m = 1; m <= k; ++m) {
                                        mvsswp(cov[ij - k + m - 1], cov[ij + m - 1]);
                                    }
                                    mvsswp(a[k-1], a[k]);
                                    mvsswp(b[k-1], b[k]);
                                    mvsswp(dl[k-1], dl[k]);
                                    int tm = infi[k-1];
                                    infi[k-1] = infi[k];
                                    infi[k] = tm;
                                    ij -= k;
                                }
                                break; // GO TO 20
                            }
                        }
                        break; // GO TO 20
                    }
                    cov[ii + j - 1] = 0.0;
                }
                // label 20
                ii += i;
                y[i-1] = 0.0;
            }
        }
    }
}

// ─── MVVLSB: Integrand subroutine ──────────────────────────────────────────
// Hot path: called for every lattice point. Keep the inner dot product as a
// plain scalar loop — vectors are tiny (avg len ~n/2 ≤ 5 for typical usage)
// and the compiler auto-vectorizes the scalar form just as well as Eigen would,
// without the Map-construction overhead per call.
inline void mvvlsb(int n, const double* w, double r,
                   const double* dl, const int* infi, const double* a,
                   const double* b, const double* cov, double* y,
                   int& nd, double& value) {
    value = 1.0;
    int infa = 0, infb = 0;
    nd = 0;
    int ij = 0;
    double ai = 0.0, bi = 0.0;

    for (int i = 1; i <= n; ++i) {
        double sum = dl[i-1];
        // Triangular mat-vec: sum += cov[ij..ij+nd-1] . y[0..nd-1]
        // Scalar loop: compiler auto-vectorizes; Eigen Map overhead not worth it
        // for the tiny lengths (avg ~2) typical in regenie (n=3..5).
        int len = std::min(i - 1, nd);
        for (int j = 0; j < len; ++j) {
            sum += cov[ij + j] * y[j];
        }
        ij += i - 1; // ij now points at diagonal of row i

        if (infi[i-1] != 0) {
            if (infa == 1) {
                ai = std::max(ai, r * a[i-1] - sum);
            } else {
                ai = r * a[i-1] - sum;
                infa = 1;
            }
        }
        if (infi[i-1] != 1) {
            if (infb == 1) {
                bi = std::min(bi, r * b[i-1] - sum);
            } else {
                bi = r * b[i-1] - sum;
                infb = 1;
            }
        }

        ij++; // past diagonal
        // Check if we should emit a variable
        if (i == n || cov[ij + nd + 1] > 0.0) {
            double di, ei;
            mvlims(ai, bi, infa + infa + infb - 1, di, ei);
            if (di >= ei) {
                value = 0.0;
                return;
            } else {
                value *= (ei - di);
                nd++;
                if (i < n) y[nd-1] = mvphnv(di + w[nd-1] * (ei - di));
                infa = 0;
                infb = 0;
            }
        }
    }
}

// ─── MVSPCL: Special cases ──────────────────────────────────────────────────
inline void mvspcl(int& nd, int nu, double* a, double* b, double* dl,
                   double* cov, int* infi, double& snu, double& vl,
                   double& er, int& inform) {
    if (inform > 0) {
        vl = 0.0;
        er = 1.0;
        return;
    }

    if (nd == 0) {
        er = 0.0;
        vl = 1.0; // Bug fix 24/03/2009
    } else if (nd == 1 && (nu < 1 || std::abs(dl[0]) == 0.0)) {
        // 1-d case for normal or central t
        vl = 1.0;
        if (infi[0] != 1) vl = mvstdt(nu, b[0] - dl[0]);
        if (infi[0] != 0) vl = vl - mvstdt(nu, a[0] - dl[0]);
        if (vl < 0.0) vl = 0.0;
        er = 2e-16;
        nd = 0;
    } else if (nd == 2 && (nu < 1 || std::abs(dl[0]) + std::abs(dl[1]) == 0.0)) {
        // 2-d case for normal or central t
        if (infi[0] != 0) a[0] -= dl[0];
        if (infi[0] != 1) b[0] -= dl[0];
        if (infi[1] != 0) a[1] -= dl[1];
        if (infi[1] != 1) b[1] -= dl[1];
        if (std::abs(cov[2]) > 0.0) {
            // 2-d nonsingular case
            double r = std::sqrt(1.0 + cov[1] * cov[1]);
            if (infi[1] != 0) a[1] /= r;
            if (infi[1] != 1) b[1] /= r;
            cov[1] /= r;
            vl = mvbvt(nu, a, b, infi, cov[1]);
            er = 1e-15;
        } else {
            // 2-d singular case
            if (infi[0] != 0) {
                if (infi[1] != 0) a[0] = std::max(a[0], a[1]);
            } else {
                if (infi[1] != 0) a[0] = a[1];
            }
            if (infi[0] != 1) {
                if (infi[1] != 1) b[0] = std::min(b[0], b[1]);
            } else {
                if (infi[1] != 1) b[0] = b[1];
            }
            if (infi[0] != infi[1]) infi[0] = 2;
            vl = 1.0;
            // Bug Fixed 28/05/2013: use b[0], a[0] directly (not subtract dl)
            if (infi[0] != 1) vl = mvstdt(nu, b[0]);
            if (infi[0] != 0) vl = vl - mvstdt(nu, a[0]);
            if (vl < 0.0) vl = 0.0;
            er = 2e-16;
        }
        nd = 0;
    } else {
        if (nu > 0) {
            snu = std::sqrt(static_cast<double>(nu));
        } else {
            nd--;
        }
    }
}

// ─── Lattice rule primes and generating vectors ─────────────────────────────
namespace lattice {

static const int P[28] = {
    31, 47, 73, 113, 173, 263, 397, 593, 907, 1361, 2053, 3079, 4621, 6947, 10427, 15641, 23473, 35221, 52837, 79259, 118891, 178349, 267523, 401287, 601943, 902933, 1354471, 2031713
};

static const int C[28][99] = {
    {12, 9, 9, 13, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3,
     12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3,
     12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3,
     12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3,
     12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 7, 3, 3, 3,
     7, 7, 7, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
     3, 3, 3, 3, 3, 3, 3, 3, 3},
    {13, 11, 17, 10, 15, 15, 15, 15, 15, 15, 22, 15, 15, 6, 6,
     6, 15, 15, 9, 13, 2, 2, 2, 13, 11, 11, 10, 15, 15, 15,
     15, 15, 15, 15, 15, 15, 6, 6, 6, 15, 15, 9, 13, 2, 2,
     2, 13, 11, 11, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 6,
     6, 6, 15, 15, 9, 13, 2, 2, 2, 13, 11, 11, 10, 10, 15,
     15, 15, 15, 15, 15, 15, 15, 6, 2, 3, 2, 3, 2, 2, 2,
     2, 2, 2, 2, 2, 2, 2, 2, 2},
    {27, 28, 10, 11, 11, 20, 11, 11, 28, 13, 13, 28, 13, 13, 13,
     14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
     14, 31, 31, 5, 5, 5, 31, 13, 11, 11, 11, 11, 11, 11, 13,
     13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14,
     14, 14, 14, 14, 14, 14, 14, 31, 31, 5, 5, 5, 11, 13, 11,
     11, 11, 11, 11, 11, 11, 13, 13, 11, 13, 5, 5, 5, 5, 14,
     13, 5, 5, 5, 5, 5, 5, 5, 5},
    {35, 27, 27, 36, 22, 29, 29, 20, 45, 5, 5, 5, 21, 21, 21,
     21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 29, 17,
     17, 17, 17, 17, 17, 17, 17, 17, 17, 23, 23, 23, 23, 23, 23,
     23, 23, 23, 23, 23, 23, 21, 27, 3, 3, 3, 24, 27, 27, 17,
     29, 29, 29, 17, 5, 5, 5, 5, 21, 21, 21, 21, 21, 21, 21,
     21, 21, 21, 21, 21, 21, 21, 21, 21, 17, 17, 17, 6, 17, 17,
     6, 3, 6, 6, 3, 3, 3, 3, 3},
    {64, 66, 28, 28, 44, 44, 55, 67, 10, 10, 10, 10, 10, 10, 38,
     38, 10, 10, 10, 10, 10, 49, 49, 49, 49, 49, 49, 49, 49, 49,
     49, 49, 49, 38, 38, 31, 4, 4, 31, 64, 4, 4, 4, 64, 45,
     45, 45, 45, 45, 45, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66,
     66, 66, 66, 66, 66, 66, 66, 66, 66, 11, 66, 66, 66, 66, 66,
     66, 66, 66, 66, 45, 11, 7, 3, 2, 2, 2, 27, 5, 3, 3,
     5, 5, 2, 2, 2, 2, 2, 2, 2},
    {111, 42, 54, 118, 20, 31, 31, 72, 17, 94, 14, 14, 11, 14, 14,
     14, 94, 10, 10, 10, 10, 14, 14, 14, 14, 14, 14, 14, 11, 11,
     11, 8, 8, 8, 8, 8, 8, 8, 18, 18, 18, 18, 18, 113, 62,
     62, 45, 45, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113,
     113, 113, 113, 113, 113, 63, 63, 53, 63, 67, 67, 67, 67, 67, 67,
     67, 67, 67, 67, 67, 67, 67, 67, 67, 51, 51, 51, 51, 51, 12,
     51, 12, 51, 5, 3, 3, 2, 2, 5},
    {163, 154, 83, 43, 82, 92, 150, 59, 76, 76, 47, 11, 11, 100, 131,
     116, 116, 116, 116, 116, 116, 138, 138, 138, 138, 138, 138, 138, 138, 138,
     101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101,
     101, 101, 101, 101, 101, 101, 116, 116, 116, 116, 116, 116, 100, 100, 100,
     100, 100, 138, 138, 138, 138, 138, 101, 101, 101, 101, 101, 101, 101, 101,
     101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 38, 38, 38, 38,
     38, 38, 38, 38, 3, 3, 3, 3, 3},
    {246, 189, 242, 102, 250, 250, 102, 250, 280, 118, 196, 118, 191, 215, 121,
     121, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 171, 171,
     171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171,
     171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171,
     171, 171, 161, 161, 161, 161, 161, 161, 161, 161, 14, 14, 14, 14, 14,
     14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 10, 10, 10,
     10, 10, 10, 103, 10, 10, 10, 10, 5},
    {347, 402, 322, 418, 215, 220, 339, 339, 339, 337, 218, 315, 315, 315, 315,
     167, 167, 167, 167, 361, 201, 124, 124, 124, 124, 124, 124, 124, 124, 124,
     124, 124, 231, 231, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90,
     90, 90, 90, 48, 48, 48, 48, 90, 90, 90, 90, 90, 90, 90, 90,
     90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90,
     243, 243, 243, 243, 243, 243, 243, 243, 243, 243, 283, 283, 283, 283, 283,
     283, 283, 283, 283, 16, 283, 16, 283, 283},
    {505, 220, 601, 644, 612, 160, 206, 206, 206, 422, 134, 518, 134, 134, 518,
     652, 382, 206, 158, 441, 179, 441, 56, 559, 559, 56, 56, 56, 56, 56,
     56, 56, 56, 56, 56, 56, 56, 56, 56, 101, 101, 56, 101, 101, 101,
     101, 101, 101, 101, 101, 193, 193, 193, 193, 193, 193, 193, 101, 101, 101,
     101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101,
     101, 101, 101, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122,
     122, 122, 122, 122, 122, 101, 101, 101, 101},
    {794, 325, 960, 528, 247, 247, 338, 366, 847, 753, 753, 236, 334, 334, 461,
     711, 652, 381, 381, 381, 652, 381, 381, 381, 381, 381, 381, 381, 226, 326,
     326, 326, 326, 326, 326, 326, 126, 326, 326, 326, 326, 326, 326, 326, 326,
     326, 326, 195, 195, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
     55, 55, 55, 55, 55, 55, 55, 55, 195, 195, 195, 195, 195, 195, 195,
     132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 387, 387, 387, 387,
     387, 387, 387, 387, 387, 387, 387, 387, 387},
    {1189, 888, 259, 1082, 725, 811, 636, 965, 497, 497, 1490, 1490, 392, 1291, 508,
     508, 1291, 1291, 508, 1291, 508, 508, 867, 867, 867, 867, 934, 867, 867, 867,
     867, 867, 867, 867, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 563, 563,
     563, 563, 1010, 1010, 1010, 208, 838, 563, 563, 563, 759, 759, 564, 759, 759,
     801, 801, 801, 801, 759, 759, 759, 759, 759, 563, 563, 563, 563, 563, 563,
     563, 563, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226,
     226, 226, 226, 226, 226, 226, 226, 226, 226},
    {1763, 1018, 1500, 432, 1332, 2203, 126, 2240, 1719, 1284, 878, 1983, 266, 266, 266,
     266, 747, 747, 127, 127, 2074, 127, 2074, 1400, 1383, 1383, 1383, 1383, 1383, 1383,
     1383, 1383, 1383, 1383, 1400, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 507, 1073, 1073,
     1073, 1073, 1990, 1990, 1990, 1990, 1990, 507, 507, 507, 507, 507, 507, 507, 507,
     507, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073,
     1073, 1073, 1073, 22, 22, 22, 22, 22, 22, 1073, 452, 452, 452, 452, 452,
     452, 318, 301, 301, 301, 301, 86, 86, 15},
    {2872, 3233, 1534, 2941, 2910, 393, 1796, 919, 446, 919, 919, 1117, 103, 103, 103,
     103, 103, 103, 103, 2311, 3117, 1101, 3117, 3117, 1101, 1101, 1101, 1101, 1101, 2503,
     2503, 2503, 2503, 2503, 2503, 2503, 2503, 429, 429, 429, 429, 429, 429, 429, 1702,
     1702, 1702, 184, 184, 184, 184, 184, 105, 105, 105, 105, 105, 105, 105, 105,
     105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105,
     105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 784, 784, 784, 784,
     784, 784, 784, 784, 784, 784, 784, 784, 784},
    {4309, 3758, 4034, 1963, 730, 642, 1502, 2246, 3834, 1511, 1102, 1102, 1522, 1522, 3427,
     3427, 3928, 915, 915, 3818, 3818, 3818, 3818, 4782, 4782, 4782, 3818, 4782, 3818, 3818,
     1327, 1327, 1327, 1327, 1327, 1327, 1327, 1387, 1387, 1387, 1387, 1387, 1387, 1387, 1387,
     1387, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 3148,
     3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148,
     3148, 3148, 1776, 1776, 1776, 3354, 3354, 3354, 925, 3354, 3354, 925, 925, 925, 925,
     925, 2133, 2133, 2133, 2133, 2133, 2133, 2133, 2133},
    {6610, 6977, 1686, 3819, 2314, 5647, 3953, 3614, 5115, 423, 423, 5408, 7426, 423, 423,
     487, 6227, 2660, 6227, 1221, 3811, 197, 4367, 351, 1281, 1221, 351, 351, 351, 7245,
     1984, 2999, 2999, 2999, 2999, 2999, 2999, 3995, 2063, 2063, 2063, 2063, 1644, 2063, 2077,
     2512, 2512, 2512, 2077, 2077, 2077, 2077, 754, 754, 754, 754, 754, 754, 754, 754,
     754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 1097, 1097, 754, 754,
     754, 754, 248, 754, 1097, 1097, 1097, 1097, 222, 222, 222, 222, 754, 1982, 1982,
     1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982},
    {9861, 3647, 4073, 2535, 3430, 9865, 2830, 9328, 4320, 5913, 10365, 8272, 3706, 6186, 7806,
     7806, 7806, 8610, 2563, 11558, 11558, 9421, 1181, 9421, 1181, 1181, 1181, 9421, 1181, 1181,
     10574, 10574, 3534, 3534, 3534, 3534, 3534, 2898, 2898, 2898, 3450, 2141, 2141, 2141, 2141,
     2141, 2141, 2141, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055,
     7055, 7055, 7055, 2831, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204,
     8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 4688, 4688,
     4688, 2831, 2831, 2831, 2831, 2831, 2831, 2831, 2831},
    {10327, 7582, 7124, 8214, 9600, 10271, 10193, 10800, 9086, 2365, 4409, 13812, 5661, 9344, 9344,
     10362, 9344, 9344, 8585, 11114, 13080, 13080, 13080, 6949, 3436, 3436, 3436, 13213, 6130, 6130,
     8159, 8159, 11595, 8159, 3436, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096,
     7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 4377, 7096, 4377, 4377, 4377, 4377, 4377,
     5410, 5410, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377,
     4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377,
     4377, 4377, 4377, 4377, 440, 440, 1199, 1199, 1199},
    {19540, 19926, 11582, 11113, 24585, 8726, 17218, 419, 4918, 4918, 4918, 15701, 17710, 4037, 4037,
     15808, 11401, 19398, 25950, 25950, 4454, 24987, 11719, 8697, 1452, 1452, 1452, 1452, 1452, 8697,
     8697, 6436, 21475, 6436, 22913, 6434, 18497, 11089, 11089, 11089, 11089, 3036, 3036, 14208, 14208,
     14208, 14208, 12906, 12906, 12906, 12906, 12906, 12906, 12906, 12906, 7614, 7614, 7614, 7614, 5021,
     5021, 5021, 5021, 5021, 5021, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145,
     10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 4544,
     4544, 4544, 4544, 4544, 4544, 8394, 8394, 8394, 8394},
    {34566, 9579, 12654, 26856, 37873, 38806, 29501, 17271, 3663, 10763, 18955, 1298, 26560, 17132, 17132,
     4753, 4753, 8713, 18624, 13082, 6791, 1122, 19363, 34695, 18770, 18770, 18770, 18770, 15628, 18770,
     18770, 18770, 18770, 33766, 20837, 20837, 20837, 20837, 20837, 20837, 6545, 6545, 6545, 6545, 6545,
     12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 30483,
     30483, 30483, 30483, 30483, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138,
     12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 9305, 11107, 11107, 11107, 11107, 11107, 11107,
     11107, 11107, 11107, 11107, 11107, 11107, 11107, 9305, 9305},
    {31929, 49367, 10982, 3527, 27066, 13226, 56010, 18911, 40574, 20767, 20767, 9686, 47603, 47603, 11736,
     11736, 41601, 12888, 32948, 30801, 44243, 53351, 53351, 16016, 35086, 35086, 32581, 2464, 2464, 49554,
     2464, 2464, 49554, 49554, 2464, 81, 27260, 10681, 2185, 2185, 2185, 2185, 2185, 2185, 2185,
     18086, 18086, 18086, 18086, 18086, 17631, 17631, 18086, 18086, 18086, 37335, 37774, 37774, 37774, 26401,
     26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 12982, 40398, 40398,
     40398, 40398, 40398, 40398, 3518, 3518, 3518, 37799, 37799, 37799, 37799, 37799, 37799, 37799, 37799,
     37799, 4721, 4721, 4721, 4721, 7067, 7067, 7067, 7067},
    {40701, 69087, 77576, 64590, 39397, 33179, 10858, 38935, 43129, 35468, 35468, 5279, 61518, 61518, 27945,
     70975, 70975, 86478, 86478, 20514, 20514, 73178, 73178, 43098, 43098, 4701, 59979, 59979, 58556, 69916,
     15170, 15170, 4832, 4832, 43064, 71685, 4832, 15170, 15170, 15170, 27679, 27679, 27679, 60826, 60826,
     6187, 6187, 4264, 4264, 4264, 4264, 4264, 45567, 32269, 32269, 32269, 32269, 62060, 62060, 62060,
     62060, 62060, 62060, 62060, 62060, 62060, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803,
     1803, 1803, 1803, 1803, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108,
     51108, 55315, 55315, 54140, 54140, 54140, 54140, 54140, 13134},
    {103650, 125480, 59978, 46875, 77172, 83021, 126904, 14541, 56299, 43636, 11655, 52680, 88549, 29804, 101894,
     113675, 48040, 113675, 34987, 48308, 97926, 5475, 49449, 6850, 62545, 62545, 9440, 33242, 9440, 33242,
     9440, 33242, 9440, 62850, 9440, 9440, 9440, 90308, 90308, 90308, 47904, 47904, 47904, 47904, 47904,
     47904, 47904, 47904, 47904, 41143, 41143, 41143, 41143, 41143, 41143, 41143, 36114, 36114, 36114, 36114,
     36114, 24997, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162,
     65162, 47650, 47650, 47650, 47650, 47650, 47650, 47650, 40586, 40586, 40586, 40586, 40586, 40586, 40586,
     38725, 38725, 38725, 38725, 88329, 88329, 88329, 88329, 88329},
    {165843, 90647, 59925, 189541, 67647, 74795, 68365, 167485, 143918, 74912, 167289, 75517, 8148, 172106, 126159,
     35867, 35867, 35867, 121694, 52171, 95354, 113969, 113969, 76304, 123709, 123709, 144615, 123709, 64958, 64958,
     32377, 193002, 193002, 25023, 40017, 141605, 189165, 189165, 141605, 189165, 189165, 141605, 141605, 141605, 189165,
     127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047,
     127047, 127047, 127047, 127047, 127047, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785,
     80822, 80822, 80822, 80822, 80822, 80822, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661,
     131661, 131661, 131661, 131661, 131661, 131661, 131661, 7114, 131661},
    {130365, 236711, 110235, 125699, 56483, 93735, 234469, 60549, 1291, 93937, 245291, 196061, 258647, 162489, 176631,
     204895, 73353, 172319, 28881, 136787, 122081, 122081, 275993, 64673, 211587, 211587, 211587, 282859, 282859, 211587,
     242821, 256865, 256865, 256865, 122203, 291915, 122203, 291915, 291915, 122203, 25639, 25639, 291803, 245397, 284047,
     245397, 245397, 245397, 245397, 245397, 245397, 245397, 94241, 66575, 66575, 217673, 217673, 217673, 217673, 217673,
     217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 210249,
     210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 94453, 94453, 94453, 94453, 94453, 94453,
     94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453},
    {333459, 375354, 102417, 383544, 292630, 41147, 374614, 48032, 435453, 281493, 358168, 114121, 346892, 238990, 317313,
     164158, 35497, 70530, 70530, 434839, 24754, 24754, 24754, 393656, 118711, 118711, 148227, 271087, 355831, 91034,
     417029, 417029, 91034, 91034, 417029, 91034, 299843, 299843, 413548, 413548, 308300, 413548, 413548, 413548, 308300,
     308300, 308300, 413548, 308300, 308300, 308300, 308300, 308300, 15311, 15311, 15311, 15311, 176255, 176255, 23613,
     23613, 23613, 23613, 23613, 23613, 172210, 204328, 204328, 204328, 204328, 121626, 121626, 121626, 121626, 121626,
     200187, 200187, 200187, 200187, 200187, 121551, 121551, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492,
     248492, 248492, 248492, 248492, 13942, 13942, 13942, 13942, 13942},
    {500884, 566009, 399251, 652979, 355008, 430235, 328722, 670680, 405585, 405585, 424646, 670180, 670180, 641587, 215580,
     59048, 633320, 81010, 20789, 389250, 389250, 638764, 638764, 389250, 389250, 398094, 80846, 147776, 147776, 296177,
     398094, 398094, 147776, 147776, 396313, 578233, 578233, 578233, 19482, 620706, 187095, 620706, 187095, 126467, 241663,
     241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 321632, 23210, 23210, 394484,
     394484, 394484, 78101, 78101, 78101, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095,
     542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 277743, 277743, 277743, 457259, 457259, 457259,
     457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259},
    {858339, 918142, 501970, 234813, 460565, 31996, 753018, 256150, 199809, 993599, 245149, 794183, 121349, 150619, 376952,
     809123, 809123, 804319, 67352, 969594, 434796, 969594, 804319, 391368, 761041, 754049, 466264, 754049, 754049, 466264,
     754049, 754049, 282852, 429907, 390017, 276645, 994856, 250142, 144595, 907454, 689648, 687580, 687580, 687580, 687580,
     978368, 687580, 552742, 105195, 942843, 768249, 307142, 307142, 307142, 307142, 880619, 880619, 880619, 880619, 880619,
     880619, 880619, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 60731, 60731,
     60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 178309, 178309, 178309, 178309, 74373, 74373,
     74373, 74373, 74373, 74373, 74373, 74373, 214965, 214965, 214965}
};

} // namespace lattice

// ─── Uniform random number generator (thread-local) ─────────────────────────
inline double mvuni() {
    thread_local std::mt19937_64 gen(std::random_device{}());
    thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}

// ─── MVKRSV: Lattice rule sum (single sample) ──────────────────────────────
// Templated on Func to allow the compiler to inline the integrand lambda,
// eliminating std::function virtual dispatch (~5-10 ns per integrand call).
template<typename Func>
inline void mvkrsv(int ndim, int kl, double* values, int prime,
                   const double* vk, int nf,
                   Func&& funsub,
                   double* x, double* r, int* pr, double* fs) {
    for (int j = 0; j < nf; ++j) values[j] = 0.0;

    // Random shifts and scramble
    for (int j = 0; j < ndim; ++j) {
        r[j] = mvuni();
        if (j < kl - 1) {
            int jp = static_cast<int>(1 + (j + 1) * r[j]);  // 1-based destination
            jp = std::min(jp, j + 1);  // clamp
            if (jp - 1 < j) {
                pr[j] = pr[jp - 1];
            }
            pr[jp - 1] = j;
        } else {
            pr[j] = j;
        }
    }

    // Precompute permuted vk once per sample so the hot k-loop has no
    // gather indirection and the inner j-loop can be auto-vectorized.
    // Algorithm spec: ndim ≤ 100, so a fixed-size stack buffer is fine.
    double vk_perm[100];
    for (int j = 0; j < ndim; ++j) vk_perm[j] = vk[pr[j]];

    // Compute lattice rule sums
    for (int k = 1; k <= prime; ++k) {
        for (int j = 0; j < ndim; ++j) {
            r[j] += vk_perm[j];
            if (r[j] >= 1.0) r[j] -= 1.0;
            x[j] = std::abs(2.0 * r[j] - 1.0);
        }
        funsub(ndim, x, nf, fs);
        for (int j = 0; j < nf; ++j) {
            values[j] += (fs[j] - values[j]) / (2.0 * k - 1.0);
        }
        for (int j = 0; j < ndim; ++j) {
            x[j] = 1.0 - x[j];
        }
        funsub(ndim, x, nf, fs);
        for (int j = 0; j < nf; ++j) {
            values[j] += (fs[j] - values[j]) / (2.0 * k);
        }
    }
}

// ─── MVKBRV: Korobov lattice rule integration ──────────────────────────────
// Templated on Func to allow inlining the integrand (see mvkrsv).
template<typename Func>
inline void mvkbrv(int ndim, int& minvls, int maxvls, int nf,
                   Func&& funsub,
                   double abseps, double releps, double* abserr,
                   double* finest, int& inform) {
    inform = 1;
    int intvls = 0;
    double varprd = 0.0;
    int sampls = MINSMP;
    int np_idx = 0;

    std::vector<double> varest(nf, 0.0);
    std::vector<double> finval(nf, 0.0);
    std::vector<double> varsqr(nf, 0.0);

    if (minvls >= 0) {
        for (int k = 0; k < nf; ++k) {
            finest[k] = 0.0;
            varest[k] = 0.0;
        }
        sampls = MINSMP;
        for (int i = std::min(ndim, 10) - 1; i < PLIM; ++i) {
            np_idx = i;
            if (minvls < 2 * sampls * lattice::P[i]) break;
        }
        sampls = std::max(MINSMP, minvls / (2 * lattice::P[np_idx]));
    }

    // Working arrays
    std::vector<double> vk(ndim);

    // Main integration loop
    while (true) {
        // Setup lattice generating vector
        vk[0] = 1.0 / lattice::P[np_idx];
        if (ndim > 1) {
            int k = 1;
            for (int i = 1; i < ndim; ++i) {
                if (i < KLIM) {
                    k = static_cast<int>(
                        std::fmod(static_cast<double>(lattice::C[np_idx][std::min(ndim - 2, KLIM - 2)]) * k,
                                  static_cast<double>(lattice::P[np_idx])));
                    vk[i] = k * vk[0];
                } else {
                    vk[i] = static_cast<int>(
                        lattice::P[np_idx] * std::pow(2.0, static_cast<double>(i - KLIM + 1) / (ndim - KLIM + 1)));
                    vk[i] = std::fmod(vk[i] / lattice::P[np_idx], 1.0);
                }
            }
        }

        for (int k = 0; k < nf; ++k) {
            finval[k] = 0.0;
            varsqr[k] = 0.0;
        }

        // Sample loop — parallelized with OpenMP
        // Each sample is an independent randomized lattice rule evaluation.
        // We parallelize across samples and reduce to compute mean + variance.
#ifdef _OPENMP
        if (sampls >= 4) {
            // Parallel path: accumulate sum and sum-of-squares across threads
            double total_sum = 0.0;
            double total_sumsq = 0.0;

            #pragma omp parallel reduction(+:total_sum, total_sumsq)
            {
                // Per-thread working arrays
                std::vector<double> t_values(nf);
                std::vector<double> t_x(ndim), t_r(ndim), t_fs(nf);
                std::vector<int> t_pr(ndim);

                #pragma omp for schedule(static)
                for (int i = 0; i < sampls; ++i) {
                    mvkrsv(ndim, KLIM, t_values.data(), lattice::P[np_idx], vk.data(),
                           nf, funsub, t_x.data(), t_r.data(), t_pr.data(), t_fs.data());
                    total_sum += t_values[0];
                    total_sumsq += t_values[0] * t_values[0];
                }
            }

            for (int k = 0; k < nf; ++k) {
                finval[k] = total_sum / sampls;
                varsqr[k] = (total_sumsq - total_sum * total_sum / sampls) / (sampls * (sampls - 1));
            }
        } else
#endif
        {
            // Sequential path (small sample count or no OpenMP)
            std::vector<double> values(nf);
            std::vector<double> x(ndim), r(ndim), fs(nf);
            std::vector<int> pr(ndim);

            for (int i = 1; i <= sampls; ++i) {
                mvkrsv(ndim, KLIM, values.data(), lattice::P[np_idx], vk.data(),
                       nf, funsub, x.data(), r.data(), pr.data(), fs.data());
                for (int k = 0; k < nf; ++k) {
                    double difint = (values[k] - finval[k]) / i;
                    finval[k] += difint;
                    varsqr[k] = (i - 2) * varsqr[k] / i + difint * difint;
                }
            }
        }

        intvls += 2 * sampls * lattice::P[np_idx];
        int kmx = 0;
        for (int k = 0; k < nf; ++k) {
            varprd = varest[k] * varsqr[k];
            finest[k] += (finval[k] - finest[k]) / (1.0 + varprd);
            if (varsqr[k] > 0.0) varest[k] = (1.0 + varprd) / varsqr[k];
            if (std::abs(finest[k]) > std::abs(finest[kmx])) kmx = k;
        }
        *abserr = 7.0 * std::sqrt(varsqr[kmx] / (1.0 + varprd)) / 2.0;

        if (*abserr <= std::max(abseps, std::abs(finest[kmx]) * releps)) {
            inform = 0;
            break;
        }

        if (np_idx < PLIM - 1) {
            np_idx++;
        } else {
            sampls = std::min(3 * sampls / 2, (maxvls - intvls) / (2 * lattice::P[np_idx]));
            sampls = std::max(MINSMP, sampls);
        }

        if (intvls + 2 * sampls * lattice::P[np_idx] > maxvls) break;
    }

    minvls = intvls;
}

// ─── Main entry point: mvtdst ───────────────────────────────────────────────
inline void mvtdst(int n, int nu, double* lower, double* upper, int* infin,
                   double* correl, double* delta, int maxpts, double abseps,
                   double releps, double* error, double* value, int* inform) {
    if (n > 1000 || n < 1) {
        *value = 0.0;
        *error = 1.0;
        *inform = 2;
        return;
    }

    // Reuse thread-local context buffers to avoid heap allocation on every call.
    // Fortran uses COMMON blocks (static global); this gives the same benefit
    // while remaining thread-safe. Vectors grow on demand, never shrink.
    static thread_local MvtContext ctx;
    MvtContext* ctxp = &ctx;  // pointer so the lambda can capture it
    int cov_n = n * (n + 1) / 2;
    if (static_cast<int>(ctx.infi.size()) < n)  ctx.infi.resize(n);
    if (static_cast<int>(ctx.a.size())    < n)  ctx.a.resize(n);
    if (static_cast<int>(ctx.b.size())    < n)  ctx.b.resize(n);
    if (static_cast<int>(ctx.dl.size())   < n)  ctx.dl.resize(n);
    if (static_cast<int>(ctx.cov.size())  < cov_n) ctx.cov.resize(cov_n);
    if (static_cast<int>(ctx.y.size())    < n)  ctx.y.resize(n);

    // Initialize: sort and compute Cholesky
    mvsort(n, lower, upper, delta, correl, infin, ctx.y.data(), true,
           ctx.nd, ctx.a.data(), ctx.b.data(), ctx.dl.data(), ctx.cov.data(),
           ctx.infi.data(), *inform);
    ctx.nu = nu;

    // Handle special cases
    double vl = 0.0, er = 0.0;
    mvspcl(ctx.nd, ctx.nu, ctx.a.data(), ctx.b.data(), ctx.dl.data(),
           ctx.cov.data(), ctx.infi.data(), ctx.snu, vl, er, *inform);

    if (*inform == 0 && ctx.nd > 0) {
        // General case: lattice rule integration
        int nd = ctx.nd;
        int ivls = 0;

        // Capture ctxp (plain pointer to thread_local) — valid for the duration of this call.
        int y_size = static_cast<int>(ctx.y.size());
        auto integrand = [ctxp, nd, y_size](int n_dim, const double* w, int nf, double* f) {
            // Each thread gets its own y scratch buffer.
            thread_local std::vector<double> tl_y;
            if (static_cast<int>(tl_y.size()) < y_size) tl_y.resize(y_size, 0.0);

            double r;
            int nd_out;
            if (ctxp->nu <= 0) {
                r = 1.0;
                mvvlsb(nd + 1, w, r, ctxp->dl.data(), ctxp->infi.data(),
                       ctxp->a.data(), ctxp->b.data(), ctxp->cov.data(),
                       tl_y.data(), nd_out, f[0]);
            } else {
                r = mvchnv(ctxp->nu, w[n_dim - 1]) / ctxp->snu;
                mvvlsb(nd, w, r, ctxp->dl.data(), ctxp->infi.data(),
                       ctxp->a.data(), ctxp->b.data(), ctxp->cov.data(),
                       tl_y.data(), nd_out, f[0]);
            }
        };

        // Effective integration dimension = nd (after mvspcl's nd-- for normal).
        // Fortran: CALL MVKBRV(ND, ...) after MVSPCL's ND=ND-1.
        // For normal: mvspcl already decremented nd by 1, so no further -1 here.
        // For t: mvspcl leaves nd unchanged; last W element is the chi-squared variate.
        int eff_ndim = nd;

        if (eff_ndim > 0) {
            double e;
            double v;
            mvkbrv(eff_ndim, ivls, maxpts, 1,
                   integrand,
                   abseps, releps, &e, &v, *inform);
            *error = e;
            *value = v;
        } else {
            *value = vl;
            *error = er;
        }
    } else {
        *value = vl;
        *error = er;
    }
}

} // namespace mvtdst

#endif // MVTDST_HPP_
