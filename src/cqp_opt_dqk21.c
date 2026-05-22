/*
 * Optimized replacement for cqp/src/dqk21.c
 *
 * Changes vs the original CQUADPACK dqk21.c:
 *  1. Tables changed from `static long double` to `static const double`.
 *     On x86-64 `long double` forces 80-bit x87 arithmetic which is scalar-
 *     only (no SSE/AVX).  Using `double` keeps everything in XMM/YMM regs
 *     and lets the compiler auto-vectorize the weight-accumulation loops.
 *  2. Function evaluations are separated into a Phase-1 pre-computation step
 *     before the weight accumulation (Phase 2).  The Phase-2 loops are then
 *     plain linear reductions over flat arrays, which the compiler can
 *     vectorize and fuse into FMA instructions with -O3 -march=native.
 *  3. pow(r, 1.5) replaced by r * sqrt(r) — avoids libm overhead.
 *  4. Intermediate resabs/resasc accumulation rewritten as separate loops
 *     over flat fval arrays — same operations, friendlier for SIMD.
 *
 * The quadrature values themselves are identical to the original (no
 * algorithmic change).  Only the implementation layout differs.
 */

#include <float.h>
#include <math.h>
#include "cquadpak.h"

/* Gauss-Kronrod 21-point nodes (positive half, index 10 = 0 = centre).
 * Values are identical to the original; type changed from long double. */
static const double XGK21[11] = {
    0.99565716302580808074,
    0.97390652851717172008,
    0.93015749135570822600,
    0.86506336668898451073,
    0.78081772658641689706,
    0.67940956829902440623,
    0.56275713466860468334,
    0.43339539412924719080,
    0.29439286270146019813,
    0.14887433898163121088,
    0.00000000000000000000
};

/* Kronrod weights for all 11 half-nodes (index 10 = centre). */
static const double WGK21[11] = {
    0.01169463886737187428,
    0.03255816230796472748,
    0.05475589657435199603,
    0.07503967481091995277,
    0.09312545458369760554,
    0.10938715880229764190,
    0.12349197626206585108,
    0.13470921731147332593,
    0.14277593857706008080,
    0.14773910490133849137,
    0.14944555400291690566
};

/* Gauss weights for the 10-point Gauss sub-rule (5 values by symmetry).
 * These correspond to XGK21 at odd indices 1,3,5,7,9. */
static const double WG10[5] = {
    0.06667134430868813759,
    0.14945134915058059315,
    0.21908636251598204400,
    0.26926671930999635509,
    0.29552422471475287017
};

double G_K21(double f(), double a, double b, double *abserr,
             double *resabs, double *resasc)
{
    const double centr  = 0.5 * (a + b);
    const double hlgth  = 0.5 * (b - a);
    const double dhlgth = fabs(hlgth);

    /* ------------------------------------------------------------------
     * Phase 1 – evaluate f at all 21 Kronrod points.
     *
     * flo[j] = f(centr - hlgth * XGK21[j])   for j = 0..9
     * fhi[j] = f(centr + hlgth * XGK21[j])   for j = 0..9
     * fc      = f(centr)                       (XGK21[10] = 0)
     * ------------------------------------------------------------------ */
    double flo[10], fhi[10];
    for (int j = 0; j < 10; j++) {
        const double absc = hlgth * XGK21[j];
        flo[j] = (*f)(centr - absc);
        fhi[j] = (*f)(centr + absc);
    }
    const double fc = (*f)(centr);

    /* ------------------------------------------------------------------
     * Phase 2 – weight accumulation.  These loops operate only on the
     * already-computed flo/fhi arrays and the constant weight tables,
     * so the compiler can vectorize them with SSE/AVX.
     *
     * fsum[j] = flo[j] + fhi[j]   (symmetric pair sum)
     * ------------------------------------------------------------------ */

    /* Pair sums — kept separate so the loops below are clean reductions. */
    double fsum[10];
    for (int j = 0; j < 10; j++)
        fsum[j] = flo[j] + fhi[j];

    /* Kronrod result: centre + 10 pair contributions. */
    double resk = fc * WGK21[10];
    for (int j = 0; j < 10; j++)
        resk += WGK21[j] * fsum[j];

    /* 10-point Gauss result: 5 pair contributions at odd Kronrod nodes. */
    double resg = 0.0;
    for (int j = 0; j < 5; j++)
        resg += WG10[j] * fsum[2 * j + 1];

    /* resabs = integral of |f|, weighted by Kronrod weights. */
    double resabs_val = WGK21[10] * fabs(fc);
    for (int j = 0; j < 10; j++)
        resabs_val += WGK21[j] * (fabs(flo[j]) + fabs(fhi[j]));

    /* resasc = integral of |f - mean|. */
    const double reskh = resk * 0.5;
    double resasc_val = WGK21[10] * fabs(fc - reskh);
    for (int j = 0; j < 10; j++)
        resasc_val += WGK21[j] * (fabs(flo[j] - reskh) + fabs(fhi[j] - reskh));

    /* ------------------------------------------------------------------
     * Scale and compute error estimate (identical to original logic).
     * ------------------------------------------------------------------ */
    const double result = resk * hlgth;
    *resabs = resabs_val * dhlgth;
    *resasc = resasc_val * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);

    if (*resasc != 0.0 && *abserr != 0.0) {
        /* min(1, (200*abserr/resasc)^1.5) — use r*sqrt(r) instead of pow */
        double r = 200.0 * (*abserr) / (*resasc);
        if (r < 1.0)
            *abserr = (*resasc) * r * sqrt(r);
        else
            *abserr = (*resasc);
    }
    if (*resabs > DBL_MIN / (50.0 * DBL_EPSILON))
        *abserr = fmax(50.0 * DBL_EPSILON * (*resabs), *abserr);

    return result;
}
