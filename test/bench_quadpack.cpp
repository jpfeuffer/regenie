/*
 * Benchmark: C (CQUADPACK) vs Fortran (QUADPACK) dqags implementations.
 *
 * Build via the CMake `bench_quadpack` target (WITH_C_QUADPACK=ON required).
 * Run from the build directory:
 *   ./bin/bench_quadpack
 *
 * Two integrands are timed:
 *  - "cheap"  : simple Gaussian exp(-x^2), dominated by G_K21 arithmetic.
 *  - "costly" : chi-squared CDF ratio, mimics the real SKATO_integral_fn
 *               cost (calls boost::math::cdf which is O(10-100) ns each).
 */

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/erf.hpp>

extern "C" {
#include "cquadpak.h"   // types only; actual call via adapter
}

/* C-linkage adapter (bench_qp_adapter.c) — avoids K&R vs C++ type mismatch. */
extern "C" double bqp_dqags_c(double fn(double), double a, double b,
                              double epsabs, double epsrel,
                              double *abserr, int *neval, int *ier);

extern "C" {
  // Fortran dqags_ (from libquad.a)
  void dqags_(double f(double*), double*, double*, double*, double*,
              double*, double*, int*, int*, int*, int*, int*, int*, double*);
}

// ── Integrands ────────────────────────────────────────────────────────────

/* Cheap integrand: Gaussian — dominated by G_K21 arithmetic overhead. */
static double cheap_fn(double x) {
    return std::exp(-x * x);
}
/* Same integrand for Fortran's pointer-to-double calling convention. */
static double cheap_fn_f(double* x) { return cheap_fn(*x); }

/* Costly integrand: chi-squared CDF ratio — mimics SKATO integrand cost.
 * Each call makes two boost::math::cdf evaluations. */
static double costly_fn(double x) {
    if (x <= 0.0) return 0.0;
    using boost::math::chi_squared;
    using boost::math::cdf;
    using boost::math::complement;
    static const chi_squared dist_a(3.0);
    static const chi_squared dist_b(7.0);
    return cdf(complement(dist_a, x)) * cdf(dist_b, x * 0.4);
}
static double costly_fn_f(double* x) { return costly_fn(*x); }

// ── Timer helper ─────────────────────────────────────────────────────────

using Clock = std::chrono::high_resolution_clock;
static double elapsed_ms(Clock::time_point t0, Clock::time_point t1) {
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// ── Run one backend repeatedly ────────────────────────────────────────────

struct BenchResult { double ms; double result; };

static BenchResult run_c(double fn_c(double), double a, double b, int N) {
    double epsabs = 1e-12, epsrel = 1e-8;
    double result = 0, abserr;
    int neval, ier;
    auto t0 = Clock::now();
    for (int i = 0; i < N; ++i)
        result = bqp_dqags_c(fn_c, a, b, epsabs, epsrel, &abserr, &neval, &ier);
    auto t1 = Clock::now();
    return {elapsed_ms(t0, t1), result};
}

static BenchResult run_fortran(double fn_f(double*), double a, double b, int N) {
    double epsabs = 1e-12, epsrel = 1e-8;
    double lo = a, hi = b;
    double result = 0, abserr;
    int neval, ier, ilimit = 500, lenw = 2000, last;
    std::vector<int>    iwork(ilimit, 0);
    std::vector<double> work(lenw, 0.0);
    auto t0 = Clock::now();
    for (int i = 0; i < N; ++i)
        dqags_(fn_f, &lo, &hi, &epsabs, &epsrel, &result, &abserr,
               &neval, &ier, &ilimit, &lenw, &last, iwork.data(), work.data());
    auto t1 = Clock::now();
    return {elapsed_ms(t0, t1), result};
}

// ── Main ──────────────────────────────────────────────────────────────────

int main() {
    printf("%-12s  %8s  %12s  %8s  %12s  %8s  %8s\n",
           "integrand", "N", "C_result", "C_ns/call",
           "F_result", "F_ns/call", "speedup");
    printf("%s\n", std::string(85, '-').c_str());

    struct Case {
        const char*      name;
        double           a, b;
        double           (*fn_c)(double);
        double           (*fn_f)(double*);
        int              N;
    } cases[] = {
        { "cheap",   0.0, 5.0,  cheap_fn,  cheap_fn_f,  200000 },
        { "costly",  0.0, 30.0, costly_fn, costly_fn_f,  20000 },
    };

    for (auto& c : cases) {
        /* Warm-up */
        run_c(c.fn_c, c.a, c.b, 200);
        run_fortran(c.fn_f, c.a, c.b, 200);

        auto rc = run_c(c.fn_c, c.a, c.b, c.N);
        auto rf = run_fortran(c.fn_f, c.a, c.b, c.N);

        double ns_c = rc.ms * 1e6 / c.N;
        double ns_f = rf.ms * 1e6 / c.N;

        printf("%-12s  %8d  %12.8f  %8.1f  %12.8f  %8.1f  %7.2fx\n",
               c.name, c.N, rc.result, ns_c, rf.result, ns_f, ns_f / ns_c);
    }

    /* Accuracy cross-check */
    {
        double lo = 0, hi = 5, epsabs = 1e-12, epsrel = 1e-8;
        double abserr;
        int neval, ier;
        double r_c = bqp_dqags_c(cheap_fn, lo, hi, epsabs, epsrel, &abserr, &neval, &ier);

        double lo2 = 0, hi2 = 5, ea = 1e-12, er = 1e-8;
        double r_f, aerr;
        int ne, ie, ilimit = 500, lenw = 2000, last;
        std::vector<int>    iw(ilimit, 0);
        std::vector<double> w(lenw, 0.0);
        dqags_(cheap_fn_f, &lo2, &hi2, &ea, &er, &r_f, &aerr,
               &ne, &ie, &ilimit, &lenw, &last, iw.data(), w.data());

        printf("\nAccuracy check (cheap integrand, a=0, b=5):\n");
        printf("  C       result = %.15g\n", r_c);
        printf("  Fortran result = %.15g\n", r_f);
        printf("  |diff|         = %.3e\n", std::fabs(r_c - r_f));
        /* Reference: integral of exp(-x^2) from 0 to 5 ≈ sqrt(pi)/2 */
        double ref = std::sqrt(M_PI) / 2.0;
        printf("  reference      = %.15g  (sqrt(pi)/2)\n", ref);
        printf("  |C - ref|      = %.3e\n", std::fabs(r_c - ref));
        printf("  |F - ref|      = %.3e\n", std::fabs(r_f - ref));
    }
    return 0;
}
