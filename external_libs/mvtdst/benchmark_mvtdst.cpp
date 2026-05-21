/*
 * Benchmark comparing C++ port vs original Fortran MVTDST.
 *
 * Build:
 *   # Build Fortran version
 *   gfortran -O2 -c mvt.f -o mvt_fortran.o
 *   gcc -O2 -c randomF77.c -o randomF77.o
 *   # Build this benchmark
 *   g++ -std=c++17 -O2 -o benchmark_mvtdst benchmark_mvtdst.cpp mvtnorm.cpp \
 *       mvt_fortran.o randomF77.o -lgfortran -lm
 *
 * Run:
 *   ./benchmark_mvtdst
 */

#include <iostream>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "mvtdst.hpp"
#include "mvtnorm.h"

// Fortran MVTDST extern
extern "C" {
    void mvtdst_(int* n, int* nu, double* lower, double* upper,
                 int* infin, double* correl, double* delta,
                 int* maxpts, double* abseps, double* releps,
                 double* error, double* value, int* inform);
}

// ─── Benchmark infrastructure ───────────────────────────────────────────────
struct BenchResult {
    double ms_total;
    double ms_per_call;
    double value;
    int iters;
};

BenchResult bench(int n_iters, std::function<double()> fn) {
    // Warmup
    double val = 0.0;
    for (int i = 0; i < 5; ++i) val = fn();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iters; ++i) val = fn();
    auto end = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    return {ms, ms / n_iters, val, n_iters};
}

void compare(const std::string& name, int iters,
             std::function<double()> cpp_fn,
             std::function<double()> fortran_fn) {
    auto cpp_r = bench(iters, cpp_fn);
    auto f_r = bench(iters, fortran_fn);

    double diff = std::abs(cpp_r.value - f_r.value);
    double speedup = f_r.ms_per_call / cpp_r.ms_per_call;

    std::cout << "─────────────────────────────────────────────" << std::endl;
    std::cout << name << " (" << iters << " iterations)" << std::endl;
    std::cout << "  C++:     " << std::fixed << std::setprecision(4) << cpp_r.ms_per_call << " ms/call"
              << "  (value=" << std::setprecision(10) << cpp_r.value << ")" << std::endl;
    std::cout << "  Fortran: " << std::fixed << std::setprecision(4) << f_r.ms_per_call << " ms/call"
              << "  (value=" << std::setprecision(10) << f_r.value << ")" << std::endl;
    std::cout << "  Speedup: " << std::setprecision(2) << speedup << "x"
              << "  |diff|=" << std::scientific << diff << std::endl;
}

int main() {
    std::cout << "╔══════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  MVTDST Benchmark: C++ port vs Fortran original     ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════════════╝" << std::endl;

    // ─── 2-D normal, rho=0.5 ────────────────────────────────────────────
    compare("2-D Normal, rho=0.5, P(X>0,Y>0)", 2000,
        []() -> double {
            int n = 2, nu = 0;
            double lower[2] = {0.0, 0.0}, upper[2] = {0.0, 0.0};
            double correl[1] = {0.5}, delta[2] = {0.0, 0.0};
            int infin[2] = {1, 1};
            int maxpts = 25000;
            double abseps = 1e-6, releps = 0.0;
            double error, value;
            int inform;
            pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        },
        []() -> double {
            int n = 2, nu = 0;
            double lower[2] = {0.0, 0.0}, upper[2] = {0.0, 0.0};
            double correl[1] = {0.5}, delta[2] = {0.0, 0.0};
            int infin[2] = {1, 1};
            int maxpts = 25000;
            double abseps = 1e-6, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    // ─── 3-D normal, equicorr 0.3 ──────────────────────────────────────
    compare("3-D Normal, equicorr 0.3, P(X>0)", 1000,
        []() -> double {
            int n = 3, nu = 0;
            double lower[3] = {0, 0, 0}, upper[3] = {0, 0, 0};
            double correl[3] = {0.3, 0.3, 0.3}, delta[3] = {0, 0, 0};
            int infin[3] = {1, 1, 1};
            int maxpts = 50000;
            double abseps = 1e-5, releps = 0.0;
            double error, value;
            int inform;
            pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        },
        []() -> double {
            int n = 3, nu = 0;
            double lower[3] = {0, 0, 0}, upper[3] = {0, 0, 0};
            double correl[3] = {0.3, 0.3, 0.3}, delta[3] = {0, 0, 0};
            int infin[3] = {1, 1, 1};
            int maxpts = 50000;
            double abseps = 1e-5, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    // ─── 5-D normal, equicorr 0.3 ──────────────────────────────────────
    compare("5-D Normal, equicorr 0.3, P(X>0)", 100,
        []() -> double {
            int n = 5, nu = 0;
            double lower[5] = {0, 0, 0, 0, 0}, upper[5] = {0, 0, 0, 0, 0};
            double correl[10] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};
            double delta[5] = {0, 0, 0, 0, 0};
            int infin[5] = {1, 1, 1, 1, 1};
            int maxpts = 100000;
            double abseps = 1e-4, releps = 0.0;
            double error, value;
            int inform;
            pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        },
        []() -> double {
            int n = 5, nu = 0;
            double lower[5] = {0, 0, 0, 0, 0}, upper[5] = {0, 0, 0, 0, 0};
            double correl[10] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};
            double delta[5] = {0, 0, 0, 0, 0};
            int infin[5] = {1, 1, 1, 1, 1};
            int maxpts = 100000;
            double abseps = 1e-4, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    // ─── 10-D normal, equicorr 0.3 ─────────────────────────────────────
    compare("10-D Normal, equicorr 0.3, P(X>0)", 20,
        []() -> double {
            int n = 10, nu = 0;
            double lower[10], upper[10], delta[10];
            int infin[10];
            for (int i = 0; i < 10; ++i) {
                lower[i] = 0.0; upper[i] = 0.0;
                delta[i] = 0.0; infin[i] = 1;
            }
            double correl[45];
            for (int i = 0; i < 45; ++i) correl[i] = 0.3;
            int maxpts = 200000;
            double abseps = 1e-3, releps = 0.0;
            double error, value;
            int inform;
            pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        },
        []() -> double {
            int n = 10, nu = 0;
            double lower[10], upper[10], delta[10];
            int infin[10];
            for (int i = 0; i < 10; ++i) {
                lower[i] = 0.0; upper[i] = 0.0;
                delta[i] = 0.0; infin[i] = 1;
            }
            double correl[45];
            for (int i = 0; i < 45; ++i) correl[i] = 0.3;
            int maxpts = 200000;
            double abseps = 1e-3, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    // ─── 3-D t with nu=5 ───────────────────────────────────────────────
    compare("3-D Student-t (nu=5), equicorr 0.5, P(X>0)", 100,
        []() -> double {
            int n = 3, nu = 5;
            double lower[3] = {0, 0, 0}, upper[3] = {0, 0, 0};
            double correl[3] = {0.5, 0.5, 0.5}, delta[3] = {0, 0, 0};
            int infin[3] = {1, 1, 1};
            int maxpts = 100000;
            double abseps = 1e-4, releps = 0.0;
            double error, value;
            int inform;
            pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        },
        []() -> double {
            int n = 3, nu = 5;
            double lower[3] = {0, 0, 0}, upper[3] = {0, 0, 0};
            double correl[3] = {0.5, 0.5, 0.5}, delta[3] = {0, 0, 0};
            int infin[3] = {1, 1, 1};
            int maxpts = 100000;
            double abseps = 1e-4, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, lower, upper, infin, correl, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    // ─── pmvnorm_complement (regenie pattern) ───────────────────────────
    compare("pmvnorm_complement 5-D (regenie pattern)", 100,
        []() -> double {
            int n = 5;
            double bound[5] = {1.0, 0.5, -0.5, 1.5, 0.0};
            double cmat[10] = {0.3, 0.2, 0.4, 0.1, 0.3, 0.5, 0.2, 0.1, 0.3, 0.4};
            double error;
            return pmvnorm_complement(n, 50000, 1e-3, bound, cmat, &error);
        },
        []() -> double {
            // Call Fortran directly to match
            int n = 5, nu = 0;
            double bound[5] = {1.0, 0.5, -0.5, 1.5, 0.0};
            double cmat[10] = {0.3, 0.2, 0.4, 0.1, 0.3, 0.5, 0.2, 0.1, 0.3, 0.4};
            int infin[5] = {1, 1, 1, 1, 1};
            double delta[5] = {0, 0, 0, 0, 0};
            int maxpts = 50000;
            double abseps = 1e-3, releps = 0.0;
            double error, value;
            int inform;
            mvtdst_(&n, &nu, bound, bound, infin, cmat, delta,
                    &maxpts, &abseps, &releps, &error, &value, &inform);
            return value;
        }
    );

    std::cout << "\n═════════════════════════════════════════════════════" << std::endl;
    std::cout << "Speedup > 1 means C++ is faster, < 1 means Fortran is faster." << std::endl;
    return 0;
}
