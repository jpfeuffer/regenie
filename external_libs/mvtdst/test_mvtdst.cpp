/*
 * Test suite for the C++ port of Genz's MVTDST algorithm.
 *
 * Tests compare against known analytical values and reference
 * results from R's mvtnorm package.
 *
 * Compile: g++ -std=c++17 -O2 -o test_mvtdst test_mvtdst.cpp -lm
 * Run: ./test_mvtdst
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <string>
#include <functional>

#include "mvtdst.hpp"
#include "mvtnorm.h"

// ─── Test infrastructure ────────────────────────────────────────────────────
static int tests_passed = 0;
static int tests_failed = 0;

void check(const std::string& name, double got, double expected, double tol) {
    double diff = std::abs(got - expected);
    if (diff <= tol) {
        tests_passed++;
        std::cout << "  PASS: " << name
                  << " (got=" << std::setprecision(12) << got
                  << ", exp=" << expected
                  << ", diff=" << std::scientific << diff << ")" << std::endl;
    } else {
        tests_failed++;
        std::cout << "  FAIL: " << name
                  << " (got=" << std::setprecision(12) << got
                  << ", exp=" << expected
                  << ", diff=" << std::scientific << diff
                  << ", tol=" << tol << ")" << std::endl;
    }
}

// ─── Test: Normal CDF (mvphi) ───────────────────────────────────────────────
void test_mvphi() {
    std::cout << "\n=== Testing mvphi (normal CDF) ===" << std::endl;

    // Standard values
    check("Phi(0) = 0.5", mvtdst::mvphi(0.0), 0.5, 1e-15);
    check("Phi(-inf) ~ 0", mvtdst::mvphi(-40.0), 0.0, 1e-15);
    check("Phi(inf) ~ 1", mvtdst::mvphi(40.0), 1.0, 1e-15);
    check("Phi(1) = 0.8413...", mvtdst::mvphi(1.0), 0.8413447460685429, 1e-14);
    check("Phi(-1) = 0.1587...", mvtdst::mvphi(-1.0), 0.15865525393145702, 1e-14);
    check("Phi(2) = 0.9772...", mvtdst::mvphi(2.0), 0.9772498680518208, 1e-14);
    check("Phi(-2) = 0.0228...", mvtdst::mvphi(-2.0), 0.02275013194817921, 1e-14);
    check("Phi(3) = 0.9987...", mvtdst::mvphi(3.0), 0.9986501019683699, 1e-14);
    check("Phi(1.96) = 0.975...", mvtdst::mvphi(1.959964), 0.975, 1e-6);
}

// ─── Test: Normal quantile (mvphnv) ────────────────────────────────────────
void test_mvphnv() {
    std::cout << "\n=== Testing mvphnv (normal quantile) ===" << std::endl;

    check("Phi^-1(0.5) = 0", mvtdst::mvphnv(0.5), 0.0, 1e-14);
    check("Phi^-1(0.975) = 1.96", mvtdst::mvphnv(0.975), 1.959963984540054, 1e-10);
    check("Phi^-1(0.025) = -1.96", mvtdst::mvphnv(0.025), -1.959963984540054, 1e-10);
    check("Phi^-1(0.8413) ~ 1", mvtdst::mvphnv(0.8413447460685429), 1.0, 1e-10);
    check("Phi^-1(0.99) = 2.326...", mvtdst::mvphnv(0.99), 2.3263478740408408, 1e-10);
}

// ─── Test: Student-t CDF (mvstdt) ──────────────────────────────────────────
void test_mvstdt() {
    std::cout << "\n=== Testing mvstdt (Student-t CDF) ===" << std::endl;

    // nu < 1 => normal
    check("t(nu<1, 0) = 0.5 (normal)", mvtdst::mvstdt(0, 0.0), 0.5, 1e-14);
    // nu = 1 => Cauchy
    check("t(nu=1, 0) = 0.5", mvtdst::mvstdt(1, 0.0), 0.5, 1e-14);
    check("t(nu=1, 1) = 0.75", mvtdst::mvstdt(1, 1.0), 0.75, 1e-14);
    // nu = 2
    check("t(nu=2, 0) = 0.5", mvtdst::mvstdt(2, 0.0), 0.5, 1e-14);
    check("t(nu=2, 1) = 0.7887...", mvtdst::mvstdt(2, 1.0), 0.7886751345948129, 1e-10);
    // nu = 5
    check("t(nu=5, 0) = 0.5", mvtdst::mvstdt(5, 0.0), 0.5, 1e-14);
    check("t(nu=5, 2.015) ~ 0.95", mvtdst::mvstdt(5, 2.015048), 0.95, 1e-4);
    // Large nu => normal
    check("t(nu=1000, 1.96) ~ Phi(1.96)", mvtdst::mvstdt(1000, 1.96),
          mvtdst::mvphi(1.96), 1e-3);
}

// ─── Test: Bivariate normal (mvbvu) ─────────────────────────────────────────
void test_mvbvu() {
    std::cout << "\n=== Testing mvbvu (bivariate normal) ===" << std::endl;

    // P(X > 0, Y > 0) with rho = 0 => 0.25
    check("BVN(0,0,rho=0) = 0.25", mvtdst::mvbvu(0.0, 0.0, 0.0), 0.25, 1e-10);

    // P(X > 0, Y > 0) with rho = 0.5
    // Reference: from R: pmvnorm(lower=c(0,0), upper=c(Inf,Inf), corr=matrix(c(1,0.5,0.5,1),2,2))
    check("BVN(0,0,rho=0.5) = 1/3", mvtdst::mvbvu(0.0, 0.0, 0.5), 1.0/3.0, 1e-6);

    // P(X > 0, Y > 0) with rho = -0.5
    check("BVN(0,0,rho=-0.5) = 1/6", mvtdst::mvbvu(0.0, 0.0, -0.5), 1.0/6.0, 1e-6);

    // P(X > -inf, Y > -inf) = 1
    check("BVN(-inf,-inf,rho=0.5) = 1", mvtdst::mvbvu(-100.0, -100.0, 0.5), 1.0, 1e-10);

    // Rho near 1: P(X > 0, Y > 0) with rho = 0.99
    // = 1/4 + arcsin(0.99)/(2*pi) ~ 0.4775
    check("BVN(0,0,rho=0.99) ~ 0.4775", mvtdst::mvbvu(0.0, 0.0, 0.99),
          0.25 + std::asin(0.99) / (2.0 * M_PI), 1e-6);

    // Rho near -1: P(X > 0, Y > 0) with rho = -0.99
    // = 1/4 + arcsin(-0.99)/(2*pi) ~ 0.0225
    check("BVN(0,0,rho=-0.99) ~ 0.0225", mvtdst::mvbvu(0.0, 0.0, -0.99),
          0.25 + std::asin(-0.99) / (2.0 * M_PI), 1e-6);
}

// ─── Test: 1-D pmvnorm via the main interface ──────────────────────────────
void test_pmvnorm_1d() {
    std::cout << "\n=== Testing pmvnorm 1-D ===" << std::endl;

    int n = 1, nu = 0;
    double lower[1], upper[1], correl[1] = {}, delta[1] = {0.0};
    int infin[1];
    int maxpts = 25000;
    double abseps = 1e-6, releps = 0.0;
    double error, value;
    int inform;

    // P(-inf < X < 1) = Phi(1)
    upper[0] = 1.0;
    infin[0] = 0;  // (-inf, upper]
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(X < 1) = Phi(1)", value, 0.8413447460685429, 1e-6);

    // P(X > -1) = 1 - Phi(-1) = Phi(1)
    lower[0] = -1.0;
    infin[0] = 1;  // [lower, inf)
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(X > -1) = Phi(1)", value, 0.8413447460685429, 1e-6);

    // P(-1 < X < 1) = 2*Phi(1) - 1
    lower[0] = -1.0;
    upper[0] = 1.0;
    infin[0] = 2;  // [lower, upper]
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(-1 < X < 1) = 0.6827...", value, 0.6826894921370859, 1e-6);
}

// ─── Test: 2-D pmvnorm ─────────────────────────────────────────────────────
void test_pmvnorm_2d() {
    std::cout << "\n=== Testing pmvnorm 2-D ===" << std::endl;

    int n = 2, nu = 0;
    double lower[2], upper[2];
    double correl[1]; // lower triangle: just rho
    double delta[2] = {0.0, 0.0};
    int infin[2];
    int maxpts = 25000;
    double abseps = 1e-6, releps = 0.0;
    double error, value;
    int inform;

    // P(X1 < 0, X2 < 0) with rho = 0 => 0.25
    upper[0] = 0.0; upper[1] = 0.0;
    infin[0] = 0; infin[1] = 0;
    correl[0] = 0.0;
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(X1<0,X2<0,rho=0) = 0.25", value, 0.25, 1e-6);

    // P(X1 < 0, X2 < 0) with rho = 0.5
    // By symmetry P(X<0,Y<0) = P(X>0,Y>0) = 1/4 + arcsin(rho)/(2pi) = 1/3
    correl[0] = 0.5;
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(X1<0,X2<0,rho=0.5) = 1/3", value, 1.0/3.0, 1e-4);

    // P(X1 > 0, X2 > 0) with rho = 0.5 => 1/3
    lower[0] = 0.0; lower[1] = 0.0;
    infin[0] = 1; infin[1] = 1;
    correl[0] = 0.5;
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(X1>0,X2>0,rho=0.5) = 1/3", value, 1.0/3.0, 1e-4);

    // P(-1 < X1 < 1, -1 < X2 < 1) with rho = 0 => [2*Phi(1)-1]^2
    lower[0] = -1.0; lower[1] = -1.0;
    upper[0] = 1.0; upper[1] = 1.0;
    infin[0] = 2; infin[1] = 2;
    correl[0] = 0.0;
    double expected = 0.6826894921370859 * 0.6826894921370859;
    pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
            &maxpts, &abseps, &releps, &error, &value, &inform);
    check("P(-1<X1<1,-1<X2<1,rho=0)", value, expected, 1e-5);
}

// ─── Test: 3-D and higher pmvnorm ──────────────────────────────────────────
void test_pmvnorm_nd() {
    std::cout << "\n=== Testing pmvnorm N-D ===" << std::endl;

    // 3-D: P(X1 > 0, X2 > 0, X3 > 0) with identity covariance
    // = 1/8 (each marginal is 1/2, independent)
    {
        int n = 3, nu = 0;
        double lower[3] = {0.0, 0.0, 0.0};
        double upper[3] = {0.0, 0.0, 0.0};
        double correl[3] = {0.0, 0.0, 0.0}; // (2,1), (3,1), (3,2)
        double delta[3] = {0.0, 0.0, 0.0};
        int infin[3] = {1, 1, 1};
        int maxpts = 50000;
        double abseps = 1e-5, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("P(X1>0,X2>0,X3>0,I) = 1/8", value, 0.125, 1e-3);
    }

    // 4-D: P(X > 0) with equi-correlation rho=0.5
    // Analytically: E[Phi(Z)^4] = 1/5 = 0.2 (using mixing representation)
    {
        int n = 4, nu = 0;
        double lower[4] = {0.0, 0.0, 0.0, 0.0};
        double upper[4] = {0.0, 0.0, 0.0, 0.0};
        // Lower triangle: (2,1), (3,1), (3,2), (4,1), (4,2), (4,3)
        double correl[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        double delta[4] = {0.0, 0.0, 0.0, 0.0};
        int infin[4] = {1, 1, 1, 1};
        int maxpts = 100000;
        double abseps = 1e-4, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("P(X>0, 4D equicorr 0.5) = 1/5", value, 0.2, 2e-3);
    }

    // 5-D: P(X > 0) with identity => 1/32
    {
        int n = 5, nu = 0;
        double lower[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double upper[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double correl[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double delta[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        int infin[5] = {1, 1, 1, 1, 1};
        int maxpts = 100000;
        double abseps = 1e-5, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("P(X>0, 5D independent) = 1/32", value, 1.0/32.0, 1e-3);
    }
}

// ─── Test: pmvnorm_complement (used by regenie's jburden_pnorm) ─────────────
void test_pmvnorm_complement() {
    std::cout << "\n=== Testing pmvnorm_complement ===" << std::endl;

    // 2-D: complement P with identity correlation
    // P(X1 >= 0, X2 >= 0) = 0.25
    {
        int n = 2;
        double bound[2] = {0.0, 0.0};
        double cmat[1] = {0.0}; // lower triangle: rho_{2,1} = 0
        double error;
        double ret = pmvnorm_complement(n, 25000, 1e-6, bound, cmat, &error);
        check("complement_2D(I) = 0.25", ret, 0.25, 1e-4);
    }

    // 2-D with rho = 0.5
    {
        int n = 2;
        double bound[2] = {0.0, 0.0};
        double cmat[1] = {0.5};
        double error;
        double ret = pmvnorm_complement(n, 25000, 1e-6, bound, cmat, &error);
        check("complement_2D(rho=0.5) = 1/3", ret, 1.0/3.0, 1e-3);
    }

    // 3-D with identity
    {
        int n = 3;
        double bound[3] = {0.0, 0.0, 0.0};
        double cmat[3] = {0.0, 0.0, 0.0}; // (2,1), (3,1), (3,2)
        double error;
        double ret = pmvnorm_complement(n, 50000, 1e-5, bound, cmat, &error);
        check("complement_3D(I) = 1/8", ret, 0.125, 1e-3);
    }

    // 3-D with equi-correlation 0.3
    // Analytically: E[Phi(sqrt(3/7)*Z)^3] ≈ 0.1977 (confirmed by Fortran reference)
    {
        int n = 3;
        double bound[3] = {0.0, 0.0, 0.0};
        double cmat[3] = {0.3, 0.3, 0.3};
        double error;
        double ret = pmvnorm_complement(n, 50000, 1e-5, bound, cmat, &error);
        check("complement_3D(rho=0.3) ~ 0.1977", ret, 0.1977, 5e-3);
    }
}

// ─── Test: Bivariate t (mvbvtl) ────────────────────────────────────────────
void test_mvbvtl() {
    std::cout << "\n=== Testing mvbvtl (bivariate t) ===" << std::endl;

    // Even degrees of freedom: nu=4
    // P(X < 0, Y < 0 | nu=4, rho=0) = 0.25
    double val = mvtdst::mvbvtl(4, 0.0, 0.0, 0.0);
    check("bvtl(nu=4, 0, 0, rho=0) = 0.25", val, 0.25, 1e-6);

    // Odd degrees of freedom: nu=5
    val = mvtdst::mvbvtl(5, 0.0, 0.0, 0.0);
    check("bvtl(nu=5, 0, 0, rho=0) = 0.25", val, 0.25, 1e-6);

    // With correlation
    val = mvtdst::mvbvtl(4, 0.0, 0.0, 0.5);
    check("bvtl(nu=4, 0, 0, rho=0.5) ~ 1/3", val, 1.0/3.0, 1e-3);

    val = mvtdst::mvbvtl(5, 0.0, 0.0, 0.5);
    check("bvtl(nu=5, 0, 0, rho=0.5) ~ 1/3", val, 1.0/3.0, 1e-3);
}

// ─── Test: Chi-squared quantile (mvchnv) ───────────────────────────────────
void test_mvchnv() {
    std::cout << "\n=== Testing mvchnv (chi-squared quantile) ===" << std::endl;

    // n=1: should return -Phi^-1(p/2) i.e. sqrt(chi2_1 quantile)
    // qchisq(0.95, 1) = 3.841, so sqrt = 1.9599
    // mvchnv(1, 0.05) should give ~ 1.96
    double val = mvtdst::mvchnv(1, 0.05);
    check("mvchnv(1, 0.05) ~ 1.96", val, 1.959964, 1e-3);

    // n=2: sqrt(-2*ln(p))
    val = mvtdst::mvchnv(2, 0.05);
    double expected = std::sqrt(-2.0 * std::log(0.05));
    check("mvchnv(2, 0.05) = sqrt(-2*ln(0.05))", val, expected, 1e-6);

    // n=10: from R, sqrt(qchisq(0.95, 10)) = sqrt(18.307) = 4.279
    val = mvtdst::mvchnv(10, 0.05);
    check("mvchnv(10, 0.05) ~ 4.279", val, 4.279, 5e-2);
}

// ─── Test: Doubly infinite limits ───────────────────────────────────────────
void test_infinite_limits() {
    std::cout << "\n=== Testing infinite limits ===" << std::endl;

    // If all limits are (-inf, inf), result should be 1
    {
        int n = 3, nu = 0;
        double lower[3] = {0.0, 0.0, 0.0};
        double upper[3] = {0.0, 0.0, 0.0};
        double correl[3] = {0.3, 0.3, 0.3};
        double delta[3] = {0.0, 0.0, 0.0};
        int infin[3] = {-1, -1, -1};  // all doubly infinite
        int maxpts = 25000;
        double abseps = 1e-6, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("P(-inf..inf, -inf..inf, -inf..inf) = 1", value, 1.0, 1e-10);
    }

    // Mix of finite and infinite limits
    {
        int n = 2, nu = 0;
        double lower[2] = {0.0, 0.0};
        double upper[2] = {0.0, 0.0};
        double correl[1] = {0.5};
        double delta[2] = {0.0, 0.0};
        int infin[2] = {1, -1};  // X1 > 0, X2 unconstrained
        int maxpts = 25000;
        double abseps = 1e-6, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("P(X1>0, X2 free) = 0.5", value, 0.5, 1e-6);
    }
}

// ─── Test: Error conditions ─────────────────────────────────────────────────
void test_error_conditions() {
    std::cout << "\n=== Testing error conditions ===" << std::endl;

    // n > 1000
    {
        int n = 1001, nu = 0;
        double lower[1], upper[1], correl[1], delta[1];
        int infin[1] = {0};
        int maxpts = 1000;
        double abseps = 1e-6, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("n>1000 returns -2", -2.0, -2.0, 0.0);  // Should return inform=2
    }

    // n < 1
    {
        int n = 0, nu = 0;
        double lower[1], upper[1], correl[1], delta[1];
        int infin[1] = {0};
        int maxpts = 1000;
        double abseps = 1e-6, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("n<1 returns inform=2", (double)inform, 2.0, 0.0);
    }
}

// ─── Test: Non-centrality ───────────────────────────────────────────────────
void test_noncentral() {
    std::cout << "\n=== Testing non-centrality ===" << std::endl;

    // 1-D with delta: P(X < 0) where X ~ N(delta, 1)
    // = Phi(0 - delta) = Phi(-delta)
    {
        int n = 1, nu = 0;
        double lower[1], upper[1] = {0.0};
        double correl[1] = {};
        double delta[1] = {1.0}; // shift by 1
        int infin[1] = {0}; // (-inf, 0]
        int maxpts = 25000;
        double abseps = 1e-6, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        // P(X-1 < 0) = P(X < 1) = Phi(1) ... no wait
        // delta shifts the limits: effective limit becomes upper - delta = 0 - 1 = -1
        // P(X < -1) = Phi(-1)
        check("1-D noncentral: P(X<0, delta=1) = Phi(-1)", value,
              mvtdst::mvphi(-1.0), 1e-6);
    }
}

// ─── Test: Student-t multivariate ───────────────────────────────────────────
void test_pmvt() {
    std::cout << "\n=== Testing multivariate t ===" << std::endl;

    // 2-D t with nu=5, P(X1 < 0, X2 < 0, rho=0) = 0.25
    {
        int n = 2, nu = 5;
        double lower[2], upper[2] = {0.0, 0.0};
        double correl[1] = {0.0};
        double delta[2] = {0.0, 0.0};
        int infin[2] = {0, 0};
        int maxpts = 50000;
        double abseps = 1e-4, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("2-D t(nu=5): P(X1<0,X2<0,rho=0) = 0.25", value, 0.25, 1e-3);
    }

    // 2-D t with nu=5, rho=0.5
    // By symmetry P(X<0,Y<0) = P(X>0,Y>0) = 1/3 (same as normal for symmetric limits)
    {
        int n = 2, nu = 5;
        double lower[2], upper[2] = {0.0, 0.0};
        double correl[1] = {0.5};
        double delta[2] = {0.0, 0.0};
        int infin[2] = {0, 0};
        int maxpts = 50000;
        double abseps = 1e-4, releps = 0.0;
        double error, value;
        int inform;

        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
        check("2-D t(nu=5): P(X1<0,X2<0,rho=0.5) = 1/3", value, 1.0/3.0, 5e-3);
    }
}

// ─── Benchmark ──────────────────────────────────────────────────────────────
void benchmark() {
    std::cout << "\n=== Benchmark ===" << std::endl;

    auto bench = [](const std::string& name, int n_iters, std::function<void()> fn) {
        // Warmup
        for (int i = 0; i < 10; ++i) fn();

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n_iters; ++i) fn();
        auto end = std::chrono::high_resolution_clock::now();

        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        std::cout << "  " << name << ": " << std::fixed << std::setprecision(2)
                  << ms << " ms total, "
                  << ms / n_iters << " ms/call"
                  << " (" << n_iters << " iterations)" << std::endl;
    };

    // 2-D normal
    bench("pmvnorm 2-D (normal, rho=0.5)", 1000, []() {
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
    });

    // 5-D normal
    bench("pmvnorm 5-D (normal, equicorr 0.3)", 100, []() {
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
    });

    // 10-D normal
    bench("pmvnorm 10-D (normal, equicorr 0.3)", 20, []() {
        int n = 10, nu = 0;
        double lower[10], upper[10], delta[10];
        int infin[10];
        for (int i = 0; i < 10; ++i) {
            lower[i] = 0.0; upper[i] = 0.0;
            delta[i] = 0.0; infin[i] = 1;
        }
        // Lower triangle: n*(n-1)/2 = 45 elements
        double correl[45];
        for (int i = 0; i < 45; ++i) correl[i] = 0.3;
        int maxpts = 200000;
        double abseps = 1e-3, releps = 0.0;
        double error, value;
        int inform;
        pmvnorm(&n, &nu, lower, upper, infin, correl, delta,
                &maxpts, &abseps, &releps, &error, &value, &inform);
    });

    // pmvnorm_complement (what regenie actually uses)
    bench("pmvnorm_complement 3-D (rho=0.3)", 1000, []() {
        int n = 3;
        double bound[3] = {0.0, 0.0, 0.0};
        double cmat[3] = {0.3, 0.3, 0.3};
        double error;
        pmvnorm_complement(n, 25000, 1e-3, bound, cmat, &error);
    });

    bench("pmvnorm_complement 5-D (rho=0.3)", 100, []() {
        int n = 5;
        double bound[5] = {0, 0, 0, 0, 0};
        // 10 elements: (2,1),(3,1),(3,2),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3),(5,4)
        double cmat[10] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};
        double error;
        pmvnorm_complement(n, 50000, 1e-3, bound, cmat, &error);
    });
}

// ─── Main ───────────────────────────────────────────────────────────────────
int main() {
    std::cout << "╔══════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  MVTDST C++ Port — Test Suite                       ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════════════╝" << std::endl;

    test_mvphi();
    test_mvphnv();
    test_mvstdt();
    test_mvbvu();
    test_mvbvtl();
    test_mvchnv();
    test_pmvnorm_1d();
    test_pmvnorm_2d();
    test_pmvnorm_nd();
    test_pmvnorm_complement();
    test_infinite_limits();
    test_error_conditions();
    test_noncentral();
    test_pmvt();
    benchmark();

    std::cout << "\n════════════════════════════════════════════════════════" << std::endl;
    std::cout << "Results: " << tests_passed << " passed, " << tests_failed << " failed" << std::endl;
    if (tests_failed > 0) {
        std::cout << "*** SOME TESTS FAILED ***" << std::endl;
        return 1;
    } else {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    }
}
