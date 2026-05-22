/*
 * C-linkage shim for bench_quadpack.cpp.
 *
 * In C, passing double(*)(double) to a K&R-style double f() parameter is
 * permitted via implicit conversion.  In C++ it is a type error.  This file
 * is compiled as C so that bench_quadpack.cpp can call these wrappers without
 * a reinterpret_cast.
 */

#include "cquadpak.h"

typedef double bench_fn_t(double);

double bqp_dqags_c(bench_fn_t *fn, double a, double b,
                   double epsabs, double epsrel,
                   double *abserr, int *neval, int *ier)
{
    return dqags(fn, a, b, epsabs, epsrel, abserr, neval, ier);
}
