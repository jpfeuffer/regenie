/* qfc_mt.h — thread-safe wrapper around Davies' exact method.
 *
 * Declares qf_mt(), a drop-in replacement for qf() with identical semantics
 * but no global state: all per-call mutable data lives on the call stack.
 * Multiple threads may call qf_mt() concurrently without synchronisation.
 *
 * The original algorithm is by R.B. Davies (Applied Statistics, 1980).
 * Thread-safety refactor for the regenie project.
 */
#ifndef QFC_MT_H
#define QFC_MT_H

#ifndef TRUE
#  define TRUE  1
#endif
#ifndef FALSE
#  define FALSE 0
#endif
#ifndef BOOL
   typedef int BOOL;
#endif
#ifndef pi
#  define pi 3.14159265358979
#endif
#ifndef log28
#  define log28 .0866
#endif

#ifdef __cplusplus
extern "C" {
#endif

double qf_mt(double* lb1, double* nc1, int* n1, int r1,
             double sigma, double c1,
             int lim1, double acc,
             double* trace, int* ifault);

#ifdef __cplusplus
}
#endif

#endif /* QFC_MT_H */
