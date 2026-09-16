/* qfc_mt.cpp — thread-safe reimplementation of Davies' exact method.
 *
 * All 16 global statics from qfc.cpp are replaced by a per-call
 * qf_state_t struct allocated on the stack of qf_mt().  Every internal
 * helper takes a qf_state_t* so multiple threads can run concurrently.
 *
 * Logic is identical to the original qf() in qfc.cpp:
 *   R.B. Davies, Algorithm AS 155, Appl. Stat. 29 (1980) 323-333.
 *
 * Global-state removal / thread-safety: regenie project.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qfc_mt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Per-call mutable state ─────────────────────────────────────────────── */

typedef struct {
  double sigsq, lmax, lmin, mean, c;
  double intl, ersm;
  int    count, r, lim, env;
  BOOL   ndtsrt, fail;
  int   *n,  *th;
  double *lb, *nc;
} qf_state_t;

/* ── Pure math helpers (read-only, thread-safe as-is) ──────────────────── */

static double exp1(double x)
{ return x < -50.0 ? 0.0 : exp(x); }

static double square(double x) { return x * x; }

static double cube(double x)   { return x * x * x; }

static double log1(double x, BOOL first)
/* if (first) log(1+x)  else  log(1+x) - x */
{
  if (fabs(x) > 0.1)
    return first ? log(1.0 + x) : (log(1.0 + x) - x);

  double s, s1, term, y, k;
  y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
  s = (first ? 2.0 : -x) * y;
  y = square(y);
  for (s1 = s + term / k; s1 != s; s1 = s + term / k)
    { k += 2.0; term *= y; s = s1; }
  return s;
}

/* ── Stateful helpers — each takes qf_state_t* ─────────────────────────── */

static void counter_mt(qf_state_t *st)
{
  st->count++;
  if (st->count > st->lim) st->env = 1;
}

static void order_mt(qf_state_t *st)
/* sort eigenvalues by absolute value into th[] */
{
  int j, k;  double lj;
  for (j = 0; j < st->r; j++) {
    lj = fabs(st->lb[j]);
    for (k = j - 1; k >= 0; k--) {
      if (lj > fabs(st->lb[st->th[k]])) st->th[k + 1] = st->th[k];
      else goto l1_order;
    }
    k = -1;
l1_order:
    st->th[k + 1] = j;
  }
  st->ndtsrt = FALSE;
}

static double errbd_mt(qf_state_t *st, double u, double *cx)
/* bound on tail probability via mgf; cutoff returned in *cx */
{
  double sum1, lj, ncj, x, y, xconst;  int j, nj;
  counter_mt(st);
  xconst = u * st->sigsq;  sum1 = u * xconst;  u = 2.0 * u;
  for (j = st->r - 1; j >= 0; j--) {
    nj = st->n[j];  lj = st->lb[j];  ncj = st->nc[j];
    x = u * lj;  y = 1.0 - x;
    xconst += lj * (ncj / y + nj) / y;
    sum1   += ncj * square(x / y) + nj * (square(x) / y + log1(-x, FALSE));
  }
  *cx = xconst;
  return exp1(-0.5 * sum1);
}

static double ctff_mt(qf_state_t *st, double accx, double *upn)
/* find ctff so that P(qf > ctff) < accx  (upn>0) or P(qf<ctff)<accx */
{
  double u1, u2, u, rb, xconst, c1, c2;
  u2 = *upn;  u1 = 0.0;  c1 = st->mean;
  rb = 2.0 * ((u2 > 0.0) ? st->lmax : st->lmin);
  for (u = u2 / (1.0 + u2 * rb); errbd_mt(st, u, &c2) > accx;
       u = u2 / (1.0 + u2 * rb))
    { u1 = u2;  c1 = c2;  u2 = 2.0 * u2; }
  for (u = (c1 - st->mean) / (c2 - st->mean); u < 0.9;
       u = (c1 - st->mean) / (c2 - st->mean)) {
    u = (u1 + u2) / 2.0;
    if (errbd_mt(st, u / (1.0 + u * rb), &xconst) > accx)
      { u1 = u;  c1 = xconst; }
    else
      { u2 = u;  c2 = xconst; }
  }
  *upn = u2;
  return c2;
}

static double truncation_mt(qf_state_t *st, double u, double tausq)
/* bound integration error due to truncation at u */
{
  double sum1, sum2, prod1, prod2, prod3, lj, ncj, x, y, err1, err2;
  int j, nj, s;
  counter_mt(st);
  sum1  = 0.0;  prod2 = 0.0;  prod3 = 0.0;  s = 0;
  sum2  = (st->sigsq + tausq) * square(u);
  prod1 = 2.0 * sum2;
  u = 2.0 * u;
  for (j = 0; j < st->r; j++) {
    lj = st->lb[j];  ncj = st->nc[j];  nj = st->n[j];
    x = square(u * lj);
    sum1 += ncj * x / (1.0 + x);
    if (x > 1.0) {
      prod2 += nj * log(x);
      prod3 += nj * log1(x, TRUE);
      s += nj;
    } else {
      prod1 += nj * log1(x, TRUE);
    }
  }
  sum1  = 0.5 * sum1;
  prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
  x = exp1(-sum1 - 0.25 * prod2) / pi;
  y = exp1(-sum1 - 0.25 * prod3) / pi;
  err1 = (s == 0)      ? 1.0 : x * 2.0 / s;
  err2 = (prod3 > 1.0) ? 2.5 * y : 1.0;
  if (err2 < err1) err1 = err2;
  x = 0.5 * sum2;
  err2 = (x <= y) ? 1.0 : y / x;
  return (err1 < err2) ? err1 : err2;
}

static void findu_mt(qf_state_t *st, double *utx, double accx)
/* find u s.t. truncation(u)<accx and truncation(u/1.2)>accx */
{
  /* divis[] is a read-only constant — static is fine even under threads */
  static const double divis[] = {2.0, 1.4, 1.2, 1.1};
  double u, ut;  int i;
  ut = *utx;  u = ut / 4.0;
  if (truncation_mt(st, u, 0.0) > accx) {
    for (u = ut; truncation_mt(st, u, 0.0) > accx; u = ut) ut *= 4.0;
  } else {
    ut = u;
    for (u /= 4.0; truncation_mt(st, u, 0.0) <= accx; u /= 4.0) ut = u;
  }
  for (i = 0; i < 4; i++) {
    u = ut / divis[i];
    if (truncation_mt(st, u, 0.0) <= accx) ut = u;
  }
  *utx = ut;
}

static void integrate_mt(qf_state_t *st, int nterm, double interv,
                         double tausq, BOOL mainx)
/* integration with nterm terms at stepsize interv */
{
  double inpi, u, sum1, sum2, sum3, x, y, z;
  int k, j, nj;
  inpi = interv / pi;
  for (k = nterm; k >= 0; k--) {
    u = (k + 0.5) * interv;
    sum1 = -2.0 * u * st->c;  sum2 = fabs(sum1);
    sum3 = -0.5 * st->sigsq * square(u);
    for (j = st->r - 1; j >= 0; j--) {
      nj = st->n[j];  x = 2.0 * st->lb[j] * u;  y = square(x);
      sum3 -= 0.25 * nj * log1(y, TRUE);
      y = st->nc[j] * x / (1.0 + y);
      z = nj * atan(x) + y;
      sum1 += z;  sum2 += fabs(z);
      sum3 -= 0.5 * x * y;
    }
    x = inpi * exp1(sum3) / u;
    if (!mainx) x *= (1.0 - exp1(-0.5 * tausq * square(u)));
    st->intl += sin(0.5 * sum1) * x;
    st->ersm += 0.5 * sum2 * x;
  }
}

static double cfe_mt(qf_state_t *st, double x)
/* coef of tausq in error when convergence factor exp1(-0.5*tausq*u^2)
   is used when df evaluated at x */
{
  double axl, axl1, axl2, sxl, sum1, lj;  int j, k, t;
  counter_mt(st);
  if (st->ndtsrt) order_mt(st);
  axl = fabs(x);  sxl = (x > 0.0) ? 1.0 : -1.0;  sum1 = 0.0;
  for (j = st->r - 1; j >= 0; j--) {
    t = st->th[j];
    if (st->lb[t] * sxl > 0.0) {
      lj   = fabs(st->lb[t]);
      axl1 = axl - lj * (st->n[t] + st->nc[t]);
      axl2 = lj / log28;
      if (axl1 > axl2) {
        axl = axl1;
      } else {
        if (axl > axl2) axl = axl2;
        sum1 = (axl - axl1) / lj;
        for (k = j - 1; k >= 0; k--)
          sum1 += (st->n[st->th[k]] + st->nc[st->th[k]]);
        goto l_cfe;
      }
    }
  }
l_cfe:
  if (sum1 > 100.0) { st->fail = TRUE; return 1.0; }
  return pow(2.0, sum1 / 4.0) / (pi * square(axl));
}

/* ── Public thread-safe entry point ────────────────────────────────────── */

double qf_mt(double* lb1, double* nc1, int* n1, int r1,
             double sigma, double c1,
             int lim1, double acc,
             double* trace, int* ifault)
/*
 * Computes the distribution function of a linear combination of non-central
 * chi-squared random variables — thread-safe variant.
 *
 * Parameters: identical to qf() in qfc.cpp.
 */
{
  qf_state_t st;
  int j, nj, nt, ntm;
  double acc1, almx, xlim, xnt, xntm;
  double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;
  double qfval;
  static const int rats[] = {1, 2, 4, 8};

  /* initialise state from arguments */
  st.r   = r1;   st.lim = lim1;  st.c  = c1;
  st.n   = n1;   st.lb  = lb1;   st.nc = nc1;
  st.env = 0;    st.count = 0;
  st.intl = 0.0; st.ersm  = 0.0;
  st.ndtsrt = TRUE;  st.fail = FALSE;
  qfval = -1.0;  acc1 = acc;
  xlim  = (double)lim1;

  for (j = 0; j < 7; j++) trace[j] = 0.0;
  *ifault = 0;

  st.th = (int*)malloc(r1 * sizeof(int));
  if (!st.th) { *ifault = 5; return qfval; }

  /* find mean, sd, max and min of lb; validate parameters */
  st.sigsq = square(sigma);  sd = st.sigsq;
  st.lmax  = 0.0;  st.lmin = 0.0;  st.mean = 0.0;
  for (j = 0; j < st.r; j++) {
    nj = n1[j];  lj = lb1[j];  ncj = nc1[j];
    if (nj < 0 || ncj < 0.0) { *ifault = 3; goto endofproc; }
    sd       += square(lj) * (2 * nj + 4.0 * ncj);
    st.mean  += lj * (nj + ncj);
    if      (st.lmax < lj) st.lmax = lj;
    else if (st.lmin > lj) st.lmin = lj;
  }
  if (sd == 0.0)
    { qfval = (c1 > 0.0) ? 1.0 : 0.0; goto endofproc; }
  if (st.lmin == 0.0 && st.lmax == 0.0 && sigma == 0.0)
    { *ifault = 3; goto endofproc; }
  sd   = sqrt(sd);
  almx = (st.lmax < -st.lmin) ? -st.lmin : st.lmax;

  /* starting values for findu, ctff */
  utx = 16.0 / sd;  up = 4.5 / sd;  un = -up;

  /* truncation point with no convergence factor */
  findu_mt(&st, &utx, 0.5 * acc1);
  if (st.env != 0) { *ifault = 4; goto endofproc; }

  /* does convergence factor help ? */
  if (st.c != 0.0 && (almx > 0.07 * sd)) {
    tausq = 0.25 * acc1 / cfe_mt(&st, st.c);
    if (st.fail) st.fail = FALSE;
    else if (truncation_mt(&st, utx, tausq) < 0.2 * acc1) {
      st.sigsq += tausq;
      findu_mt(&st, &utx, 0.25 * acc1);
      trace[5] = sqrt(tausq);
    }
  }
  if (st.env != 0) { *ifault = 4; goto endofproc; }
  trace[4] = utx;  acc1 = 0.5 * acc1;

  /* find RANGE of distribution, quit if outside this */
l1:
  d1 = ctff_mt(&st, acc1, &up) - st.c;
  if (d1 < 0.0) { qfval = 1.0; goto endofproc; }
  if (st.env != 0) { *ifault = 4; goto endofproc; }
  d2 = st.c - ctff_mt(&st, acc1, &un);
  if (d2 < 0.0) { qfval = 0.0; goto endofproc; }
  if (st.env != 0) { *ifault = 4; goto endofproc; }

  /* find integration interval */
  intv  = 2.0 * pi / ((d1 > d2) ? d1 : d2);
  xnt   = utx / intv;  xntm = 3.0 / sqrt(acc1);
  if (xnt > xntm * 1.5) {
    /* parameters for auxiliary integration */
    if (xntm > xlim) { *ifault = 1; goto endofproc; }
    ntm   = (int)floor(xntm + 0.5);
    intv1 = utx / ntm;  x = 2.0 * pi / intv1;
    if (x <= fabs(st.c)) goto l2;
    /* convergence factor */
    tausq = 0.33 * acc1 / (1.1 * (cfe_mt(&st, st.c - x) + cfe_mt(&st, st.c + x)));
    if (st.env != 0) { *ifault = 4; goto endofproc; }
    if (st.fail) goto l2;
    acc1 = 0.67 * acc1;
    /* auxiliary integration */
    integrate_mt(&st, ntm, intv1, tausq, FALSE);
    if (st.env != 0) { *ifault = 4; goto endofproc; }
    xlim -= xntm;  st.sigsq += tausq;
    trace[2]++;  trace[1] += ntm + 1;
    findu_mt(&st, &utx, 0.25 * acc1);  acc1 = 0.75 * acc1;
    if (st.env != 0) { *ifault = 4; goto endofproc; }
    goto l1;
  }

  /* main integration */
l2:
  trace[3] = intv;
  if (xnt > xlim) { *ifault = 1; goto endofproc; }
  nt = (int)floor(xnt + 0.5);
  integrate_mt(&st, nt, intv, 0.0, TRUE);
  if (st.env != 0) { *ifault = 4; goto endofproc; }
  trace[2]++;  trace[1] += nt + 1;
  qfval    = 0.5 - st.intl;
  trace[0] = st.ersm;

  /* test for round-off error (allow for radix 8 or 16 machines) */
  up = st.ersm;  x = up + acc / 10.0;
  for (j = 0; j < 4; j++)
    if (rats[j] * x == rats[j] * up) *ifault = 2;

endofproc:
  free(st.th);
  trace[6] = (double)st.count;
  return qfval;
}

#ifdef __cplusplus
}
#endif
