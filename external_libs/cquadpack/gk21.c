/* 21-point Gauss-Kronrod rule, hand-ported to C from the public-domain
 * QUADPACK Fortran77 routine DQK21 (Piessens & de Doncker, 1983), whose
 * source ships in this repo at external_libs/quadpack/dqk21.f. Ported to
 * complete external_libs/cquadpack (dqags.c calls G_K21 but no C
 * implementation of it was otherwise available).
 */

#include <float.h>
#include <math.h>
#include "cquadpak.h"

double G_K21(double f(), double a, double b, double *abserr,
             double *resabs, double *resasc)
{
  /* abscissae and weights of the 21-point kronrod rule; xgk(2),xgk(4),...
   * are abscissae of the 10-point gauss rule, the remaining are those
   * optimally added to the 10-point gauss rule to give the 21-point rule. */
  static const double xgk[11] = {
    0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
    0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
    0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
    0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
    0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
    0.000000000000000000000000000000000
  };
  static const double wgk[11] = {
    0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821
  };
  static const double wg[5] = {
    0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
    0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
    0.295524224714752870173892994651338
  };

  double fv1[10], fv2[10];
  double centr  = 0.5 * (a + b);
  double hlgth  = 0.5 * (b - a);
  double dhlgth = fabs(hlgth);

  double fc   = f(centr);
  double resg = 0.0;
  double resk = wgk[10] * fc;
  *resabs = fabs(resk);

  int j;
  for (j = 0; j < 5; j++) {
    int jtw = 2 * j + 1; /* xgk(2),xgk(4),...,xgk(10) -> 0-based 1,3,5,7,9 */
    double absc  = hlgth * xgk[jtw];
    double fval1 = f(centr - absc);
    double fval2 = f(centr + absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    double fsum = fval1 + fval2;
    resg += wg[j] * fsum;
    resk += wgk[jtw] * fsum;
    *resabs += wgk[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 5; j++) {
    int jtwm1 = 2 * j; /* xgk(1),xgk(3),...,xgk(9) -> 0-based 0,2,4,6,8 */
    double absc  = hlgth * xgk[jtwm1];
    double fval1 = f(centr - absc);
    double fval2 = f(centr + absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    double fsum = fval1 + fval2;
    resk += wgk[jtwm1] * fsum;
    *resabs += wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
  }

  double reskh = resk * 0.5;
  *resasc = wgk[10] * fabs(fc - reskh);
  for (j = 0; j < 10; j++)
    *resasc += wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh));

  double result = resk * hlgth;
  *resabs *= dhlgth;
  *resasc *= dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if (*resasc != 0.0 && *abserr != 0.0)
    *abserr = *resasc * min(1.0, pow(200.0 * (*abserr) / (*resasc), 1.5));
  if (*resabs > uflow / (50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs), *abserr);

  return result;
}
