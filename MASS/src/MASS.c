/*
 *  MASS/MASS.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 */

#include <stdio.h>
#include <math.h>
#include <S.h>

#if defined(SPLUS_VERSION) && SPLUS_VERSION >= 4000 && SPLUS_VERSION < 5000
#  include <newredef.h>
#endif

#include "verS.h"

#define min9(a,b) (a < b?a:b)
#define max9(a,b) (a > b?a:b)
#define abs9(a) (a > 0?a:-a)

/* -----------------------------------------------------------------
 *  Former sammon.c
 */

void
VR_sammon(double *dd, Sint *nn, Sint *kd, double *Y, Sint *niter,
	  double *stress, Sint *trace, double *aa, double *tol)
{
  int     i, j, k, m, n = *nn, nd = *kd;
  double *xu, *xv, *e1, *e2;
  double  dpj, dq, dr, dt;
  double  xd, xx;
  double  e, epast, eprev, tot, d, d1, ee, magic = *aa;

  xu = Calloc(nd * n, double);
  xv = Calloc(nd, double);
  e1 = Calloc(nd, double);
  e2 = Calloc(nd, double);

  epast = eprev = 1.0;
  magic = magic;

  /* Error in distances */
  e = tot = 0.0;
  for (j = 1; j < n; j++)
    for (k = 0; k < j; k++) {
      d = dd[k * n + j];
      if (d <= 0.0)
	PROBLEM "%s\n", "some distance is zero or negative"
	  RECOVER(NULL_ENTRY);
      tot += d;
      d1 = 0.0;
      for (m = 0; m < nd; m++) {
	xd = Y[j + m * n] - Y[k + m * n];
	d1 += xd * xd;
      }
      ee = d - sqrt(d1);
      e += (ee * ee / d);
    }
  e /= tot;
  if (*trace) {
    printf("Initial stress        : %7.5f\n", e);
    fflush(stdout);
  }
  epast = eprev = e;
  
  /* Iterate */
  for (i = 1; i <= *niter; i++) {
  CORRECT:
    for (j = 0; j < n; j++) {
      for (m = 0; m < nd; m++)
	e1[m] = e2[m] = 0.0;
      for (k = 0; k < n; k++) {
	if (j == k)
	  continue;
	d1 = 0.0;
	for (m = 0; m < nd; m++) {
	  xd = Y[j + m * n] - Y[k + m * n];
	  d1 += xd * xd;
	  xv[m] = xd;
	}
	dpj = sqrt(d1);

	/* Calculate derivatives */
	dt = dd[k * n + j];
	dq = dt - dpj;
	dr = dt * dpj;
	for (m = 0; m < nd; m++) {
	  e1[m] += xv[m] * dq / dr;
	  e2[m] += (dq - xv[m] * xv[m] * (1.0 + dq / dpj) / dpj) / dr;
	}
      }
      /* Correction */
      for (m = 0; m < nd; m++)
	xu[j + m * n] = Y[j + m * n] + magic * e1[m] / fabs(e2[m]);
    }

    /* Error in distances */
    e = 0.0;
    for (j = 1; j < n; j++)
      for (k = 0; k < j; k++) {
	d = dd[k * n + j];
	d1 = 0.0;
	for (m = 0; m < nd; m++) {
	  xd = xu[j + m * n] - xu[k + m * n];
	  d1 += xd * xd;
	}
	ee = d - sqrt(d1);
	e += (ee * ee / d);
      }
    e /= tot;
    if (e > eprev) {
      e = eprev;
      magic = magic * 0.2;
      if (magic > 1.0e-3) goto CORRECT;
      if (*trace) {
	printf("stress after %3d iters: %7.5f\n", i - 1, e);
	fflush(stdout);
      }
      break;
    }
    magic *= 1.5;
    if (magic > 0.5) magic = 0.5;
    eprev = e;

    /* Move the centroid to origin and update */
    for (m = 0; m < nd; m++) {
      xx = 0.0;
      for (j = 0; j < n; j++)
	xx += xu[j + m * n];
      xx /= n;
      for (j = 0; j < n; j++)
	Y[j + m * n] = xu[j + m * n] - xx;
    }

    if (i % 10 == 0) {
      if (*trace) {
	printf("stress after %3d iters: %7.5f, magic = %5.3f\n", i, e, magic);
	fflush(stdout);
      }
      if (e > epast - *tol)
	break;
      epast = e;
    }
  }
  *stress = e;
  Free(xu); Free(xv); Free(e1); Free(e2); 
}

/*
 * ----------------------------------------------------------
 *  Former isoMDS.c

    C code for mds S-Plus library, which implements Kruskal's MDS.
    (c) B.D. Ripley, May 1995.
 *
 */

static Sint *ord;		/* ranks of dissimilarities */
static Sint *ord2;		/* inverse ordering (which one is
				   rank i?) */
static Sint n;			/* number of  dissimilarities */
static Sint nr;			/* number of data points */
static Sint nc;			/* # cols of  fitted configuration */
static int  dimx;		/* Size of configuration array */
static double *x;		/* configuration */
static double *d;		/* dissimilarities */
static double *y;		/* fitted distances (in rank of d
				   order) */
static double *yf;		/* isotonic regression fitted values
				   (ditto) */

void    VR_mds_fn(double *, double *, Sint *, double *, Sint *,
		  double *, Sint *, Sint *, double *, Sint *);
static void vmmin(int, double *, double *, int, int);
static void errmsg(char *);

/*
 *  Download the data.
 */
void
VR_mds_init_data(Sint *pn, Sint *pc, Sint *pr, Sint *orde,
		 Sint *ordee, double *xx)
{
  int     i;
  n = *pn;
  nr = *pr;
  nc = *pc;
  dimx = nr * nc;
  ord = Calloc(n, Sint);
  ord2 = Calloc(n, Sint);
  x = Calloc(dimx, double);
  d = Calloc(n, double);
  y = Calloc(n, double);
  yf = Calloc(n, double);
  if (!ord | !ord2 | !x | !d | !y | !yf) errmsg("Not enough space");
  for (i = 0; i < n; i++) ord[i] = orde[i];
  for (i = 0; i < n; i++) ord2[i] = ordee[i];
  for (i = 0; i < dimx; i++) x[i] = xx[i];
}

void
VR_mds_unload()
{
  Free(ord); Free(ord2); Free(x); Free(d); Free(y); Free(yf);
}


static void
calc_dist(double *x)
{
  int     r1, r2, c, index;
  double  tmp, tmp1;
  index = 0;
  for (r1 = 0; r1 < nr; r1++)
    for (r2 = r1 + 1; r2 < nr; r2++) {
      tmp = 0.0;
      for (c = 0; c < nc; c++) {
	tmp1 = x[r1 + c * nr] - x[r2 + c * nr];
	tmp += tmp1 * tmp1;
      }
      d[index++] = sqrt(tmp);
    }
  for (index = 0; index < n; index++) y[index] = d[ord[index]];
}

static double fminfn(double *x)
{
  double  ssq;
  Sint    do_derivatives = 0;
  calc_dist(x);
  VR_mds_fn(y, yf, &n, &ssq, ord2, x, &nr, &nc, 0, &do_derivatives);
  return (ssq);
}

static void fmingr(double *x, double *der)
{
  double  ssq;
  Sint    do_derivatives = 1;
  calc_dist(x);
  VR_mds_fn(y, yf, &n, &ssq, ord2, x, &nr, &nc, der, &do_derivatives);
}

void
VR_mds_dovm(double *val, Sint *maxit, Sint *trace, double *xx)
{
  int     i;
  vmmin(dimx, x, val, (int) *maxit, (int) *trace);
  for (i = 0; i < dimx; i++) xx[i] = x[i];
}

/*
Does isotonic regression.
*/

void
VR_mds_fn(double *y, double *yf, Sint *pn, double *pssq, Sint *pd,
	  double *x, Sint *pr, Sint *pncol, double *der, 
	  Sint *do_derivatives)
{
  int     n = *pn, i, ip, known, u, s, r = *pr, ncol = *pncol, k;
  double  tmp, ssq, *yc, slope, tstar, sstar;

  yc = Calloc((n + 1), double);
  yc[0] = 0.0;
  tmp = 0.0;
  for (i = 0; i < n; i++) {
    tmp += y[i];
    yc[i + 1] = tmp;
  }
  known = 0;
  do {
    slope = 1.0e+200;
    for (i = known + 1; i <= n; i++) {
      tmp = (yc[i] - yc[known]) / (i - known);
      if (tmp < slope) {
	slope = tmp;
	ip = i;
      }
    }
    for (i = known; i < ip; i++)
      yf[i] = (yc[ip] - yc[known]) / (ip - known);
  } while ((known = ip) < n);

  sstar = 0.0;
  tstar = 0.0;
  for (i = 0; i < n; i++) {
    tmp = y[i] - yf[i];
    sstar += tmp * tmp;
    tstar += y[i] * y[i];
  }
  ssq = 100 * sqrt(sstar / tstar);
  *pssq = ssq;
  Free(yc);
  if (!(*do_derivatives)) return;
  /* get derivatives */
  for (u = 0; u < r; u++) {
    for (i = 0; i < ncol; i++) {
      tmp = 0.0;
      for (s = 0; s < r; s++) {
	if (s > u) k = r * u - u * (u + 1) / 2 + s - u;
	else if (s < u) k = r * s - s * (s + 1) / 2 + u - s;
	if (s != u) {
	  k = pd[k - 1];
	  tmp += ((y[k] - yf[k]) / sstar
		  - y[k] / tstar) * (x[u + r * i] - x[s + r * i]) / y[k];
	}
      }
      der[u + i * r] = tmp * ssq;
    }
  }
}

/*  From here on, code borrowed from nnet library  */


static void errmsg(char *string)
{
  PROBLEM "%s\n", string RECOVER(NULL_ENTRY);
}

static double *vect(int n)
{
  double *v;
  v = (double *) Calloc(n, double);
  if (!v) errmsg("allocation failure in vect()");
  return v;
}

static void free_vect(double *v)
{
  Free(v);
}

static double **Lmatrix(int n)
{
  int     i;
  double **m;

  m = (double **) Calloc(n, double *);
  if (!m) errmsg("fail1 in Lmatrix()");

  for (i = 0; i < n; i++) {
    m[i] = (double *) Calloc((i + 1), double);
    if (!m[i]) errmsg("fail2 in Lmatrix()");
  }
  return m;
}

static void free_Lmatrix(double **m, int n)
{
  int     i;
  for (i = n - 1; i >= 0; i--) Free(m[i]);
  Free(m);
}

typedef unsigned char Boolean;

#define false 0

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0
#define abstol 		1.0e-2
#define reltol 		1.0e-3
#define REPORT		5


/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */

static void
vmmin(int n, double *b, double *Fmin, int maxit, int trace)
{
  Boolean accpoint, enough;
  double *g, *t, *X, *c, **B;
  int     count, funcount, gradcount;
  double  f, gradproj;
  int     i, j, ilast, iter = 0;
  double  s, steplength;
  double  D1, D2;

  g = vect(n);
  t = vect(n);
  X = vect(n);
  c = vect(n);
  B = Lmatrix(n);
  f = fminfn(b);
  if (trace) {
    printf("initial  value %f \n", f);
    fflush(stdout);
  }
  {
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(b, g);
    iter++;
    ilast = gradcount;

    do {
      if (ilast == gradcount) {
	for (i = 0; i < n; i++) {
	  for (j = 0; j < i; j++)
	    B[i][j] = 0.0;
	  B[i][i] = 1.0;
	}
      }

      for (i = 0; i < n; i++) {
	X[i] = b[i];
	c[i] = g[i];
      }
      gradproj = 0.0;
      for (i = 0; i < n; i++) {
	s = 0.0;
	for (j = 0; j <= i; j++) s -= B[i][j] * g[j];
	for (j = i + 1; j < n; j++) s -= B[j][i] * g[j];
	t[i] = s;
	gradproj += s * g[i];
      }

      if (gradproj < 0.0) {	/* search direction is downhill */
	steplength = 1.0;
	accpoint = false;
	do {
	  count = 0;
	  for (i = 0; i < n; i++) {
	    b[i] = X[i] + steplength * t[i];
	    if (reltest + X[i] == reltest + b[i])	/* no change */
	      count++;
	  }
	  if (count < n) {
	    f = fminfn(b);
	    funcount++;
	    accpoint = (f <= *Fmin + gradproj * steplength * acctol);

	    if (!accpoint) {
	      steplength *= stepredn;
	    }
	  }
	} while (!(count == n || accpoint));
	enough = (f > abstol) && (f < (1.0 - reltol) * (*Fmin));
	/* stop if value if small or if relative change is low */
	if (!enough) count = n;
	if (count < n) {	/* making progress */
	  *Fmin = f;
	  fmingr(b, g);
	  gradcount++;
	  iter++;
	  D1 = 0.0;
	  for (i = 0; i < n; i++) {
	    t[i] = steplength * t[i];
	    c[i] = g[i] - c[i];
	    D1 += t[i] * c[i];
	  }
	  if (D1 > 0) {
	    D2 = 0.0;
	    for (i = 0; i < n; i++) {
	      s = 0.0;
	      for (j = 0; j <= i; j++)
		s += B[i][j] * c[j];
	      for (j = i + 1; j < n; j++)
		s += B[j][i] * c[j];
	      X[i] = s;
	      D2 += s * c[i];
	    }
	    D2 = 1.0 + D2 / D1;
	    for (i = 0; i < n; i++) {
	      for (j = 0; j <= i; j++)
		B[i][j] += (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
	    }
	  }
	  else {		/* D1 < 0 */
	    ilast = gradcount;
	  }
	}
	else {		/* no progress */
	  if (ilast < gradcount) {
	    count = 0;
	    ilast = gradcount;
	  }
	}
      }
      else {			/* uphill search */
	count = 0;
	if (ilast == gradcount) count = n;
	else ilast = gradcount;
	/* Resets unless has just been reset */
      }
      if (iter % REPORT == 0 && trace) {
	printf("iter%4d value %f\n", iter, f);
	fflush(stdout);
      } if (iter >= maxit)
	break;
      if (gradcount - ilast > 2 * n)
	ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
  }
  if (trace) {
    printf("final  value %f \n", *Fmin);
    if (iter < maxit)
      printf("converged\n");
    else
      printf("stopped after %i iterations\n", iter);
  }
  free_vect(g);
  free_vect(t);
  free_vect(X);
  free_vect(c);
  free_Lmatrix(B, n);
}

/* -----------------------------------------------------------------
 *  Former ucv.c
 */

#if !defined(PI)		/* it is currently defined in
				   S_tokens.h */
#define PI 3.14159265
#endif
#define DELMAX 1000
/* Avoid slow and possibly error-producing underflows by cutting off at
   plus/minus sqrt(DELMAX) std deviations */
/* Formulae (6.67) and (6.69) of Scott (1992), the latter corrected. */

void
VR_ucv_bin(Sint *n, Sint *nb, singl *d, Sint *x, singl *h, singl *u)
{
  int     i, nn = *n, nbin = *nb;
  singl   delta, hh = (*h) / 4, sum, term;

  sum = 0.0;
  for (i = 0; i < nbin; i++) {
    delta = i * (*d) / hh;
    delta *= delta;
    if (delta >= DELMAX) break;
    term = exp(-delta / 4) - sqrt(8.0) * exp(-delta / 2);
    sum += term * x[i];
  }
  *u = 1 / (2 * nn * hh * sqrt(PI)) + sum / (nn * nn * hh * sqrt(PI));
}

void
VR_bcv_bin(Sint *n, Sint *nb, singl *d, Sint *x, singl *h, singl *u)
{
  int     i, nn = *n, nbin = *nb;
  singl   delta, hh = (*h) / 4, sum, term;

  sum = 0.0;
  for (i = 0; i < nbin; i++) {
    delta = i * (*d) / hh;
    delta *= delta;
    if (delta >= DELMAX) break;
    term = exp(-delta / 4) * (delta * delta - 12 * delta + 12);
    sum += term * x[i];
  }
  *u = 1 / (2 * nn * hh * sqrt(PI)) + sum / (64 * nn * nn * hh * sqrt(PI));
}


void
VR_phi4_bin(Sint *n, Sint *nb, singl *d, Sint *x, singl *h, singl *u)
{
  int     i, nn = *n, nbin = *nb;
  singl   delta, sum, term;

  sum = 0.0;
  for (i = 0; i < nbin; i++) {
    delta = i * (*d) / (*h);
    delta *= delta;
    if (delta >= DELMAX) break;
    term = exp(-delta / 2) * (delta * delta - 6 * delta + 3);
    sum += term * x[i];
  }
  sum = 2 * sum + nn * 3;	/* add in diagonal */
  *u = sum / (nn * (nn - 1) * pow(*h, 5.0) * sqrt(2 * PI));
}

void
VR_phi6_bin(Sint *n, Sint *nb, singl *d, Sint *x, singl *h, singl *u)
{
  int     i, nn = *n, nbin = *nb;
  singl   delta, sum, term;

  sum = 0.0;
  for (i = 0; i < nbin; i++) {
    delta = i * (*d) / (*h);
    delta *= delta;
    if (delta >= DELMAX) break;
    term = exp(-delta / 2) *
      (delta * delta * delta - 15 * delta * delta + 45 * delta - 15);
    sum += term * x[i];
  }
  sum = 2 * sum - 15 * nn;	/* add in diagonal */
  *u = sum / (nn * (nn - 1) * pow(*h, 7.0) * sqrt(2 * PI));
}

void
VR_den_bin(Sint *n, Sint *nb, singl *d, singl *x, Sint *cnt)
{
  int     i, j, ii, jj, iij, nn = *n;
  singl   xmin, xmax, rang, dd;

  for (i = 0; i < *nb; i++) cnt[i] = 0;
  xmin = xmax = x[0];
  for (i = 1; i < nn; i++) {
    xmin = min9(xmin, x[i]);
    xmax = max9(xmax, x[i]);
  }
  rang = (xmax - xmin) * 1.01;
  *d = dd = rang / (*nb);
  for (i = 1; i < nn; i++) {
    ii = x[i] / dd;
    for (j = 0; j < i; j++) {
      jj = x[j] / dd;
      iij = abs9((ii - jj));
      cnt[iij]++;
    }
  }
}
