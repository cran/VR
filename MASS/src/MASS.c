/*
 *  MASS/MASS.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-2004
 */

#include <R.h>
#include <R_ext/Applic.h>

#ifndef max
#  define max(a,b) ((a) > (b) ? (a) : (b))
#  define min(a,b) ((a) < (b) ? (a) : (b))
#endif

#define abs9(a) (a > 0 ? a:-a)

/* -----------------------------------------------------------------
 *  Former sammon.c
 */

void
VR_sammon(double *dd, Sint *nn, Sint *kd, double *Y, Sint *niter,
	  double *stress, Sint *trace, double *aa, double *tol)
{
    int   i, j, k, m, n = *nn, nd = *kd;
    double *xu, *xv, *e1, *e2;
    double dpj, dq, dr, dt;
    double xd, xx;
    double e, epast, eprev, tot, d, d1, ee, magic = *aa;

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
	    if (ISNAN(d)) continue;
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
    if (*trace) Rprintf("Initial stress        : %7.5f\n", e);
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
		dt = dd[k * n + j];
		if (ISNAN(dt)) continue;
		d1 = 0.0;
		for (m = 0; m < nd; m++) {
		    xd = Y[j + m * n] - Y[k + m * n];
		    d1 += xd * xd;
		    xv[m] = xd;
		}
		dpj = sqrt(d1);

		/* Calculate derivatives */
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
		if (ISNAN(d)) continue;
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
		Rprintf("stress after %3d iters: %7.5f\n", i - 1, e);
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
		Rprintf("stress after %3d iters: %7.5f, magic = %5.3f\n", i, e, magic);
	    }
	    if (e > epast - *tol)
		break;
	    epast = e;
	}
    }
    *stress = e;
    Free(xu);
    Free(xv);
    Free(e1);
    Free(e2);
}

/*
 * ----------------------------------------------------------
 *  Former isoMDS.c

    C code for mds S-Plus library, which implements Kruskal's MDS.
    (c) B.D. Ripley, May 1995.
 *
 */

static Sint *ord;		/* ranks of dissimilarities */
static Sint *ord2;		/* inverse ordering (which one is rank i?) */
static Sint n;			/* number of  dissimilarities */
static Sint nr;			/* number of data points */
static Sint nc;			/* # cols of  fitted configuration */
static int dimx;		/* Size of configuration array */
static double *x;		/* configuration */
static double *d;		/* dissimilarities */
static double *y;		/* fitted distances (in rank of d order) */
static double *yf;		/* isotonic regression fitted values (ditto) */

void
VR_mds_fn(double *, double *, Sint *, double *, Sint *,
	  double *, Sint *, Sint *, double *, Sint *);

/*
 *  Download the data.
 */
void
VR_mds_init_data(Sint *pn, Sint *pc, Sint *pr, Sint *orde,
		 Sint *ordee, double *xx)
{
    int   i;

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
    int   r1, r2, c, index;
    double tmp, tmp1;

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
    for (index = 0; index < n; index++)
	y[index] = d[ord[index]];
}

static double
fminfn(int nn, double *x, void *dummy)
{
    double ssq;
    Sint  do_derivatives = 0;

    calc_dist(x);
    VR_mds_fn(y, yf, &n, &ssq, ord2, x, &nr, &nc, 0, &do_derivatives);
    return (ssq);
}

static void
fmingr(int nn, double *x, double *der, void *dummy)
{
    double ssq;
    Sint  do_derivatives = 1;

    calc_dist(x);
    VR_mds_fn(y, yf, &n, &ssq, ord2, x, &nr, &nc, der, &do_derivatives);
}

#define abstol 		1.0e-2
#define REPORT		5

void
VR_mds_dovm(double *val, Sint *maxit, Sint *trace, double *xx, double *tol)
{
    int   i, ifail, fncount, grcount, *mask;

    mask = (int *) R_alloc(dimx, sizeof(int));
    for (i = 0; i < dimx; i++) mask[i] = 1;
    vmmin(dimx, x, val, fminfn, fmingr, (int) *maxit, (int) *trace, mask,
	  abstol, *tol, REPORT, NULL, &fncount, &grcount, &ifail);
    for (i = 0; i < dimx; i++)
	xx[i] = x[i];
}

/*
 *  Does isotonic regression.
 */

void
VR_mds_fn(double *y, double *yf, Sint *pn, double *pssq, Sint *pd,
	  double *x, Sint *pr, Sint *pncol, double *der,
	  Sint *do_derivatives)
{
    int   n = *pn, i, ip=0, known, u, s, r = *pr, ncol = *pncol, k=0;
    double tmp, ssq, *yc, slope, tstar, sstar;

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
		if (s > u)
		    k = r * u - u * (u + 1) / 2 + s - u;
		else if (s < u)
		    k = r * s - s * (s + 1) / 2 + u - s;
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

/* -----------------------------------------------------------------
 *  Former ucv.c
 */

#if !defined(M_PI)
#  define M_PI 3.141592653589793238462643383280
#endif
#define DELMAX 1000
/* Avoid slow and possibly error-producing underflows by cutting off at
   plus/minus sqrt(DELMAX) std deviations */
/* Formulae (6.67) and (6.69) of Scott (1992), the latter corrected. */

void
VR_ucv_bin(Sint *n, Sint *nb, Sfloat *d, Sint *x, Sfloat *h, Sfloat *u)
{
    int   i, nn = *n, nbin = *nb;
    Sfloat delta, hh = (*h) / 4, sum, term;

    sum = 0.0;
    for (i = 0; i < nbin; i++) {
	delta = i * (*d) / hh;
	delta *= delta;
	if (delta >= DELMAX) break;
	term = exp(-delta / 4) - sqrt(8.0) * exp(-delta / 2);
	sum += term * x[i];
    }
    *u = 1 / (2 * nn * hh * sqrt(M_PI)) + sum / (nn * nn * hh * sqrt(M_PI));
}

void
VR_bcv_bin(Sint *n, Sint *nb, Sfloat *d, Sint *x, Sfloat *h, Sfloat *u)
{
    int   i, nn = *n, nbin = *nb;
    Sfloat delta, hh = (*h) / 4, sum, term;

    sum = 0.0;
    for (i = 0; i < nbin; i++) {
	delta = i * (*d) / hh;
	delta *= delta;
	if (delta >= DELMAX) break;
	term = exp(-delta / 4) * (delta * delta - 12 * delta + 12);
	sum += term * x[i];
    }
    *u = 1 / (2 * nn * hh * sqrt(M_PI)) + sum / (64 * nn * nn * hh * sqrt(M_PI));
}


void
VR_phi4_bin(Sint *n, Sint *nb, Sfloat *d, Sint *x, Sfloat *h, Sfloat *u)
{
    int   i, nn = *n, nbin = *nb;
    Sfloat delta, sum, term;

    sum = 0.0;
    for (i = 0; i < nbin; i++) {
	delta = i * (*d) / (*h);
	delta *= delta;
	if (delta >= DELMAX) break;
	term = exp(-delta / 2) * (delta * delta - 6 * delta + 3);
	sum += term * x[i];
    }
    sum = 2 * sum + nn * 3;	/* add in diagonal */
    *u = sum / (nn * (nn - 1) * pow(*h, 5.0) * sqrt(2 * M_PI));
}

void
VR_phi6_bin(Sint *n, Sint *nb, Sfloat *d, Sint *x, Sfloat *h, Sfloat *u)
{
    int   i, nn = *n, nbin = *nb;
    Sfloat delta, sum, term;

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
    *u = sum / (nn * (nn - 1) * pow(*h, 7.0) * sqrt(2 * M_PI));
}

void
VR_den_bin(Sint *n, Sint *nb, Sfloat *d, Sfloat *x, Sint *cnt)
{
    int   i, j, ii, jj, iij, nn = *n;
    Sfloat xmin, xmax, rang, dd;

    for (i = 0; i < *nb; i++) cnt[i] = 0;
    xmin = xmax = x[0];
    for (i = 1; i < nn; i++) {
	xmin = min(xmin, x[i]);
	xmax = max(xmax, x[i]);
    }
    rang = (xmax - xmin) * 1.01;
    *d = dd = rang / (*nb);
    for (i = 1; i < nn; i++) {
	ii = (int) (x[i] / dd);
	for (j = 0; j < i; j++) {
	    jj = (int) (x[j] / dd);
	    iij = abs9((ii - jj));
	    cnt[iij]++;
	}
    }
}

/*
#include "R_ext/Rdynload.h"

R_CMethodDef CEntries[] = {
    {"VR_bcv_bin", (DL_FUNC) &VR_bcv_bin, 6},
    {"VR_den_bin", (DL_FUNC) &VR_den_bin, 5},
    {"VR_mds_dovm", (DL_FUNC) &VR_mds_dovm, 5},
    {"VR_mds_fn", (DL_FUNC) &VR_mds_fn, 10},
    {"VR_mds_init_data", (DL_FUNC) &VR_mds_init_data, 6},
    {"VR_mds_unload", (DL_FUNC) &VR_mds_unload, 0},
    {"VR_phi4_bin", (DL_FUNC) &VR_phi4_bin, 6},
    {"VR_phi6_bin", (DL_FUNC) &VR_phi6_bin, 6},
    {"VR_sammon", (DL_FUNC) &VR_sammon, 9},
    {"VR_ucv_bin", (DL_FUNC) &VR_ucv_bin, 6},
    {NULL, NULL, 0}
};

void R_init_MASS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
*/
