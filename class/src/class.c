/*
 *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 */

#include <S.h>
/* for DOUBLE_XMAX and PROBLEM */

#define EPS 1e-4		/* relative test of equality of distances */

#include "verS.h"

static void
errmsg(char *string)
{
    PROBLEM "%s\n", string RECOVER(NULL_ENTRY);
}

void
VR_knn1(Sint *pntr, Sint *pnte, Sint *p, double *train, Sint *class,
	double *test, Sint *res, Sint *votes, Sint *nc, double *dists)
{
    int   npat, index, i, j, k, ntr = *pntr, nte = *pnte, nind, ntie, *ind;
    double dm, dist, tmp;

    S_EVALUATOR
	RANDIN;
    ind = Calloc(ntr, int);
    for (npat = 0; npat < nte; npat++) {
	dm = DOUBLE_XMAX;
	for (j = 0; j < ntr; j++) {
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = test[npat + k * nte] - train[j + k * ntr];
		dist += tmp * tmp;
	    }
	    if (dist <= dm * (1 + EPS)) {
		if (dist < dm * (1 - EPS))
		    nind = 0;
		else
		    nind++;
		dm = dist;
		ind[nind] = j;
	    }
	}
	for (i = 0; i < *nc; i++)
	    votes[i] = 0;
/* nind is the number of tied minima, minus one */
	if (nind == 0)
	    index = class[ind[0]];
	else {
	    for (i = 0; i <= nind; i++)
		votes[class[ind[i]]]++;
	    j = votes[0];
/*
     This uses `reservoir sampling' to choose amongst ties at random
     on a single pass.

 */
	    index = 0;
	    ntie = 1;
	    for (i = 1; i <= *nc; i++)
		if (votes[i] > j) {
		    ntie = 1;
		    index = i;
		    j = votes[i];
		} else if (votes[i] == j) {
		    if (++ntie * UNIF < 1.0)
			index = i;
		}
	}
	res[npat] = index;
	dists[npat] = dm;
    }
    RANDOUT;
    Free(ind);
}


#define MAX_TIES 1000
/* Not worth doing this dynamically -- limits k + # ties + fence, in fact */


void
VR_knn(Sint *kin, Sint *lin, Sint *pntr, Sint *pnte, Sint *p, 
       double *train, Sint *class, double *test, Sint *res, double *pr, 
       Sint *votes, Sint *nc, Sint *cv, Sint *use_all)
{
    int   i, index, j, k, k1, kinit = *kin, kn, l = *lin, mm, npat, ntie,
          ntr = *pntr, nte = *pnte, extras;
    int   pos[MAX_TIES];
    int   j1, j2, needed, t;
    double dist, tmp, nndist[MAX_TIES];

    S_EVALUATOR
	RANDIN;
/*
    Use a `fence' in the (k+1)st position to avoid special cases.
    Simple insertion sort will suffice since k will be small.
 */

    for (npat = 0; npat < nte; npat++) {
	kn = kinit;
	for (k = 0; k < kn; k++)
	    nndist[k] = 0.99 * DOUBLE_XMAX;
	for (j = 0; j < ntr; j++) {
	    if ((*cv > 0) && (j == npat))
		continue;
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = test[npat + k * nte] - train[j + k * ntr];
		dist += tmp * tmp;
	    }
/* Use `fuzz' since distance computed could depend on order of coordinates */
	    if (dist <= nndist[kinit - 1] * (1 + EPS))
		for (k = 0; k <= kn; k++)
		    if (dist < nndist[k]) {
			for (k1 = kn; k1 > k; k1--) {
			    nndist[k1] = nndist[k1 - 1];
			    pos[k1] = pos[k1 - 1];
			}
			nndist[k] = dist;
			pos[k] = j;
/* Keep an extra distance if the largest current one ties with current kth */
			if (nndist[kn] <= nndist[kinit - 1])
			    if (++kn == MAX_TIES - 1)
				errmsg("too many ties in knn");
			break;
		    }
	    nndist[kn] = 0.99 * DOUBLE_XMAX;
	}

	for (j = 0; j <= *nc; j++)
	    votes[j] = 0;
	if (*use_all) {
	    for (j = 0; j < kinit; j++)
		votes[class[pos[j]]]++;
	    extras = 0;
	    for (j = kinit; j < kn; j++) {
		if (nndist[j] > nndist[kinit - 1] * (1 + EPS))
		    break;
		extras++;
		votes[class[pos[j]]]++;
	    }
	} else {
	    for (j = 0; j < kinit; j++) {
		if (nndist[j] >= nndist[kinit - 1] * (1 - EPS))
		    break;
		votes[class[pos[j]]]++;
	    }
/* Use reservoir sampling to choose amongst the ties */
	    j1 = j;
	    needed = kinit - j1;
	    t = needed;
	    for (j = j1; j < kn; j++) {
		if (nndist[j] > nndist[kinit - 1] * (1 + EPS))
		    break;
		if (++t * UNIF < needed) {
		    j2 = j1 + UNIF * needed;
		    class[pos[j2]] = class[pos[j]];
		}
	    }
	    extras = 0;
	    for (j = j1; j < kinit; j++)
		votes[class[pos[j]]]++;
	}
/* Use reservoir sampling to choose amongst the ties */
	ntie = 1;
	if (l > 0)
	    mm = l - 1 + extras;
	else
	    mm = 0;
	index = 0;
	for (i = 1; i <= *nc; i++)
	    if (votes[i] > mm) {
		ntie = 1;
		index = i;
		mm = votes[i];
	    } else if (votes[i] == mm && votes[i] >= l) {
		if (++ntie * UNIF < 1.0)
		    index = i;
	    }
	res[npat] = index;
	pr[npat] = (double) mm / (kinit + extras);
    }
    RANDOUT;
}


#define min9(a,b) ((a < b)?a:b)

void
VR_olvq(double *alpha, Sint *pn, Sint *p, double *x, Sint *cl, 
	Sint *pncodes, double *xc, Sint *clc, Sint *niter, 
	Sint *iters)
{
    int   index, iter, j, k, n = *pn, ncodes = *pncodes, npat, s;
    double *al;
    double dist, dm, tmp;

    al = Calloc(ncodes, double);
    for (j = 0; j < ncodes; j++)
	al[j] = *alpha;
    for (iter = 0; iter < *niter; iter++) {
	npat = iters[iter];
	dm = DOUBLE_XMAX;
	for (j = 0; j < ncodes; j++) {
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = x[npat + k * n] - xc[j + k * ncodes];
		dist += tmp * tmp;
	    }
	    if (dist < dm) {
		dm = dist;
		index = j;
	    }
	}
	s = 2 * (clc[index] == cl[npat]) - 1;
	for (k = 0; k < *p; k++)
	    xc[index + k * ncodes] += s * al[index] *
		(x[npat + k * n] - xc[index + k * ncodes]);
	al[index] = min9(*alpha, al[index] / (1 + s * al[index]));
    }
    Free(al);
}

void
VR_lvq1(double *alpha, Sint *pn, Sint *p, double *x, Sint *cl, 
	Sint *pncodes, double *xc, Sint *clc, Sint *niter, 
	Sint *iters)
{
    int   index, iter, j, k, n = *pn, ncodes = *pncodes, npat, s;
    double alpha_t;
    double dist, dm, tmp;

    for (iter = 0; iter < *niter; iter++) {
	npat = iters[iter];
	alpha_t = *alpha * (*niter - iter) / (double) *niter;
	dm = DOUBLE_XMAX;
	for (j = 0; j < ncodes; j++) {
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = x[npat + k * n] - xc[j + k * ncodes];
		dist += tmp * tmp;
	    }
	    if (dist < dm) {
		dm = dist;
		index = j;
	    }
	}
	s = 2 * (clc[index] == cl[npat]) - 1;
	for (k = 0; k < *p; k++)
	    xc[index + k * ncodes] += s * alpha_t *
		(x[npat + k * n] - xc[index + k * ncodes]);
    }
}

void
VR_lvq2(double *alpha, double *win, Sint *pn, Sint *p, double *x, 
	Sint *cl, Sint *pncodes, double *xc, Sint *clc, 
	Sint *niter, Sint *iters)
{
    int   index, iter, j, k, n = *pn, ncodes = *pncodes, nindex, npat, ntmp;
    double alpha_t;
    double dist, dm, ndm, tmp;

    for (iter = 0; iter < *niter; iter++) {
	npat = iters[iter];
	alpha_t = *alpha * (*niter - iter) / (double) *niter;
	ndm = dm = DOUBLE_XMAX;

	/* Find two nearest codebook vectors */
	for (j = 0; j < ncodes; j++) {
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = x[npat + k * n] - xc[j + k * ncodes];
		dist += tmp * tmp;
	    }
	    if (dist < dm) {
		ndm = dm;
		nindex = index;
		dm = dist;
		index = j;
	    } else if (dist < ndm) {
		ndm = dist;
		nindex = j;
	    }
	}
	if (clc[index] != clc[nindex]) {
	    if (((clc[index] == cl[npat]) || (clc[nindex] == cl[npat]))
		&& dm / ndm > (1 - *win) / (1 + *win)) {
		if (clc[nindex] == cl[npat]) {
		    ntmp = index;
		    index = nindex;
		    nindex = ntmp;
		}
		for (k = 0; k < *p; k++) {
		    xc[index + k * ncodes] += alpha_t *
			(x[npat + k * n] - xc[index + k * ncodes]);
		    xc[nindex + k * ncodes] -= alpha_t *
			(x[npat + k * n] - xc[nindex + k * ncodes]);
		}
	    }
	}
    }
}

void
VR_lvq3(double *alpha, double *win, double *epsilon, Sint *pn, Sint *p,
	double *x, Sint *cl, Sint *pncodes, double *xc, Sint *clc,
	Sint *niter, Sint *iters)
{
    int   index, iter, j, k, n = *pn, ncodes = *pncodes, nindex, npat, ntmp;
    double alpha_t;
    double dist, dm, ndm, tmp;

    for (iter = 0; iter < *niter; iter++) {
	npat = iters[iter];
	alpha_t = *alpha * (*niter - iter) / (double) *niter;
	ndm = dm = DOUBLE_XMAX;
	/* Find two nearest codebook vectors */
	for (j = 0; j < ncodes; j++) {
	    dist = 0.0;
	    for (k = 0; k < *p; k++) {
		tmp = x[npat + k * n] - xc[j + k * ncodes];
		dist += tmp * tmp;
	    }
	    if (dist < dm) {
		ndm = dm;
		nindex = index;
		dm = dist;
		index = j;
	    } else if (dist < ndm) {
		ndm = dist;
		nindex = j;
	    }
	}
	if (clc[index] != clc[nindex]) {
	    if (((clc[index] == cl[npat]) || (clc[nindex] == cl[npat]))
		&& dm / ndm > (1 - *win) / (1 + *win)) {
		if (clc[nindex] == cl[npat]) {
		    ntmp = index;
		    index = nindex;
		    nindex = ntmp;
		}
		for (k = 0; k < *p; k++) {
		    xc[index + k * ncodes] += alpha_t *
			(x[npat + k * n] - xc[index + k * ncodes]);
		    xc[nindex + k * ncodes] -= alpha_t *
			(x[npat + k * n] - xc[nindex + k * ncodes]);
		}
	    }
	} else if (clc[index] == cl[npat]) {
	    for (k = 0; k < *p; k++) {
		xc[index + k * ncodes] += *epsilon * alpha_t *
		    (x[npat + k * n] - xc[index + k * ncodes]);
		xc[nindex + k * ncodes] += *epsilon * alpha_t *
		    (x[npat + k * n] - xc[nindex + k * ncodes]);
	    }
	}
    }
}
