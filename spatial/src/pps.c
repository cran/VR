/*
 *  spatial/pps.c by W. N. Venables and B. D. Ripley.  Copyright (C) 1994-9
 */

#include <S.h>
#include <math.h>
#define min9(a,b) ((a < b)?a:b)
#define max9(a,b) ((a > b)?a:b)
#include "verS.h"

static singl xl0, yl0, xu0, yu0;
static singl pi = 3.14159265;


static void
errmsg(char *string)
{
  PROBLEM "%s\n", string RECOVER(NULL_ENTRY);
}

static void
testinit(void)
{
  if ((xu0 == xl0) || (yu0 == yl0))
    errmsg("not initialized -- use ppregion");
}


void
VR_ppset(singl *area)
{
  xl0 = *area++;
  xu0 = *area++;
  yl0 = *area++;
  yu0 = *area;
  testinit();
}

void
VR_ppget(singl *xx)
{
  singl *z = xx;
  testinit();
  *z++ = xl0; *z++ = xu0; *z++ = yl0; *z++ = yu0;
}

static singl 
edge(singl x, singl y, singl a)
{
  singl   b, c, c1, c2, r[6], w;
  int     i;
  w = x - xl0;
  if (w > y - yl0) w = y - yl0;
  if (w > xu0 - x) w = xu0 - x;
  if (w > yu0 - y) w = yu0 - y;
  if (a <= w) return (0.5);
  r[4] = r[0] = x - xl0;
  r[5] = r[1] = yu0 - y;
  r[2] = xu0 - x;
  r[3] = y - yl0;
  b = 0.0;
  for (i = 1; i <= 4; i++)
    if (r[i] < a) {
      if (r[i] == 0.0)
	b += pi;
      else {
	c = acos(r[i] / a);
	c1 = atan(r[i - 1] / r[i]);
	c2 = atan(r[i + 1] / r[i]);
	b += min9(c, c1);
	b += min9(c, c2);
      }
    }
  if (b < 6.28) return (1.0 / (2.0 - b / pi));
  return (0.0);
}

void
VR_sp_pp2(singl *x, singl *y, Sint *npt, Sint *k,
	  singl *h, singl *dmin, singl *lm, singl *fs)
{
  int     n = *npt, kk = *k, k1, i, j, ib;
  singl   ax, ay, xi, yi, sarea, g, dm, alm;
  double  a, x1, y1, rr, fss = *fs, fs1, s1;
  testinit();
  ax = xu0 - xl0;
  ay = yu0 - yl0;
  sarea = sqrt(ax * ay);
  dm = fss;
  g = 2.0 / (n * n);
  fs1 = min9(fss, 0.5 * sqrt(ax * ax + ay * ay));
  s1 = kk / fss;
  k1 = floor(s1 * fs1 + 1e-3);
  *k = k1;
  rr = fs1 * fs1;
  for (i = 0; i < kk; i++) h[i] = 0.0;
  for (i = 1; i < n; i++) {
    xi = x[i];
    yi = y[i];
    for (j = 0; j < i; j++) {
      x1 = x[j] - xi;
      y1 = y[j] - yi;
      a = x1 * x1 + y1 * y1;
      if (a < rr) {
	a = sqrt(a);
	dm = min9(a, dm);
	ib = floor(s1 * a);
	if (ib < k1)
	  h[ib] += g * (edge(xi, yi, a) + edge(x[j], y[j], a));
      }
    }
  }
  a = 0.0;
  alm = 0.0;
  for (i = 0; i < k1; i++) {
    a += h[i];
    h[i] = sqrt(a / pi) * sarea;
    alm = max9(alm, fabs(h[i] - (i + 1) / s1));
  }
  *dmin = dm;
  *lm = alm;
}

void
VR_pdata(Sint *npt, singl *x, singl *y)
{
  S_EVALUATOR
  int     i;
  singl   ax, ay;
  testinit();
  ax = xu0 - xl0;
  ay = yu0 - yl0;
  RANDIN;
  for (i = 0; i < *npt; i++) {
    x[i] = xl0 + ax * UNIF;
    y[i] = yl0 + ay * UNIF;
  }
  RANDOUT;
}


void
VR_simpat(Sint *npt, singl *x, singl *y, singl *c,
	  singl *r, Sint *init)
{
  S_EVALUATOR
  int     i, id, j, mm, n = *npt;
  singl   cc, rr, ax, ay, d, x1, y1, u;
  testinit();
  cc = *c;
  if (cc >= 1.0) {
    VR_pdata(npt, x, y);
    return;
  }
  RANDIN;
  ax = xu0 - xl0;
  ay = yu0 - yl0;
  rr = (*r) * (*r);
  mm = 4 * n;
  if (*init > 0) mm = 10 * mm;
  for (i = 1; i <= mm; i++) {
    id = floor(n *UNIF);
    x[id] = x[0];
    y[id] = y[0];
    do {
      x[0] = xl0 + ax *UNIF;
      y[0] = yl0 + ay *UNIF;
      u =UNIF;
      d = 1.0;
      for (j = 1; j < n; j++) {
	x1 = x[j] - x[0];
	y1 = y[j] - y[0];
	if (x1 * x1 + y1 * y1 < rr) {
	  d *= cc;
	  if (d < u) continue;
	}
      }
    }
    while (d < u);
  }
  RANDOUT;
}

void
VR_simmat(Sint *npt, singl *x, singl *y, singl *r)
{
  S_EVALUATOR
  int     i, icnt, j, n = *npt;
  singl   x1, y1, rr, ax, ay;
  testinit();
  RANDIN;
  ax = xu0 - xl0;
  ay = yu0 - yl0;
  rr = (*r) * (*r);
  for (i = 0; i < n; i++) {
    do {
      icnt = 0;
      x[i] = xl0 + ax *UNIF;
      y[i] = yl0 + ay *UNIF;
      if (i > 0)
	for (j = 0; j < i; j++) {
	  x1 = x[i] - x[j];
	  y1 = y[i] - y[j];
	  if (x1 * x1 + y1 * y1 < rr) {
	    icnt = 1;
	    break;
	  }
	}
    }
    while (icnt > 0);
  }
  RANDOUT;
}

void
VR_plike(singl *x, singl *y, Sint *npt, singl *c,
	 singl *r, Sint *ng, singl *target, singl *res)
{
  singl   ar, rr, suma = 0, sumb = 0, xi, yi, x1, y1, c1, cc = *c;
  int     ic, i1, i2, j, n = *npt, g = *ng;
  testinit();
  ar = (*r);
  rr = ar * ar;
  if(cc <= 0) {
    *res = - *target;
    return;
  }
  for (i1 = 0; i1 < g; i1++) {
    xi = xl0 + ar + (xu0 - xl0 - 2 * ar) * i1 / (g - 1);
    for (i2 = 0; i2 < g; i2++) {
      yi = yl0 + ar + (yu0 - yl0 - 2 * ar) * i2 / (g - 1);
      ic = 0;
      for (j = 0; j < n; j++) {
	x1 = x[j] - xi;
	y1 = y[j] - yi;
	if (x1 * x1 + y1 * y1 < rr) ic += 1;
      }
      c1 = (ic > 0) ? pow(cc, (double) ic) : 1;
      suma += ic * c1;
      sumb += c1;
    }
  }
  *res = suma / sumb - *target;
}
