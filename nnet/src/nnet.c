/*  nnet/nnet.c by W. N. Venables and B. D. Ripley  Copyright (C) 1992-9
 *
 * weights are stored in order of their destination unit.
 * the array Conn gives the source unit for the weight (0 = bias unit)
 * the array Nconn gives the number of first weight connecting to each unit,
 * so the weights connecting to unit i are Nconn[i] ... Nconn[i+1] - 1.
 *
 */

#include <S.h>
#include <stdio.h>
#include <math.h>
#if defined(SPLUS_VERSION) && SPLUS_VERSION >= 4000 && SPLUS_VERSION < 5000
#include <newredef.h>
#endif

#ifdef USING_R
typedef double Sdata;
#define DataByCopy 1  /* use this to copy data into C, e.g. for R */
#else
typedef float Sdata;  /* type of data, weights and Hessian in caller */
#endif
#include "verS.h"


static double *vect(int n);
static void free_vect(double *v);
static double **matrix(int nrh, int nch);
static void free_matrix(double **m, int nrh, int nch);
static void vmmin(int n, double *b, double *Fmin, int maxit, int trace,
		  Sint *mask, double abstol, double reltol);
static void fpass(Sdata *input, Sdata *goal, Sdata wx, int nr);
static void Build_Net(int ninputs, int nhidden, int noutputs);
static Sdata *svect(int n);
static void free_svect(Sdata *v);


static int Epoch;
static double *Decay;
static double TotalError;

static int Nunits;
static int Ninputs;
static int FirstHidden;
static int FirstOutput;
static int Noutputs;
static int NSunits;		/* number of sigmoid units */
static int Nweights;
static int Entropy;
static int Linout;
static int Softmax;
static int Censored;

static double *Outputs;
static double *ErrorSums;
static double *Errors;
static int *Nconn;
static int *Conn;
static double *wts;
static double *Slopes;
static double *Probs;

static int NTrain;
static int NTest;
static Sdata *TrainIn;
static Sdata *TrainOut;
static Sdata *Weights;

static Sdata *toutputs;

static void
errmsg(char *string)
{
  PROBLEM "%s\n", string RECOVER(NULL_ENTRY);
}

void
VR_set_net(Sint *n, Sint *nconn, Sint *conn,
	   double *decay, Sint *nsunits, Sint *entropy,
	   Sint *softmax, Sint *censored)
{
  int     i;

  Build_Net((int) n[0], (int) n[1], (int) n[2]);
  for (i = 0; i <= Nunits; i++) Nconn[i] = nconn[i];
  Nweights = Nconn[Nunits];
  Conn = Calloc(Nweights, int);
  wts = Calloc(Nweights, double);
  Slopes = Calloc(Nweights, double);
  Probs = Calloc(Nweights, double);
  Decay = Calloc(Nweights, double);
  for (i = 0; i < Nweights; i++) Conn[i] = conn[i];
  Epoch = 0;
  for (i = 0; i < Nweights; i++) Decay[i] = decay[i];
  TotalError = 0.0;
  NSunits = *nsunits;
  Entropy = *entropy;
  Linout = (NSunits < Nunits);
  Softmax = *softmax;
  Censored = *censored;
}

void 
VR_unset_net()
{
  Free(Conn); Free(wts); Free(Slopes); Free(Probs); Free(Decay);
  Free(Nconn); Free(Outputs); Free(ErrorSums); Free(Errors); Free(toutputs);
}

void
VR_set_train(Sint *ntr, Sdata *train, Sdata *weights)
{
  int     i;

  NTrain = *ntr;
#ifdef DataByCopy
  TrainIn = svect(NTrain * Ninputs);
  TrainOut = svect(NTrain * Noutputs);
  Weights = svect(NTrain);
  for (i = 0; i < NTrain*Ninputs; i++) TrainIn[i] = *train++;
  for (i = 0; i < NTrain*Noutputs; i++) TrainOut[i] = *train++;
  for (i = 0; i < NTrain; i++) Weights[i] = *weights++;
#else
  TrainIn = train;
  TrainOut = train + Ninputs * NTrain;
  Weights = weights;
#endif
}

void 
VR_unset_train()
{
#ifdef DataByCopy
  free_svect(TrainIn); free_svect(TrainOut); free_svect(Weights);
#endif

}

void
VR_nntest(Sint *ntest, Sdata *test, Sdata *result, double *inwts)
{
  int     i, j;

  for (i = 0; i < Nweights; i++) wts[i] = inwts[i];
  NTest = *ntest;
  if (Nweights == 0) errmsg("No model set");

  for (i = 0; i < Noutputs; i++) toutputs[i] = 0.5;
  for (j = 0; j < NTest; j++) {
    fpass(test + j, toutputs, 1.0, NTest);
    if (Softmax)
      for (i = 0; i < Noutputs; i++) 
	result[j + NTest*i] = Probs[FirstOutput + i];
    else
      for (i = 0; i < Noutputs; i++)
	result[j + NTest*i] = Outputs[FirstOutput + i];
  }
}


static void
Build_Net(int ninputs, int nhidden, int noutputs)
{
  Nunits = 1 + ninputs + nhidden + noutputs;
  Nconn = Calloc(Nunits+1, int);
  Outputs = Calloc(Nunits, double);
  ErrorSums = Calloc(Nunits, double);
  Errors = Calloc(Nunits, double);
  toutputs = Calloc(Nunits, Sdata);
  Ninputs = ninputs;
  FirstHidden = 1 + ninputs;
  FirstOutput = 1 + ninputs + nhidden;
  Noutputs = noutputs;
  Outputs[0] = 1.0;
}

static double
sigmoid(double sum)
{
  if (sum < -15.0) return (0.0);
  else if (sum > 15.0) return (1.0);
  else return (1.0 / (1.0 + exp(-sum)));
}
#define EPS 1.0E-80

static double
E(double y, double t)
{
  double  dif, sum = 0;
  if (Entropy) {
    if (t > 0) sum -= t * log((y + EPS) / t);
    if (t < 1) sum -= (1 - t) * log((1 - y + EPS) / (1 - t));
  } else {
    dif = y - t;
    sum = dif * dif;
  }
  return (sum);
}


static void
fpass(Sdata *input, Sdata *goal, Sdata wx, int nr)
{
  int     i, j;
  double  sum, t, thisError;

  for (i = 0; i < Ninputs; i++) Outputs[i + 1] = input[i*nr];

  for (j = FirstHidden; j < Nunits; j++) {
    sum = 0.0;
    for (i = Nconn[j]; i < Nconn[j + 1]; i++)
      sum += Outputs[Conn[i]] * wts[i];
    if (j < NSunits) sum = sigmoid(sum);
    Outputs[j] = sum;
  }
  if (Softmax) {
    sum = 0.0;
    /* avoid overflows by re-normalizing */
    t = Outputs[FirstOutput];
    for (i = FirstOutput + 1; i < Nunits; i++)
      if (Outputs[i] > t) t = Outputs[i];
    for (i = FirstOutput; i < Nunits; i++) {
      Probs[i] = exp(Outputs[i] - t);
      sum += Probs[i];
    }
    thisError = 0.0;
    for (i = FirstOutput; i < Nunits; i++) {
      Probs[i] = Probs[i] / sum;
      t = goal[i - FirstOutput];
      if(Censored) {
	if(t == 1) thisError += Probs[i];
      } else if (t > 0) {
	if(Probs[i] > 0) TotalError -= wx * t * log(Probs[i]);
	else TotalError += wx * 1000;
      }
    }
    if (Censored) {
      if (thisError > 0) TotalError -= wx * log(thisError);
      else TotalError += wx * 1000;   
    }
  } else
    for (i = FirstOutput; i < Nunits; i++)
      TotalError += wx * E(Outputs[i], goal[i - FirstOutput]);
}

static double
sigmoid_prime(double value)
{
  return (value * (1.0 - value));
}

static double
sigmoid_prime_prime(double value)
{
  return (value * (1.0 - value) * (1.0 - 2.0 * value));
}

static void
bpass(Sdata *goal, Sdata wx)
{
  int     i, j, cix;
  double  sum, denom;

  if (Softmax) {
    if (Censored) {
      denom = 0.0;
      for (i = FirstOutput; i < Nunits; i++)
	if (goal[i - FirstOutput] == 1) denom += Probs[i];
      for (i = FirstOutput; i < Nunits; i++) {
	ErrorSums[i] = Probs[i];
	if (goal[i - FirstOutput] == 1) ErrorSums[i] -= Probs[i]/denom;
      }
    } else {
      sum = 0.0;
      for (i = FirstOutput; i < Nunits; i++)
	sum += goal[i - FirstOutput];
      for (i = FirstOutput; i < Nunits; i++)
	ErrorSums[i] = sum * Probs[i] - goal[i - FirstOutput];
    }
  } else if (Entropy)
    for (i = FirstOutput; i < Nunits; i++)
      ErrorSums[i] = Outputs[i] - goal[i - FirstOutput];
  else
    for (i = FirstOutput; i < Nunits; i++) {
      ErrorSums[i] = 2 * (Outputs[i] - goal[i - FirstOutput]);
      if (i < NSunits)
	ErrorSums[i] *= sigmoid_prime(Outputs[i]);
    }
  for (i = FirstHidden; i < FirstOutput; i++) ErrorSums[i] = 0.0;

  for (j = Nunits - 1; j >= FirstHidden; j--) {
    Errors[j] = ErrorSums[j];
    if (j < FirstOutput)
      Errors[j] *= sigmoid_prime(Outputs[j]);
    for (i = Nconn[j]; i < Nconn[j + 1]; i++) {
      cix = Conn[i];
      ErrorSums[cix] += Errors[j] * wts[i];
      Slopes[i] += wx * Errors[j] * Outputs[cix];
    }
  }
}

void
VR_dfunc(double *p, double *df, double *fp)
{
  int     i, j;
  double  sum1;

  for (i = 0; i < Nweights; i++) wts[i] = p[i];
  for (j = 0; j < Nweights; j++) Slopes[j] = 2 * Decay[j] * wts[j];
  TotalError = 0.0;
  for (i = 0; i < NTrain; i++) {
    for (j = 0; j < Noutputs; j++) toutputs[j] = TrainOut[i + NTrain*j];
    fpass(TrainIn + i, toutputs, Weights[i], NTrain);
    bpass(toutputs, Weights[i]);
  }
  sum1 = 0.0;
  for (i = 0; i < Nweights; i++) sum1 += Decay[i] * p[i] * p[i];
  *fp = TotalError + sum1;
  for (j = 0; j < Nweights; j++) df[j] = Slopes[j];
  Epoch++;
}

static double
fminfn(double *p)
{
  int     i, j;
  double  sum1;

  for (i = 0; i < Nweights; i++) wts[i] = p[i];
  TotalError = 0.0;
  for (i = 0; i < NTrain; i++) {
    for (j = 0; j < Noutputs; j++) toutputs[j] = TrainOut[i + NTrain*j];
    fpass(TrainIn + i, toutputs, Weights[i], NTrain);
  }
  
  sum1 = 0.0;
  for (i = 0; i < Nweights; i++) sum1 += Decay[i] * p[i] * p[i];
  Epoch++;
  return (TotalError + sum1);
}

static void
fmingr(double *p, double *df)
{
  int     i, j;

  for (i = 0; i < Nweights; i++) wts[i] = p[i];
  for (j = 0; j < Nweights; j++) Slopes[j] = 2 * Decay[j] * wts[j];
  TotalError = 0.0;
  for (i = 0; i < NTrain; i++) {
    for (j = 0; j < Noutputs; j++) toutputs[j] = TrainOut[i + NTrain*j];
    fpass(TrainIn + i, toutputs, Weights[i], NTrain);
    bpass(toutputs, Weights[i]);
  }
  for (j = 0; j < Nweights; j++) df[j] = Slopes[j];
  Epoch++;
}


static double *vect(int n)
{
  double  *v;

  v = Calloc(n, double);
  if (!v) errmsg("allocation failure in vect()");
  return v;
}

static void free_vect(double *v)
{
   Free(v);
}

static double **matrix(int nrh, int nch)
{
  int     i;
  double **m;

  m = Calloc((nrh + 1), double *);
  if (!m) errmsg("allocation failure 1 in matrix()");

  for (i = 0; i <= nrh; i++) {
    m[i] = Calloc((nch + 1), double);
    if (!m[i]) errmsg("allocation failure 2 in matrix()");
  }
  return m;
}

static void free_matrix(double **m, int nrh, int nch)
{
  int     i;

  for (i = nrh; i >= 0; i--) Free(m[i]);
  Free(m);
}

static double **Lmatrix(int n)
{
  int     i;
  double **m;

  m = Calloc(n, double *);
  if (!m) errmsg("fail1 in Lmatrix()");

  for (i = 0; i < n; i++) {
    m[i] = Calloc((i + 1), double);
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

#ifdef DataByCopy
static Sdata *svect(int n)
{
  Sdata  *v;

  v = Calloc(n, Sdata);
  if (!v) errmsg("allocation failure in svect()");
  return v;
}

static void free_svect(Sdata *v)
{
   Free(v);
}
#endif

void
VR_dovm(Sint *Nw, double *wts, double *Fmin,
	Sint *maxit, Sint *trace, Sint *mask,
	double *abstol, double *reltol)
{
  vmmin((int) *Nw, wts, Fmin, (int) *maxit, (int) *trace, mask,
	*abstol, *reltol);
}


typedef unsigned char Boolean;

#define false 0

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0
#define REPORT		10


/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */

static void
vmmin(int n0, double *b, double *Fmin, int maxit, int trace, Sint *mask, 
      double abstol, double reltol)
{
  Boolean accpoint, enough;
  double  *g, *t, *X, *c, **B;
  int     count, funcount, gradcount;
  double  f, gradproj;
  int     i, j, ilast, iter = 0;
  double  s, steplength;
  double  D1, D2;
  int     n, *l;

  l = Calloc(n0, int);
  if (!l) errmsg("allocation failure in dovm()");
  n = 0;
  for (i = 0; i < n0; i++)
    if (mask[i]) l[n++] = i;

  g = vect(n0);
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
	X[i] = b[l[i]];
	c[i] = g[l[i]];
      }
      gradproj = 0.0;
      for (i = 0; i < n; i++) {
	s = 0.0;
	for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	t[i] = s;
	gradproj += s * g[l[i]];
      }

      if (gradproj < 0.0) {	/* search direction is downhill */
	steplength = 1.0;
	accpoint = false;
	do {
	  count = 0;
	  for (i = 0; i < n; i++) {
	    b[l[i]] = X[i] + steplength * t[i];
	    if (reltest + X[i] == reltest + b[l[i]])	/* no change */
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
	    c[i] = g[l[i]] - c[i];
	    D1 += t[i] * c[i];
	  }
	  if (D1 > 0) {
	    D2 = 0.0;
	    for (i = 0; i < n; i++) {
	      s = 0.0;
	      for (j = 0; j <= i; j++) s += B[i][j] * c[j];
	      for (j = i + 1; j < n; j++) s += B[j][i] * c[j];
	      X[i] = s;
	      D2 += s * c[i];
	    }
	    D2 = 1.0 + D2 / D1;
	    for (i = 0; i < n; i++) {
	      for (j = 0; j <= i; j++)
		B[i][j] += (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
	    }
	  } else {		/* D1 < 0 */
	    ilast = gradcount;
	  }
	} else {		/* no progress */
	  if (ilast < gradcount) {
	    count = 0;
	    ilast = gradcount;
	  }
	}
      } else {			/* uphill search */
	count = 0;
	if (ilast == gradcount) count = n;
	else ilast = gradcount;
	/* Resets unless has just been reset */
      }
      if (iter % REPORT == 0 && trace) {
	printf("iter%4d value %f\n", iter, f);
	fflush(stdout);
      }
      if (iter >= maxit) break;
      if (gradcount - ilast > 2 * n)
	ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
  }
  if (trace) {
    printf("final  value %f \n", *Fmin);
    if (iter < maxit) printf("converged\n");
    else printf("stopped after %i iterations\n", iter);
  }
  free_vect(g);
  free_vect(t);
  free_vect(X);
  free_vect(c);
  free_Lmatrix(B, n);
  Free(l);
}


static double **H, *h, *h1, **w;

static void
pHessian(Sdata *input, Sdata *goal, Sdata wx, int nr)
{
  int     i, to1, to2, from1, from2, j, j1, j2, first1, first2;
  double  out, s, sum1, sum2, t, tmp, tot, P;

  fpass(input, goal, 1.0, nr);
  bpass(goal, 1.0);

  /* Formulae from Ripley (1996, p.152) */
  if (Softmax) {
    for (i = 0; i < Nunits; i++) {
      sum1 = 0.0; sum2 = 0.0; tot = 0.0; P = 0.0;
      for (j = FirstOutput; j < Nunits; j++) {
	sum1 += w[i][j] * Probs[j];
	t = goal[j - FirstOutput];
	P += t * Probs[j];
	sum2 += w[i][j] * Probs[j] * t;
	tot += t;
      }
      h[i] = sum1; h1[i] = sum2/P;
    }
    if(Censored) tot = 1;
    for (to1 = 0; to1 < Nunits; to1++)
      for (j1 = Nconn[to1]; j1 < Nconn[to1 + 1]; j1++) {
	from1 = Conn[j1];
	first1 = (to1 < FirstOutput);
	for (to2 = 0; to2 < Nunits; to2++)
	  for (j2 = Nconn[to2]; j2 < Nconn[to2 + 1]; j2++)
	    if (j2 <= j1) {
	      from2 = Conn[j2];
	      first2 = (to2 < FirstOutput);
	      if ((!first1) && (!first2)) {  /* both -> output */
		if(Censored) {
		  tmp = -Probs[to1] * Probs[to2] * 
		    (1 - goal[to1 - FirstOutput] * goal[to2 - FirstOutput]/P/P);
		  if (to1 == to2) tmp += Probs[to1] * (1 - goal[to1 - FirstOutput]/P);
		  H[j1][j2] += wx * (tmp * Outputs[from1] * Outputs[from2]);	  
		} else {
		    tmp = -Probs[to1] * Probs[to2];
		    if (to1 == to2) tmp += Probs[to1];
		    H[j1][j2] += wx * tot * (tmp * Outputs[from1] * Outputs[from2]);
		}
	      } else if (first1 && first2) {  /* both -> hidden */
		sum1 = sum2 = 0.0;
		for (i = FirstOutput; i < Nunits; i++) {
		  sum1 += Errors[i] * w[to1][i];
		  tmp = w[to1][i] * w[to2][i] * Probs[i];
		  if(Censored) tmp *= (1 - goal[i - FirstOutput]/P);
		  sum2 += tmp;
		}
		if(Censored) {
		  sum2 += - h[to1] * h[to2] + h1[to1] * h1[to2];
		  s = sigmoid_prime(Outputs[to1]) * sigmoid_prime(Outputs[to2])
		    * sum2;
		} else  {
		  sum2 -= h[to1] * h[to2];
		  s = sigmoid_prime(Outputs[to1]) * sigmoid_prime(Outputs[to2])
		    * tot * sum2;
		}
		if (to1 == to2)
		  s += sigmoid_prime_prime(Outputs[to1]) * sum1;
		H[j1][j2] += wx * (s * Outputs[from1] * Outputs[from2]);
	      } else {  /* one -> hidden, one -> output */
		if (to1 < to2) {
		  tmp = w[to1][to2] - h[to1];
		  if(Censored) tmp += goal[to2 - FirstOutput]/P *
				 (h1[to1] - w[to1][to2]);
		  H[j1][j2] += wx * (Outputs[from1] * sigmoid_prime(Outputs[to1])
				     * (Outputs[from2] * Probs[to2] * tmp * tot
					+ ((to1 == from2) ? Errors[to2] : 0)));
		} else {
		  tmp = w[to2][to1] - h[to2];
		  if(Censored) tmp += goal[to1 - FirstOutput]/P *
				 (h1[to2] - w[to2][to1]);
		  H[j1][j2] += wx * (Outputs[from2] * sigmoid_prime(Outputs[to2])
				     * (Outputs[from1] * Probs[to1] * tmp * tot
					+ ((to2 == from1) ? Errors[to1] : 0)));
		}
	      }
	    }
      }
  } else {  /* Not softmax */
    for (i = FirstOutput; i < Nunits; i++) {
      out = Outputs[i];
      s = sigmoid_prime(out);
      t = goal[i - FirstOutput];
      if (Linout) h[i] = 2;
      else if (Entropy) h[i] = out * (1 - out);
      else h[i] = sigmoid_prime_prime(out) * 2 * (out - t) + 2 * s * s;
    }
    for (to1 = 0; to1 < Nunits; to1++)
      for (j1 = Nconn[to1]; j1 < Nconn[to1 + 1]; j1++) {
	from1 = Conn[j1];
	first1 = (to1 < FirstOutput);
	for (to2 = 0; to2 < Nunits; to2++)
	  for (j2 = Nconn[to2]; j2 < Nconn[to2 + 1]; j2++)
	    if (j2 <= j1) {
	      from2 = Conn[j2];
	      first2 = (to2 < FirstOutput);
	      if ((!first1) && (!first2)) {  /* both -> output */
		if (to1 == to2)
		  H[j1][j2] += wx * (h[to1] * Outputs[from1] * Outputs[from2]);
	      }
	      else if (first1 && first2) {  /* both -> hidden */
		sum1 = sum2 = 0.0;
		for (i = FirstOutput; i < Nunits; i++) {
		  sum1 += Errors[i] * w[to1][i];
		  sum2 += w[to1][i] * w[to2][i] * h[i];
		}
		s = sigmoid_prime(Outputs[to1]) * sigmoid_prime(Outputs[to2])
		  * sum2;
		if (to1 == to2)
		  s += sigmoid_prime_prime(Outputs[to1]) * sum1;
		H[j1][j2] += wx * (s * Outputs[from1] * Outputs[from2]);
	      }
	      else {  /* one -> hidden, one -> output */
		if (to1 < to2) {
		  H[j1][j2] += wx * (Outputs[from1] * sigmoid_prime(Outputs[to1])
				     * (Outputs[from2] * w[to1][to2] * h[to2]
					+ ((to1 == from2) ? Errors[to2] : 0)));
		}
		else {
		  H[j1][j2] += wx * (Outputs[from2] * sigmoid_prime(Outputs[to2])
				     * (Outputs[from1] * w[to2][to1] * h[to1]
					+ ((to2 == from1) ? Errors[to1] : 0)));
		}
	      }
	    }
      }
  }
}

#define max9(a,b) a>b?a:b
#define min9(a,b) a<b?a:b

void
VR_nnHessian(double *inwts, Sdata *Hess)
{
  int     i, j;

  for (i = 0; i < Nweights; i++) wts[i] = inwts[i];
  H = Lmatrix(Nweights);
  h = vect(Nunits); h1 = vect(Nunits);
  w = matrix(Nunits, Nunits);
  for (i = 0; i < Nweights; i++)
    for (j = 0; j <= i; j++) H[i][j] = 0.0;
  for (j = FirstOutput; j < Nunits; j++)
    for (i = FirstHidden; i < FirstOutput; i++) w[i][j] = 0.0;
  for (j = FirstOutput; j < Nunits; j++)
    for (i = Nconn[j]; i < Nconn[j + 1]; i++) w[Conn[i]][j] = wts[i];
  for (i = 0; i < NTrain; i++) {
    for (j = 0; j < Noutputs; j++) toutputs[j] = TrainOut[i + NTrain*j];
    pHessian(TrainIn + i, toutputs, Weights[i], NTrain);
  }
  for (i = 0; i < Nweights; i++) H[i][i] += 2 * Decay[i];

  for (i = 0; i < Nweights; i++)
    for (j = 0; j < Nweights; j++) *Hess++ = H[max9(i, j)][min9(i, j)];
  free_Lmatrix(H, Nweights);
  free_vect(h);free_vect(h1);
  free_matrix(w, Nunits, Nunits);
}
static int p, q;

static int
Zcompar(const Sdata *a, const Sdata *b)
{
  int     i;
  for (i = 0; i < p; i++)
    if (a[i] != b[i]) return (2 * (a[i] > b[i]) - 1);
  return (0);
}

/* Z is transposed, so (p+q) x n */

void
VR_summ2(Sint *n0, Sint *p0, Sint *q0, Sdata *Z, Sint *na)
{
  int     n = *n0, m;
  int     i, j, k, l;
  p = *p0;
  q = *q0;
  m = p + q;
  qsort(Z, n, m * sizeof(Sdata), Zcompar);
  j = 0;
  for (i = 1; i < n; i++) {
    k = -1;
    for (l = 0; l < p; l++)
      if (Z[l + i * m] != Z[l + (i - 1) * m]) {
	k = l;
	break;
      }
    if (k >= 0) {
      j++;
      for (l = 0; l < m; l++) Z[l + j * m] = Z[l + i * m];
    }
    else
      for (l = p; l < m; l++) Z[l + j * m] += Z[l + i * m];
  }
  *na = j + 1;
}

/* find maximum column: designed for probabilities. Uses
   reservoir sampling to break ties at random */

void    VR_max_col(double *matrix, int *pnr, int *nc, int *maxes)
{
  int     r, c, m, nr = *pnr, ntie;
  double  a, b;

  S_EVALUATOR
  RANDIN;
  for (r = 0; r < nr; r++) {
    m = 0;
    ntie = 1;
    a = matrix[r];
    for (c = 1; c < *nc; c++) {
      if ((b = matrix[r + c * nr]) >= a) {
	if (fabs(a - b) < 1e-5) {
	  ntie++;
	  if (ntie * UNIF < 1.0) m = c;
	} else {
	  ntie = 1;
	  a = b;
	  m = c;
	}
      }
    }
    maxes[r] = m + 1;
  }
  RANDOUT;
}
