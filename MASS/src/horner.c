static double horner(double x, double *b, int n)
{
    int i;
    double p = b[n];
    for(i = n-1; i >= 0; i--)
	p = b[i] + x * p;
    return p;
}

void
poly(long int *m, double *p, double *x, long int *n, double *b)
{
    long i;
    for (i = 0; i < *m; i++)
	p[i] = horner(x[i], b, (int)*n);
}
