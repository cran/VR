#if( defined(SPLUS_VERSION) && SPLUS_VERSION > 5100 )
#  define RANDIN  seed_in((long *)NULL, S_evaluator)
#  define RANDOUT seed_out((long *)NULL, S_evaluator)
#  define UNIF unif_rand(S_evaluator)
#else
#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()
#  define Salloc(n, t) (t *)S_alloc(n, sizeof(t))
#  define S_EVALUATOR
#endif

#ifdef USING_R
  typedef double Sfloat;
  typedef int Sint;
# define SINT_MAX INT_MAX
#else
  typedef float Sfloat;
  typedef long Sint;
# define SINT_MAX LONG_MAX
#endif

