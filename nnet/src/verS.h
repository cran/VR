#if( defined(SPLUS_VERSION) && SPLUS_VERSION >= 5000 )
#  define RANDIN  seedin((long *)NULL, S_evaluator)
#  define RANDOUT seedout((long *)NULL, S_evaluator)
#  define UNIF unif_rand(S_evaluator)
#else
#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()
#  define Salloc(n, t) (t *)S_alloc(n, sizeof(t))
#  define S_EVALUATOR
#endif

#ifdef USING_R
typedef int Sint;
typedef double singl;
#else
typedef long Sint;
typedef float singl;
#endif

