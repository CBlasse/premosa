/*****************************************************************************************\
*                                                                                         *
*  Distribution generator data abstraction                                                *
*     One can create a distribution generator for a number of parameterized distribution  *
*     types and then generate events with that distribution                               *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  January 2007                                                                  *
*                                                                                         *
*  (c) June 19, '09, Dr. Gene Myers and Howard Hughes Medical Institute                   *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mylib.h"
#include "utilities.h"
#include "cdf.h"

#undef DEBUG

#define SIZEOF(x) ((int) sizeof(x))

/* Windows doesn't have drand48 in its native library and drand48 isn't thread-safe.
   So we built our own version that generates *exactly* the same value sequence given
   any specific seed state.
*/

static pthread_mutex_t Drand_Mutex = PTHREAD_MUTEX_INITIALIZER;

static uint64 myrand48_a = 0x5deece66dull;
static uint64 myrand48_c = 0xbull;

double erand(uint64 *state)
{ uint64  temp;
  double *tempd = ((double *) (&temp));

  *state = ((*state) * myrand48_a + myrand48_c) & 0xffffffffffffull;

  temp = 0x3ff0000000000000ull | (*state << 4);

  return (*tempd - 1.0);
}

static uint64 myrand_state = 0x1234abcd330eull;

double drand()
{ double rez;
  pthread_mutex_lock(&Drand_Mutex);
  rez = erand(&myrand_state);
  pthread_mutex_unlock(&Drand_Mutex);
  return (rez);
}

uint64 eseed(uint32 seedval)
{ return ((((uint64) seedval) << 16) | 0x330e); }

void dseed(uint32 seedval)
{ pthread_mutex_lock(&Drand_Mutex);
  myrand_state = eseed(seedval);
  pthread_mutex_unlock(&Drand_Mutex);
}

typedef struct
  { int    rcnt;
    uint64 state;
  } Generator;

static Generator *newgen(string routine)
{ Generator *gen = (Generator *) Guarded_Malloc(sizeof(Generator),routine);
  gen->rcnt  = 1;
  gen->state = 0x1234abcd330eull;
  return (gen);
}

/* BASIC UTILITIES:
     Declaration of internal implementation, memory manager, and bin_search
*/

typedef enum { NORMAL, BINOMIAL, POISSON, BERNOUILLI,
               GEOMETRIC, EXPONENTIAL, FAIRCOIN, UNIFORM } ctype;

typedef struct
  { ctype         kind;    /* type of distribution                                             */
    Generator    *grand;   /* CDF's random number generator                                    */
    double        parm1;   /* mean if NORMAL, p if GEOMETRIC, a if EXPONENTIAL, low if UNIFORM */
    double        parm2;   /* stdev if NORMAL, hgh-low if UNIFORM                              */
    int           parm3;   /* n if FAIRCOIN                                                    */
                           /* For BINOMIAL, POISSON, and BERNOUILLI:                           */
    int           low;     /*   Lowest significant domain element                              */
    int           spn;     /*   Span of significant domain elements                            */
    double       *tab;     /*   Explicit table of cdf of significant els                       */
  } Cdf;

static inline int cdf_tsize(Cdf *cdf)
{ return (SIZEOF(double)*(cdf->spn+1)); }

MANAGER -iofkc CDF(Cdf) tab:tsize

CDF *Copy_CDF(CDF *cdf)
{ Cdf *c = copy_cdf((Cdf *) cdf);
  c->grand = newgen("Copy_CDF");
  c->grand->state = ((Cdf *) cdf)->grand->state;
  return ((CDF *) c);
}

void Free_CDF(CDF *cdf)
{ Cdf *c = (Cdf *) cdf;
  if (c->grand->rcnt > 1)
    c->grand->rcnt -= 1;
  else
    free(c->grand);
  free_cdf(c);
}

void Kill_CDF(CDF *cdf)
{ Cdf *c = (Cdf *) cdf;
  if (c->grand->rcnt > 1)
    c->grand->rcnt -= 1;
  else
    free(c->grand);
  kill_cdf(c);
}

static void init_unorm();

CDF *Read_CDF(FILE *input)
{ Cdf *cdf = read_cdf(input);
  if (cdf != NULL)
    { cdf->grand = newgen("Read_CDF");
      cdf->grand->rcnt = 1;
      fread(&(cdf->grand->state),sizeof(uint64),1,input);
      if (cdf->kind == NORMAL)
        init_unorm();
    }
  return ((CDF *) cdf);
}

void Write_CDF(CDF *cdf, FILE *output)
{ Cdf *c = (Cdf *) cdf;
  write_cdf(c,output);
  fwrite(&(c->grand->state),sizeof(uint64),1,output);
}

static int bin_search(int len, double *tab, double y)
{ int l, m, r;

/* Searches tab[0..len] for min { r : y < tab[r] }.
       Assumes y < 1, tab[0] = 0 and tab[len] = 1.
       So returned index is in [1,len].               */

  l = 0;
  r = len;
  while (l < r)
    { m = (l+r) >> 1;
      if (y < tab[m])
        r = m; 
      else
        l = m+1;
    }
  return (r);
}

/* NORMAL DISTRIBUTION ROUTINES:
     void   init_unorm();
     double sample_unorm();
     CDF   *Normal_CDF(double mean, double stdev)
*/

#define UNORM_LEN 60000
#define UNORM_MAX   6.0

static double unorm_table[UNORM_LEN+1];  /* Upper half of cdf of N(0,1) */
static double unorm_scale;
static int    unorm_defined = 0;

static pthread_mutex_t Normal_Mutex = PTHREAD_MUTEX_INITIALIZER;

static void init_unorm()
{ double del, sum, x;
  int    i;

  pthread_mutex_lock(&Normal_Mutex);
  if (unorm_defined)
    { pthread_mutex_unlock(&Normal_Mutex);
      return;
    }
  unorm_defined = 1;

  unorm_scale = del = UNORM_MAX / UNORM_LEN;

		/* Integrate pdf, x >= 0 half only. */
  sum = 0;
  for (i = 0; i < UNORM_LEN; i++)
    { x = i * del;
      unorm_table[i] = sum;
      sum += exp(-.5*x*x) * del; 
    }
  unorm_table[UNORM_LEN] = sum;

		/* Normalize cdf */
  sum *= 2.;
  for (i = 0; i < UNORM_LEN; i++)
    unorm_table[i] /= sum;
  unorm_table[UNORM_LEN] = 1.;

#ifdef DEBUG
  printf("Truncated tail is < %g\n",
          exp(-.5*UNORM_MAX*UNORM_MAX)/(sum*(1.-exp(-UNORM_MAX))) );
  printf("Diff between last two entries is %g\n",.5-unorm_table[UNORM_LEN-1]);

  printf("\n  CDF:\n");
  for (i = 0; i <= UNORM_LEN; i += 100)
    printf("%6.2f: %10.9f\n",i*del,unorm_table[i]);
#endif

  pthread_mutex_unlock(&Normal_Mutex);
}

static double sample_unorm(double x)
{ double y;
  int    f;

		/* Map [0,1) random var to upper-half of cdf */
  if (x >= .5)
    y = x-.5;
  else
    y = .5-x;
		/* Bin. search upper-half cdf */
  f = bin_search(UNORM_LEN,unorm_table,y);
#ifdef DEBUG
  printf("Normal search %g -> %g -> %d",x,y,f); 
#endif

		/* Linear interpolate between table points */
  y = (f - (unorm_table[f]-y) / (unorm_table[f] - unorm_table[f-1]) ) * unorm_scale; 

		/* Map upper-half var back to full range */
  if (x < .5) y = -y;
#ifdef DEBUG
  printf(" -> %g\n",y);
#endif

  return (y);
}

CDF *G(Normal_CDF)(double mean, double stdev)
{ Cdf *bin;

  bin = new_cdf(SIZEOF(double),"Normal_CDF");
  bin->kind  = NORMAL;
  bin->grand = newgen("Normal_CDF");
  bin->parm1 = mean;
  bin->parm2 = stdev;
  bin->spn   = 0;
  init_unorm();
  return ((CDF *) bin);
}


/* BINOMIAL DISTRIBUTION ROUTINE:
      CDF *Binomial_CDF(int n, double p);  n >= 1 & p in [0,1]
      Table uses 400 * 3.2 ^ log_10(n) bytes
*/

CDF *G(Binomial_CDF)(int n, double p)
{ double *tab;
  double  pek, nxt, sum;
  double  pm1, var;
  int     i, k, c;
  int     low, hgh, spn;
  Cdf    *bin;

		/* Compute p-value at k = pn */
  pm1 = 1.-p;
  var = p*pm1;

  k = (int) (p*n);
  if (p <= .5)
    { pek = 1.;
      c = i = k-1;
      while (i >= 2*k-n || c >= 0)
        { if (pek < n && c >= 0)
            { pek *= (n-c) / (c+1.);
              c -= 1;
            }
          else if (i >= 0)
            { pek *= var;
              i -= 1;
            }
          else
            { pek *= pm1;
              i -= 1;
            }
        }
    }
  else
    { if (k < n) k += 1;
      pek = 1.;
      c = i = n-k-1;
      while (i >= n-2*k || c >= 0)
        { if (pek < n && c >= 0)
            { pek *= (n-c) / (c+1.);
              c -= 1;
            }
          else if (i >= 0)
            { pek *= var;
              i -= 1;
            }
          else
            { pek *= p;
              i -= 1;
            }
        }
    }

#ifdef DEBUG
  printf("Binomial(%d,%d):\n",n,p);
  printf("  Peak at %d : %g\n",k,pek);
#endif

		/* Compute length of non-zero tails on either side of peak */
  nxt = pek;
  for (i = k-1; i >= 0; i--)
    { nxt *= ((i+1)*pm1)/((n-i)*p);  
      if (nxt < 1e-50) break;
    }
  low = i;

  nxt = pek;
  for (i = k+1; i <= n; i++)
    { nxt *= (((n-i)+1)*p)/(i*pm1);  
      if (nxt < 1e-50) break;
    }
  hgh = i-1;
  spn = hgh - low;

#ifdef DEBUG
  printf("  Span is [%d,%d] of length %d\n",low,hgh,spn+1);
#endif

		/* Allocate binomial data structure */

  bin = new_cdf((spn+1)*SIZEOF(double),"Binomial_CDF");
  bin->kind  = BINOMIAL;
  bin->grand = newgen("Binomial_CDF");
  bin->low   = low;
  bin->spn   = spn;
  tab        = bin->tab;

		/* Compute tails again, this time storing them */
  tab[k-low] += pek;

  nxt = pek;
  for (i = k-1; i > low; i--)
    tab[i-low] = nxt *= ((i+1)*pm1)/((n-i)*p);  

  nxt = pek;
  for (i = k+1; i <= hgh; i++)
    tab[i-low] = nxt *= (((n-i)+1)*p)/(i*pm1);  

#ifdef DEBUG
  if (low >= 0)
    { sum = ((low+1)*pm1)/((n-low)*p);
      printf("  Low tail truncated < %g (r = %g)\n",tab[1]/(1.-sum),sum);
    }
  if (hgh < n)
    { sum = ((n-hgh)*p)/((hgh+1)*pm1);
      printf("  Hgh tail truncated < %g (r = %g)\n",tab[spn]/(1.-sum),sum);
    }
  printf("  PDF:\n");
  for (i = 1; i <= spn; i++)
    printf("    %4d: %10.9f\n",low+i,tab[i]);
#endif

		/* Compute cdf and normalize sum */
  sum = 0.;
  for (i = 0; i < spn; i++)
    { tab[i] = sum;
      sum += tab[i+1];
    }
  for (i = 1; i < spn; i++)
    tab[i] /= sum;
  tab[spn] = 1.;

#ifdef DEBUG
  printf("  CDF:\n");
  for (i = 0; i <= spn; i++)
    printf("    %4d: %10.9f\n",low+i,tab[i]);
#endif
 
  return ((CDF *) bin);
}

/* POISSON DISTRIBUTION ROUTINES:
      CDF *Poisson_CDF(double a);  a > 0
      Table uses 960 * 3 ^ log_10(a) bytes
*/

CDF *G(Poisson_CDF)(double a)
{ double  *tab;
  double   pek, nxt, sum;
  double   nap;
  int      i, k, c;
  int      low, hgh, spn;
  Cdf     *pois;

  nap = exp(1.);
  c = i = k = (int) a;
  pek = 1.;
  while (i > 0 || c > 0)
    if (pek < a && c > 0)
      { pek *= a;
        c -= 1;
      }
    else
      { pek /= (nap*i);
        i -= 1;
      }

#ifdef DEBUG
  printf("Poisson(%g):\n",a);
  printf("  Peak at %d : %g\n",k,pek);
#endif

		/* Compute length of non-zero tails on either side of peak */
  nxt = pek;
  for (i = k-1; i >= 0; i--)
    { nxt *= (i+1.)/a;  
      if (nxt < 1e-50) break;
    }
  low = i;

  nxt = pek;
  for (i = k+1; 1; i++)
    { nxt *= a/i;  
      if (nxt < 1e-50) break;
    }
  hgh = i-1;
  spn = hgh - low;

#ifdef DEBUG
  printf("  Span is [%d,%d] of length %d\n",low,hgh,spn+1);
#endif

		/* Allocate Poisson data structure */

  pois = new_cdf((spn+1)*SIZEOF(double),"Poisson_CDF");
  pois->kind  = POISSON;
  pois->grand = newgen("Poisson_CDF");
  pois->low   = low;
  pois->spn   = spn;
  tab         = pois->tab;

		/* Compute tails again, this time storing them */
  tab[k-low] = pek;

  nxt = pek;
  for (i = k-1; i > low; i--)
    tab[i-low] = nxt *= (i+1.)/a;  

  nxt = pek;
  for (i = k+1; i <= hgh; i++)
    tab[i-low] = nxt *= a/i;  

#ifdef DEBUG
  if (low >= 0)
    { sum = (low+1)/a;
      printf("  Low tail truncated < %g (r = %g)\n",tab[1]/(1.-sum),sum);
    }
  sum = a/(hgh+1.);
  printf("  Hgh tail truncated < %g (r = %g)\n",tab[spn]/(1.-sum),sum);
  printf("  PDF:\n");
  for (i = 1; i <= spn; i++)
    printf("    %4d: %10.9f\n",low+i,tab[i]);
#endif

		/* Compute cdf and normalize sum */
  sum = 0.;
  for (i = 0; i < spn; i++)
    { tab[i] = sum;
      sum += tab[i+1];
    }
  for (i = 1; i < spn; i++)
    tab[i] /= sum;
  tab[spn] = 1.;

#ifdef DEBUG
  printf("  CDF:\n");
  for (i = 0; i <= spn; i++)
    printf("    %4d: %10.9f\n",low+i,tab[i]);
#endif
 
  return ((CDF *) pois);
}

CDF *G(Bernouilli_CDF)(int n, double *weight)
{ int     i;
  double  sum;
  double *tab;
  Cdf    *bin;

  bin = new_cdf((n+1)*SIZEOF(double),"Bernouilli_CDF");
  bin->kind  = BERNOUILLI;
  bin->grand = newgen("Bernouilli_CDF");
  bin->low   = 0;
  bin->spn   = n;

  tab = bin->tab;
  sum = 0.;
  for (i = 0; i < n; i++)
    { tab[i] = sum;
      sum += weight[i];
    }

  for (i = 1; i < n; i++)
    tab[i] /= sum;
  tab[n] = 1.;

#ifdef DEBUG
  printf("  CDF:\n");
  for (i = 0; i <= n; i++)
    printf("    %4d: %10.9f\n",i,tab[i]);
#endif

  return ((CDF *) bin);
}

CDF *G(Geometric_CDF)(double p)
{ Cdf  *bin;

  bin = new_cdf(SIZEOF(double),"Geometric_CDF");
  bin->kind  = GEOMETRIC;
  bin->grand = newgen("Geometric_CDF");
  bin->parm1 = p;
  bin->spn   = 0;
  return ((CDF *) bin);
}

CDF *G(Exponential_CDF)(double a)
{ Cdf  *bin;

  bin = new_cdf(SIZEOF(double),"Exponential_CDF");
  bin->kind  = EXPONENTIAL;
  bin->grand = newgen("Exponential_CDF");
  bin->parm1 = a;
  bin->spn   = 0;
  return ((CDF *) bin);
}

CDF *G(Uniform_CDF)(double low, double hgh)
{ Cdf  *bin;

  bin = new_cdf(SIZEOF(double),"Uniform_CDF");
  bin->kind  = UNIFORM;
  bin->grand = newgen("Uniform_CDF");
  bin->parm1 = low;
  bin->parm2 = hgh - low;
  bin->spn   = 0;
  return ((CDF *) bin);
}

CDF *G(FairCoin_CDF)(int n)
{ Cdf  *bin;

  bin = new_cdf(SIZEOF(double),"FairCoin_CDF");
  bin->kind  = FAIRCOIN;
  bin->grand = newgen("FairCoin_CDF");
  bin->parm3 = n;
  bin->spn   = 0;
  return ((CDF *) bin);
}

void Seed_CDF(CDF *cdf, uint32 seedval)
{ ((Cdf *) cdf)->grand->state = eseed(seedval); } 

void Link_CDF(CDF *ref, CDF *sub)
{ Cdf *r = (Cdf *) ref;
  Cdf *s = (Cdf *) sub;
  if (s->grand->rcnt != 1)
    { fprintf(stderr,"CDF is already linked (Link_CDF)\n");
      exit (1);
    }
  free(s->grand);
  s->grand = r->grand;
  r->grand->rcnt += 1;
}

void Unlink_CDF(CDF *cdf)
{ Cdf       *c = (Cdf *) cdf;
  Generator *g = c->grand;
  if (g->rcnt > 1)
    { g->rcnt -= 1;
      c->grand = newgen("Unlink_CDF");
      c->grand->state = g->state;
    }
} 

double Sample_CDF(CDF *cdf)
{ int    f;
  double x;
  Cdf   *bin;

  bin = (Cdf *) cdf;

  x = erand(&(bin->grand->state));
  switch (bin->kind)
  { case FAIRCOIN:
      f = (int) (x * bin->parm3); 
      return (f+1.);
    case UNIFORM:
      return ( bin->parm1 + bin->parm2*x );
    case EXPONENTIAL:
      return ( - log(1.-x) / bin->parm1 );
    case GEOMETRIC:
      f = (int) (x = log(1.-x)/log(1.-bin->parm1));
      if (f < x) f += 1;
      x = f;
      return (x);
    case BINOMIAL:
    case BERNOUILLI:
    case POISSON:
      f = bin_search(bin->spn,bin->tab,x);
#ifdef DEBUG
      printf("CDF search %g --> %d\n",x,f+bin->low);
#endif
      return ((double) (bin->low + f));
    case NORMAL:
      return (bin->parm1 + bin->parm2 * sample_unorm(x));
  }
  return (0.);
}
