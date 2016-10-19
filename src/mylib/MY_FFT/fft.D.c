/*****************************************************************************************\
*                                                                                         *
*  Double Precision FFT library                                                           *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  June 2009                                                                     *
*                                                                                         *
*  (c) June 19, '09, Dr. Gene Myers and Howard Hughes Medical Institute                   *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fft.D.h"

   //   Only needed if compiled independently of fft.F.c, define ALONE

  /*  Important optimizations:

      * sort+fft in the first log_4 n passes
      * radix 4 sort in all (but possibly one) of the remaining passes
      * use a table of roots of unity for all powers up to 4096
      * segment all phases with a span less than the size of the L2_CACHE
      * use mmx registers in the critical loops
      * size is limited as want to use a table for all the sort+fft passes
          ==> 8096 * 1024 = 8M max for complex arrays, 16M for real arrays

      For multi-dimensional arrays, we also:

      * cache each line on the fly doing the 1st sort+fft step as data is move to the cache.
      * currrently only wrote SSE3 code for the case were each dim can use the table of
           roots of unity ==> 4K max for any dimension, (8K for 1st dim of a real array)
  */

//  Critical limiting constants:

#define L2_CACHE     262144        //  Problems bigger than this get cross cut
#define NINETY         1024        //  Ninety degree ticks (problems > 4*NINETY are not tabled)
#define MAX_REAL_DIM    128        //  Maximum dimensionality of a real multi-dim. fft

#ifndef FFT_PROLOGUE

#define FFT_PROLOGUE

#define ONE_EIGHTY   (2*NINETY)

#define TPI 6.28318530717959       //  2*Pi

typedef long long int64; //  64-bit ints needed for vars over an multi-dimensional
                         //    array's size

static void print_index(int ndim, int *dims, int k)      //  print multi-d index
{ if (ndim > 1)
    { print_index(ndim-1,dims+1,k/dims[0]);
      printf(",");
    }
  printf("%d",k%dims[0]);
}

#ifdef _MSC_VER

#include <windows.h>

#pragma warning( disable:4996 )   // Turn off deprecation warnings on Windows

//  WINDOW pthreads

  //  Mutex macros

typedef SRWLOCK pthread_mutex_t;

#define PTHREAD_MUTEX_INITIALIZER RTL_SRWLOCK_INIT

#define pthread_mutex_lock(m)  AcquireSRWLockExclusive(m)

#define pthread_mutex_unlock(m) ReleaseSRWLockExclusive(m)

  //  Condition variable macros

typedef CONDITION_VARIABLE pthread_cond_t;

#define PTHREAD_COND_INITIALIZER RTL_CONDITION_VARIABLE_INIT

#define pthread_cond_signal(c) WakeConditionVariable(c)

#define pthread_cond_broadcast(c) WakeAllConditionVariable(c)

#define pthread_cond_wait(c,m) SleepConditionVariableSRW(c,m,INFINITE,0)

  //  Current thread ID

typedef int pthread_id;

static pthread_id pthread_tag()
{ return (GetCurrentThreadId()); }

static int pthread_is_this(pthread_id id)
{ return (GetCurrentThreadId() == id); }

#else   //  Small extension to pthreads!

#include <pthread.h>

typedef pthread_t pthread_id;

#define pthread_tag() pthread_self()

static inline int pthread_is_this(pthread_id id)
{ return (pthread_equal(pthread_self(),id)); }

#endif

  //  If an FFT routine encounters an exception it trys to grab the "error resource" with 
  //    grab_message and if successfull then sets the global strings FFT_Estring/Esource 
  //    before returning NULL.  A user can subsequently get these strings with FFT_Error_String()
  //    and FFT_Error_Source() and release the "error resource" with a call to FFT_Error_Release().

static pthread_mutex_t FFT_Err_Mutex = PTHREAD_MUTEX_INITIALIZER;

static int        FFT_Error = 0;
static char       FFT_Estring[75];
static char       FFT_Esource[25];
static pthread_id FFT_Ethread;

static int grab_message()
{ pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error) return (1);
  FFT_Error   = 1;
  FFT_Ethread = pthread_tag();
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (0);
}

static int my_message()
{ int mine;
  pthread_mutex_lock(&FFT_Err_Mutex);
  mine = (FFT_Error && pthread_is_this(FFT_Ethread));
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (mine);
}

/* CORINNA
char *FFT_Error_Source()
{ pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error && pthread_is_this(FFT_Ethread))
    return (FFT_Esource);
  else
    return (NULL);
  pthread_mutex_unlock(&FFT_Err_Mutex);
}

char *FFT_Error_String()
{ pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error && pthread_is_this(FFT_Ethread))
    return (FFT_Estring);
  else
    return (NULL);
  pthread_mutex_unlock(&FFT_Err_Mutex);
}

void FFT_Error_Release()
{ pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error && pthread_is_this(FFT_Ethread))
    FFT_Error = 0;
  pthread_mutex_unlock(&FFT_Err_Mutex);
}

  //  All the FFT routines assume n is a power of 2.  As a convenience to help you pad
  //    your vectors to this size, Next_Power_Of_2f returns the smallest power of 2 greater
  //    than or equal to m.

int Next_Power_Of_2(int m)
{ int n;

  n = 2;
  while (n < m)
    n <<= 1;
  return (n);
}
 CORINNA */

#endif  //  FFT_PROLOGUE

/*********************************************************************************************\
 *                                                                                            *
 *  PROLOGUE AND UTILITIES:                                                                   *
 *     sse3 macros, powers of 2 and error strings, print routines, roots of unity tables      *
 *                                                                                            *
\*********************************************************************************************/

  //  Includes and macros for intel intrinsics if present (-msse3 option)
  //     A "__m128d" contains 2 double values that are operated on by intrinsic functions.
  //     In most cases we encode a complex # a+ib in an __m128d as [a,b], we call this
  //     normal form.  For a twiddle factor x+iy, however, we occasionally encode them
  //     into two __m128d's, one containing [x,x] (the real part) and the other containing
  //     [y,y].  We call this twiddle form.  It is easy(er) to do complex multiply of a normal
  //     form times a twiddle form.

#ifdef __SSE3__

#include <pmmintrin.h>

#define LOADP(x)    _mm_load_pd((double *) (x))      //  load x[0], x[1]
#define LOADUP(x)   _mm_load1_pd((double *) (x))     //  load x[0], x[0]
#define STOREP(x,y) _mm_store_pd((double *) (x),y)   //  store m128d y into x[0],x[1]
#define SETP(r,i)   _mm_set_pd(i,r)                  //  load i, r

#define SHUFFLE(x)  _mm_shuffle_pd(x,x,1)            //  flip the quad words of m128d x
#define ADDP(x,y)   _mm_add_pd(x,y)                  //  add, sub, mul, div, xor quad words in ||
#define SUBP(x,y)   _mm_sub_pd(x,y)
#define MULP(x,y)   _mm_mul_pd(x,y)
#define XORP(x,y)   _mm_xor_pd(x,y)
#define ADDSUB(x,y) _mm_addsub_pd(x,y)               //  x0-y0, x1+y1

  //  Multiply d in normal form and (r,i) in twiddle form placing the result in p in normal form.

#define CMULTIPLY(p,d,r,i)              \
 p = MULP(r,d);				\
 p = ADDSUB(p,MULP(i,SHUFFLE(d)));

  //  Multiply a in normal form and (r,i) in twiddle form placing the result back in (r,i).

#define NEXT_ROOT(a,r,i)        	\
  CMULTIPLY(r,a,r,i); 		      	\
  i = _mm_unpackhi_pd(r,r);     	\
  r = _mm_movedup_pd(r);

  //   Handy __m128d constants

static __m128d Conjugate_d = { +0., -0. };   //  Used for conjugating (xor)
static __m128d FlipReal_d  = { -0., +0. };   //  Used for anti-conjugating (xor)
static __m128d Negate_d    = { -0., -0. };   //  Used for negating (xor)
static __m128d Half_d      = {  .5,  .5 };   //  Used for the real-valued fft's
static __m128d Unity_d     = {  1.,  0. };   //  Unity

#endif

  //  Print complex and real-valued 1d and multi-d arrays with a title

void Print_Complex_1d(int n, Complex_d *array, char *title)
{ int i;

  printf("\n%s:\n",title);
  for (i = 0; i < n; i++)
    { printf("   %4d",i);
      if (fabs(array[i].real) < 1e-12)
        printf(": 0 + i ");
      else
        printf(": %.5g + i ",array[i].real);
      if (fabs(array[i].imag) < 1e-12)
        printf("0\n");
      else
        printf("%.5g\n",array[i].imag);
    }
}

void Print_Real_1d(int n, double *array, char *title)
{ int i;

  printf("\n%s:\n",title);
  for (i = 0; i < n; i++)
    { printf("   %4d",i);
      if (fabs(array[i]) < 1e-12)
        printf(": 0\n");
      else
        printf(": %.5g\n",array[i]);
    }
}

void Print_Complex_nd(int ndim, int *dims, Complex_d *array, char *title)
{ int i, d;

  for (d = 1; d < ndim; d++)
    dims[d] *= dims[d-1];

  printf("\n%s:\n",title);
  for (i = 0; i < dims[ndim-1]; i++)
    { printf("   ");
      print_index(ndim,dims,i);
      if (fabs(array[i].real) < 1e-13)
        printf(": 0 + i ");
      else
        printf(": %.5g + i ",array[i].real);
      if (fabs(array[i].imag) < 1e-13)
        printf("0\n");
      else
        printf("%.5g\n",array[i].imag);
    }

  for (d = ndim-1; d > 0; d--)
    dims[d] /= dims[d-1];
}

void Print_Real_nd(int ndim, int *dims, double *array, char *title)
{ int i, d;

  for (d = 1; d < ndim; d++)
    dims[d] *= dims[d-1];

  printf("\n%s:\n",title);
  for (i = 0; i < dims[ndim-1]; i++)
    { printf("   ");
      print_index(ndim,dims,i);
      if (fabs(array[i]) < 1e-13)
        printf(": 0\n");
      else
        printf(": %.5g\n",array[i]);
    }

  for (d = ndim-1; d > 0; d--)
    dims[d] /= dims[d-1];
}

  //  Table of roots of unity Root[k] = (e^(Aki),e^(2Aki),e^(3Aki)) for A = TPI/(4*NINETY)
  //  Coot is the conjugate of Root (useful for inversion).

typedef struct
  { Complex_d pow1;
    Complex_d pow2;
    Complex_d pow3;
  } Trig_Block_d;

static Trig_Block_d Root_d[NINETY];
static Trig_Block_d Coot_d[NINETY];
static int          Trig_Firstime_d = 1;

static pthread_mutex_t FFT_D_MUTEX = PTHREAD_MUTEX_INITIALIZER;

static void init_trig_tables_d()
{ int    i;
  double inc, ang, mul;

  inc = TPI / (4*NINETY);
  ang = 0.;
  for (i = 0; i < NINETY; i++)
    { mul = ang;
      Coot_d[i].pow1.real = Root_d[i].pow1.real = cos(mul);
      Coot_d[i].pow1.imag = - (Root_d[i].pow1.imag = sin(mul));
      mul += ang;
      Coot_d[i].pow2.real = Root_d[i].pow2.real = cos(mul);
      Coot_d[i].pow2.imag = - (Root_d[i].pow2.imag = sin(mul));
      mul += ang;
      Coot_d[i].pow3.real = Root_d[i].pow3.real = cos(mul);
      Coot_d[i].pow3.imag = - (Root_d[i].pow3.imag = sin(mul));
      ang += inc;
    }
}


/*********************************************************************************************\
 *                                                                                            *
 *  1-D COMPLEX FFT:                                                                          *
 *       radix2[_table]_d, radix4[_table]_d, first2_table_d, sort2_table_d, argcheck_1d       *
 *                                                                                            *
\*********************************************************************************************/

  //  Radix 2 fft of span s (using twiddle table):  h = 2s and h | n.
  //    For all j s.t. j mod h < s
  //        D[0] = D[0] + u*D[1]
  //        D[1] = D[0] - u*D[1]
  //    where
  //      D[x] = data[j+x*s], u = d(w)^(j mod s), w_h = e^(-TPI/h)i
  //          and d(x) = x  if table = Root_d
  //                   = x* if table = Coot_d

static void radix2_table_d(Complex_d *data, int n, int s, int h, Trig_Block_d *table)
{ int a, v;

  a = NINETY/s;
  for (v = 0; v < n; v += h)
    { Complex_d *j0 = data+v;
      Complex_d *j1 = j0+s;
      Complex_d *je = j1;

      Trig_Block_d *t = table;

      while (j0 != je)
#ifdef __SSE3__
        { __m128d p, d;
          __m128d ui, ur;

          ur = LOADUP(&(t->pow2.real));
          ui = LOADUP(&(t->pow2.imag));

          d = LOADP(j1);
          CMULTIPLY(p,d,ur,ui);

          d = LOADP(j0);
          STOREP(j1,SUBP(d,p));
          STOREP(j0,ADDP(d,p));

#else
        { double dr, di;
          double pr, pi;
          double ur, ui;

          ur = t->pow2.real;
          ui = t->pow2.imag;

          dr = j1->real;
          di = j1->imag;

          pr = ur * dr - ui * di;
          pi = ur * di + ui * dr;

          dr = j0->real;
          di = j0->imag;

          j1->real = dr - pr;
          j1->imag = di - pi;
          j0->real = dr + pr;
          j0->imag = di + pi;
#endif

          j0 += 1;
          j1 += 1;
          t  += a;
        }
    }
}

  //  Radix 2 fft of span s (not using twiddle table):  h = 2s and h | n.

static void radix2_d(Complex_d *data, int n, int s, int h, double direct)
{ int    v;
  double theta = direct/h;

#ifdef __SSE3__
  __m128d ure = SETP(1.,1.);
  __m128d uri = SETP(0.,0.);

  __m128d ang = SETP(cos(theta),sin(theta));
#else
  double cos0  = cos(theta);
  double sin0  = sin(theta);
#endif

  for (v = 0; v < n; v += h)
    { Complex_d *j0 = data+v;
      Complex_d *j1 = j0+s;
      Complex_d *je = j1;

#ifdef __SSE3__
      __m128d    ur, ui;

      ur = ure;    // ur,ui holds w_h^u in twiddle form after iteration u of the loop
      ui = uri;

      while (j0 != je)
        { __m128d d, p;

          d = LOADP(j1);
          CMULTIPLY(p,d,ur,ui);

          d = LOADP(j0);
          STOREP(j1,SUBP(d,p));
          STOREP(j0,ADDP(d,p));

          j0 += 1;
          j1 += 1;

          NEXT_ROOT(ang,ur,ui);
        }
#else
      double ur, ui;

      ur = 1.0;      // ur + i*ui holds the w_h^u after iteration u of the loop
      ui = 0.0;

      while (j0 != je)
        { double dr, di;
          double pr, pi;

          dr = j1->real;
          di = j1->imag;

          pr = ur * dr - ui * di;
          pi = ur * di + ui * dr;

          dr = j0->real;
          di = j0->imag;

          j1->real = dr - pr;
          j1->imag = di - pi;
          j0->real = dr + pr;
          j0->imag = di + pi;

          j0 += 1;
          j1 += 1;

          pr = ur * cos0 - ui * sin0;
          ui = ur * sin0 + ui * cos0;
          ur = pr;
        }
#endif
    }
}

  //  Radix 4 fft of span s (using twiddle table):  h = 4s and h | n.
  //    For all j s.t. j mod h < s
  //        D[0] = D[0] +   u*D[2] + u^2*D[1] +   u^3*D[3]
  //        D[1] = D[0] + p*u*D[2] - u^2*D[1] - p*u^3*D[3]
  //        D[2] = D[0] -   u*D[2] + u^2*D[1] -   u^3*D[3]
  //        D[3] = D[0] - p*u*D[2] - u^2*D[1] + p*u^3*D[3]
  //    where
  //      D[x] = data[j+x*s], u = d(w)^(j mod s), w = e^(-TPI/h)i, p = d(i)
  //          and d(x) = x  if table = Root_d
  //                   = x* if table = Coot_d

static void radix4_table_d(Complex_d *data, int n, int s, int h, Trig_Block_d *table)
{ int v, a;
  int sign = (table == Coot_d);

#ifdef __SSE3__
  __m128d neg;

  if (sign)
    neg = LOADP(&FlipReal_d);
  else
    neg = LOADP(&Conjugate_d);
#endif

  a = NINETY/s; 
  for (v = 0; v < n; v += h)
    { Complex_d *j0 = data + v;
      Complex_d *j1 = j0+s;
      Complex_d *j2 = j1+s;
      Complex_d *j3 = j2+s;
      Complex_d *je = j1;

      Trig_Block_d *t = table;

      while (j0 != je)

#ifdef __SSE3__
        { __m128d vr, vi;
          __m128d d, e;
          __m128d t0, t1, t2, t3;

          vr = LOADUP(&(t->pow1.real));
          vi = LOADUP(&(t->pow1.imag));

          d  = LOADP(j2);
          CMULTIPLY(t2,d,vr,vi);

          vr = LOADUP(&(t->pow2.real));
          vi = LOADUP(&(t->pow2.imag));

          d  = LOADP(j1);
          CMULTIPLY(t1,d,vr,vi);

          vr = LOADUP(&(t->pow3.real));
          vi = LOADUP(&(t->pow3.imag));

          d  = LOADP(j3);
          CMULTIPLY(t3,d,vr,vi);

          t0 = LOADP(j0);

          d  = ADDP(t0,t1);
          e  = ADDP(t2,t3);

          STOREP(j0,ADDP(d,e));
          STOREP(j2,SUBP(d,e));

          d  = SUBP(t0,t1);
          e  = XORP(SHUFFLE(SUBP(t3,t2)),neg);

          STOREP(j1,ADDP(d,e));
          STOREP(j3,SUBP(d,e));

#else
        { double r0, r1, r2, r3;
          double i0, i1, i2, i3;
          double dr, di;
          double vr, vi;

          vr = t->pow1.real;
          vi = t->pow1.imag;

          dr = j2->real;
          di = j2->imag;

          r2 = vr * dr - vi * di;
          i2 = vr * di + vi * dr;

          vr = t->pow2.real;
          vi = t->pow2.imag;

          dr = j1->real;
          di = j1->imag;

          r1 = vr * dr - vi * di;
          i1 = vr * di + vi * dr;

          vr = t->pow3.real;
          vi = t->pow3.imag;

          dr = j3->real;
          di = j3->imag;

          r3 = vr * dr - vi * di;
          i3 = vr * di + vi * dr;

          r0 = j0->real;
          i0 = j0->imag;

          dr = r0 + r1;
          di = i0 + i1;
          vr = r2 + r3;
          vi = i2 + i3;

          j0->real = dr + vr;
          j0->imag = di + vi;
          j2->real = dr - vr;
          j2->imag = di - vi;

          dr = r0 - r1;
          di = i0 - i1;
          if (sign)
            { vr = i2 - i3;
              vi = r3 - r2;
            }
          else
            { vr = i3 - i2;
              vi = r2 - r3;
            }

          j1->real = dr + vr;
          j1->imag = di + vi;
          j3->real = dr - vr;
          j3->imag = di - vi;
#endif

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;
          t  += a;
        }
    }
}

  //  Radix 4 fft of span s (not using twiddle table):  h = 4s and h | n.

static void radix4_d(Complex_d *data, int n, int s, int h, double direct)
{ int    v;
  double theta = direct/h;
  double cos0  = cos(theta);
  double sin0  = sin(theta);
  int    sign  = (direct < 0);

#ifdef __SSE3__
  __m128d angr = LOADUP(&cos0);
  __m128d angi = LOADUP(&sin0);
  __m128d neg;

  if (sign < 0.)
    neg = LOADP(&FlipReal_d);
  else
    neg = LOADP(&Conjugate_d);
#endif

  for (v = 0; v < n; v += h)
    { Complex_d *j0 = data + v;
      Complex_d *j1 = j0+s;
      Complex_d *j2 = j1+s;
      Complex_d *j3 = j2+s;

#ifdef __SSE3__
      Complex_d *je = j1;
      __m128d    uc = Unity_d;  // uc holds w_h^u in normal form
                                //    after iteration u of the loop
      while (j0 != je)
        { __m128d d, e;
          __m128d t0, t1, t2, t3;
          __m128d ur, ui;

          ui = _mm_unpackhi_pd(uc,uc);   //  ui,ur = twiddle form of uc
          ur = _mm_movedup_pd(uc);

          d  = LOADP(j2);
          CMULTIPLY(t2,d,ur,ui);

          NEXT_ROOT(uc,ur,ui);           //  ui,ur = twiddle form of uc^2

          d  = LOADP(j1);
          CMULTIPLY(t1,d,ur,ui);

          NEXT_ROOT(uc,ur,ui);           //  ui,ur = twiddle form of uc^3

          d  = LOADP(j3);
          CMULTIPLY(t3,d,ur,ui);

          d  = MULP(angr,uc);
          uc = ADDSUB(d,MULP(angi,SHUFFLE(uc)));   //   uc *= ang (= [w_h^2])

          t0 = LOADP(j0);

          d  = ADDP(t0,t1);
          e  = ADDP(t2,t3);

          STOREP(j0,ADDP(d,e));
          STOREP(j2,SUBP(d,e));

          d  = SUBP(t0,t1);
          e  = XORP(SHUFFLE(SUBP(t3,t2)),neg);

          STOREP(j1,ADDP(d,e));
          STOREP(j3,SUBP(d,e));

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;
        }
    }

#else
      int    u;
      double ur, ui;

      ur = 1.0;         // ur + i*ui is w_h^u after iteration u of the loop
      ui = 0.0;

      for (u = 0; u < s; u++)
        { double r0, r1, r2, r3;
          double i0, i1, i2, i3;
          double dr, di;
          double vr, vi;

          dr = j2->real;
          di = j2->imag;

          r2 = ur * dr - ui * di;
          i2 = ur * di + ui * dr;

          vr = ur * ur - ui * ui;
          vi = 2. * ur * ui;

          dr = j1->real;
          di = j1->imag;

          r1 = vr * dr - vi * di;
          i1 = vr * di + vi * dr;

          dr = vr * ur - vi * ui;
          vi = vr * ui + vi * ur;
          vr = dr;

          dr = j3->real;
          di = j3->imag;

          r3 = vr * dr - vi * di;
          i3 = vr * di + vi * dr;

          r0 = j0->real;
          i0 = j0->imag;

          dr = r0 + r1;
          di = i0 + i1;
          vr = r2 + r3;
          vi = i2 + i3;

          j0->real = dr + vr;
          j0->imag = di + vi;
          j2->real = dr - vr;
          j2->imag = di - vi;

          dr = r0 - r1;
          di = i0 - i1;
          if (sign)
            { vr = i2 - i3;
              vi = r3 - r2;
            }
          else
            { vr = i3 - i2;
              vi = r2 - r3;
            }

          j1->real = dr + vr;
          j1->imag = di + vi;
          j3->real = dr - vr;
          j3->imag = di - vi;

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;

          dr = ur * cos0 - ui * sin0; 
          ui = ur * sin0 + ui * cos0;
          ur = dr;
        }
    }

#endif
}

  //  Radix 2 auto-sort (combination of fft stage and sorting pass): s < h and s*h = n/2
  //    For all j s.t. j mod 2h < h and j mod 2s < s
  //        D[0,0] = D[0,0] + u*D[1,0]
  //        D[0,1] = D[0,0] - u*D[1,0]
  //        D[1,0] = D[0,1] + u*D[1,1]
  //        D[1,1] = D[0,1] - u*D[1,1]
  //    where
  //      D[x,y] = data[j+x*h+y*s], u = d(w)^(j mod s), w = e^(-TPI/h)i,
  //          and d(x) = x  if table = Root_d
  //                   = x* if table = Coot_d

static void sort2_table_d(Complex_d *data, int n, int s, int h, Trig_Block_d *table)
{ int    h2, s2;
  int    v, w, a;

  h2 = h*2;
  s2 = s*2;
  a  = NINETY/s;
  for (v = 0; v < n; v += h2)
   for (w = 0; w < h; w += s2)
    { Complex_d *j0 = data + (v+w);
      Complex_d *j1 = j0+s;
      Complex_d *j2 = j0+h;
      Complex_d *j3 = j2+s;
      Complex_d *je = j1;

      Trig_Block_d *t = table;

      while (j0 != je)

#ifdef __SSE3__
        { __m128d p, q, d;
          __m128d ui, ur;

          ur = LOADUP(&(t->pow2.real));
          ui = LOADUP(&(t->pow2.imag));

          d = LOADP(j2);
          CMULTIPLY(p,d,ur,ui);

          d = LOADP(j3);
          CMULTIPLY(q,d,ur,ui);

          d = LOADP(j1);
          STOREP(j3,SUBP(d,q));
          STOREP(j2,ADDP(d,q));

          d = LOADP(j0);
          STOREP(j1,SUBP(d,p));
          STOREP(j0,ADDP(d,p));

#else
        { double dr, di;
          double pr, pi;
          double qr, qi;
          double ur, ui;

          ur = t->pow2.real;
          ui = t->pow2.imag;

          dr = j2->real;
          di = j2->imag;

          pr = ur * dr - ui * di;
          pi = ur * di + ui * dr;

          dr = j3->real;
          di = j3->imag;

          qr = ur * dr - ui * di;
          qi = ur * di + ui * dr;

          dr = j1->real;
          di = j1->imag;

          j3->real = dr - qr;
          j3->imag = di - qi;
          j2->real = dr + qr;
          j2->imag = di + qi;

          dr = j0->real;
          di = j0->imag;

          j1->real = dr - pr;
          j1->imag = di - pi;
          j0->real = dr + pr;
          j0->imag = di + pi;
#endif

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;
          t  += a;
        }
    }
}

  //  Radix 2 auto-sort specialized for s=1, h = n/2

static void first2_table_d(Complex_d *data, int n)
{ Complex_d *j0 = data;
  Complex_d *j1 = j0+1;
  Complex_d *j2 = j0+n/2;
  Complex_d *j3 = j2+1;
  Complex_d *je = j2;

  while (j0 != je)
#ifdef __SSE3__
    { __m128d p, d, q;

      q = LOADP(j2);
      p = LOADP(j3);
      d = LOADP(j1);
      STOREP(j3,SUBP(d,p));
      STOREP(j2,ADDP(d,p));

      d = LOADP(j0);
      STOREP(j1,SUBP(d,q));
      STOREP(j0,ADDP(d,q));
#else
    { double dr, di;
      double pr, pi;
      double qr, qi;

      pr = j2->real;
      pi = j2->imag;
      qr = j3->real;
      qi = j3->imag;
      dr = j1->real;
      di = j1->imag;
      j3->real = dr - qr;
      j3->imag = di - qi;
      j2->real = dr + qr;
      j2->imag = di + qi;

      dr = j0->real;
      di = j0->imag;
      j1->real = dr - pr;
      j1->imag = di - pi;
      j0->real = dr + pr;
      j0->imag = di + pi;
#endif
      j0 += 2;
      j1 += 2;
      j2 += 2;
      j3 += 2;
    }
}

  //  Check that n is a power of 2 and that m does not exceed the max problem size

static int argcheck_1d(int n, int m, char *source)
{ int s;

  s = n;
  while (s >= 2)
    { if ((s & 0x1) != 0)
        { if (grab_message()) return (1);
          strcpy(FFT_Esource,source);
          sprintf(FFT_Estring,"n = %d is not a power of 2",n);
          return (1);
        }
      s >>= 1;
    }

  if (m > NINETY*NINETY*8)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"n is larger than hardcoded maximum of %dM",
                          ((n/m)*(NINETY*NINETY))>>17);
      return (1);
    }

  return (0);
}

  //  Basic 1-dimenstional FFT-algorithm.  The FFT is performed in-place within 'data' and
  //    for convenience a pointer to data is returned by FFT_1d.  If invert is non-zero then
  //    the inverse Fourier Transform is performed.

Complex_d *FFT_1d(int n, Complex_d *data, int invert)
{ double        direct;
  Trig_Block_d *table;

  pthread_mutex_lock(&FFT_D_MUTEX);    //  Set up trig tables
  if (Trig_Firstime_d)
    { Trig_Firstime_d = 0;
      init_trig_tables_d();
    }
  pthread_mutex_unlock(&FFT_D_MUTEX);

  if (invert)            //  Establish direction
    { direct = -TPI;
      table  = Coot_d;
    }
  else
    { direct =  TPI;
      table  = Root_d;
    }

  if (argcheck_1d(n,n,"FFT_1d"))    //  Check that n is a power of 2 and not too large
    return (NULL);

  if (n <= 2)           //  For n = 1 nothing to do, for n = 2 a tabled radix 2 possibly
    { if (n == 1)       //    followed by normalization (if invert on)
        return (data);
      radix2_table_d(data,2,1,2,table);
    }

  else                //  For n > 2:
    { int i, s, h;    //    Auto sort passes for the first half, radix-4 and maybe one radix-2
      int ss, hh, nn; //    for the remainder

      if (n < L2_CACHE/2)  //  If the h-span of a run phases fits in the L2-cache then do
        nn = n;            //     them all in a cross-cut of the passes in order to gain
      else                 //     L2 cache coherence.  Hugely important for big, big n.
        nn = L2_CACHE/2;

      first2_table_d(data,n);
      for (s = 2, h = n/4; s < h && L2_CACHE/4 < h; s <<= 1, h >>= 1)
        sort2_table_d(data,n,s,h,table);

      for (i = 0; i < n; i += L2_CACHE/2)
        { for (ss = s, hh = h; ss < hh; ss <<= 1, hh >>= 1)
            sort2_table_d(data+i,nn,ss,hh,table);

          for (hh = ss<<2; hh <= nn; ss <<= 2, hh <<= 2)
            if (ss > NINETY)
              radix4_d(data+i,nn,ss,hh,direct);
            else
              radix4_table_d(data+i,nn,ss,hh,table);
        }

      for (s = ss, h = hh; h <= n; s <<= 2, h <<= 2)
        if (s > NINETY)
          radix4_d(data,n,s,h,direct);
        else
          radix4_table_d(data,n,s,h,table);

      if (s < n)
        { h = (s << 1);
          if (s > NINETY)
            radix2_d(data,n,s,h,direct);
          else
            radix2_table_d(data,n,s,h,table);
        }
    }

  if (invert)               //  Normalize the inverse transform if invert is on
    { Complex_d *d;
      Complex_d *e = data+n;
      double     v = 1./n;

#ifdef __SSE3__
      __m128d p = LOADUP(&v);

      for (d = data; d < e; d++)
        STOREP(d,MULP(p,LOADP(d)));
#else
      for (d = data; d < e; d++)
        { d->real *= v;
          d->imag *= v;
        }
#endif
    }

  return (data);
}


/*********************************************************************************************\
 *                                                                                            *
 *  1-D REAL FFT & INVERSE:                                                                   *
 *                                                                                            *
\*********************************************************************************************/

  //  FFT-algorithms optimized for the case of a real-valued time-domain data.
  //
  //    The forward transform, Real_FFT_1d, takes a double array of length n, and in-place produces 
  //    a Complex_d array C of length n/2 that is the first half of the conjugate symmetric FFT
  //    of the real data with the exception that F_(n/2) (which is real) is tucked into C[0].imag
  //    (which works as F_0 is also real).  Again, the pointer returned by Real_FFT_1d is really the
  //    same as rdata, the FFT is performed in-place. 

Complex_d *Real_FFT_1d(int n, double *rdata)
{ Complex_d *rfft = (Complex_d *) rdata;
  int n2 = n/2;
  int n4 = n/4;
 
  if (argcheck_1d(n,n2,"Real_FFT_1d"))   //  Check that n is a power of 2 and not too large
    return (NULL);

  if (n == 1) return (rfft);

  FFT_1d(n2,rfft,0);

  { double zr = rfft[0].real;   //  Special case for the 0^th and n/2^th elements
    double zi = rfft[0].imag;   //     whose x'forms are real 

    rfft[0].real = zr + zi;     //  Put F_0 in .real and F_n2 in .imag of 0^th element
    rfft[0].imag = zr - zi;
  }

  if (n == 2) return (rfft);

  if (n4 > NINETY)                 //  Code when n too big for tabled twiddles
    { double theta = TPI / n;
      double cos0  = cos(theta);
      double sin0  = sin(theta);

#ifdef __SSE3__
      __m128d ang = SETP(cos0,sin0);
      __m128d ur  = SETP(cos0,cos0);
      __m128d ui  = SETP(sin0,sin0);
#else
      double ur   = cos0;
      double ui   = sin0;
#endif
      Complex_d *ck = rfft + 1;
      Complex_d *cj = rfft + (n2-1);

      while (ck < cj)
#ifndef __SSE3__
        { double kr = ck->real;    //  Forward xform:  Ck = rfft[k]  Cj = rfft[M-k]
          double ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) - i(w^k)(Ck-Cj*))
          double jr = cj->real;    //    Cj = .5 ((Cj+Ck*) + i(w^k)*(Cj-Ck*))
          double ji = cj->imag;    //  where w = e^(-TPI/n)i and k = ck-rfft

          double f0r = .5 * (kr + jr);
          double f0i = .5 * (ki - ji);
          double f1r = .5 * (ji + ki);
          double f1i = .5 * (jr - kr);

          kr = ur*f1r - ui*f1i;    // (ur,ui) is the n'th root of unity to the power k
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r + kr;
          ck->imag =   f0i + ki;
          cj->real =   f0r - kr;
          cj->imag = -(f0i - ki);

          kr = ur * cos0 - ui * sin0;
          ui = ur * sin0 + ui * cos0;
          ur = kr;
#else
        { __m128d fk, fj;
          __m128d f0, f1;

          fk = LOADP(ck);
          fj = LOADP(cj);

          f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
          f1 = MULP(Half_d,ADDSUB(fj,fk));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));  //  [ur,ui] twiddle form of w_n^u
          STOREP(ck,ADDP(f0,fk));                         //     where u = ck-rfft
          STOREP(cj,XORP(Conjugate_d,SUBP(f0,fk)));

          NEXT_ROOT(ang,ur,ui);
#endif
          ck += 1;
          cj -= 1;
        }
    }

  else                               //  Code when the twiddles can come from a table
    { int           a = NINETY/n4;
      Trig_Block_d *t = Root_d + a;

      Complex_d *ck = rfft + 1;
      Complex_d *cj = rfft + (n2-1);

      while (ck < cj)
#ifndef __SSE3__
        { double kr = ck->real;
          double ki = ck->imag;
          double jr = cj->real;
          double ji = cj->imag;

          double f0r = .5 * (kr + jr);
          double f0i = .5 * (ki - ji);
          double f1r = .5 * (ji + ki);
          double f1i = .5 * (jr - kr);

          double ur = t->pow1.real;
          double ui = t->pow1.imag;

          kr = ur*f1r - ui*f1i;
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r + kr;
          ck->imag =   f0i + ki;
          cj->real =   f0r - kr;
          cj->imag = -(f0i - ki);
#else
        { __m128d fk, fj;
          __m128d f0, f1;
          __m128d ui, ur;

          fk = LOADP(ck);
          fj = LOADP(cj);

          f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
          f1 = MULP(Half_d,ADDSUB(fj,fk));

          ur = LOADUP(&(t->pow1.real));
          ui = LOADUP(&(t->pow1.imag));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
          STOREP(ck,ADDP(f0,fk));
          STOREP(cj,XORP(Conjugate_d,SUBP(f0,fk)));
#endif
          ck += 1;
          cj -= 1;
          t  += a;
        }
    }

  return (rfft);
} 


  //    The inverse transform, Real_FFT_Inverse_1d, takes a complex half-matrix as produced by
  //    Real_FFT_1d, and produces a real-valued result *in-place*.  That is, the pointer returned is
  //    exactly rfft (coerced to be double *).  Note carefully that n is the length of the
  //    resulting real array and is twice the length of rfft.

double *Real_FFT_Inverse_1d(int n, Complex_d *rfft)
{ int n2 = n/2;
  int n4 = n/4;

  if (argcheck_1d(n,n2,"Real_FFT_Inverse_1d"))   //  Check that n is a power of 2 and not too large
    return (NULL);

  pthread_mutex_lock(&FFT_D_MUTEX);    //  Set up trig tables
  if (Trig_Firstime_d)
    { Trig_Firstime_d = 0;
      init_trig_tables_d();
    }
  pthread_mutex_unlock(&FFT_D_MUTEX);

  if (n == 1) return ((double *) rfft);

  { double zr = rfft[0].real;   //  Special case for reconstituting 0th element that depends on
    double zi = rfft[0].imag;   //    the 0^th and n/2^th value of the x'form packed in 0^th el.

    rfft[0].real = .5*(zr + zi);
    rfft[0].imag = .5*(zr - zi);
  }

  if (n == 2) return ((double *) rfft);

  if (n4 > NINETY)                   // Code when n too big for tabled twiddles
    { double theta = - TPI / n;
      double cos0  = cos(theta);
      double sin0  = sin(theta);

#ifdef __SSE3__
      __m128d ang = SETP(cos0,sin0);
      __m128d ur  = _mm_movedup_pd(ang);
      __m128d ui  = _mm_unpackhi_pd(ang,ang); 
#else
      double ur    = cos0;
      double ui    = sin0;
#endif

      Complex_d *ck = rfft + 1;
      Complex_d *cj = rfft + (n2-1);

      while (ck < cj)
#ifndef __SSE3__
        { double kr = ck->real;    //  Inverse xform:  Ck = rfft[k]  Cj = rfft[M-k]
          double ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) + i(w^-k)*(Ck-Cj*))
          double jr = cj->real;    //    Cj = .5 ((Cj+Ck*) - i(w^-k)(Cj-Ck*))
          double ji = cj->imag;    //  where w = e^(-TPI/n)i and k = ck-rfft

          double f0r = .5 * (kr + jr);
          double f0i = .5 * (ki - ji);
          double f1r = .5 * (ji + ki);
          double f1i = .5 * (jr - kr);

          kr = ur*f1r - ui*f1i;    // (ur,ui) is the n'th root of unity to the power -k
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r - kr;
          ck->imag =   f0i - ki;
          cj->real =   f0r + kr;
          cj->imag = -(f0i + ki);

          kr = ur * cos0 - ui * sin0;
          ui = ur * sin0 + ui * cos0;
          ur = kr;
#else
        { __m128d fk, fj;
          __m128d f0, f1;

          fk = LOADP(ck);
          fj = LOADP(cj);

          f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
          f1 = MULP(Half_d,ADDSUB(fj,fk));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));  //  [ur,ui] twiddle form of w_n^u*
          STOREP(ck,SUBP(f0,fk));                         //     where u = ck-rfft
          STOREP(cj,XORP(Conjugate_d,ADDP(f0,fk)));

          NEXT_ROOT(ang,ur,ui);
#endif
          ck += 1;
          cj -= 1;
        }
    }

  else                               //  Code when the twiddles can come from a table
    { int           a = NINETY/n4;
      Trig_Block_d *t = Coot_d + a;

      Complex_d *ck = rfft + 1;
      Complex_d *cj = rfft + (n2-1);

      while (ck < cj)
#ifndef __SSE3__
        { double kr = ck->real;
          double ki = ck->imag;
          double jr = cj->real;
          double ji = cj->imag;

          double f0r = .5 * (kr + jr);
          double f0i = .5 * (ki - ji);
          double f1r = .5 * (ji + ki);
          double f1i = .5 * (jr - kr);

          double ur = t->pow1.real;
          double ui = t->pow1.imag;

          kr = ur*f1r - ui*f1i;
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r - kr;
          ck->imag =   f0i - ki;
          cj->real =   f0r + kr;
          cj->imag = -(f0i + ki);
#else
        { __m128d fk, fj;
          __m128d f0, f1;
          __m128d ur, ui;

          fk = LOADP(ck);
          fj = LOADP(cj);

          f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
          f1 = MULP(Half_d,ADDSUB(fj,fk));

          ur = LOADUP(&(t->pow1.real));
          ui = LOADUP(&(t->pow1.imag));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
          STOREP(ck,SUBP(f0,fk));
          STOREP(cj,XORP(Conjugate_d,ADDP(f0,fk)));
#endif
          ck += 1;
          cj -= 1;
          t  += a;
        }
    }

  return ((double *) FFT_1d(n2,rfft,1));
} 


/*********************************************************************************************\
 *                                                                                            *
 *  1-D CONVOLUTION & CORRELATION                                                             *
 *                                                                                            *
\*********************************************************************************************/

  //  Complex_Convolution_1d performs the term-wise complex multiplication of fft1 and fft2
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //        FFT_1d(n,Complex_Convolution_1d(n,FFT_1d(n,data1,0),FFT_1d(n,data2,0)),1);

static Complex_d *complex_convolution_1d(int64 n, Complex_d *fft1, Complex_d *fft2)
{ Complex_d *s = fft1;
  Complex_d *f = fft2;
  Complex_d *h = fft2 + n;

  while (f < h)
#ifndef __SSE3__
    { double kr = s->real;
      double ki = s->imag;
      double jr = f->real;
      double ji = f->imag;

      s->real = kr*jr - ki*ji;
      s->imag = ki*jr + kr*ji;
#else
    { __m128d kr = LOADUP(&(s->real));
      __m128d ki = LOADUP(&(s->imag));
      __m128d  j = LOADP(f);
      __m128d  p;

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);
#endif
      f += 1;
      s += 1;
    }

  return (fft1);
}

Complex_d *Complex_Convolution_1d(int n, Complex_d *fft1, Complex_d *fft2)
{ return (complex_convolution_1d((int64) n, fft1, fft2)); }

  //  Complex_Correlation_1d performs the term-wise complex multiplication of fft1 and fft2
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //        FFT_1d(n,Complex_Correlation_1d(n,FFT_1d(n,data1,0),FFT_1d(n,data2,0)),1);

static Complex_d *complex_correlation_1d(int64 n, Complex_d *fft1, Complex_d *fft2)
{ Complex_d *s = fft1;
  Complex_d *f = fft2;
  Complex_d *h = fft2 + n;
#ifdef __SSE3__
  __m128d    conj = LOADP(&Conjugate_d);
#endif

  while (f < h)
#ifndef __SSE3__
    { double kr = s->real;
      double ki = s->imag;
      double jr = f->real;
      double ji = f->imag;

      s->real = kr*jr + ki*ji;
      s->imag = ki*jr - kr*ji;

#else
    { __m128d kr = LOADUP(&(s->real));
      __m128d ki = LOADUP(&(s->imag));
      __m128d  j = XORP(LOADP(f),conj);
      __m128d  p;

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);
#endif
      s += 1;
      f += 1;
    }

  return (fft1);
}

Complex_d *Complex_Correlation_1d(int n, Complex_d *fft1, Complex_d *fft2)
{ return (complex_correlation_1d((int64) n, fft1, fft2)); }

  //  Real_Convolution_1d effects a term-wise multiplication of the fft of two vectors encoded
  //    as packed half-sized complex vectors of length ** n/2 ** as produced by Real_FFT_1d.  The
  //    result is produced in-place in rfft1 and a double pointer to it is returned.
  //    The code works just fine if rfft2 is the same as array rfft1.  So a real convolution of
  //    of two real arrays data1 and data2, in-place within data 1, is accomplised with the code:
  //
  //  Real_FFT_Inverse_1d(n,Real_Convolution_1d(n,Real_FFT_1d(n,data1,0),Real_FFT_1d(n,data2,0)),1);

Complex_d *Real_Convolution_1d(int n, Complex_d *rfft1, Complex_d *rfft2)
{ Complex_d *f = rfft2;
  Complex_d *s = rfft1;
  Complex_d *g = rfft2 + n/2;

  //  Basically all one is doing is element-wise multiplying the FFTs of rfft1 and rfft2
  //    with the subtlety that the 0th element codes 2 real valued terms instead of 1.

  if (n == 1)
    { s->real *= f->real;
      return (rfft1);
    }

  if ((n & 0x1) != 0)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"Real_Convolution_1d");
      sprintf(FFT_Estring,"n = %d must be even",n);
      return (NULL);
    }

  s->real *= f->real;
  s->imag *= f->imag;

  s += 1;
  f += 1;

  while (f < g)
#ifdef __SSE3__
    { __m128d kr = LOADUP(&(s->real));
      __m128d ki = LOADUP(&(s->imag));
      __m128d  j = LOADP(f);
      __m128d  p;

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);
#else
    { double kr = s->real;
      double ki = s->imag;
      double jr = f->real;
      double ji = f->imag;

      s->real = kr*jr - ki*ji;
      s->imag = ki*jr + kr*ji;
#endif
      s += 1;
      f += 1;
    }

  return (rfft1);
}

  //  Real_Correlation_1d effects a term-wise multiplication of the fft of a vector and the
  //    conjugate of another, both encoded as packed half-sized complex vectors of length ** n/2 **
  //    as produced by Real_FFT_1d.  The result is produced in-place in rfft1 and a double
  //    pointer to it is returned.  The code works just fine if rfft2 is the same as array rfft1.
  //    So a real convolution of two real arrays data1 and data2, in-place within data 1, is
  //    accomplised with the code:
  //
  //  Real_FFT_Inverse_1d(n,Real_Correlation_1d(n,Real_FFT_1d(n,data1,0),Real_FFT_1d(n,data2,0)),1);

Complex_d *Real_Correlation_1d(int n, Complex_d *rfft1, Complex_d *rfft2)
{ Complex_d *f = rfft2;
  Complex_d *s = rfft1;
  Complex_d *g = rfft2 + n/2;
#ifdef __SSE3__
  __m128d conj = LOADP(&Conjugate_d);
#endif

  //  Basically all one is doing is element-wise multiplying the FFT of rfft1 and the element-wise
  //    conjugate of the FFT of rfft2 with the subtlety that the 0th element codes 2 real valued
  //    terms instead of 1.

  if (n == 1)
    { s->real *= f->real;
      return (rfft1);
    }

  if ((n & 0x1) != 0)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"Real_Correlation_1d");
      sprintf(FFT_Estring,"n must be even");
      return (NULL);
    }

  s->real *= f->real;
  s->imag *= f->imag;

  s += 1;
  f += 1;

  while (f < g)
#ifdef __SSE3__
    { __m128d kr = LOADUP(&(s->real));
      __m128d ki = LOADUP(&(s->imag));
      __m128d  j = XORP(LOADP(f),conj);
      __m128d  p;

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);
#else
    { double kr = s->real;
      double ki = s->imag;
      double jr = f->real;
      double ji = f->imag;

      s->real = kr*jr + ki*ji;
      s->imag = ki*jr - kr*ji;
#endif
      s += 1;
      f += 1;
    }

  return (rfft1);
}


/*********************************************************************************************\
 *                                                                                            *
 *  NORMALIZATION                                                                             *
 *                                                                                            *
\*********************************************************************************************/

  //  For detecting an "optimal" correlation overlap, it has been found that finding the maximum
  //    "normalized" correlation score in the spatial domain gives the desired overlap.  The
  //    normalized correlation of interval R versus S is
  //
  //                        SUM_(x,y)_in_(R,S) (f_x-m_R)*(g_x-m_S)/((|R|-1)*s_R*s_S)
  //
  //    where m_R is the mean value of f over R, m_S is the mean value of g over S, and s_R and
  //    s_S are the standard deviations over R and S, respectively.  Note that this is tricky
  //    because the ranges R and S are different for each correlation term.  However, it can
  //    still be done in linear time by incrementally computing the necessary sums as R and S
  //    are progressive intervals of f and g.  Specifically, the sum above equals
  //
  //         (Correlation(f_R,g_S) - n*F*G) / ((n-1)*((F2/n-(F/n)^2)*(G2/n-(G/n)^2))^.5)
  //
  //    where n = |R| = |S|, F = SUM_x_in_R f_x, F2 = SUM_x_in_R f_x^2, etc.
  //
  //  Normalize_1d normalizes the scores in the spatial correlation vector cdata that
  //    is presumed to have been produced by applying FFT routines to copies of idata
  //    and tdata that had been 0-padded to be of length cdim.  Note that idim+tdim <= cdim or
  //    the vectors were insufficiently padded to give every overlap correlation.  The
  //    normalization takes effect directly on cdata and for convenience a pointer
  //    to it is returned.

double *Normalize_1d(int idim, double *idata, int tdim, double *tdata, int cdim, double *cdata)
{ double IA, ID;
  double TA, TD;

  int    b, e;
  int    f, c;
  int    w, t;

  if (idim > cdim || tdim > cdim)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"Normalize_1d");
      sprintf(FFT_Estring,"Operand dimensions > correlation dimension");
      return (NULL);
    }

  //  At position b, let I = [b,e] intersect [0,idim-1] = the el.s in idata for con. value at b
  //                 and T = [f,c] intersect [0,tdim-1] = the el.s in tdata for con. value at b
  //     where [x,y] = [0,y] union [x,cdim-1] if x > y
  //
  //  Then IA = SUM_x_in_I idata[x], ID = SUM_x_in_I idata[x]^2
  //   and TA = SUM_x_in_T tdata[x], TD = SUM_x_in_T tdata[x]^2
  //   and  w = |I| = |T|

  w = (idim + tdim) - cdim;
  if (w < 0)
    w = 0;

  IA = ID = TA = TD = 0.;
  c = tdim;
  for (b = 0; b < w; b++)
    { double s = idata[b];
      IA += s;
      ID += s*s;
      s = tdata[--c];
      TA += s;
      TD += s*s;
    }

  c = cdim-1;
  if (idim == cdim)
    { f = c;
      b = 0;
    }
  else
    { f = c-idim;
      b = idim;
    }
  e = b + tdim;
  if (e >= cdim)
    e -= cdim;

  t = b;
  do
    { { double ia = IA/w;
        double ta = TA/w;
        double dn = sqrt((ID/w - ia*ia)*(TD/w - ta*ta));
        if (dn == 0.)
          cdata[b] = 0.;
        else
          cdata[b] = ((cdata[b] - ia*TA) / ((w-1)*dn));
      }

      if (b < idim)
        { double s = idata[b];
          IA -= s;
          ID -= s*s;
          w  -= 1;
        }
      if (++b >= cdim)
        b = 0;

      if (e < idim)
        { double s = idata[e];
          IA += s;
          ID += s*s;
          w  += 1;
        }
      if (++e >= cdim)
        e = 0;

      if (c < tdim)
        { double s = tdata[c];
          TA -= s;
          TD -= s*s;
        }
      if (c-- <= 0)
        c = cdim-1;

      if (f < tdim)
        { double s = tdata[f];
          TA += s;
          TD += s*s;
        }
      if (f-- <= 0)
        f = cdim-1;
    }
  while (b != t);

  return (cdata);
}


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL COMPLEX FFT:                                                            *
 *       init_multi_fft_d, fft_a_dim_d                                                        *
 *                                                                                            *
\*********************************************************************************************/

static int64 init_multi_fft_d(int d, int *dims, int *dmax)
{ int64 volume;
  int   max;

  { int n, i, s;           //  Check that n is a power of 2 and not too big
                           //    for every dimension
    volume = 1;
    max = 0;
    for (i = 0; i < d; i++)
      { n = dims[i];

        s = n;
        while (s >= 2)
          { if ((s & 0x1) != 0)
              { if (grab_message()) return (0);
                sprintf(FFT_Esource,"FFT_nd");
                sprintf(FFT_Estring,"dims[%d] = %d is not a power of 2",d,n);
                return (0);
              }
            s >>= 1;
          }

        if (n > NINETY*4)
          { if (grab_message()) return (0);
            sprintf(FFT_Esource,"FFT_nd");
            sprintf(FFT_Estring,"dims[%d] is larger than hardcoded maximum of %dK",d,(NINETY)>>8);
            return (0);
          }

        if (n > max)
          max = n;
        volume *= n;
      }
  }

  pthread_mutex_lock(&FFT_D_MUTEX);  //  Initialize the twiddle tables
  if (Trig_Firstime_d)
    { Trig_Firstime_d = 0;
      init_trig_tables_d();
    }
  pthread_mutex_unlock(&FFT_D_MUTEX);

  *dmax = max;
  return (volume);
}

  //  Perform an FFT on a given dimension, other than the first, of the matrix 'data'.  The
  //  matrix contains 'size' elements, there are 'n' > 1 elements in a row of the given dimension,
  //  and those elements are 'step' elements apart in the matrix.  Invert signals the direction
  //  of the transform.

  //  The strategy is to move elements from a row to the supplied cache, performing the
  //  first step of the auto-sort as you do so.  Once the transform of a row is
  //  complete than it is moved back into the matrix, being normalized at the same time
  //  if invert is on.

static void fft_a_dim_d(int64 size, int64 step, int n, int invert,
                        Complex_d *data, Complex_d *cache)
{ int n2      = n/2;
  int n4      = n/4;
  int64 ostep = step*n;
  int64 hstep = step*n2;

  int64 i, j;
  int64 jtop;

  double        direct;
  Trig_Block_d *table;

  if (invert)             //  Establish direction
    { direct = -TPI;
      table  = Coot_d;
    }
  else
    { direct =  TPI;
      table  = Root_d;
    }

  //  Annoying, but there is no point caching rows of 2 dimensions, so it is handled
  //    as a special case below where a radix 2 FFT of span 1 is performed in place

  if (n == 2)
    { for (j = 0; j < size; j += ostep)
        { Complex_d *s0 = data + j;
          Complex_d *s2 = s0 + hstep;
          Complex_d *se = s2;

          while (s0 != se)
#ifdef __SSE3__
            { __m128d p, q;

              p = LOADP(s0);
              q = LOADP(s2);
              STOREP(s0,ADDP(p,q));
              STOREP(s2,SUBP(p,q));
#else
            { double kr = s0->real;
              double ki = s0->imag;
              double jr = s2->real;
              double ji = s2->imag;

              s0->real = kr + jr;
              s0->imag = ki + ji;
              s2->real = kr - jr;
              s2->imag = ki - ji;

#endif
              s0 += 1;
              s2 += 1;
            }
        }
      if (invert)
        { Complex_d *s = data;
          Complex_d *e = data + size;

          while (s != e)
#ifdef __SSE3__
            { __m128d q = MULP(Half_d,LOADP(s));
              STOREP(s,q);
#else
            { s->real *= .5;
              s->imag *= .5;
#endif
              s += 1;
            }
        }
      return;
    }

  //  NB: n >= 4 if you get here.

  for (j = 0; j < size; j += ostep)
    { jtop = j+step;

      for (i = j; i < jtop; i++)  //  Move and auto sort the row starting at element i into cache
        {
          { Complex_d *s0 = data + i;
            Complex_d *s2 = s0 + hstep;
            Complex_d *j0 = cache;
            Complex_d *j2 = j0 + n2;
            Complex_d *je = j2;

            while (j0 != je)
#ifdef __SSE3__
              { __m128d p, q;

                p = LOADP(s0);
                s0 += step;
                q = LOADP(s2);
                s2 += step;

                STOREP(j0,ADDP(p,q));
                j0 += 1;
                STOREP(j0,SUBP(p,q));
                j0 += 1;

                p = LOADP(s0);
                s0 += step;
                q = LOADP(s2);
                s2 += step;

                STOREP(j2,ADDP(p,q));
                j2 += 1;
                STOREP(j2,SUBP(p,q));
                j2 += 1;
              }
#else
              { double pr, pi;
                double dr, di;

                pr = s2->real;
                pi = s2->imag;
                s2 += step;

                dr = s0->real;
                di = s0->imag;
                s0 += step;

                j0->real = dr + pr;
                j0->imag = di + pi;
                j0 += 1;
                j0->real = dr - pr;
                j0->imag = di - pi;
                j0 += 1;

                pr = s2->real;
                pi = s2->imag;
                s2 += step;

                dr = s0->real;
                di = s0->imag;
                s0 += step;

                j2->real = dr + pr;
                j2->imag = di + pi;
                j2 += 1;
                j2->real = dr - pr;
                j2->imag = di - pi;
                j2 += 1;
              }
#endif
          }

          { int s, h;   //  Sort passes for the first half, radix-4 and maybe one radix-2 pass
                        //    for the remainder
            for (s = 2, h = n4; s < h; s <<= 1, h >>= 1)
              sort2_table_d(cache,n,s,h,table);

            for (h = (s << 2); h <= n; s <<= 2, h <<= 2)
              radix4_table_d(cache,n,s,h,table);

            if (s < n)
              { h = (s << 1);
                radix2_table_d(cache,n,s,h,table);
              }
          }

          { Complex_d *s = data+i;    //  Return the cached row to the matrix
            Complex_d *d = cache;     //    and normalize if invert is on
            Complex_d *e = d + n;

            if (invert)
              { double v = 1./n;
#ifdef __SSE3__
                __m128d p = LOADUP(&v);
                __m128d q;

                while (d != e)
                  { q = MULP(p,LOADP(d));
                    STOREP(s,q);
                    s += step;
                    d += 1;
                  }
#else
                while (d != e)
                  { s->real = d->real * v;
                    s->imag = d->imag * v;
                    s += step;
                    d += 1;
                  }
#endif
              }
            else
              { while (d != e)
                  { *s = *d++;
                    s += step;
                  }
              }
          }
        }
    }
}

  //  Basic n-dimenstional FFT-algorithm.  The FFT is performed in-place within 'data' and
  //    for convenience a pointer to data is returned by FFT_nd.  If invert is non-zero then
  //    the inverse Fourier Transform is performed.

Complex_d *FFT_nd(int ndim, int *dims, Complex_d *data, int invert)
{ int64     size;
  int64     step;
  int64     k;
  int       d;
  int       dmax;
  Complex_d *cache;
  Complex_d  Cache[4096];

  size = init_multi_fft_d(ndim,dims,&dmax);
  if (size <= 0)
    return (NULL);

  if (dmax > 4096)
    cache = (Complex_d *) malloc(sizeof(Complex_d)*((size_t) dmax));
  else
    cache = Cache;
  if (cache == NULL)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"FFT_nd");
      sprintf(FFT_Estring,"Out of memory");
      return (NULL);
    }

  step = dims[0];      // Do innermost rows in place with 1-d routines (if the dimesnion is > 1!)
  if (step > 1)
    for (k = 0; k < size; k += step)
      FFT_1d((int) step,data+k,invert);

  for (d = 1; d < ndim; d++)    // Do rows of each higher dimension (if the dimension is > 1!)
    if (dims[d] > 1)
      { fft_a_dim_d(size,step,dims[d],invert,data,cache);
        step *= dims[d];
      }

  if (dmax > 4096)
    free(cache);

  return (data);
}


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL REAL FFT & INVERSE: argcheck_nd                                         *
 *                                                                                            *
\*********************************************************************************************/

  //  Common argument checks for real multi-dim. FFT routines

static int argcheck_nd(int ndim, int *dims, char *source)
{ if (ndim > MAX_REAL_DIM)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"Array cannot have more than %d dimensions",MAX_REAL_DIM);
      return (1);
    }
  if (dims[0]%4 != 0 && dims[0] != 2 && dims[0] != 1)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"dims[0] = %d is not a power of 2",dims[0]);
      return (1);
    }
  if (dims[0] > NINETY*8)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"dims[0] is larger than hardcoded maximum of %dK",(NINETY)>>7);
      return (1);
    }
  return (0);
}

  //  Multi-dimensional FFT-algorithms optimized for the case of a real-valued time-domain data.
  //
  //    The forward transform, Real_FFT_nd, takes a dobule array of size s, and *in-place* produces
  //    a Complex_d array c of size s/2 that is the conjugate symmetric FFT of the real data for
  //    the first half of the lowest dimension [0..M-1] where M = dims[0]/2.  In fact, some M-terms
  //    are essential and tucked into 0-terms that are redundant, so the complex array is not
  //    directly interpretable as the lower half of the FFT, but rather is an encoding of all the
  //    terms not inferable by symmetry and hence sufficient for computing correlations and
  //    convolutions.
  //
  //    The pointer C returned by FFT is really the same as rdata, the FFT is performed in-place.
  //    To fully document the encoding in C, if a = (ik,...,i1) then a* = (nk-ik mod nk, ...,
  //    n1-i1 mod n1 ) where k = ndims-1 and nx = dims[x].  Then the values of C encode values
  //    of the real-valued FFT F as follows:
  //
  //           C[a,0] = F[a,0]                      if a < a*   (F[a,M]  = F[a*,M]*)
  //                  = F[a*,M]                     if a > a*   (F[a*,0] = F[a,0]*)
  //                  = F[a,0].real + i F[a,M].real if a = a*   (F[a,0] and F[a,M] are real!)
  //           C[a,x] = F[a,x]                      x > 0
  //
  //    where a < a* is with respect to lexicographical order
  //
  //    For the improbable but still possible case that you ask for the real FFT of a matrix
  //    that has a run of lower dimensions that are 1, then the lower dimensions are effectively
  //    ignored in computing a representation, albeit the dimensions of the result are not
  //    alterred (at least from the perspective of the user)

Complex_d *Real_FFT_nd(int in_ndim, int *in_dims, double *rdata)
{ Complex_d *rfft = (Complex_d *) rdata;
  int        M, A, *dims;
  int        ndim;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int64      size, nsub, offs;

  if (argcheck_nd(in_ndim,in_dims,"Real_FFT_nd"))
    return (NULL);

  { int i;      //  If lowest dims are 1, then effectively ignore them (by modifying ndim and dims)
                //                     in what follows
    for (i = 0; i < in_ndim; i++)
      if (in_dims[i] > 1)
        break;
    if (i == in_ndim)
      return (rfft);
    dims = in_dims+i;
    ndim = in_ndim-i;
  }

  M = dims[0] / 2;
  A = NINETY / M;

  dims[0] = M;                  //  Compute H, the fft of rdata as if it were a complex array
  if (FFT_nd(ndim,dims,rfft,0) == NULL)
    { if (my_message())
        sprintf(FFT_Esource,"Real_FFT_nd");
      dims[0] = M * 2;
      return (NULL);
    }
  dims[0] = M * 2;

  //  Given H the fft of rfft, the FFT of rdata for (a,x) in [0,nk-1] x ... [0,n1-1] x [0,M] is:
  //
  //      F[a,x] = .5*(H[a,x] + H[a*,x*]*) - .5*i*w^x*(H[a,x] - H[a*,x*]*)
  //
  //  where if a = (ik,...,i1) then a* = (nk-ik mod nk, ..., n1-i1 mod n1 ) where k = ndims-1,
  //  nx = dims[x], x* = M-x mod M, w = e^(-TPI/(2*M))i, and M = dims[0]/2.
  //
  //  Notice that we only have a complex array that has M elements in the first dimension as
  //  opposed to the M+1 needed to record F[a,x] for x in [0,M] (as opposed to [0,M-1]).
  //  Fortunately there are symmetries that allow us to pack all the non-redundant values into
  //  half the space as given by entity C[a,x] defined at the header of the routine.  To explain
  //  that a bit further:
  //
  //  The basic symmetry of a real-valued FFT is that F[a,x] = F[a*,2M-x]*.  We put F[a,x] into
  //  C[a,x] for all a and x in [1,M-1].  By the symmetry we thus have all non-redundant values
  //  save those in F[a,0] and F[a,M].  We pack all the non-redundant values from these two slices
  //  into the single slice C[a,0].  The good news is that F[a,M] = F[a*,M]* and F[a,0] = F[a*,0]*
  //  so as long as a != a* if suffices to keep F[a,0] and F[a*,M] for a < a*.  So we place F[a,0]
  //  into C[a,0] and F[a*,M] into C[a*,0] for a < a*.  Now the only values not recorded are F[a,0]
  //  and F[a,M] when a = a* and the only space for them is in C[a,0].  The good news is that both
  //  of these values of any given a = a* are real (F[a,0] = F[a*,0]* = F[a,0]* ==> F[a,0] real)!
  //  So we pack F[a,0] in the real part of C[a,0] and F[a,M] into the imaginary part.
  //
  //  The tricky coding then is that we must know a* for any given a as the algorithm proceeds.
  //  A little algebra reveals that if the index of [a,0] = [ik,...,i1,0] is A then the index
  //  of A* = [a*,0] is:
  //
  //         A* = Nk - (A + SUM_(ix=0) Nx)
  //
  //             where Nx = SUM_j=1^x ( PI_i=1^j nj ) * M

  { int d;                     //  Get total size of array and set up dimensional terms needed
                               //    to compute complement indices
    size = M;                  //  a = (nk-1,...,n1-1) = (cntr[d])_d=1,k
    nsub = 0;                  //  base[d] = Nd
    offs = 0;                  //  offs = SUM_(ix=0) Nx for current (cntr[d])
    for (d = 1; d < ndim; d++)
      { size   *= dims[d];
        nsub   += size;
        base[d] = size;
        cntr[d] = dims[d]-1;
        if (dims[d] == 1)
          offs += size;
      }
  }

  { int64 a, b;

    for (a = 0; a < size; a += M)      //  For each index A = [a,0] for a in [0,size/M]
      { { int d;

          for (d = 1; d < ndim; d++)   //  Compute b = a* (O(1) amortized time per iteration!)
            if (dims[d] > 1)           //  Advance counters for a and update offs accordingly
              { cntr[d] += 1;
                if (cntr[d] == 1)
                  { offs -= base[d];
                    break;
                  }
                else if (cntr[d] >= dims[d])
                  { cntr[d] = 0;
                    offs += base[d];
                  }
                else
                  break;
              }
          b = nsub - (offs + a);       //  Apply the formula above
        }

        { Trig_Block_d *t  = Root_d + A;
          Complex_d    *ck, *cj, *ce;

          if (a < b)                   //  Handle [a,0] where a < a*
            { ck = rfft + a;
              cj = rfft + b;

              { double kr = ck->real;    //  Forward xform for 0 row:  k = [a,0]  j = [b,0]
                double ki = ck->imag;    //    Ck = F[a,0] = .5 ((Ck+Cj*) - i(Ck-Cj*))
                double jr = cj->real;    //    Cj = F[a,M] = .5 ((Ck+Cj*) + i(Ck-Cj*))
                double ji = cj->imag;    //  Can do this as F[b,0] = F[a*,0] and F[b,M] = F[a*,M]

                double f0r = .5 * (kr + jr);
                double f0i = .5 * (ki - ji);
                double f1r = .5 * (ji + ki);
                double f1i = .5 * (jr - kr);

                ck->real =   f0r + f1r;
                ck->imag =   f0i + f1i;
                cj->real =   f0r - f1r;
                cj->imag = -(f0i - f1i);
              }

              ce  = ck + M;            //  Let x run from 1 to M-1 in the loops below
              cj += (M-1);
              ck += 1;
            }
          else if (a == b)            //  Handle [a,0] where a = a*
            { ck = rfft + a;

              { double kr = ck->real;   //  Ck.real = F_0,a (is real) = Ck.real + Ck.imag
                double ki = ck->imag;   //  Ck.imag = F_M,a (is real) = Ck.real - Ck.imag

                ck->real = kr + ki;
                ck->imag = kr - ki;
              }

              ce  = ck + M/2;         //  a = a* so let x run from 1 to M/2-1 in the loop below
              cj  = ck + (M-1);
              ck += 1;
            }
          else
            continue;

          //  Handle [a,x] where x > 0 and a <= a* in the remainder

          while (ck < ce)
#ifndef __SSE3__
            { double kr = ck->real;    //  Forward xform:  k = [x,a]  j = [M-x,b]
              double ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) - i(w^k)(Ck-Cj*))
              double jr = cj->real;    //    Cj = .5 ((Cj+Ck*) + i(w^k)*(Cj-Ck*))
              double ji = cj->imag;

              double f0r = .5 * (kr + jr);
              double f0i = .5 * (ki - ji);
              double f1r = .5 * (ji + ki);
              double f1i = .5 * (jr - kr);

              double ur = t->pow2.real;
              double ui = t->pow2.imag;

              kr = ur*f1r - ui*f1i;
              ki = ur*f1i + ui*f1r;

              ck->real =   f0r + kr;
              ck->imag =   f0i + ki;
              cj->real =   f0r - kr;
              cj->imag = -(f0i - ki);
#else
            { __m128d fk, fj;
              __m128d f0, f1;
              __m128d ui, ur;

              fk = LOADP(ck);
              fj = LOADP(cj);

              f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
              f1 = MULP(Half_d,ADDSUB(fj,fk));

              ur = LOADUP(&(t->pow2.real));
              ui = LOADUP(&(t->pow2.imag));

              fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
              STOREP(ck,ADDP(f0,fk));
              STOREP(cj,XORP(Conjugate_d,SUBP(f0,fk)));
#endif
              ck += 1;
              cj -= 1;
              t  += A;
            }
        }
      }
  }

  return (rfft);
} 

  //    The inverse transform, Real_FFT_Inverse_nd, takes a complex half-matrix as produced by
  //    Real_FFT_nd, and produces a real-valued result *in-place*.  That is, the pointer returned is
  //    exactly rfft (coerced to be double *).  Note carefully that dims[0] is the length of the
  //    0th dimension of the result double array and twice that of the 0th dimension of rfft.

double *Real_FFT_Inverse_nd(int in_ndim, int *in_dims, Complex_d *rfft)
{ int   M, A, *dims;
  int   ndim;
  int64 base[MAX_REAL_DIM];
  int   cntr[MAX_REAL_DIM];
  int64 size, nsub, offs;

  if (argcheck_nd(in_ndim,in_dims,"Real_FFT_Inverse_nd"))
    return (NULL);

  { int i;

    for (i = 0; i < in_ndim; i++)
      if (in_dims[i] > 1)
        break;
    if (i == in_ndim)
      return ((double *) rfft);
    dims = in_dims+i;
    ndim = in_ndim-i;
  }

  //  Comments are sparse because the code is a mirror image of Real_FFT.  First one takes
  //    the packed representation C of the FFT F, and inverts F to get back to the complex
  //    matrix H that then has the inverse FFT applied to it.  The essential relationship
  //    for this inversion is that:
  //
  //      H[a,x] = .5*(F[a,x] + F[a*,x*]*) + .5*i*w^-x*(F[a,x] - F[a*,x*]*)
  //
  //    Note that this is nearly identical to the forward transform save for a couple of
  //    sign inversions (i.e. between the two major parts and the power of w).

  M = dims[0] / 2;
  A = NINETY / M;

  { int d;

    size = M;
    nsub = 0;
    offs = 0;
    for (d = 1; d < ndim; d++)
      { size   *= dims[d];
        nsub   += size;
        base[d] = size;
        cntr[d] = dims[d]-1;
        if (dims[d] == 1)
          offs += size;
      }
  }

  { int64 a, b;

    for (a = 0; a < size; a += M)
      { { int d;

          for (d = 1; d < ndim; d++)
            if (dims[d] > 1)
              { cntr[d] += 1;
                if (cntr[d] == 1)
                  { offs -= base[d];
                    break;
                  }
                else if (cntr[d] >= dims[d])
                  { cntr[d] = 0;
                    offs += base[d];
                  }
                else
                  break;
              }
          b = nsub - (offs + a);
        }

        { Trig_Block_d *t  = Coot_d + A;
          Complex_d    *ck, *cj, *ce;

          if (a < b)
            { ck = rfft + a;
              cj = rfft + b;

              { double kr = ck->real;
                double ki = ck->imag;
                double jr = cj->real;
                double ji = cj->imag;

                double f0r = .5 * (kr + jr);
                double f0i = .5 * (ki - ji);
                double f1r = .5 * (ji + ki);
                double f1i = .5 * (jr - kr);

                ck->real =   f0r - f1r;
                ck->imag =   f0i - f1i;
                cj->real =   f0r + f1r;
                cj->imag = -(f0i + f1i);
              }

              ce  = ck + M;
              ck += 1;
              cj += (M-1);
            }
          else if (a == b)
            { ck = rfft + a;

              { double kr = ck->real;
                double ki = ck->imag;

                ck->real = .5*(kr + ki);
                ck->imag = .5*(kr - ki);
              }

              ce  = ck + M/2;
              cj  = ck + (M-1);
              ck += 1;
            }
          else
            continue;

          while (ck < ce)
#ifndef __SSE3__
            { double kr = ck->real;
              double ki = ck->imag;
              double jr = cj->real;
              double ji = cj->imag;

              double f0r = .5 * (kr + jr);
              double f0i = .5 * (ki - ji);
              double f1r = .5 * (ji + ki);
              double f1i = .5 * (jr - kr);

              double ur = t->pow2.real;
              double ui = t->pow2.imag;

              kr = ur*f1r - ui*f1i;
              ki = ur*f1i + ui*f1r;

              ck->real =   f0r - kr;
              ck->imag =   f0i - ki;
              cj->real =   f0r + kr;
              cj->imag = -(f0i + ki);
#else
            { __m128d fk, fj;
              __m128d f0, f1;
              __m128d ui, ur;

              fk = LOADP(ck);
              fj = LOADP(cj);

              f0 = MULP(Half_d,ADDSUB(fk,XORP(fj,Negate_d)));
              f1 = MULP(Half_d,ADDSUB(fj,fk));

              ur = LOADUP(&(t->pow2.real));
              ui = LOADUP(&(t->pow2.imag));

              fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
              STOREP(ck,SUBP(f0,fk));
              STOREP(cj,XORP(Conjugate_d,ADDP(f0,fk)));
#endif
              ck += 1;
              cj -= 1;
              t  += A;
            }
        }
      }
  }

  dims[0] = M;
  if (FFT_nd(ndim,dims,rfft,1) == NULL)
    { if (my_message())
        sprintf(FFT_Esource,"Real_FFT_Inverse_nd");
      dims[0] = M * 2;
      return (NULL);
    }
  dims[0] = M * 2;

  return ((double *) rfft);
} 


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL CONVOLUTION & CORRELATION: argcheck_cd                                  *
 *                                                                                            *
\*********************************************************************************************/

  //  Complex_Convolution/Correlation_nd performs a convolution/correlation in the spectral domain
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //      FFT_nd(k,dims,Complex_Convolution_nd(k,dims,FFT_nd(n,data1,0),FFT_nd(k,dims,data2,0)),1);

Complex_d *Complex_Convolution_nd(int ndim, int *dims, Complex_d *fft1, Complex_d *fft2)
{ int   d;
  int64 size;

  size = 1;
  for (d = 0; d < ndim; d++)
    size *= dims[d];
  return (complex_convolution_1d(size,fft1,fft2));
}

Complex_d *Complex_Correlation_nd(int ndim, int *dims, Complex_d *fft1, Complex_d *fft2)
{ int   d;
  int64 size;

  size = 1;
  for (d = 0; d < ndim; d++)
    size   *= dims[d];
  return (complex_correlation_1d(size,fft1,fft2));
}

  //  Common argument checks for real multi-dim. convolution & correlation routines

static int argcheck_cd(int ndim, int *dims, char *source)
{ if (ndim > MAX_REAL_DIM)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"Array cannot have more than %d dimensions",MAX_REAL_DIM);
      return (1);
    }
  if (dims[0]%2 != 0 && dims[0] != 1)
    { if (grab_message()) return (1);
      strcpy(FFT_Esource,source);
      sprintf(FFT_Estring,"dims[0] = %d must be even",dims[0]);
      return (1);
    }
  return (0);
}

  //  Real_Convolution/Correlation_nd performs a convolution/correlation in the spectral domain
  //    of the packed, half-space complex encodings of the fft of real arrays.  Note carefully
  //    that the smallest dimension of rfft1 and rfft2 is dims[0]/2, and the result, produced
  //    in place in rfft1 is returned as a pointer to a double array whose smallest
  //    dimension is dims[0].  rfft1 and rfft2 can be the same array.  Only Real_Convolution_nd
  //    is commented, the correlation routine is identical save that one multiplies the conjugate
  //    of the second array's elements.

Complex_d *Real_Convolution_nd(int in_ndim, int *in_dims, Complex_d *rfft1, Complex_d *rfft2)
{ int64      size;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int       *dims;
  int        ndim;
  Complex_d *s, *f, *g, *h;

  if (argcheck_cd(in_ndim,in_dims,"Real_Convolution_nd"))
    return (NULL);

  { int i;       //  If lowest dims are 1, then effectively ignore them (by modifying ndim and dims)
                 //                     in what follows
    for (i = 0; i < in_ndim; i++)
      if (in_dims[i] > 1)
        break;
    if (i == in_ndim)
      { rfft1->real *= rfft2->real;
        return (rfft1);
      }
    dims = in_dims+i;
    ndim = in_ndim-i;
  }

  //  rfft1 and rfft2 are packed encodings as described in the comments for Real_FFT_nd.
  //  To compute the packed encoding of the convolution, it suffices to multiply all terms
  //  pairwise.  For elements C[a,0] for which a = a* we have to be careful as really two
  //  real valued terms are encoded in the one complex element, where as for all other elements
  //  of C a straight complex multiply is in order.  Note that there are 2^k values of a for
  //  which a = a* and they are a in {0,nk/2} x {0,nk-1/2} x ... x {0,n1/2}!  We use a counter

  { int d;               //  { cntr[d] }_d encodes the current special value that will be next
                         //     i.e., a = PI_d=1^k base[d] * cntr[d]
    size = dims[0]/2;    //        where base[d] = dims[d]/2 * PI_j=1^d-1 dims[j] * dims[0]/2
    for (d = 1; d < ndim; d++)
      { base[d] = (dims[d]/2) * size;
        cntr[d] = 0;
        size *= dims[d];
      }
  }

  g = rfft2;
  h = rfft2 + size;
  s = rfft1;
  f = rfft2;
  do
    { { int d;         //  At a 2 real-valued element a (pointed at by g&f in rfft2 and s in rfft1)
                       //  Advance g to the next special element (or the end of rfft2 if last one)
        for (d = 1; d < ndim; d++)
          if (dims[d] > 1)
            { if (cntr[d] == 1)
                { g -= base[d];
                  cntr[d]  = 0;
                }
              else
                { g += base[d];
                  cntr[d]  = 1;
                  break;
                }
            }
          if (d >= ndim)
            g = h;
      }

      s->real *= f->real;      //  Handle the special element
      s->imag *= f->imag;

      s += 1;
      f += 1;

      while (f < g)            //  Then do term-wise products up to the next one (in g)
#ifdef __SSE3__
        { __m128d kr = LOADUP(&(s->real));
          __m128d ki = LOADUP(&(s->imag));
          __m128d  j = LOADP(f);
          __m128d  p;

          CMULTIPLY(p,j,kr,ki);
          STOREP(s,p);
#else
        { double kr = s->real;
          double ki = s->imag;
          double jr = f->real;
          double ji = f->imag;

          s->real = kr*jr - ki*ji;
          s->imag = ki*jr + kr*ji;
#endif
          f += 1;
          s += 1;
        }
    }
  while (g < h);

  return (rfft1);
}

Complex_d *Real_Correlation_nd(int in_ndim, int *in_dims, Complex_d *rfft1, Complex_d *rfft2)
{ int64      size;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int       *dims;
  int        ndim;
  Complex_d *s, *f, *g, *h;
#ifdef __SSE3__
  __m128d    conj = LOADP(&Conjugate_d);
#endif

  if (argcheck_cd(in_ndim,in_dims,"Real_Correlation_nd"))
    return (NULL);

  { int i;

    for (i = 0; i < in_ndim; i++)
      if (in_dims[i] > 1)
        break;
    if (i == in_ndim)
      { rfft1->real *= rfft2->real;
        return (rfft1);
      }
    dims = in_dims+i;
    ndim = in_ndim-i;
  }

  { int d;

    size = dims[0]/2;
    for (d = 1; d < ndim; d++)
      { base[d] = (dims[d]/2) * size;
        cntr[d] = 0;
        size *= dims[d];
      }
  }

  g = rfft2;
  h = rfft2 + size;
  s = rfft1;
  f = rfft2;
  do
    { { int d;

        for (d = 1; d <= ndim; d++)
          if (dims[d] > 1)
            { if (cntr[d] == 1)
                { g -= base[d];
                  cntr[d]  = 0;
                }
              else
                { g += base[d];
                  cntr[d]  = 1;
                  break;
                }
            }
        if (d >= ndim)
          g = h;
      }

      s->real *= f->real;
      s->imag *= f->imag;

      s += 1;
      f += 1;

      while (f < g)
#ifdef __SSE3__
        { __m128d kr = LOADUP(&(s->real));
          __m128d ki = LOADUP(&(s->imag));
          __m128d  j = XORP(LOADP(f),conj);
          __m128d  p;

          CMULTIPLY(p,j,kr,ki);
          STOREP(s,p);
#else
        { double kr = s->real;
          double ki = s->imag;
          double jr = f->real;
          double ji = f->imag;

          s->real = kr*jr + ki*ji;
          s->imag = ki*jr - kr*ji;
#endif
          f += 1;
          s += 1;
        }
    }
  while (g < h);

  return (rfft1);
}


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL NORMALIZATION: normalize[_0]                                            *
 *                                                                                            *
\*********************************************************************************************/

  //  Normalize_nd normalizes the scores in the ndim-dimensional spatial correlation
  //    matrix cdata that is presumed to have been produced by applying FFT routines to copies
  //    of idata and tdata that had been 0-padded to be of dimensions cdims.  Note that
  //    idims[x]+tdims[x] <= cdims[x] or the arrays were insufficiently padded to give every
  //    overlap correlation.  The normalization takes effect directly on cdata and for
  //    convenience a pointer to it is returned.
  //
  //  Please read the comments prefacing Normalize_1d to understand what sums
  //    need to be produced and how they are used.  In this scenario R and S are multi-dimensional
  //    and develop progressively in each dimension.  Sums are accumulated at each dimensional
  //    level k from n-1 down to 1 recursively.  At level k, normalize_d(ms,idx,area) will be
  //    called with every value of idx representing the partial coordinate (i_n-1,...,i_k)
  //    where n = ndims and i_x is in [0,cdims[x]-1] foreach x in [k,n-1].
  //    idx = |(i_n-1,...,i_k)|_cdims specifically has the value:
  //
  //          ((i_n-1 * cdims[n-2] + i_n-2) * cdims[n-3] + ... i_k+1) * cdims[k] + i_k
  //
  //    which is in the range from 0 to (PI_x=k^n-1 cdims[x])-1.  When called with a particular
  //    partial coordinate idx the following are true:
  //
  //      ms->ndim   = k
  //      ms->im_len = PI_x=0^k-1 idims[x]
  //      ms->tm_len = PI_x=0^k-1 tdims[x]
  //      area       = PI_x=k^n-1 I_x = PI_x=k^n-1 T_x
  //
  //    where I_x = [b_x,e_x] intersect [0,idims[x]] and
  //          T_x = [f_x,c_x] intersect [0,tdims[x]] as in the 1D case
  //
  //    where b_x = i_x, e_x = (b_x + idims[x]) mod cdims[x] and
  //          f_x = -i_x mod cdims[x], c_x = (f_x + tdims[x]) mod cdims[x]
  //
  //      for u in I_n-1 x ... x I_k and v in [0,ms->im_len-1]
  //          ms->im_ave[v] = SUM_u idata[|u|_idims*ms->im_len + v]
  //          ms->im_std[v] = SUM_u idata[|u|_idims*ms->im_len + v] ^ 2
  //
  //      for u in T_n-1 x ... x T_k and v in [0,ms->tm_len-1]
  //          ms->tm_ave[v] = SUM_u tdata[|u|_tdims*ms->tm_len + v]
  //          ms->tm_std[v] = SUM_u tdata[|u|_tdims*ms->tm_len + v] ^ 2

typedef struct
  { int     ndim;
    int64   im_len,  tm_len;
    double *im_ave, *tm_ave;
    double *im_std, *tm_std;
  } Partial_Sums_d;

typedef struct
  { int    *idims;
    int    *tdims;
    int    *cdims;
    double *cdata;
  } Partial_Args_d;

  //  The three routines below each compute progressive overlap sums accross a given dimension
  //    each being tuned to a different scenario.  Normalize_nd starts the recursion
  //    and computes tables of sums and sums of squares.  normalize_d computes intermediate sums
  //    and recurses downward to the next level.  normalize_0_d handles the deepest level of the
  //    recursion where terms of data are normalized.

static void normalize_0_d(Partial_Args_d *ag, Partial_Sums_d *ms, int64 idx, int64 area)
{ double IA, ID;
  double TA, TD;

  double *isums, *idevs;
  double *tsums, *tdevs;
  int64   b, e, c, f;
  int64   w, t;

  int id = ag->idims[0];
  int td = ag->tdims[0];
  int cd = ag->cdims[0];

  double *data = ag->cdata + idx*cd;

  isums = ms->im_ave;
  idevs = ms->im_std;
  tsums = ms->tm_ave;
  tdevs = ms->tm_std;

  { double *ai = isums;
    double *di = idevs;
    double *at = tsums + td;
    double *dt = tdevs + td;

    w = (id + td) - cd;
    if (w < 0)
      w = 0;

    IA = ID = TA = TD = 0.;

    for (b = 0; b < w; b++)
      { IA += *ai++;
        ID += *di++;
        TA += *--at;
        TD += *--dt;
      }

    w *= area;
  }

  c = cd-1;
  if (id == cd)
    { f = c;
      b = 0;
    }
  else
    { f = c-id;
      b = id;
    }
  e = b + td;
  if (e >= cd)
    e -= cd;

  t = b;
  do
    {
      { double ia = IA/w;    //  normalize data[x] (vs. recursing to deeper levels)
        double ta = TA/w;
        double dn = sqrt((ID/w - ia*ia)*(TD/w - ta*ta));
#ifdef DEBUG_NORMALIZE
        printf("      (%lld)%lld: %g %llu %g %g %g %g\n",idx*cd+b,b,data[b],w,IA,TA,ID,TD);
#endif
        if (dn == 0.)
          data[b] = 0.;
        else
          data[b] = ((data[b] - ia*TA) / ((w-1)*dn));
      }

      if (b < id)
        { IA -= isums[b];
          ID -= idevs[b];
          w -= area;
        }
      if (++b >= cd)
        b = 0;

      if (e < id)
        { IA += isums[e];
          ID += idevs[e];
          w += area;
        }
      if (++e >= cd)
        e = 0;

      if (c < td)
        { TA -= tsums[c];
          TD -= tdevs[c];
        }
      if (c-- <= 0)
        c = cd-1;

      if (f < td)
        { TA += tsums[f];
          TD += tdevs[f];
        }
      if (f-- <= 0)
        f = cd-1;
    }
  while (b != t);
}

static void normalize_d(Partial_Args_d *ag, Partial_Sums_d *ms, int64 idx, int64 area)
{ Partial_Sums_d ps;

  int64   in,  tn;
  double *IA, *ID;
  double *TA, *TD;

  double *isums, *idevs;
  double *tsums, *tdevs;
  int64   b, e, c, f;
  int64   w, i, t;

  int ndim = ms->ndim-1;

  int id = ag->idims[ndim];
  int td = ag->tdims[ndim];
  int cd = ag->cdims[ndim];

  ps.ndim   = ndim;                 //  Prepare parameter packet ps and array sum ptrs.
  ps.im_len = in = ms->im_len/id;
  ps.tm_len = tn = ms->tm_len/td;

  isums = ms->im_ave;
  idevs = ms->im_std;
  tsums = ms->tm_ave;
  tdevs = ms->tm_std;

  ps.im_ave = IA = tdevs + ms->tm_len;
  ps.im_std = ID = IA + in;
  ps.tm_ave = TA = ID + in;
  ps.tm_std = TD = TA + tn;

  { double *ai = isums;
    double *di = idevs;
    double *at = tsums + td*tn;
    double *dt = tdevs + td*tn;

    w = (id + td) - cd;
    if (w < 0)
      w = 0;

    for (e = 0; e < in; e++)
      IA[e] = ID[e] = 0.;
    for (e = 0; e < tn; e++)
      TA[e] = TD[e] = 0.;

    for (b = 0; b < w; b++)
      { for (e = 0; e < in; e++)
          { IA[e] += *ai++;
            ID[e] += *di++;
          }
        for (e = tn-1; e >= 0; e--)
          { TA[e] += *--at;
            TD[e] += *--dt;
          }
      }

    w *= area;
  }

  c = cd-1;
  if (id == cd)
    { f = c;
      b = 0;
    }
  else
    { f = c-id;
      b = id;
    }
  e = b + td;
  if (e >= cd)
    e -= cd;

  t = b;
  do
    {
#ifdef DEBUG_NORMALIZE
      { int uu;
        printf("  %lld: %llu:",b,w);
        for (uu = 0; uu < in; uu++)
          printf(" %g",IA[uu]);
        printf(" :");
        for (uu = 0; uu < tn; uu++)
          printf(" %g",TA[uu]);
        printf("\n");
      }
#endif

      if (ndim > 1)
        normalize_d( ag, &ps, idx*cd + b, w );
      else
        normalize_0_d( ag, &ps, idx*cd + b, w );

      if (b < id)
        { double *ai = isums + b*in;
          double *ad = idevs + b*in;
          for (i = 0; i < in; i++)
            { IA[i] -= *ai++;
              ID[i] -= *ad++;
            }
          w -= area;
        }
      if (++b >= cd)
        b = 0;

      if (e < id)
        { double *ai = isums + e*in;
          double *ad = idevs + e*in;
          for (i = 0; i < in; i++)
            { IA[i] += *ai++;
              ID[i] += *ad++;
            }
          w += area;
        }
      if (++e >= cd)
        e = 0;

      if (c < td)
        { double *ai = tsums + c*tn;
          double *ad = tdevs + c*tn;
          for (i = 0; i < tn; i++)
            { TA[i] -= *ai++;
              TD[i] -= *ad++;
            }
        }
      if (c-- <= 0)
        c = cd-1;

      if (f < td)
        { double *ai = tsums + f*tn;
          double *ad = tdevs + f*tn;
          for (i = 0; i < tn; i++)
            { TA[i] += *ai++;
              TD[i] += *ad++;
            }
        }
      if (f-- <= 0)
        f = cd-1;
    }
  while (b != t);
}

double *Normalize_nd(int ndim, int *idims, double *idata,
                               int *tdims, double *tdata,
                               int *cdims, double *cdata)
{ double *cumulative;
  int64   isize, tsize, space;

  if (ndim == 1)
    return (Normalize_1d(idims[0],idata,tdims[0],tdata,cdims[0],cdata));

  //  At each level k, 2*(PI_x=0^k-1 idims[x] + PI_x=0^k-1 tdims[x]) floating point
  //    values are needed for storing partial sums.  Compute the total required over
  //    all ndim-1 levels of the recursion and allocate a big enough vector.

  { int d;

    isize = 1;
    tsize = 1;
    space = 0;
    for (d = 0; d < ndim; d++)
      { if (idims[d] > cdims[d] || tdims[d] > cdims[d])
          { if (grab_message()) return (NULL);
            sprintf(FFT_Esource,"Normalize_nd");
            sprintf(FFT_Estring,"Operand %d-dimension > correlation dimension",d);
            return (NULL);
          }
        space += isize + tsize;
        isize *= idims[d];
        tsize *= tdims[d];
      }

    cumulative = (double *) malloc(sizeof(double)*2*((size_t) space));
    if (cumulative == NULL)
      { if (grab_message()) return (NULL);
        sprintf(FFT_Esource,"Normalize_nd");
        sprintf(FFT_Estring,"Out of memory");
        return (NULL);
      }
  }

  ndim -= 1;

  { Partial_Sums_d ps;
    Partial_Args_d ag;

    int64   in,  tn;
    double *IA, *ID;
    double *TA, *TD;

    int64  b, e, c, f;
    int64  w, i, t;

    int id = idims[ndim];
    int td = tdims[ndim];
    int cd = cdims[ndim];

    ag.idims = idims;     //  Establish global parameters
    ag.tdims = tdims;
    ag.cdims = cdims;
    ag.cdata = cdata;

    ps.ndim   = ndim;
    ps.im_len = in = isize/id;
    ps.tm_len = tn = tsize/td;

    ps.im_ave = IA = cumulative;
    ps.im_std = ID = IA + in;
    ps.tm_ave = TA = ID + in;
    ps.tm_std = TD = TA + tn;

    { double *ai = idata;
      double *at = tdata + tsize;

      w = (id + td) - cd;
      if (w < 0)
        w = 0;

      for (i = 0; i < in; i++)
        IA[i] = ID[i] = 0.;
      for (i = 0; i < tn; i++)
        TA[i] = TD[i] = 0.;

      for (b = 0; b < w; b++)
        { for (i = 0; i < in; i++)
            { double s = *ai++;
              IA[i] += s;
              ID[i] += s*s;
            }
          for (i = tn-1; i >= 0; i--)
            { double s = *--at;
              TA[i] += s;
              TD[i] += s*s;
            }
        }
    }

    c = cd-1;
    if (id == cd)
      { f = c;
        b = 0;
      }
    else
      { f = c-id;
        b = id;
      }
    e = b + td;
    if (e >= cd)
      e -= cd;

    t = b;
    do
      {
#ifdef DEBUG_NORMALIZE
        { int uu;
          printf("  %lld: %llu:",b,w);
          for (uu = 0; uu < in; uu++)
            printf(" %g",IA[uu]);
          printf(" :");
          for (uu = 0; uu < tn; uu++)
            printf(" %g",TA[uu]);
          printf("\n");
        }
#endif

        if (ndim > 1)
          normalize_d( &ag, &ps, b, w );
        else
          normalize_0_d( &ag, &ps, b, w );

        if (b < id)
          { double *ai = idata + b*in;
            for (i = 0; i < in; i++)
              { double s = *ai++;
                IA[i] -= s;
                ID[i] -= s*s;
              }
            w -= 1;
          }
        if (++b >= cd)
          b = 0;
  
        if (e < id)
          { double *ai = idata + e*in;
            for (i = 0; i < in; i++)
              { double s = *ai++;
                IA[i] += s;
                ID[i] += s*s;
              }
            w += 1;
          }
        if (++e >= cd)
          e = 0;
  
        if (c < td)
          { double *ai = tdata + c*tn;
            for (i = 0; i < tn; i++)
              { double s = *ai++;
                TA[i] -= s;
                TD[i] -= s*s;
              }
          }
        if (c-- <= 0)
          c = cd-1;
  
        if (f < td)
          { double *ai = tdata + f*tn;
            for (i = 0; i < tn; i++)
              { double s = *ai++;
                TA[i] += s;
                TD[i] += s*s;
              }
          }
        if (f-- <= 0)
          f = cd-1;
      }
    while (b != t);
  }

  free(cumulative);

  return (cdata);
}

//  Must undefine all macros in order that concatenating fft.F.c and
//    fft.D.c into one file works.

#undef LOADP
#undef LOADUP
#undef SETP

#undef STOREP

#undef SHUFFLE
#undef ADDP
#undef SUBP
#undef MULP
#undef XORP
#undef ADDSUB
#undef CMULTIPLY
#undef NEXT_ROOT

#undef L2_CACHE     
#undef NINETY         
#undef MAX_REAL_DIM    
