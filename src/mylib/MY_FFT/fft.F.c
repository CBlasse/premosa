/*****************************************************************************************\
*                                                                                         *
*  Float Point FFT library                                                                *
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

#include "fft.F.h"

#undef DEBUG_NORMALIZE

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
      * for the SSE3 variant solve 2 line fft's (adjacent in the 1st dimension) at a time
           as this better uses the MMX registers.
      * currrently only wrote SSE3 code for the case were each dim can use the table of
           roots of unity ==> 4K max for any dimension, (8K for 1st dim of a real array)
  */

//  Critical limiting constants:

#define L2_CACHE     262144        //  Problems bigger than this get cross cut
#define NINETY         1024        //  Ninety degree ticks (problems > 4*NINETY are not tabled)
#define MAX_REAL_DIM    128        //  Maximum dimensionality of a real multi-dim. fft

//  Common code between fft.F.c and fft.D.c, use FFT_PROLOGUE to insure they get declared
//    only once in a file that is the concatenation (in either order) of the two files.

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
  printf("%u",k%dims[0]);
}

#ifdef _MSC_VER

#pragma warning( disable:4996 )   //  Turn of deprecation warnings on Windows

#include <windows.h>

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
{ int ret;
  pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error)
    ret = 1;
  else
    { FFT_Error   = 1;
      FFT_Ethread = pthread_tag();
      ret = 0;
    }
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (ret);
}

static int my_message()
{ int mine;
  pthread_mutex_lock(&FFT_Err_Mutex);
  mine = (FFT_Error && pthread_is_this(FFT_Ethread));
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (mine);
}

char *FFT_Error_Source()
{ char *ret;
  pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error && pthread_is_this(FFT_Ethread))
    ret = FFT_Esource;
  else
    ret = NULL;
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (ret);
}

char *FFT_Error_String()
{ char *ret;
  pthread_mutex_lock(&FFT_Err_Mutex);
  if (FFT_Error && pthread_is_this(FFT_Ethread))
    ret = FFT_Estring;
  else
    ret = NULL;
  pthread_mutex_unlock(&FFT_Err_Mutex);
  return (ret);
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

#endif  //  FFT_PROLOGUE


/*********************************************************************************************\
 *                                                                                            *
 *  PROLOGUE AND UTILITIES:                                                                   *
 *     sse3 macros, powers of 2 and error strings, print routines, roots of unity tables      *
 *                                                                                            *
\*********************************************************************************************/

  //  Includes and macros for intel intrinsics if present (-msse3 option)
  //     A "__m128" contains 4 floating point values that are operated on by intrinsic functions.
  //     In most cases we encode 2 complex #s a+ib and c+id in an __m128 as [a,b,c,d], we call this
  //     normal form.  For two twiddle factors x+iy and u+iv, however, we occasionally encode them
  //     into two __m128's, one containing [x,x,u,u] (the real parts) and the other containing
  //     [y,y,v,v].  We call this twiddle form.  It is easy(er) to do complex multiply of a normal
  //     form times a twiddle form.

#ifdef __SSE3__

#include <pmmintrin.h>

#define LOADP(x)      _mm_load_ps((float *) (x))      //  load x[0], x[1], x[2], x[3]
#define LOADHI(u,x)   _mm_loadl_pi(u,(__m64 *) (x))   //  load x into u[0], u[1]
#define LOADLO(u,x)   _mm_loadh_pi(u,(__m64 *) (x))   //  load x into u[2], u[3]

#define LOADUP(x)     _mm_load1_ps((float *) (x))     //  load x, x, x, x
#define SETP(r,i,s,j) _mm_set_ps(j,s,i,r)             //  load (r,i) and (s,j)

#define STOREP(x,u)   _mm_store_ps((float *) (x),u)   //  store m128 u into x[0],x[1], x[2], x[3]
#define STOREHI(x,u)  _mm_storel_pi((__m64 *) (x),u)  //  store u[0], u[1] into x
#define STORELO(x,u)  _mm_storeh_pi((__m64 *) (x),u)  //  store u[2], u[3] into x

#define SHUFFLE(x)    _mm_shuffle_ps(x,x,177)         //  flip the quad words of m128 x
#define ADDP(x,y)     _mm_add_ps(x,y)                 //  add, sub, mul, div, xor quad words in ||
#define SUBP(x,y)     _mm_sub_ps(x,y)
#define MULP(x,y)     _mm_mul_ps(x,y)
#define XORP(x,y)     _mm_xor_ps(x,y)
#define ADDSUB(x,y)   _mm_addsub_ps(x,y)              //  x0-y0, x1+y1, x2-y2, x3+y3

  //  Multiply d in normal form and (r,i) in twiddle form placing the result in p in normal form.

#define CMULTIPLY(p,d,r,i)              \
 p = MULP(r,d);				\
 p = ADDSUB(p,MULP(i,SHUFFLE(d)));

  //  Multiply a in normal form and (r,i) in twiddle form placing the result back in (r,i).

#define NEXT_ROOT(a,r,i)        	\
  CMULTIPLY(r,a,r,i); 		      	\
  i = _mm_movehdup_ps(r);		\
  r = _mm_moveldup_ps(r);

  //   Handy __m128 constants

static __m128 Conjugate = { +0., -0., +0., -0. };   //  Used for conjugating (xor)
static __m128 FlipReal  = { -0., +0., -0., +0. };   //  Used for anti-conjugating (xor)
static __m128 Negate    = { -0., -0., -0., -0. };   //  Used for negating (xor)
static __m128 Half      = {  .5,  .5, .5, .5 };     //  Used for the real-valued fft's

#endif

  //  Print complex and real-valued 1d and multi-d arrays with a title

void Print_Complex_1f(int n, Complex_f *array, char *title)
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

void Print_Real_1f(int n, float *array, char *title)
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

void Print_Complex_nf(int ndim, int *dims, Complex_f *array, char *title)
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

void Print_Real_nf(int ndim, int *dims, float *array, char *title)
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
  { Complex_f pow1;
    Complex_f pow2;
    Complex_f pow3;
  } Trig_Block;

static Trig_Block Root[NINETY];
static Trig_Block Coot[NINETY];
static int        Trig_Firstime = 1;

static pthread_mutex_t FFT_F_MUTEX = PTHREAD_MUTEX_INITIALIZER;

static void init_trig_tables()
{ int   i;
  float    inc, ang, mul;

  inc = (float) (TPI / (4*NINETY));
  ang = 0.;
  for (i = 0; i < NINETY; i++)
    { mul = ang;
      Coot[i].pow1.real = Root[i].pow1.real = (float) cos(mul);
      Coot[i].pow1.imag = - (Root[i].pow1.imag = (float) sin(mul));
      mul += ang;
      Coot[i].pow2.real = Root[i].pow2.real = (float) cos(mul);
      Coot[i].pow2.imag = - (Root[i].pow2.imag = (float) sin(mul));
      mul += ang;
      Coot[i].pow3.real = Root[i].pow3.real = (float) cos(mul);
      Coot[i].pow3.imag = - (Root[i].pow3.imag = (float) sin(mul));
      ang += inc;
    }
}


/*********************************************************************************************\
 *                                                                                            *
 *  1-D COMPLEX FFT:                                                                          *
 *       radix2[_table], radix4[_table], first2_table, sort2_table, argcheck_1f               *
 *                                                                                            *
\*********************************************************************************************/

  //  Radix 2 fft of span s (using twiddle table):  h = 2s and h | n.
  //    For all j s.t. j mod h < s
  //        D[0] = D[0] + u*D[1]
  //        D[1] = D[0] - u*D[1]
  //    where
  //      D[x] = data[j+x*s], u = d(w)^(j mod s), w_h = e^(-TPI/h)i
  //          and d(x) = x  if table = Root
  //                   = x* if table = Coot

static void radix2_table(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int a, v;

  a = NINETY/s;

#ifdef __SSE3__
  if (n > 2)
    for (v = 0; v < n; v += h)
      { Complex_f *j0 = data+v;
        Complex_f *j1 = j0+s;
        Complex_f *je = j1;

        Trig_Block *t = table;

        __m128 ur = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall

        while (j0 != je)
          { __m128 p, d;
            __m128 ui;

            ur = LOADHI(ur,&(t->pow2));
            t += a;
            ur = LOADLO(ur,&(t->pow2));
            t += a;
            ui = _mm_movehdup_ps(ur);
            ur = _mm_moveldup_ps(ur);

            d = LOADP(j1);
            CMULTIPLY(p,d,ur,ui);

            d = LOADP(j0);
            STOREP(j1,SUBP(d,p));
            STOREP(j0,ADDP(d,p));

            j0 += 2;
            j1 += 2;
          }
      }
  else

#endif

    for (v = 0; v < n; v += h)
      { Complex_f *j0 = data+v;
        Complex_f *j1 = j0+s;
        Complex_f *je = j1;

        Trig_Block *t = table;

        while (j0 != je)
          { float dr, di;
            float pr, pi;
            float ur, ui;

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

            j0 += 1;
            j1 += 1;
            t  += a;
          }
    }
}

  //  Radix 2 fft of span s (not using twiddle table):  h = 2s and h | n.

static void radix2(Complex_f *data, int n, int s, int h, float direct)
{ int    v;
  float  theta = direct/h;
  float  cos0  = (float) cos(theta);
  float  sin0  = (float) sin(theta);

#ifdef __SSE3__
  __m128 ure = SETP(1.,1.,cos0,cos0);
  __m128 umi = SETP(0.,0.,sin0,sin0);

  cos0 = (float) cos(2*theta);
  sin0 = (float) sin(2*theta);

  __m128 ang = SETP(cos0,sin0,cos0,sin0);
#endif

  for (v = 0; v < n; v += h)
    { Complex_f *j0 = data+v;
      Complex_f *j1 = j0+s;
      Complex_f *je = j1;

#ifdef __SSE3__
      __m128   ur, ui;

      ur = ure;    // ur,ui holds w_h^(2u) and w_h^(2u+1) in twiddle form
      ui = umi;    //    after iteration u of the loop

      while (j0 != je)
        { __m128 d, p;

          d = LOADP(j1);
          CMULTIPLY(p,d,ur,ui);

          d = LOADP(j0);
          STOREP(j1,SUBP(d,p));
          STOREP(j0,ADDP(d,p));

          j0 += 2;
          j1 += 2;

          NEXT_ROOT(ang,ur,ui);
        }
#else
      float ur, ui;

      ur = 1.0;    // ur + i*ui holds the w_h^u after iteration u of the loop
      ui = 0.0;

      while (j0 != je)
        { float dr, di;
          float pr, pi;

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
  //          and d(x) = x  if table = Root
  //                   = x* if table = Coot

static void radix4_table(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int v, a;
  int sign = (table == Coot);

#ifdef __SSE3__
  __m128 neg;

  if (sign)
    neg = LOADP(&FlipReal);
  else
    neg = LOADP(&Conjugate);
#endif

  a = NINETY/s; 
  for (v = 0; v < n; v += h)
    { Complex_f *j0 = data + v;
      Complex_f *j1 = j0+s;
      Complex_f *j2 = j1+s;
      Complex_f *j3 = j2+s;
      Complex_f *je = j1;

      Trig_Block *t = table;

#ifdef __SSE3__

      Trig_Block *x = table + a;

      __m128 vr = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall

      while (j0 != je)
        { __m128 vi;
          __m128 d, e;
          __m128 t0, t1, t2, t3;

          vr = LOADHI(vr,&(t->pow1));
          vr = LOADLO(vr,&(x->pow1));
          vi = _mm_movehdup_ps(vr);
          vr = _mm_moveldup_ps(vr);

          d  = LOADP(j2);
          CMULTIPLY(t2,d,vr,vi);

          vr = LOADHI(vr,&(t->pow2));
          vr = LOADLO(vr,&(x->pow2));
          vi = _mm_movehdup_ps(vr);
          vr = _mm_moveldup_ps(vr);

          d  = LOADP(j1);
          CMULTIPLY(t1,d,vr,vi);

          vr = LOADHI(vr,&(t->pow3));
          vr = LOADLO(vr,&(x->pow3));
          vi = _mm_movehdup_ps(vr);
          vr = _mm_moveldup_ps(vr);

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

          j0 += 2;
          j1 += 2;
          j2 += 2;
          j3 += 2;
          t = x+a;
          x = t+a;
        }

#else
      while (j0 != je)
        { float r0, r1, r2, r3;
          float i0, i1, i2, i3;
          float dr, di;
          float vr, vi;

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

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;
          t  += a;
        }
#endif
    }
}

  //  Radix 4 fft of span s (not using twiddle table):  h = 4s and h | n.

static void radix4(Complex_f *data, int n, int s, int h, float direct)
{ int     v;
  float   theta = direct/h;
  float   cos0  = (float) cos(theta);
  float   sin0  = (float) sin(theta);
  int     sign  = (direct < 0);

#ifdef __SSE3__
  __m128 uin  = SETP(1.,0.,cos0,sin0);

  cos0 = (float) cos(2.*theta);
  sin0 = (float) sin(2.*theta);

  __m128 angr = LOADUP(&cos0);
  __m128 angi = LOADUP(&sin0);

  __m128 neg;

  if (sign)
    neg = LOADP(&FlipReal);
  else
    neg = LOADP(&Conjugate);
#endif

  for (v = 0; v < n; v += h)
    { Complex_f *j0 = data + v;
      Complex_f *j1 = j0+s;
      Complex_f *j2 = j1+s;
      Complex_f *j3 = j2+s;
      Complex_f *je = j1;

#ifdef __SSE3__
      __m128   uc = uin;   // uc holds w_h^(2u) and w_h^(2u+1) in normal form
                           //    after iteration u of the loop
      while (j0 != je)
        { __m128 d, e;
          __m128 t0, t1, t2, t3;
          __m128 ui, ur;

          ui = _mm_movehdup_ps(uc);  //  ui,ur = twiddle form of uc
          ur = _mm_moveldup_ps(uc);

          d  = LOADP(j2);
          CMULTIPLY(t2,d,ur,ui);

          NEXT_ROOT(uc,ur,ui);       //   ui,ui = twiddle form of uc^2

          d  = LOADP(j1);
          CMULTIPLY(t1,d,ur,ui);

          NEXT_ROOT(uc,ur,ui);       //   ui,ui = twiddle form of uc^3

          d  = LOADP(j3);
          CMULTIPLY(t3,d,ur,ui);

          d  = MULP(angr,uc);
          uc = ADDSUB(d,MULP(angi,SHUFFLE(uc)));    //  uc *= ang (= [w_h^2,w_h^2])

          t0 = LOADP(j0);

          d  = ADDP(t0,t1);
          e  = ADDP(t2,t3);

          STOREP(j0,ADDP(d,e));
          STOREP(j2,SUBP(d,e));

          d  = SUBP(t0,t1);
          e  = XORP(SHUFFLE(SUBP(t3,t2)),neg);

          STOREP(j1,ADDP(d,e));
          STOREP(j3,SUBP(d,e));

          j0 += 2;
          j1 += 2;
          j2 += 2;
          j3 += 2;
        }
#else
      float ur, ui;

      ur = 1.0;    // ur + i*ui holds the w_h^u after iteration u of the loop
      ui = 0.0;

      while (j0 != je)
        { float r0, r1, r2, r3;
          float i0, i1, i2, i3;
          float dr, di;
          float vr, vi;

          dr = j2->real;
          di = j2->imag;

          r2 = ur * dr - ui * di;
          i2 = ur * di + ui * dr;

          vr = ur * ur - ui * ui;
          vi = (float) (2. * ur * ui);

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
#endif
    }
}

  //  Radix 2 auto-sort (combination of fft stage and sorting pass): s < h and s*h = n/2
  //    For all j s.t. j mod 2h < h and j mod 2s < s
  //        D[0,0] = D[0,0] + u*D[1,0]
  //        D[0,1] = D[0,0] - u*D[1,0]
  //        D[1,0] = D[0,1] + u*D[1,1]
  //        D[1,1] = D[0,1] - u*D[1,1]
  //    where
  //      D[x,y] = data[j+x*h+y*s], u = d(w)^(j mod s), w = e^(-TPI/h)i,
  //          and d(x) = x  if table = Root
  //                   = x* if table = Coot

static void sort2_table(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int    h2, s2;
  int    v, w, a;

  h2 = h*2;
  s2 = s*2;
  a  = NINETY/s;
  for (v = 0; v < n; v += h2)
   for (w = 0; w < h; w += s2)
    { Complex_f *j0 = data + (v+w);
      Complex_f *j1 = j0+s;
      Complex_f *j2 = j0+h;
      Complex_f *j3 = j2+s;
      Complex_f *je = j1;

      Trig_Block *t = table;

#ifdef __SSE3__
      __m128 ur = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall

      while (j0 != je)
        { __m128 p, q, d;
          __m128 ui;

          ur = LOADHI(ur,&(t->pow2));
          t += a;
          ur = LOADLO(ur,&(t->pow2));
          t += a;
          ui = _mm_movehdup_ps(ur);
          ur = _mm_moveldup_ps(ur);

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

          j0 += 2;
          j1 += 2;
          j2 += 2;
          j3 += 2;
        }
#else
      while (j0 != je)
        { float dr, di;
          float pr, pi;
          float qr, qi;
          float ur, ui;

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

          j0 += 1;
          j1 += 1;
          j2 += 1;
          j3 += 1;
          t  += a;
        }
#endif
    }
}

  //  Radix 2 auto-sort specialized for s=1, h = n/2

static void first2_table(Complex_f *data, int n)
{ Complex_f *j0 = data;
  Complex_f *j2 = j0 + n/2;
  Complex_f *je = j2;

#ifdef __SSE3__

  while (j0 != je)
    { __m128 p, q, d;

      p = LOADP(j0);
      q = LOADP(j2);

      d = ADDP(p,q);
      q = SUBP(p,q);
      p = d;

      p = _mm_movelh_ps(p,q);
      q = _mm_movehl_ps(q,d);

      STOREP(j0,p);
      STOREP(j2,q);

      j0 += 2;
      j2 += 2;
    }
#else

  Complex_f *j1 = j0+1;
  Complex_f *j3 = j2+1;

  while (j0 != je)
    { float dr, di;
      float pr, pi;
      float qr, qi;

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

      j0 += 2;
      j1 += 2;
      j2 += 2;
      j3 += 2;
    }
#endif
}

  //  Check that n is a power of 2 and that m does not exceed the max problem size

static int argcheck_1f(int n, int m, char *source)
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
  //    for convenience a pointer to data is returned by FFT_1f.  If invert is non-zero then
  //    the inverse Fourier Transform is performed.

Complex_f *FFT_1f(int n, Complex_f *data, int invert)
{ float       direct;
  Trig_Block *table;

  pthread_mutex_lock(&FFT_F_MUTEX);    //  Set up trig tables
  if (Trig_Firstime)
    { Trig_Firstime = 0;
      init_trig_tables();
    }
  pthread_mutex_unlock(&FFT_F_MUTEX);

  if (invert)           //  Establish direction
    { direct = (float) (-TPI);
      table  = Coot;
    }
  else
    { direct = (float)  TPI;
      table  = Root;
    }

  if (argcheck_1f(n,n,"FFT_1f"))    //  Check that n is a power of 2 and not too large
    return (NULL);

  if (n <= 2)           //  For n = 1 nothing to do, for n = 2 a tabled radix 2 possibly
    { if (n == 1)       //    followed by normalization (if invert on)
        return (data);
      radix2_table(data,2,1,2,table);
    }

  else                  //  For n > 2:
    { int i, s, h;      //    Auto sort passes for the first half, radix-4 and maybe one radix-2
      int ss, hh, nn;   //    for the remainder

      if (n < L2_CACHE)    //  If the h-span of a run phases fits in the L2-cache then do
        nn = n;            //     them all in a cross-cut of the passes in order to gain
      else                 //     L2 cache coherence.  Hugely important for big, big n.
        nn = L2_CACHE;

      first2_table(data,n);
      for (s = 2, h = n/4; s < h && L2_CACHE/2 < h; s <<= 1, h >>= 1)
        sort2_table(data,n,s,h,table);

      for (i = 0; i < n; i += L2_CACHE)
        { for (ss = s, hh = h; ss < hh; ss <<= 1, hh >>= 1)
            sort2_table(data+i,nn,ss,hh,table);

          for (hh = ss<<2; hh <= nn; ss <<= 2, hh <<= 2)
            if (ss > NINETY)
              radix4(data+i,nn,ss,hh,direct);
            else
              radix4_table(data+i,nn,ss,hh,table);
        }

      for (s = ss, h = hh; h <= n; s <<= 2, h <<= 2)
        if (s > NINETY)
          radix4(data,n,s,h,direct);
        else
          radix4_table(data,n,s,h,table);

      if (s < n)
        { h = (s << 1);
          if (s > NINETY)
            radix2(data,n,s,h,direct);
          else
            radix2_table(data,n,s,h,table);
        }
    }

  if (invert)                  //  Normalize the inverse transform if invert is on
    { Complex_f *d;
      Complex_f *e = data+n;
      float      v = (float) (1./n);

#ifdef __SSE3__
      __m128 p = LOADUP(&v);

      for (d = data; d < e; d += 2)
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
  //    The forward transform, Real_FFT_1f, takes a float array of length n, and in-place produces 
  //    a Complex_f array C of length n/2 that is the first half of the conjugate symmetric FFT
  //    of the real data with the exception that F_(n/2) (which is real) is tucked into C[0].imag
  //    (which works as F_0 is also real).  Again, the pointer returned by Real_FFT_1f is really the
  //    same as rdata, the FFT is performed in-place. 

Complex_f *Real_FFT_1f(int n, float *rdata)
{ Complex_f *rfft = (Complex_f *) rdata;
  int n2 = n/2;
  int n4 = n/4;

  if (argcheck_1f(n,n2,"Real_FFT_1f")) //  Check that n is a power of 2 and not too large
    return (NULL);

  if (n == 1) return (rfft);

  FFT_1f(n2,rfft,0);

  { float zr = rfft[0].real;   //  Special case for the 0^th and n/2^th elements
    float zi = rfft[0].imag;   //     whose x'forms are real 

    rfft[0].real = (float) (zr + zi);     //  Put F_0 in .real and F_n2 in .imag of 0^th element
    rfft[0].imag = (float) (zr - zi);
  }

  if (n == 2) return (rfft);

  if (n4 > NINETY)               //  Code when n too big for tabled twiddles
    { float theta = (float) (TPI / n);
      float cos0  = (float) cos(theta);
      float sin0  = (float) sin(theta);
      float vr    = cos0;
      float vi    = sin0;

      Complex_f *ck = rfft + 1;
      Complex_f *cj = rfft + (n2-1);

#ifdef __SSE3__
      float cos20 = (float) cos(2.*theta);   //  For SSE3 1st loop will be for k = 2+3.
      float sin20 = (float) sin(2.*theta);   //    So setup twiddles accordingly
      float cos30 = (float) cos(3.*theta);
      float sin30 = (float) sin(3.*theta);

      __m128 ang = SETP(cos20,sin20,cos20,sin20);
      __m128 ur  = SETP(cos20,cos20,cos30,cos30);
      __m128 ui  = SETP(sin20,sin20,sin30,sin30);
      __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
#else
      while (ck < cj)    //  Sneaky: if SSE3 then do the loop once to take care of k=1 case
#endif

        { float kr = ck->real;    //  Forward xform:  Ck = rfft[k]  Cj = rfft[M-k]
          float ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) - i(w^k)(Ck-Cj*))
          float jr = cj->real;    //    Cj = .5 ((Cj+Ck*) + i(w^k)*(Cj-Ck*))
          float ji = cj->imag;    //  where w = e^(-TPI/n)i and k = ck-rfft

          float f0r = (float) (.5 * (kr + jr));
          float f0i = (float) (.5 * (ki - ji));
          float f1r = (float) (.5 * (ji + ki));
          float f1i = (float) (.5 * (jr - kr));

          kr = vr*f1r - vi*f1i;    // (vr,vi) is the n'th root of unity to the power k
          ki = vr*f1i + vi*f1r;

          ck->real =   f0r + kr;
          ck->imag =   f0i + ki;
          cj->real =   f0r - kr;
          cj->imag = -(f0i - ki);

          kr = vr * cos0 - vi * sin0;
          vi = vr * sin0 + vi * cos0;
          vr = kr;

          ck += 1;
          cj -= 1;
        }

#ifdef __SSE3__               //  Sneaky: if not SSE3 then ck >= cj at this point, otherise
      while (ck < cj)         //    ck-rfft is 2!  Remember doing 2 at a time.
        { __m128 fk;
          __m128 f0, f1;

          fk = LOADP(ck);
          fj = LOADHI(fj,cj);
          fj = LOADLO(fj,cj-1);

          f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
          f1 = MULP(Half,ADDSUB(fj,fk));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));   //  [ur,ui] twiddle form of w_n^u and
          STOREP(ck,ADDP(f0,fk));                          //     w_n^(u+1) where u = ck-rfft

          fj = XORP(Conjugate,SUBP(f0,fk));
          STOREHI(cj,fj);
          STORELO(cj-1,fj);

          NEXT_ROOT(ang,ur,ui);

          ck += 2;
          cj -= 2;
        }
#endif
    }

  else                               //  Code when the twiddles can come from a table
    { int         a = NINETY/n4;
      Trig_Block *t = Root + a;

      Complex_f *ck = rfft + 1;
      Complex_f *cj = rfft + (n2-1);

#ifndef __SSE3__
      while (ck < cj)                //  Handle all if not SSE3, handle k = 1 if SSE3
#endif
        { float kr = ck->real;
          float ki = ck->imag;
          float jr = cj->real;
          float ji = cj->imag;

          float f0r = (float) (.5 * (kr + jr));
          float f0i = (float) (.5 * (ki - ji));
          float f1r = (float) (.5 * (ji + ki));
          float f1i = (float) (.5 * (jr - kr));

          float ur = t->pow1.real;
          float ui = t->pow1.imag;

          kr = ur*f1r - ui*f1i;
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r + kr;
          ck->imag =   f0i + ki;
          cj->real =   f0r - kr;
          cj->imag = -(f0i - ki);

          ck += 1;
          cj -= 1;
          t  += a;
        }

#ifdef __SSE3__
      __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
      __m128 ur  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall

      while (ck < cj)              //  No rep's if not SSE3, reps starting from k=2 if SSE3
        { __m128 fk;
          __m128 f0, f1;
          __m128 ui;

          fk = LOADP(ck);
          fj = LOADHI(fj,cj);
          fj = LOADLO(fj,cj-1);

          f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
          f1 = MULP(Half,ADDSUB(fj,fk));

          ur = LOADHI(ur,&(t->pow1));
          t += a;
          ur = LOADLO(ur,&(t->pow1));
          t += a;
          ui = _mm_movehdup_ps(ur);
          ur = _mm_moveldup_ps(ur);

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
          STOREP(ck,ADDP(f0,fk));

          fj = XORP(Conjugate,SUBP(f0,fk));
          STOREHI(cj,fj);
          STORELO(cj-1,fj);

          ck += 2;
          cj -= 2;
        }
#endif
    }

  return (rfft);
} 


  //    The inverse transform, Real_FFT_Inverse_1f, takes a complex half-matrix as produced by
  //    Real_FFT_1f, and produces a real-valued result *in-place*.  That is, the pointer returned is
  //    exactly rfft (coerced to be float *).  Note carefully that n is the length of the
  //    resulting real array and is twice the length of rfft.

float *Real_FFT_Inverse_1f(int n, Complex_f *rfft)
{ int n2 = n/2;
  int n4 = n/4;
 
  if (argcheck_1f(n,n2,"Real_FFT_Inverse_1f")) //  Check that n is a power of 2 and not too large
    return (NULL);

  pthread_mutex_lock(&FFT_F_MUTEX);     //  Set up trig tables
  if (Trig_Firstime)
    { Trig_Firstime = 0;
      init_trig_tables();
    }
  pthread_mutex_unlock(&FFT_F_MUTEX);
 
  if (n == 1) return ((float *) rfft);

  { float zr = rfft[0].real;   //  Special case for reconstituting 0^th element that depends on
    float zi = rfft[0].imag;   //    the 0^th and n/2^th value of the x'form packed in 0^th el.

    rfft[0].real = (float) (.5*(zr + zi));
    rfft[0].imag = (float) (.5*(zr - zi));
  }

  if (n == 2) return ((float *) rfft);

  if (n4 > NINETY)               //  Code when n too big for tabled twiddles
    { float theta = (float) (- TPI / n);
      float cos0  = (float) cos(theta);
      float sin0  = (float) sin(theta);
      float vr    = cos0;
      float vi    = sin0;

      Complex_f *ck = rfft + 1;
      Complex_f *cj = rfft + (n2-1);

#ifdef __SSE3__
      float cos20 = (float) cos(2.*theta);   //  For SSE3 1st loop will be for k = 2+3.
      float sin20 = (float) sin(2.*theta);   //    So setup twiddles accordingly
      float cos30 = (float) cos(3.*theta);
      float sin30 = (float) sin(3.*theta);

      __m128 ang = SETP(cos20,sin20,cos20,sin20);
      __m128 ur  = SETP(cos20,cos20,cos30,cos30);
      __m128 ui  = SETP(sin20,sin20,sin30,sin30);
      __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
#else
      while (ck < cj)    //  Sneaky: if SSE3 then do the loop once to take care of k=1 case
#endif

        { float kr = ck->real;    //  Inverse xform:  Ck = rfft[k]  Cj = rfft[M-k]
          float ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) + i(w^-k)*(Ck-Cj*))
          float jr = cj->real;    //    Cj = .5 ((Cj+Ck*) - i(w^-k)(Cj-Ck*))
          float ji = cj->imag;    //  where w = e^(-TPI/n)i and k = ck-rfft

          float f0r = (float) (.5 * (kr + jr));
          float f0i = (float) (.5 * (ki - ji));
          float f1r = (float) (.5 * (ji + ki));
          float f1i = (float) (.5 * (jr - kr));

          kr = vr*f1r - vi*f1i;    // (vr,vi) is the n'th root of unity to the power -k
          ki = vr*f1i + vi*f1r;

          ck->real =   f0r - kr;
          ck->imag =   f0i - ki;
          cj->real =   f0r + kr;
          cj->imag = -(f0i + ki);

          { float x;

            x  = vr * cos0 - vi * sin0;
            vi = vr * sin0 + vi * cos0;
            vr = x;
          }

          ck += 1;
          cj -= 1;
        }

#ifdef __SSE3__               //  Sneaky: if not SSE3 then ck >= cj at this point, otherise
      while (ck < cj)         //    ck-rfft is 2!  Remember doing 2 at a time.
        { __m128 fk;
          __m128 f0, f1;

          fk = LOADP(ck);
          fj = LOADHI(fj,cj);
          fj = LOADLO(fj,cj-1);

          f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
          f1 = MULP(Half,ADDSUB(fj,fk));

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));   //  [ur,ui] twiddle form of w_n^u* and
          STOREP(ck,SUBP(f0,fk));                          //     w_n^(u+1)* where u = ck-rfft

          fj = XORP(Conjugate,ADDP(f0,fk));
          STOREHI(cj,fj);
          STORELO(cj-1,fj);

          NEXT_ROOT(ang,ur,ui);

          ck += 2;
          cj -= 2;
        }
#endif
    }

  else                               //  Code when the twiddles can come from a table
    { int         a = NINETY/n4;
      Trig_Block *t = Coot + a;

      Complex_f *ck = rfft + 1;
      Complex_f *cj = rfft + (n2-1);

#ifndef __SSE3__
      while (ck < cj)                //  Handle all if not SSE3, handle k = 1 if SSE3
#endif
        { float kr = ck->real;
          float ki = ck->imag;
          float jr = cj->real;
          float ji = cj->imag;

          float f0r = (float) (.5 * (kr + jr));
          float f0i = (float) (.5 * (ki - ji));
          float f1r = (float) (.5 * (ji + ki));
          float f1i = (float) (.5 * (jr - kr));

          float ur = t->pow1.real;
          float ui = t->pow1.imag;

          kr = ur*f1r - ui*f1i;
          ki = ur*f1i + ui*f1r;

          ck->real =   f0r - kr;
          ck->imag =   f0i - ki;
          cj->real =   f0r + kr;
          cj->imag = -(f0i + ki);

          ck += 1;
          cj -= 1;
          t  += a;
        }

#ifdef __SSE3__
      __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
      __m128 ur  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall

      while (ck < cj)              //  No rep's if not SSE3, reps starting from k=2 if SSE3
        { __m128 fk;
          __m128 f0, f1;
          __m128 ui;

          fk = LOADP(ck);
          fj = LOADHI(fj,cj);
          fj = LOADLO(fj,cj-1);

          f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
          f1 = MULP(Half,ADDSUB(fj,fk));

          ur = LOADHI(ur,&(t->pow1));
          t += a;
          ur = LOADLO(ur,&(t->pow1));
          t += a;
          ui = _mm_movehdup_ps(ur);
          ur = _mm_moveldup_ps(ur);

          fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
          STOREP(ck,SUBP(f0,fk));

          fj = XORP(Conjugate,ADDP(f0,fk));
          STOREHI(cj,fj);
          STORELO(cj-1,fj);

          ck += 2;
          cj -= 2;
        }
#endif
    }

  return ((float *) FFT_1f(n2,rfft,1));
} 


/*********************************************************************************************\
 *                                                                                            *
 *  1-D CONVOLUTION & CORRELATION                                                             *
 *                                                                                            *
\*********************************************************************************************/

  //  Complex_Convolution_1f performs the term-wise complex multiplication of fft1 and fft2
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //        FFT_1f(n,Complex_Convolution_1f(n,FFT_1f(n,data1,0),FFT_1f(n,data2,0)),1);

static Complex_f *complex_convolution_1f(int64 n, Complex_f *fft1, Complex_f *fft2)
{ Complex_f *s = fft1;
  Complex_f *f = fft2;
  Complex_f *h = fft2 + n;

#ifdef __SSE3__
  if (n == 1)
#else
  while (f < h)
#endif
    { float kr = s->real;
      float ki = s->imag;
      float jr = f->real;
      float ji = f->imag;

      s->real = kr*jr - ki*ji;
      s->imag = ki*jr + kr*ji;

      s += 1;
      f += 1;
    }

#ifdef __SSE3__
  while (f < h)
    { __m128 kr = LOADP(s);
      __m128 j  = LOADP(f);
      __m128 ki, p;

      ki = _mm_movehdup_ps(kr);
      kr = _mm_moveldup_ps(kr);

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);

      f += 2;
      s += 2;
    }
#endif

  return (fft1);
}

Complex_f *Complex_Convolution_1f(int n, Complex_f *fft1, Complex_f *fft2)
{ return (complex_convolution_1f((int64) n, fft1, fft2)); }

  //  Complex_Correlation_1f performs the term-wise complex multiplication of fft1 and fft2
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //        FFT_1f(n,Complex_Correlation_1f(n,FFT_1f(n,data1,0),FFT_1f(n,data2,0)),1);

static Complex_f *complex_correlation_1f(int64 n, Complex_f *fft1, Complex_f *fft2)
{ Complex_f *s = fft1;
  Complex_f *f = fft2;
  Complex_f *h = fft2 + n;
#ifdef __SSE3__
  __m128     conj = LOADP(&Conjugate);
#endif

#ifdef __SSE3__
  if (n == 1)
#else
  while (f < h)
#endif
    { float kr = s->real;
      float ki = s->imag;
      float jr = f->real;
      float ji = f->imag;

      s->real = kr*jr + ki*ji;
      s->imag = ki*jr - kr*ji;

      s += 1;
      f += 1;
    }

#ifdef __SSE3__
  while (f < h)
    { __m128 kr = LOADP(s);
      __m128 j  = XORP(LOADP(f),conj);
      __m128 ki, p;

      ki = _mm_movehdup_ps(kr);
      kr = _mm_moveldup_ps(kr);

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);

      f += 2;
      s += 2;
    }
#endif

  return (fft1);
}

Complex_f *Complex_Correlation_1f(int n, Complex_f *fft1, Complex_f *fft2)
{ return (complex_correlation_1f((int64) n, fft1, fft2)); }

  //  Real_Convolution_1f effects a term-wise multiplication of the fft of two vectors encoded
  //    as packed half-sized complex vectors of length ** n/2 ** as produced by Real_FFT_1f.  The
  //    result is produced in-place in rfft1 and a floating point pointer to it is returned.
  //    The code works just fine if rfft2 is the same as array rfft1.  So a real convolution of
  //    of two real arrays data1 and data2, in-place within data 1, is accomplised with the code:
  //
  //  Real_FFT_Inverse_1f(n,Real_Convolution_1f(n,Real_FFT_1f(n,data1,0),Real_FFT_1f(n,data2,0)),1);

Complex_f *Real_Convolution_1f(int n, Complex_f *rfft1, Complex_f *rfft2)
{ Complex_f *s = rfft1;
  Complex_f *f = rfft2;
  Complex_f *g = rfft2 + n/2;

  //  Basically all one is doing is element-wise multiplying the FFTs of rfft1 and rfft2
  //    with the subtlety that the 0th element codes 2 real valued terms instead of 1.

  if (n == 1)
    { s->real *= f->real;
      return (rfft1);
    }

  if ((n & 0x1) != 0)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"Real_Convolution_1f");
      sprintf(FFT_Estring,"n = %d must be even",n);
      return (NULL);
    }

  s->real *= f->real;
  s->imag *= f->imag;

  s += 1;
  f += 1;

#ifndef __SSE3__
  while (f < g)
#else
  if (f < g)
#endif
    { float kr = s->real;
      float ki = s->imag;
      float jr = f->real;
      float ji = f->imag;

      s->real = kr*jr - ki*ji;
      s->imag = ki*jr + kr*ji;

      f += 1;
      s += 1;
    }
#ifdef __SSE3__
  while (f < g)
    { __m128 kr = LOADP(s);
      __m128 j  = LOADP(f);
      __m128 ki, p;

      ki = _mm_movehdup_ps(kr);
      kr = _mm_moveldup_ps(kr);

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);

      f += 2;
      s += 2;
    }
#endif

  return (rfft1);
}

  //  Real_Correlation_1f effects a term-wise multiplication of the fft of a vector and the
  //    conjugate of another, both encoded as packed half-sized complex vectors of length ** n/2 **
  //    as produced by Real_FFT_1f.  The result is produced in-place in rfft1 and a floating point
  //    pointer to it is returned.  The code works just fine if rfft2 is the same as array rfft1.
  //    So a real convolution of two real arrays data1 and data2, in-place within data 1, is
  //    accomplised with the code:
  //
  //  Real_FFT_Inverse_1f(n,Real_Correlation_1f(n,Real_FFT_1f(n,data1,0),Real_FFT_1f(n,data2,0)),1);

Complex_f *Real_Correlation_1f(int n, Complex_f *rfft1, Complex_f *rfft2)
{ Complex_f *s = rfft1;
  Complex_f *f = rfft2;
  Complex_f *g = rfft2 + n/2;
#ifdef __SSE3__
  __m128     conj = LOADP(&Conjugate);
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
      sprintf(FFT_Esource,"Real_Correlation_1f");
      sprintf(FFT_Estring,"n must be even");
      return (NULL);
    }

  s->real *= f->real;
  s->imag *= f->imag;

  s += 1;
  f += 1;

#ifndef __SSE3__
  while (f < g)
#else
  if (f < g)
#endif
    { float kr = s->real;
      float ki = s->imag;
      float jr = f->real;
      float ji = f->imag;

      s->real = kr*jr + ki*ji;
      s->imag = ki*jr - kr*ji;

      f += 1;
      s += 1;
    }
#ifdef __SSE3__
  while (f < g)
    { __m128 kr = LOADP(s);
      __m128 j  = XORP(LOADP(f),conj);
      __m128 ki ,p;

      ki = _mm_movehdup_ps(kr);
      kr = _mm_moveldup_ps(kr);

      CMULTIPLY(p,j,kr,ki);
      STOREP(s,p);

      f += 2;
      s += 2;
    }
#endif

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
  //  Normalize_1f normalizes the scores in the spatial correlation vector cdata that
  //    is presumed to have been produced by applying FFT routines to copies of idata
  //    and tdata that had been 0-padded to be of length cdim.  Note that idim+tdim <= cdim or
  //    the vectors were insufficiently padded to give every overlap correlation.  The
  //    normalization takes effect directly on cdata and for convenience a pointer
  //    to it is returned.

float *Normalize_1f(int idim, float *idata, int tdim, float *tdata, int cdim, float *cdata)
{ float IA, ID;
  float TA, TD;

  int    b, e;
  int    f, c;
  int    w, t;

  if (idim > cdim || tdim > cdim)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"Normalize_1f");
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
    { float s = idata[b];
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
          cdata[b] = (float) ((cdata[b] - ia*TA) / ((w-1)*dn));
      }

      if (b < idim)
        { float s = idata[b];
          IA -= s;
          ID -= s*s;
          w  -= 1;
        }
      if (++b >= cdim)
        b = 0;

      if (e < idim)
        { float s = idata[e];
          IA += s;
          ID += s*s;
          w  += 1;
        }
      if (++e >= cdim)
        e = 0;

      if (c < tdim)
        { float s = tdata[c];
          TA -= s;
          TD -= s*s;
        }
      if (c-- <= 0)
        c = cdim-1;

      if (f < tdim)
        { float s = tdata[f];
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
 *       init_multi_fft, fft_a_dim, (radix2|radix4|sort2)_array                               *
 *                                                                                            *
\*********************************************************************************************/

static int64 init_multi_fft(int d, int *dims, int *dmax)
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
                sprintf(FFT_Esource,"FFT_nf");
                sprintf(FFT_Estring,"dims[%d] = %d is not a power of 2",d,n);
                return (0);
              }
            s >>= 1;
          }

        if (n > NINETY*4)
          { if (grab_message()) return (0);
            sprintf(FFT_Esource,"FFT_nf");
            sprintf(FFT_Estring,"dims[%d] is larger than hardcoded maximum of %dK",d,(NINETY)>>8);
            return (0);
          }

        if (n > max)
          max = n;
        volume *= n;
      }
  }

  pthread_mutex_lock(&FFT_F_MUTEX);      //  Initialize the twiddle tables
  if (Trig_Firstime)
    { Trig_Firstime = 0;
      init_trig_tables();
    }
  pthread_mutex_unlock(&FFT_F_MUTEX);

  *dmax = max;
  return (volume);
}

#ifdef __SSE3__

  // With the SSE optimizations on, the line cache contains 2 rows that are
  //   interleaved, i.e., if row 1 is [ a0, a1, a2, ...] and row 2 is [ b0, b1, b2, ...]
  //   then the cache contains [ a0, b0, a1, b1, a2, b2, ...].  So we have to have distinct
  //   codes to perform the FFT passes over the cache as the 1-dimensional versions process
  //   a single row 2-elements at a time.

  // Radix 2 FFT of span s/2 with twiddle tables.  Same function as radix2_table save
  //   that s and h are twice as large and 2 problems are performed in parallel.

static void radix2_array(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int v;
  int a;

  a = ONE_EIGHTY/s;
  for (v = 0; v < n; v += h)
    { Complex_f *j0 = data+v;
      Complex_f *j1 = j0+s;
      Complex_f *je = j1;

      Trig_Block *t = table;

      while (j0 != je)
        { __m128 p, d;
          __m128 ui, ur;

          ur = LOADUP(&(t->pow2.real));
          ui = LOADUP(&(t->pow2.imag));

          d = LOADP(j1);
          CMULTIPLY(p,d,ur,ui);

          d = LOADP(j0);
          STOREP(j1,SUBP(d,p));
          STOREP(j0,ADDP(d,p));

          j0 += 2;
          j1 += 2;
          t  += a;
        }
    }
}

  // Radix 4 FFT of span s/2 with twiddle tables.  Same function as radix4_table save
  //   that s and h are twice as large and 2 problems are performed in parallel.

static void radix4_array(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int v, a;

  __m128 neg;

  if (table == Coot)
    neg = LOADP(&FlipReal);
  else
    neg = LOADP(&Conjugate);

  a = ONE_EIGHTY/s;
  for (v = 0; v < n; v += h)
    { Complex_f *j0 = data + v;
      Complex_f *j1 = j0+s;
      Complex_f *j2 = j1+s;
      Complex_f *j3 = j2+s;
      Complex_f *je = j1;

      Trig_Block *t = table;

      while (j0 != je)
        { __m128 vr, vi;
          __m128 d, e;
          __m128 t0, t1, t2, t3;

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

          j0 += 2;
          j1 += 2;
          j2 += 2;
          j3 += 2;
          t  += a;
        }
    }
}

  // Radix 2 auto-sort.  Same function as sort2_table save that s and h are twice
  //    as large and 2 problems are performed in parallel.

static void sort2_array(Complex_f *data, int n, int s, int h, Trig_Block *table)
{ int    h2, s2;
  int    v, w, a;

  h2 = h*2;
  s2 = s*2;

  a = ONE_EIGHTY/s;
  for (v = 0; v < n; v += h2)
   for (w = 0; w < h; w += s2)
    { Complex_f *j0 = data + (v+w);
      Complex_f *j1 = j0+s;
      Complex_f *j2 = j0+h;
      Complex_f *j3 = j2+s;
      Complex_f *je = j1;

      Trig_Block *t = table;

      while (j0 != je)
        { __m128 p, q, d;
          __m128 ui, ur;

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

          j0 += 2;
          j1 += 2;
          j2 += 2;
          j3 += 2;
          t  += a;
        }
    }
}

#endif

  //  Perform an FFT on a given dimension, other than the first, of the matrix 'data'.  The
  //  matrix contains 'size' elements, there are 'n' > 1 elements in a row of the given dimension,
  //  and those elements are 'step' elements apart in the matrix.  Invert signals the direction
  //  of the transform.

  //  The strategy is to move elements from a row to the supplied cache, performing the
  //  first step of the auto-sort as you do so.  The SSE3 version operates on two rows
  //  that are adjacent in the first dimension.  Once the transform of a row (or two) is
  //  complete than it is moved back into the matrix, being normalized at the same time
  //  if invert is on.

static void fft_a_dim(int64 size, int64 step, int n, int invert,
                      Complex_f *data, Complex_f *cache)
{ int   n2    = n/2;
  int   n4    = n/4;
  int64 ostep = step*n;
  int64 hstep = step*n2;

  int64 i, j; 
  int64 jtop;

  float       direct;
  Trig_Block *table;

  if (invert)          //   Establish direction
    { direct = (float) (-TPI);
      table  = Coot;
    }
  else
    { direct = (float)  TPI;
      table  = Root;
    }

  //  Annoying, but there is no point caching rows of 2 dimensions, so it is handled
  //    as a special case below where a radix 2 FFT of span 1 is performed in place

  if (n == 2)
    { for (j = 0; j < size; j += ostep)
        { Complex_f *s0 = data + j;
          Complex_f *s2 = s0 + hstep;
          Complex_f *se = s2;

#ifdef __SSE3__
          if (step >= 2)         //  Sneaky: if the step size isn't at least 2 then there is
                                 //     only one row to cache and you can't use the SSE3 strategy
            while (s0 != se)
              { __m128 p, q;

                p = LOADP(s0);
                q = LOADP(s2);
                STOREP(s0,ADDP(p,q));
                STOREP(s2,SUBP(p,q));
                s0 += 2;
                s2 += 2;
              }

          else
#endif
            while (s0 != se)
              { float kr = s0->real;
                float ki = s0->imag;
                float jr = s2->real;
                float ji = s2->imag;

                s0->real = kr + jr;
                s0->imag = ki + ji;
                s2->real = kr - jr;
                s2->imag = ki - ji;

                s0 += 1;
                s2 += 1;
              }
        }
      if (invert)
        { Complex_f *s = data;
          Complex_f *e = data + size;

          while (s != e)
#ifdef __SSE3__
            { __m128 q = MULP(Half,LOADP(s));
              STOREP(s,q);
              s += 2;
            }
#else
            { s->real *= .5;
              s->imag *= .5;
              s += 1;
            }
#endif
        }
      return;
    }

  //  NB: n >= 4 if you get here.

#ifdef __SSE3__

  if (step >= 2)     //     Note carefully that if step is 1, then you cannot use the
                     //       SSE3 strategy of solving 2 rows at a time.
    for (j = 0; j < size; j += ostep)
      { jtop = j+step;

        for (i = j; i < jtop; i += 2)     //  Move and auto sort 2 rows starting at elements i and
          {                               //    i+1 into cache
            { Complex_f *s0 = data + i;
              Complex_f *s2 = s0 + hstep;
              Complex_f *j0 = cache;
              Complex_f *j2 = j0 + n;
              Complex_f *je = j2;

              while (j0 != je)
                { __m128 p, q;

                  p = LOADP(s0);
                  s0 += step;
                  q = LOADP(s2);
                  s2 += step;

                  STOREP(j0,ADDP(p,q));
                  j0 += 2;
                  STOREP(j0,SUBP(p,q));
                  j0 += 2;

                  p = LOADP(s0);
                  s0 += step;
                  q = LOADP(s2);
                  s2 += step;

                  STOREP(j2,ADDP(p,q));
                  j2 += 2;
                  STOREP(j2,SUBP(p,q));
                  j2 += 2;
                }
            }

            { int s, h;     // Sort passes for the first half, radix-4 and maybe one radix-2 pass
              int nn = 2*n; //   for the remainder.  Remember s and h, twice there normal values.

              for (s = 4, h = n/2; s < h; s <<= 1, h >>= 1)
                sort2_array(cache,nn,s,h,table);

              for (h = (s<<2); h <= nn; s <<= 2, h <<= 2)
                radix4_array(cache,nn,s,h,table);

              if (s < nn)
                { h = (s << 1);
                  radix2_array(cache,nn,s,h,table);
                }
            }

            { Complex_f *s = data + i;  //  Return the cached rows to the matrix, normalizing
              Complex_f *d = cache;     //    if invert is on
              Complex_f *e = d + 2*n;

              if (invert)
                { float  v  = (float) (1./n);
                  __m128 p = LOADUP(&v);
                  __m128 q;

                  while (d != e)
                    { q = MULP(p,LOADP(d));
                      STOREP(s,q);
                      s += step;
                      d += 2;
                    }
                }
              else
                { while (d != e)
                    { s[0] = *d++;
                      s[1] = *d++;
                      s += step;
                    }
                }
            }
          }
      }

  else

#endif

    for (j = 0; j < size; j += ostep)
      { jtop = j+step;

        for (i = j; i < jtop; i++)  //  Move and auto sort the row starting at element i into cache
          {
            { Complex_f *s0 = data + i;
              Complex_f *s2 = s0 + hstep;
              Complex_f *j0 = cache;
              Complex_f *j2 = j0 + n2;
              Complex_f *je = j2;

              while (j0 != je)
                { float pr, pi;
                  float dr, di;

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
            }

            { int s, h;   //  Sort passes for the first half, radix-4 and maybe one radix-2 pass
                          //    for the remainder
              for (s = 2, h = n4; s < h; s <<= 1, h >>= 1)
                sort2_table(cache,n,s,h,table);

              for (h = (s << 2); h <= n; s <<= 2, h <<= 2)
                radix4_table(cache,n,s,h,table);

              if (s < n)
                { h = (s << 1);
                  radix2_table(cache,n,s,h,table);
                }
            }

            { Complex_f *s = data + i;    //  Return the cached row to the matrix
              Complex_f *d = cache;       //    and normalize if invert is on
              Complex_f *e = d + n;

              if (invert)
                { float v = (float) (1./n);

                  while (d != e)
                    { s->real = d->real * v;
                      s->imag = d->imag * v;
                      s += step;
                      d += 1;
                    }
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
  //    for convenience a pointer to data is returned by FFT_nf.  If invert is non-zero then
  //    the inverse Fourier Transform is performed.

Complex_f *FFT_nf(int ndim, int *dims, Complex_f *data, int invert)
{ int64      size;
  int64      step;
  int64      k;
  int        d;
  int        dmax;
  Complex_f *cache;        //  One or two lines of a matrix are cached here.
  Complex_f  Cache[4096];  //  One if SSE3 off, and two if SSE3 optimizations are on

  size = init_multi_fft(ndim,dims,&dmax);
  if (size <= 0)
    return (NULL);

#ifdef __SSE3__
  if (2*dmax > 4096)
    cache = (Complex_f *) malloc(2*sizeof(Complex_f)*((size_t) dmax));
  else
    cache = Cache;
#else
  if (dmax > 4096)
    cache = (Complex_f *) malloc(sizeof(Complex_f)*((size_t) dmax));
  else
    cache = Cache;
#endif
  if (cache == NULL)
    { if (grab_message()) return (NULL);
      sprintf(FFT_Esource,"FFT_nf");
      sprintf(FFT_Estring,"Out of memory");
      return (NULL);
    }

  step = dims[0];      // Do innermost rows in place with 1-d routines (if the dimesnion is > 1!)
  if (step > 1)
    for (k = 0; k < size; k += step)
      FFT_1f((int) step,data+k,invert);

  for (d = 1; d < ndim; d++)    // Do rows of each higher dimension (if the dimension is > 1!)
    if (dims[d] > 1)
      { fft_a_dim(size,step,dims[d],invert,data,cache);
        step *= dims[d];
      }

#ifdef __SSE3__
  if (2*dmax > 4096)
    free(cache);
#else
  if (dmax > 4096)
    free(cache);
#endif

  return (data);
}


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL REAL FFT & INVERSE: argcheck_nf                                         *
 *                                                                                            *
\*********************************************************************************************/

  //  Common argument checks for real multi-dim. FFT routines

static int argcheck_nf(int ndim, int *dims, char *source)
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
  //    The forward transform, Real_FFT_nf, takes a float array of size s, and *in-place* produces
  //    a Complex_f array c of size s/2 that is the conjugate symmetric FFT of the real data for
  //    the first half of the lowest dimension [0..M-1] where M = dims[0]/2.  In fact, some M-terms
  //    are essential and tucked into 0-terms that are redundant, so the complex array is not
  //    directly interpretable as the lower half of the FFT, but rather is an encoding of all the
  //    terms not inferable by symmetry and hence sufficient for computing correlations and
  //    convolutions.
  //
  //    The pointer C returned by FFT is really the same as rdata, the FFT is performed in-place.
  //    To fully document the encoding in C, if a = (ik,...,i1) then a* = (nk-ik mod nk, ...,
  //    n1-i1 mod n1 ) where k = ndim-1 and nx = dims[x].  Then the values of C encode values
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

Complex_f *Real_FFT_nf(int in_ndim, int *in_dims, float *rdata)
{ Complex_f *rfft = (Complex_f *) rdata;
  int        M, A, *dims;
  int        ndim;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int64      size, nsub, offs;

  if (argcheck_nf(in_ndim,in_dims,"Real_FFT_nf"))
    return (NULL);

  { int i;       //  If lowest dims are 1, then effectively ignore them (by modifying ndim and dims)
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
  if (FFT_nf(ndim,dims,rfft,0) == NULL)
    { if (my_message())
        sprintf(FFT_Esource,"Real_FFT_nf");
      dims[0] = M * 2;
      return (NULL);
    }
  dims[0] = M * 2;


  //  Given H the fft of rfft, the FFT of rdata for (a,x) in [0,nk-1] x ... [0,n1-1] x [0,M] is:
  //
  //      F[a,x] = .5*(H[a,x] + H[a*,x*]*) - .5*i*w^x*(H[a,x] - H[a*,x*]*)
  //
  //  where if a = (ik,...,i1) then a* = (nk-ik mod nk, ..., n1-i1 mod n1 ) where k = ndim-1,
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
                    offs   += base[d];
                  }
                else
                  break;
              }
          b = nsub - (offs + a);        //  Apply the formula above
        }

        { Trig_Block *t  = Root + A;
          Complex_f  *ck, *cj, *ce;

#ifdef __SSE3__
          __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
          __m128 ur  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
#endif

          if (a < b)                  //  Handle [a,0] where a < a*
            { ck = rfft + a;
              cj = rfft + b;

              { float kr = ck->real;    //  Forward xform for 0 row:  k = [a,0]  j = [b,0]
                float ki = ck->imag;    //    Ck = F[a,0] = .5 ((Ck+Cj*) - i(Ck-Cj*))
                float jr = cj->real;    //    Cj = F[a,M] = .5 ((Ck+Cj*) + i(Ck-Cj*))
                float ji = cj->imag;    //  Can do this as F[b,0] = F[a*,0] and F[b,M] = F[a*,M]

                float f0r = (float) (.5 * (kr + jr));
                float f0i = (float) (.5 * (ki - ji));
                float f1r = (float) (.5 * (ji + ki));
                float f1i = (float) (.5 * (jr - kr));

                ck->real =   f0r + f1r;
                ck->imag =   f0i + f1i;
                cj->real =   f0r - f1r;
                cj->imag = -(f0i - f1i);
              }

              ce  = ck + M;            //  Let x run from 1 to M-1 in the loops beow
              cj += (M-1);
              ck += 1;
            }
          else if (a == b)             //  Handle [a,0] where a = a*
            { ck = rfft + a;

              { float kr = ck->real;     //  Ck.real = F_0,a (is real) = Ck.real + Ck.imag
                float ki = ck->imag;     //  Ck.imag = F_M,a (is real) = Ck.real - Ck.imag

                ck->real = kr + ki;
                ck->imag = kr - ki;
              }

              ce  = ck + M/2;          //  a = a* so let x run from 1 to M/2-1 in the loops below
              cj  = ck + (M-1);
              ck += 1;
            }
          else
            continue;

          //  Handle [a,x] where x > 0 and a <= a* in the remainder

#ifdef __SSE3__
          if (ck < ce)     //  Sneaky: Do first rep of SSE3 unoptimized to get to even count
#else
          while (ck < ce)
#endif

            { float kr = ck->real;    //  Forward xform:  k = [a,x]  j = [b,M-x]
              float ki = ck->imag;    //    Ck = .5 ((Ck+Cj*) - i(w^k)(Ck-Cj*))
              float jr = cj->real;    //    Cj = .5 ((Cj+Ck*) + i(w^k)*(Cj-Ck*))
              float ji = cj->imag;

              float f0r = (float) (.5 * (kr + jr));
              float f0i = (float) (.5 * (ki - ji));
              float f1r = (float) (.5 * (ji + ki));
              float f1i = (float) (.5 * (jr - kr));

              float ur = t->pow2.real;
              float ui = t->pow2.imag;

              kr = ur*f1r - ui*f1i;
              ki = ur*f1i + ui*f1r;

              ck->real =   f0r + kr;
              ck->imag =   f0i + ki;
              cj->real =   f0r - kr;
              cj->imag = -(f0i - ki);

              ck += 1;
              cj -= 1;
              t  += A;
            }

#ifdef __SSE3__
          while (ck < ce)
            { __m128 fk;
              __m128 f0, f1;
              __m128 ui;

              fk = LOADP(ck);
              fj = LOADHI(fj,cj);
              fj = LOADLO(fj,cj-1);

              f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
              f1 = MULP(Half,ADDSUB(fj,fk));

              ur = LOADHI(ur,&(t->pow2));
              t += A;
              ur = LOADLO(ur,&(t->pow2));
              t += A;
              ui = _mm_movehdup_ps(ur);
              ur = _mm_moveldup_ps(ur);

              fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
              STOREP(ck,ADDP(f0,fk));

              fj = XORP(Conjugate,SUBP(f0,fk));
              STOREHI(cj,fj);
              STORELO(cj-1,fj);

              ck += 2;
              cj -= 2;
            }
#endif
        }
      }
  }

  return (rfft);
} 

  //    The inverse transform, Real_FFT_Inverse_nf, takes a packed complex half-matrix as produced
  //    by Real_FFT_nf, and produces a real-valued result *in-place*.  That is, the pointer returned
  //    is exactly rfft (coerced to be double *).  Note carefully that dims[0] is the length of the
  //    0th dimension of the result double array and twice that of the 0th dimension of rfft.

float *Real_FFT_Inverse_nf(int in_ndim, int *in_dims, Complex_f *rfft)
{ int   M, A, *dims;
  int   ndim;
  int64 base[MAX_REAL_DIM];
  int   cntr[MAX_REAL_DIM];
  int64 size, nsub, offs;

  if (argcheck_nf(in_ndim,in_dims,"Real_FFT_Inverse_nf"))
    return (NULL);

  { int i;

    for (i = 0; i < in_ndim; i++)
      if (in_dims[i] > 1)
        break;
    if (i == in_ndim)
      return ((float *) rfft);
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
                    offs   += base[d];
                  }
                else
                  break;
              }
          b = nsub - (offs + a);
        }

        { Trig_Block *t = Coot + A;
          Complex_f  *ck, *cj, *ce;

#ifdef __SSE3__
          __m128 fj  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
          __m128 ur  = { 0., 0., 0., 0. };    //  Necessary only to shut up -Wall
#endif

          if (a < b)
            { ck = rfft + a;
              cj = rfft + b;

              { float kr = ck->real;
                float ki = ck->imag;
                float jr = cj->real;
                float ji = cj->imag;

                float f0r = (float) (.5 * (kr + jr));
                float f0i = (float) (.5 * (ki - ji));
                float f1r = (float) (.5 * (ji + ki));
                float f1i = (float) (.5 * (jr - kr));

                ck->real =   f0r - f1r;
                ck->imag =   f0i - f1i;
                cj->real =   f0r + f1r;
                cj->imag = -(f0i + f1i);
              }

              ce  = ck + M;
              cj += (M-1);
              ck += 1;
            }
          else if (a == b)
            { ck = rfft + a;

              { float kr = ck->real;
                float ki = ck->imag;

                ck->real = (float) (.5*(kr + ki));
                ck->imag = (float) (.5*(kr - ki));
              }

              ce  = ck + M/2;
              cj  = ck + (M-1);
              ck += 1;
            }
          else
            continue;

#ifdef __SSE3__
          if (ck < ce)
#else
          while (ck < ce)
#endif
            { float kr = ck->real;
              float ki = ck->imag;
              float jr = cj->real;
              float ji = cj->imag;

              float f0r = (float) (.5 * (kr + jr));
              float f0i = (float) (.5 * (ki - ji));
              float f1r = (float) (.5 * (ji + ki));
              float f1i = (float) (.5 * (jr - kr));

              float ur = t->pow2.real;
              float ui = t->pow2.imag;

              kr = ur*f1r - ui*f1i;
              ki = ur*f1i + ui*f1r;

              ck->real =   f0r - kr;
              ck->imag =   f0i - ki;
              cj->real =   f0r + kr;
              cj->imag = -(f0i + ki);

              ck += 1;
              cj -= 1;
              t  += A;
            }

#ifdef __SSE3__
          while (ck < ce)
            { __m128 fk;
              __m128 f0, f1;
              __m128 ui;

              fk = LOADP(ck);
              fj = LOADHI(fj,cj);
              fj = LOADLO(fj,cj-1);

              f0 = MULP(Half,ADDSUB(fk,XORP(fj,Negate)));
              f1 = MULP(Half,ADDSUB(fj,fk));

              ur = LOADHI(ur,&(t->pow2));
              t += A;
              ur = LOADLO(ur,&(t->pow2));
              t += A;
              ui = _mm_movehdup_ps(ur);
              ur = _mm_moveldup_ps(ur);

              fk = ADDSUB(MULP(ur,SHUFFLE(f1)),MULP(ui,f1));
              STOREP(ck,SUBP(f0,fk));

              fj = XORP(Conjugate,ADDP(f0,fk));
              STOREHI(cj,fj);
              STORELO(cj-1,fj);

              ck += 2;
              cj -= 2;
            }
#endif
        }
      }
  }

  dims[0] = M;
  if (FFT_nf(ndim,dims,rfft,1) == NULL)
    { if (my_message())
        sprintf(FFT_Esource,"Real_FFT_Inverse_nf");
      dims[0] = M * 2;
      return (NULL);
    }
  dims[0] = M * 2;

  return ((float *) rfft);
} 


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL CONVOLUTION & CORRELATION: argcheck_cf                                  *
 *                                                                                            *
\*********************************************************************************************/

  //  Complex_Convolution/Correlation_nf performs a convolution/correlation in the spectral domain
  //    in-place within fft1 and returns a pointer to fft1.  The code works fine if fft1 and
  //    fft2 are the same array.  So a complex convolution of data1 and data2, in-place within
  //    data 1, is accomplised with the code:
  //
  //      FFT_nf(k,dims,Complex_Convolution_nf(k,dims,FFT_nf(n,data1,0),FFT_nf(k,dims,data2,0)),1);

Complex_f *Complex_Convolution_nf(int ndim, int *dims, Complex_f *fft1, Complex_f *fft2)
{ int   d;
  int64 size;

  size = 1;
  for (d = 0; d < ndim; d++)
    size *= dims[d];
  return (complex_convolution_1f(size,fft1,fft2));
}

Complex_f *Complex_Correlation_nf(int ndim, int *dims, Complex_f *fft1, Complex_f *fft2)
{ int   d;
  int64 size;

  size = 1;
  for (d = 0; d < ndim; d++)
    size *= dims[d];
  return (complex_correlation_1f(size,fft1,fft2));
}

  //  Common argument checks for real multi-dim. convolution & correlation routines

static int argcheck_cf(int ndim, int *dims, char *source)
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

  //  Real_Convolution/Correlation_nf performs a convolution/correlation in the spectral domain
  //    of the packed, half-space complex encodings of the fft of real arrays.  Note carefully
  //    that the smallest dimension of rfft1 and rfft2 is dims[0]/2, and the result, produced
  //    in place in rfft1 is returned as a pointer to a floating point array whose smallest
  //    dimension is dims[0].  rfft1 and rfft2 can be the same array.  Only Real_Convolution_nf
  //    is commented, the correlation routine is identical save that one multiplies the conjugate
  //    of the second array's elements.

Complex_f *Real_Convolution_nf(int in_ndim, int *in_dims, Complex_f *rfft1, Complex_f *rfft2)
{ int64      size;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int       *dims;
  int        ndim;
  Complex_f *s, *f, *g, *h;

  if (argcheck_cf(in_ndim,in_dims,"Real_Convolution_nf"))
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

  //  rfft1 and rfft2 are packed encodings as described in the comments for Real_FFT_nf.
  //  To compute the packed encoding of the convolution, it suffices to multiply all terms
  //  pairwise.  For elements C[a,0] for which a = a* we have to be careful as really two
  //  real valued terms are encoded in the one complex element, where as for all other elements
  //  of C a straight complex multiply is in order.  Note that there are 2^k values of a for
  //  which a = a* and they are a in {0,nk/2} x {0,nk-1/2} x ... x {0,n1/2}!  We use a counter
  //  scheme to click through the special values of a.

  { int d;               //  { cntr[d] }_d encodes the current special value that will be next
                         //     i.e., a = PI_d=1^k base[d] * cntr[d]
    size = dims[0]/2;    //        where base[d] = dims[d]/2 * PI_j=1^d-1 dims[j] * dims[0]/2
    for (d = 1; d < ndim; d++)
      { base[d] = (dims[d]/2) * size;
        cntr[d] = 0;
        size   *= dims[d];
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

      { float kr, ki;      //  On a special element (its even BTW), handle it
        float jr, ji;

        s->real *= f->real;
        s->imag *= f->imag;

        s += 1;
        f += 1;

#ifndef __SSE3__           //  Then do term-wise products up to the next one (in g)
        while (f < g)
#else
        if (f < g)         //  Sneaky: Do first rep of SSE3 unoptimized to get to even offset
#endif
          { kr = s->real;
            ki = s->imag;
            jr = f->real;
            ji = f->imag;

            s->real = kr*jr - ki*ji;
            s->imag = ki*jr + kr*ji;

            s += 1;
            f += 1;
          }
      }

#ifdef __SSE3__
      while (f < g)
        { __m128 kr = LOADP(s);
          __m128 j  = LOADP(f);
          __m128 ki, p;

          ki = _mm_movehdup_ps(kr);
          kr = _mm_moveldup_ps(kr);

          CMULTIPLY(p,j,kr,ki);
          STOREP(s,p);

          f += 2;
          s += 2;
        }
#endif
    }
  while (g < h);

  return (rfft1);
}

Complex_f *Real_Correlation_nf(int in_ndim, int *in_dims, Complex_f *rfft1, Complex_f *rfft2)
{ int64      size;
  int64      base[MAX_REAL_DIM];
  int        cntr[MAX_REAL_DIM];
  int       *dims;
  int        ndim;
  Complex_f *s, *f, *g, *h;
#ifdef __SSE3__
  __m128     conj = LOADP(&Conjugate);
#endif

  if (argcheck_cf(in_ndim,in_dims,"Real_Correlation_nf"))
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

      { float kr, ki;
        float jr, ji;

        s->real *= f->real;
        s->imag *= f->imag;

        s += 1;
        f += 1;

#ifndef __SSE3__
        while (f < g)
#else
        if (f < g)
#endif
          { kr = s->real;
            ki = s->imag;
            jr = f->real;
            ji = f->imag;

            s->real = kr*jr + ki*ji;
            s->imag = ki*jr - kr*ji;

            s += 1;
            f += 1;
          }
      }

#ifdef __SSE3__
      while (f < g)
        { __m128 kr = LOADP(s);
          __m128 j  = XORP(LOADP(f),conj);
          __m128 ki, p;

          ki = _mm_movehdup_ps(kr);
          kr = _mm_moveldup_ps(kr);

          CMULTIPLY(p,j,kr,ki);
          STOREP(s,p);

          f += 2;
          s += 2;
        }
#endif
    }
  while (g < h);

  return (rfft1);
}


/*********************************************************************************************\
 *                                                                                            *
 *  MULTI-DIMENSIONAL NORMALIZATION: normalize[_0]                                            *
 *                                                                                            *
\*********************************************************************************************/

  //  Normalize_nf normalizes the scores in the ndim-dimensional spatial correlation
  //    matrix cdata that is presumed to have been produced by applying FFT routines to copies
  //    of idata and tdata that had been 0-padded to be of dimensions cdims.  Note that
  //    this implies that idims[x] <= cdims[x] and tdims[x] <= cdims[x].  The normalization
  //    takes effect directly on cdata and for convenience a pointer to it is returned.
  //
  //  Please read the comments prefacing Normalize_1f to understand what sums
  //    need to be produced and how they are used.  In this scenario R and S are multi-dimensional
  //    and develop progressively in each dimension.  Sums are accumulated at each dimensional
  //    level k from n-1 down to 0 recursively.  At level k, normalize(ms,idx,area) will be
  //    called with every value of idx representing the partial coordinate (i_n-1,...,i_k)
  //    where n = ndims and i_x is in [0,cdims[x]-1] for each x in [k,n-1].
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
    float  *im_ave, *tm_ave;
    float  *im_std, *tm_std;
  } Partial_Sums;

typedef struct
  { int   *idims;
    int   *tdims;
    int   *cdims;
    float *cdata;
  } Partial_Args;

  //  The three routines below each compute progressive overlap sums accross a given dimension
  //    each being tuned to a different scenario.  Normalize_nf starts the recursion
  //    and computes tables of sums and sums of squares.  normalize computes intermediate sums
  //    and recurses downward to the next level.  normalize_0 handles the deepest level of the
  //    recursion where terms of data are normalized.

static void normalize_0(Partial_Args *ag, Partial_Sums *ms, int64 idx, int64 area)
{ float IA, ID;
  float TA, TD;

  float *isums, *idevs;
  float *tsums, *tdevs;
  int64  b, e, c, f;
  int64  w, t;

  int id = ag->idims[0];
  int td = ag->tdims[0];
  int cd = ag->cdims[0];

  float *data = ag->cdata + idx*cd;

  isums = ms->im_ave;
  idevs = ms->im_std;
  tsums = ms->tm_ave;
  tdevs = ms->tm_std;

  { float *ai = isums;
    float *di = idevs;
    float *at = tsums + td;
    float *dt = tdevs + td;

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
          data[b] = (float) ((data[b] - ia*TA) / ((w-1)*dn));
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

static void normalize(Partial_Args *ag, Partial_Sums *ms, int64 idx, int64 area)
{ Partial_Sums ps;

  int64  in,  tn;
  float *IA, *ID;
  float *TA, *TD;

  float *isums, *idevs;
  float *tsums, *tdevs;
  int64  b, e, c, f;
  int64  w, i, t;

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

  { float *ai = isums;
    float *di = idevs;
    float *at = tsums + td*tn;
    float *dt = tdevs + td*tn;

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
        normalize( ag, &ps, idx*cd + b, w );
      else
        normalize_0( ag, &ps, idx*cd + b, w );

      if (b < id)
        { float *ai = isums + b*in;
          float *ad = idevs + b*in;
          for (i = 0; i < in; i++)
            { IA[i] -= *ai++;
              ID[i] -= *ad++;
            }
          w -= area;
        }
      if (++b >= cd)
        b = 0;

      if (e < id)
        { float *ai = isums + e*in;
          float *ad = idevs + e*in;
          for (i = 0; i < in; i++)
            { IA[i] += *ai++;
              ID[i] += *ad++;
            }
          w += area;
        }
      if (++e >= cd)
        e = 0;

      if (c < td)
        { float *ai = tsums + c*tn;
          float *ad = tdevs + c*tn;
          for (i = 0; i < tn; i++)
            { TA[i] -= *ai++;
              TD[i] -= *ad++;
            }
        }
      if (c-- <= 0)
        c = cd-1;

      if (f < td)
        { float *ai = tsums + f*tn;
          float *ad = tdevs + f*tn;
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

float *Normalize_nf(int ndim, int *idims, float *idata,
                              int *tdims, float *tdata,
                              int *cdims, float *cdata)
{ float *cumulative;
  int64  isize, tsize, space;

  if (ndim == 1)
    return (Normalize_1f(idims[0],idata,tdims[0],tdata,cdims[0],cdata));

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
            sprintf(FFT_Esource,"Normalize_nf");
            sprintf(FFT_Estring,"Operand %d-dimension > correlation dimension",d);
            return (NULL);
          }
        space += isize + tsize;
        isize *= idims[d];
        tsize *= tdims[d];
      }

    cumulative = (float *) malloc(sizeof(float)*2*((size_t) space));
    if (cumulative == NULL)
      { if (grab_message()) return (NULL);
        sprintf(FFT_Esource,"Normalize_nf");
        sprintf(FFT_Estring,"Out of memory");
        return (NULL);
      }
  }

  ndim -= 1;

  { Partial_Sums ps;
    Partial_Args ag;

    int64  in,  tn;
    float *IA, *ID;
    float *TA, *TD;

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

    { float *ai = idata;
      float *at = tdata + tsize;

      w = (id + td) - cd;
      if (w < 0)
        w = 0;

      for (i = 0; i < in; i++)
        IA[i] = ID[i] = 0.;
      for (i = 0; i < tn; i++)
        TA[i] = TD[i] = 0.;

      for (b = 0; b < w; b++)
        { for (i = 0; i < in; i++)
            { float s = *ai++;
              IA[i] += s;
              ID[i] += s*s;
            }
          for (i = tn-1; i >= 0; i--)
            { float s = *--at;
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
          normalize( &ag, &ps, b, w );
        else
          normalize_0( &ag, &ps, b, w );

        if (b < id)
          { float *ai = idata + b*in;
            for (i = 0; i < in; i++)
              { float s = *ai++;
                IA[i] -= s;
                ID[i] -= s*s;
              }
            w -= 1;
          }
        if (++b >= cd)
          b = 0;
  
        if (e < id)
          { float *ai = idata + e*in;
            for (i = 0; i < in; i++)
              { float s = *ai++;
                IA[i] += s;
                ID[i] += s*s;
              }
            w += 1;
          }
        if (++e >= cd)
          e = 0;
  
        if (c < td)
          { float *ai = tdata + c*tn;
            for (i = 0; i < tn; i++)
              { float s = *ai++;
                TA[i] -= s;
                TD[i] -= s*s;
              }
          }
        if (c-- <= 0)
          c = cd-1;
  
        if (f < td)
          { float *ai = tdata + f*tn;
            for (i = 0; i < tn; i++)
              { float s = *ai++;
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
#undef LOADHI
#undef LOADLO

#undef LOADUP
#undef SETP

#undef STOREP
#undef STOREHI
#undef STORELO

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
