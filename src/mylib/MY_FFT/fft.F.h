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

#ifndef _FLOAT_FFTS

#define _FLOAT_FFTS

#ifdef __cplusplus
extern "C" {
#endif

typedef struct      //  Complex number data type
  { float  real;
    float  imag;
  } Complex_f;

//  ALL POINTERS TO COMPLEX AND FLOATING DATA MUST BE 16-BYTE ALIGNED IF COMPILED
//     WITH THE SSE3 OPTIMIZATIONS

#ifndef _DOUBLE_FFTS

//  All the FFT routines assume n is a power of 2.  As a convenience to help you pad
//    your vectors to this size, Next_Power_Of_2 returns the smallest power of 2 greater
//    than or equal to m.

int Next_Power_Of_2(int m);

//  If an FFT routine returns an error value, then a text description of the string can be
//    accessed by calling FFT_Error_String() and the name of the routine that threw
//    the error can be accessed by calling FFT_Error_Source().  If there has been no
//    error or the error does not belong to the current thread, then the routines return NULL.
//  NB: To make this mechanism thread-safe, only one error source and message is available
//    at any time.  That is, the first thread to report an error grabs the "error resource"
//    and doesn't let it go for another thread to use until FFT_Error_Release is called.  While
//    the thread has the error resource it can print the error source and message at its leisure.
//    While quite limiting, I figure this is enough for debug purposes which is exactly the scenario
//    involved as all error messages reported by this module are fatal and indicate that the
//    caller has committed a programming error (unlike, say, not being able to open a file).

char *FFT_Error_Source();
char *FFT_Error_String();
void  FFT_Error_Release();

#endif

/*****************************************************************************************\
*                                                                                         *
*  1-dimensional FFT algorithms                                                           *
*                                                                                         *
\*****************************************************************************************/

//  Basic 1-dimenstional FFT-algorithm.  The FFT is performed in-place within 'data' and
//    for convenience a pointer to data is returned by FFT_1f.  If invert is non-zero then
//    the inverse Fourier Transform is performed.  NULL is returned if an error was detected
//    (see list at the end of this header).

Complex_f *FFT_1f(int n, Complex_f *data, int invert);

//  FFT-algorithms optimized for the case of a real-valued time-domain data.
//
//    The forward transform, Real_FFT_1f, takes a float array of length n, and *in-place* produces
//    a Complex_f array c of length n/2 that is the first half of the conjugate symmetric FFT
//    of the real data with the exception that F_(n/2) (which is real) is tucked into c[0].imag
//    (which works as F_0 is also real).  Again, the pointer returned by FFT_1f is equal to
//    rdata, the FFT is performed in-place.
//
//    The inverse transform, Real_FFT_Inverse_1f, takes a complex half-matrix as produced by
//    Real_FFT_1f, and produces a real-valued result *in-place*.  That is, the pointer returned is
//    exactly rfft (coerced to be float *), the real and imaginary parts of rfft *must be
//    contiguous in memory.  Note carefully that n is the length of the resulting real array
//    and is twice the length of rfft.
//
//    Both routines return NULL if an error was detected (see list at end of this header).

Complex_f *Real_FFT_1f(int n, float *rdata);

float     *Real_FFT_Inverse_1f(int n, Complex_f *rfft);

//  FFT-algorithms for convolution and correlation in the *spectral* domain.  Convolution of two
//    ffts in the spectral domain is simply a matter of taking the term-wise product of the elements
//    of the two ffts.  Likewise correlation requires taking the term-wise product of the elements
//    of the first fft and the conjugate of the elements of the second fft.  So a convolution
//    or correlation in the *time* domain is effected by the template:
//
//        FFT_Inverse( Convolution/Correlation( FFT(data1), FFT(data2)) )
//
//  Complex_Convolution/Correlation_1f performs the spectral domain convolution of fft1 and fft2
//    in-place within fft1 and returns a pointer to fft1.  fft1 and fft2 can be the same array.
//  Real_Convolution/Correlation_1f performs the same operation but on the *encoding* of two packed
//    half-size complex vectors rfft1 and rfft2 of length ** n/2 ** as produced by Real_FFT_1f.
//    These routines return NULL if an error was detected (see list at the end of this header).
//
//  As an example, performing the convolution of two real vectors d1 and d2 of length n,
//    is accomplished with the code:
//
//        Real_FFT_Inverse_1f(n,Real_Convolution_1f(n,Real_FFT_1f(n,d1,0),Real_FFT_1f(n,d2,0)),1);
//
//  The convolution is in d1 upon completion, and d2 contains the Real FFT of its original value.

Complex_f *Complex_Convolution_1f(int n, Complex_f *fft1, Complex_f *fft2);
Complex_f *Complex_Correlation_1f(int n, Complex_f *fft1, Complex_f *fft2);

Complex_f *Real_Convolution_1f(int n, Complex_f *rfft1, Complex_f *rfft2);
Complex_f *Real_Correlation_1f(int n, Complex_f *rfft1, Complex_f *rfft2);

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
//    are progressive intervals of f and g of the same length.
//  Normalize_1f normalizes the scores in the spatial correlation vector cdata
//    that is presumed to have been produced by applying FFT routines to copies of idata
//    and tdata that had been 0-padded to be of length cdim.  Note that idim+tdim <= cdim or
//    the vectors were insufficiently padded to give every overlap correlation.  The normalization
//    takes effect directly on cdata and for convenience a pointer to it is returned.
//  NULL is returned if an error was detected (see list at the end of this header).

float *Normalize_1f(int idim, float *idata,
                    int tdim, float *tdata,
                    int cdim, float *cdata);

//  Utility routines to print complex and real 1-dimensional arrays with a title

void Print_Real_1f(int n, float *data, char *title);
void Print_Complex_1f(int n, Complex_f *data, char *title);


/*****************************************************************************************\
*                                                                                         *
*  Multi-dimensional FFT algorithms                                                       *
*                                                                                         *
\*****************************************************************************************/

//  Basic n-dimensional FFT-algorithm.  The FFT is performed in-place within 'data' and
//    for convenience a pointer to data is returned by FFT_nf.  If invert is non-zero then
//    the inverse Fourier Transform is performed.  NULL is returned if an error was detected
//    (see list at the end of this header).

Complex_f *FFT_nf(int ndim, int *dims, Complex_f *data, int invert);

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
//    The pointer C returned by Real_FFT_nf is really the same as rdata, the FFT is performed
//    in-place. To fully document the encoding in C, if a = (i_K,...,i_1) then
//    a* = (N_K-i_K mod N_K, ..., N_1-i_1 mod N_1) where K = ndims-1 and N_x = dims[x].  Then
//    the values of C encode value of the real-valued FFT F as follows:
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
//
//    The inverse transform, Real_FFT_Inverse_nf, takes a complex half-matrix as produced by
//    Real_FFT_nf, and produces a real-valued result *in-place*.  That is, the pointer returned is
//    exactly rfft (coerced to be float *).  Note carefully that dims[0] is the length of the
//    0'th dimension of the result float array and twice that of the 0th dimension of rfft.
//
//    Both routines return NULL if an error was detected (see list at end of this header).

Complex_f *Real_FFT_nf(int ndim, int *dims, float *rdata);

float     *Real_FFT_Inverse_nf(int ndim, int *dims, Complex_f *rfft);

//  Complex_Convolution/Correlation_nf performs the spectral domain convolution of fft1 and fft2
//    in-place within fft1 and returns a pointer to fft1.  fft1 and fft2 can be the same array.
//  Real_Convolution/Correlation_nf performs the same operation but on the *encoding* of two packed
//    half-size complex vectors rfft1 and rfft2 whose first dimension is ** dims[0]/2 ** as
//    produced by Real_FFT_nf.  These routines return NULL if an error was detected (see list
//    at the end of this header).
//
//  As an example, performing the convolution of two real arrays d1 and d2 of dimension k,
//    is accomplished with the code:
//
//        Real_FFT_Inverse_nf(k,dim,Real_Convolution_nf(k,dim,Real_FFT_nf(k,dim,d1,0),
//                                                            Real_FFT_nf(k,dim,d2,0)),1);
//
//  The convolution is in d1 upon completion, and d2 contains the Real FFT of its original value.

Complex_f *Complex_Convolution_nf(int ndim, int *dims, Complex_f *fft1, Complex_f *fft2);
Complex_f *Complex_Correlation_nf(int ndim, int *dims, Complex_f *fft1, Complex_f *fft2);

Complex_f *Real_Convolution_nf(int ndim, int *dims, Complex_f *rfft1, Complex_f *rfft2);
Complex_f *Real_Correlation_nf(int ndim, int *dims, Complex_f *rfft1, Complex_f *rfft2);

//  Normalize_nf normalizes the scores in the ndim-dimensional spatial correlation
//    matrix cdata that is presumed to have been produced by applying FFT routines to copies
//    of idata and tdata that had been 0-padded to be of dimensions cdims.  Note that
//    idims[x]+tdims[x] <= cdims[x] or the arrays were insufficiently padded to give every
//    overlap correlation.  The normalization takes effect directly on cdata and for
//    convenience a pointer to it is returned.
//  This routine uses O(I/in + T/tn) auxilliary space where I and T are the sizes of idata
//    and tdata and in and tn are the lengths of the ndim-1^st dimensions thereof.  This
//    space is allocated and freed with each call in order to be re-entrant.
//  NULL is returned if an error was detected (see list at end of this header).

float *Normalize_nf(int ndim, int *idims, float *idata,
                              int *tdims, float *tdata,
                              int *cdims, float *cdata);

//  Utility routines to print complex and real multi-dimensional arrays with a title

void Print_Real_nf(int ndim, int *dims, float *data, char *title);
void Print_Complex_nf(int ndim, int *dims, Complex_f *data, char *title);

#ifdef __cplusplus
}
#endif

#endif

/******** LIST OF ERRORS & COMPILED CONSTANTS *****************************************

ERROR MESSAGES:  Any comments are in [] brackets on a seperate line.  For some limit
  errors changing one or more defined constants can remove the error, in which case
  a comment detailing the parameter(s) is given.

        FFT_1f
        Real_FFT_1f
        Real_FFT_Inverse_1f
            n = ? is not a power of 2
            n is larger than hardcoded maximum of ?M
                [maximum is a 8*NINETY^2 (8M) save twice that (16M) for a real FFT]

        Real_Convolution_1f
        Real_Correlation_1f
            n = ? must be even

        Normalize_1f
            Sum of operand dimensions > correlation dimension

        FFT_nf
            dims[?] = ? is not a power of 2
            dims[?] is larger than hardcorded maximum of ?K
                 [maximum is 4*NINETY (4K)]
            Out of memory
                 [needs working storage <= max(dims[i])*8 bytes,
                  uses 32Kb of stack for each thread if multi-threaded]

        Real_FFT_nf
        Real_FFT_Inverse_nf
            dims[?] = ? is not a power of 2
            dims[?] is larger than hardcoded maximum of ?K
                 [maximum is 4*NINETY (4K) save twice that (8K) for dims[0])]
            Out of memory
                 [needs working storage <= max(dims[i])*8 bytes,
                  uses 32Kb of stack for each thread if multi-threaded]
            Array cannot have more than MAX_REAL_DIM(128) dimensions
                 [maximum is MAX_REAL_DIM (128)]

        Real_Convolution_nf
        Real_Correlation_nf
            dims[0] = ? must be even
            Array cannot have more than MAX_REAL_DIM(128) dimensions
                 [maximum is MAX_REAL_DIM (128)]

        Normalize_nf
            Operand %d-dimension > correlation dimension"
            Out of memory
                 [needs working storage = O(size(I)/dims_I[ndim-1] + size(T)/dims_T[ndim-1])]

COMPILED CONSTANTS: you can change these if you wish.  Indeed L2_CACHE and NINETY should be
  tuned for a particular L1 and L2 cache size.  The values below are set for a 2.33 GHz Intel
  Core 2 Duo.  I suggest when you compile for a new machine that you try multiples or fractions of
  2 of the set values until those giving the best performance are found.  It is unlikely you will
  ever need to change MAX_REAL_DIM.

       L2_CACHE     262144   //  Problems bigger than this get cross cut
       NINETY         1024   //  Ninety degree ticks (problems > 4*NINETY are not tabled)
       MAX_REAL_DIM    128   //  Maximum dimensionality of a real multi-dim. fft

*/
