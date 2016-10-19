#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fft.F.h"
#include "fft.D.h"

#define QSIZE 16
#define GSIZE 8388608

#define HSIZE  4
#define DSIZE  4096

int main(int argc, char *argv[])
{ void *block;

  (void) argc;
  (void) argv;

  if (GSIZE > DSIZE*DSIZE)
    block = (void *) malloc(sizeof(double)*GSIZE*4);
  else
    block = (void *) malloc(sizeof(double)*DSIZE*DSIZE*4);

  { Complex_f *data = (Complex_f *) block;
    int        i, k;

    printf("\n1-DIMENSIONAL FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { data[i].real = (float) i;
            data[i].imag = (float) (QSIZE-i);
          }

        Print_Complex_1f(QSIZE,data,"In");

        FFT_1f(QSIZE,data,0);

        Print_Complex_1f(QSIZE,data,"FFT");

        FFT_1f(QSIZE,data,1);

        Print_Complex_1f(QSIZE,data,"FFT_1(FFT)");
      }
  }

  { Complex_f *data = (Complex_f *) block;
    int        i, k;

    printf("\n1-DIMENSIONAL FFT TEST:\n\n");

    for (k = 1; k <= GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          { data[i].real = (float) ((i % 256) + 1);
            data[i].imag = (float) (((k-i) % 256) + 1);
          }

        if (FFT_1f(k,data,0) == NULL) goto error;
        if (FFT_1f(k,data,1) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { dev = fabs(data[i].real - ((i%256)+1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= ((i%256)+1.);
              if (dev > mxdev) mxdev = dev;
              avg += dev;
              dev = fabs(data[i].imag - (((k-i) % 256) + 1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= ((k-i) % 256) + 1.;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/(2.*k));
        }
      }
  }

  { float     *data = (float *) block;
    int        i, k;
    Complex_f *fft;

    printf("\n1-DIMENSIONAL *REAL* FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          data[i] = (float) i;

        Print_Real_1f(QSIZE,data,"Input");

        fft = Real_FFT_1f(QSIZE,data);

        Print_Complex_1f(QSIZE/2,fft,"Real_FFT");

        Real_FFT_Inverse_1f(QSIZE,fft);

        Print_Real_1f(QSIZE,data,"F_1(F(data))");
      }
  }

  { float     *data = (float *) block;
    int        i, k;
    Complex_f *fft;

    printf("\n1-DIMENSIONAL *REAL* FFT TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = (float) ((i%256)+1);

        if ((fft = Real_FFT_1f(k,data)) == NULL) goto error;
        if (Real_FFT_Inverse_1f(k,fft) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { dev = fabs(data[i] - (i%256+1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= (i%256+1);
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { float *rata = (float *) block;
    float *fata = rata + QSIZE;
    int          i, k;

    printf("\n1-DIMENSIONAL REAL CONVOLUTION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { rata[i] = (float) (i%2+1);
            fata[i] = (float) (i%2+1);
          }

        Print_Real_1f(QSIZE,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_1f(QSIZE,Real_Convolution_1f(QSIZE,Real_FFT_1f(QSIZE,fata),
                                                            Real_FFT_1f(QSIZE,rata)));

        Print_Real_1f(QSIZE,fata,"Convolution");
      }
  }

  { float *rata = (float *) block;
    float *fata = rata + QSIZE;
    int          i, k;

    printf("\n1-DIMENSIONAL REAL CORRELATION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { rata[i] = (float) (i%2+1);
            fata[i] = (float) (i%2+1);
          }

        Print_Real_1f(QSIZE,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_1f(QSIZE,Real_Correlation_1f(QSIZE,Real_FFT_1f(QSIZE,fata),
                                                            Real_FFT_1f(QSIZE,rata)));

        Print_Real_1f(QSIZE,fata,"Correlation");
      }
  }

  { float *data = (float *) block;
    float *rata = data + 2*GSIZE;
    int    i, k;

    printf("\n1-DIMENSIONAL REAL CONVOLUTION TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = rata[i] = (float) (i%2+1);

        if (Real_FFT_1f(k,data) == NULL) goto error;
        if (Real_FFT_1f(k,rata) == NULL) goto error;
        if (Real_Convolution_1f(k,(Complex_f *) data,(Complex_f *) rata) == NULL) goto error;
        if (Real_FFT_Inverse_1f(k,(Complex_f *) data) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { if (k == 1)
                x = 1;
              else
                x = 2*k + ((i+1)%2)*(k/2);

              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { float *data = (float *) block;
    float *rata = data + 2*GSIZE;
    int    i, k;

    printf("\n1-DIMENSIONAL REAL CORRELATION TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = rata[i] = (float) (i%2+1);

        if (Real_FFT_1f(k,data) == NULL) goto error;
        if (Real_FFT_1f(k,rata) == NULL) goto error;
        if (Real_Correlation_1f(k,(Complex_f *) data,(Complex_f *) rata) == NULL) goto error;
        if (Real_FFT_Inverse_1f(k,(Complex_f *) data) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { if (k == 1)
                x = 1;
              else
                x = 2*k + ((i+1)%2)*(k/2);

              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { Complex_f   *data = (Complex_f *) block;
    int          dims[2];
    int          i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)

      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              data[d].real = (float) d;
              data[d].imag = (float) (base1-d);
            }

        Print_Complex_nf(2,dims,data,"Input");

        FFT_nf(2,dims,data,0);

        Print_Complex_nf(2,dims,data,"FFT");

        FFT_nf(2,dims,data,1);

        Print_Complex_nf(2,dims,data,"FFT_1(FFT)");
      }
  }

  { Complex_f   *data = (Complex_f *) block;
    int          dims[2];
    int          h, k, i, j, d;

    printf("\nMULTI-DIMENSIONAL FFT TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= DSIZE; k *= 2)
      { dims[0] = k;
        dims[1] = h;

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*k+j;
              data[d].real = (float) ((d % 256) + 1);
              data[d].imag = (float) (abs((k-d) % 256) + 1);
            }

        if (FFT_nf(2,dims,data,0) == NULL) goto error;
        if (FFT_nf(2,dims,data,1) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < dims[1]; i++)
            for (j = 0; j < dims[0]; j++)
              { d = i*k+j;
                x = (d % 256) + 1;
                dev = fabs(data[d].real - x);
                if (dev > mxdiff) mxdiff = dev;
                dev /= x;
                if (dev > mxdev) mxdev = dev;
                avg += dev;
                x = abs((k-d) % 256) + 1;
                dev = fabs(data[d].imag - x);
                if (dev > mxdiff) mxdiff = dev;
                dev /= x;
                if (dev > mxdev) mxdev = dev;
                avg += dev;
              }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { float        *data = (float *) block;
    int           dims[2];
    Complex_f    *fft;
    int           k, i;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL *REAL* FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < dims[0]*dims[1]; i++)
          data[i] = (float) i;

        Print_Real_nf(2,dims,data,"Input");

        fft = Real_FFT_nf(2,dims,data);

        dims[0] /= 2;
        Print_Complex_nf(2,dims,fft,"Real_FFT");
        dims[0] *= 2;

        Real_FFT_Inverse_nf(2,dims,fft);

        Print_Real_nf(2,dims,data,"F_1(F(data))");
      }
  }

  { float        *data = (float *) block;
    int           dims[2];
    Complex_f    *fft;
    int           h, k, i;

    printf("\nMULTI-DIMENSIONAL *REAL* FFT TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= 2*DSIZE; k *= 2)
      { dims[0] = k;
        dims[1] = h;

        for (i = 0; i < k*h; i++)
          data[i] = (float) ((i % 256) + 1);

        if ((fft = Real_FFT_nf(2,dims,data)) == NULL) goto error;
        if (Real_FFT_Inverse_nf(2,dims,fft) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k*h; i++)
            { x = (i % 256) + 1;
              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { float        *rata = (float *) block;
    float        *fata = rata + HSIZE*HSIZE;
    int           dims[2];
    int           i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL REAL CONVOLUTION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              rata[d] = (float) (d % 2 + 1);
              fata[d] = (float) (d % 2 + 1);
            }

        Print_Real_nf(2,dims,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_nf(2,dims,
          Real_Convolution_nf(2,dims,Real_FFT_nf(2,dims,fata),Real_FFT_nf(2,dims,rata)));

        Print_Real_nf(2,dims,fata,"Convolution");
      }
  }

  { float        *rata = (float *) block;
    float        *fata = rata + HSIZE*HSIZE;
    int           dims[2];
    int           i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL REAL CORRELATION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              rata[d] = (float) (d % 2 + 1);
              fata[d] = (float) (d % 2 + 1);
            }

        Print_Real_nf(2,dims,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_nf(2,dims,
          Real_Correlation_nf(2,dims,Real_FFT_nf(2,dims,fata),Real_FFT_nf(2,dims,rata)));

        Print_Real_nf(2,dims,fata,"Correlation");
      }
  }

  { float        *rata = (float *) block;
    float        *fata = rata + 2*DSIZE*DSIZE;
    int           dims[2];
    int           h, k, i;

    printf("\nMULTI-DIMENSIONAL REAL CONVOLUTION TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= 2*DSIZE; k *= 2)
      { int size;

        dims[0] = k;
        dims[1] = h;
        size    = k*h;

        for (i = 0; i < size; i++)
          { rata[i] = (float) ((i % 2) + 1);
            fata[i] = (float) ((i % 2) + 1);
          }

        if (Real_FFT_nf(2,dims,fata) == NULL) goto error;
        if (Real_FFT_nf(2,dims,rata) == NULL) goto error;
        if (Real_Convolution_nf(2,dims,(Complex_f *) fata,(Complex_f *) rata) == NULL) goto error;
        if (Real_FFT_Inverse_nf(2,dims,(Complex_f *) fata) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < size; i++)
            { if (size == 1)
                x = 1;
              else
                { x = 4 + (i+1)%2;
                  x = x * k * h / 2;
                }
              dev = fabs(fata[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { float       *rata = (float *) block;
    float       *fata = rata + 2*DSIZE*DSIZE;
    int          dims[2];
    int          h, k, i;

    printf("\nMULTI-DIMENSIONAL REAL CORRELATION TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= 2*DSIZE; k *= 2)
      { int size;

        dims[0] = k;
        dims[1] = h;
        size    = k*h;

        for (i = 0; i < size; i++)
          { rata[i] = (float) ((i % 2) + 1);
            fata[i] = (float) ((i % 2) + 1);
          }

        if (Real_FFT_nf(2,dims,fata) == NULL) goto error;
        if (Real_FFT_nf(2,dims,rata) == NULL) goto error;
        if (Real_Correlation_nf(2,dims,(Complex_f *) fata,(Complex_f *) rata) == NULL) goto error;
        if (Real_FFT_Inverse_nf(2,dims,(Complex_f *) fata) == NULL) goto error;

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < size; i++)
            { if (size == 1)
                x = 1;
              else
                { x = 4 + (i+1)%2;
                  x = x * k * h / 2;
                }
              dev = fabs(fata[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { Complex_d *data = (Complex_d *) block;
    int        i, k;

    printf("\n1-DIMENSIONAL FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { data[i].real = i;
            data[i].imag = QSIZE-i;
          }

        Print_Complex_1d(QSIZE,data,"In");

        FFT_1d(QSIZE,data,0);

        Print_Complex_1d(QSIZE,data,"FFT");

        FFT_1d(QSIZE,data,1);

        Print_Complex_1d(QSIZE,data,"FFT_1(FFT)");
      }
  }

  { Complex_d *data = (Complex_d *) block;
    int        i, k;

    printf("\n1-DIMENSIONAL FFT TEST:\n\n");

    for (k = 1; k <= GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          { data[i].real = (i % 256) + 1;
            data[i].imag = ((k-i) % 256) + 1;
          }

        FFT_1d(k,data,0);
        FFT_1d(k,data,1);

        { double avg, dev, mxdev, mxdiff;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { dev = fabs(data[i].real - ((i%256)+1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= ((i%256)+1.);
              if (dev > mxdev) mxdev = dev;
              avg += dev;
              dev = fabs(data[i].imag - (((k-i) % 256) + 1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= ((k-i) % 256) + 1.;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/(2.*k));
        }
      }
  }

  { double    *data = (double *) block;
    int        i, k;
    Complex_d *fft;

    printf("\n1-DIMENSIONAL *REAL* FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          data[i] = i;

        Print_Real_1d(QSIZE,data,"Input");

        fft = Real_FFT_1d(QSIZE,data);

        Print_Complex_1d(QSIZE/2,fft,"Real_FFT");

        Real_FFT_Inverse_1d(QSIZE,fft);

        Print_Real_1d(QSIZE,data,"F_1(F(data))");
      }
  }

  { double    *data = (double *) block;
    int        i, k;
    Complex_d *fft;

    printf("\n1-DIMENSIONAL *REAL* FFT TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = (i%256)+1;

        fft = Real_FFT_1d(k,data);
        Real_FFT_Inverse_1d(k,fft);

        { double avg, dev, mxdev, mxdiff;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { dev = fabs(data[i] - (i%256+1));
              if (dev > mxdiff) mxdiff = dev;
              dev /= (i%256+1);
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { double *rata = (double *) block;
    double *fata = rata + QSIZE;
    int     i, k;

    printf("\n1-DIMENSIONAL REAL CONVOLUTION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { rata[i] = i%2+1;
            fata[i] = i%2+1;
          }

        Print_Real_1d(QSIZE,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_1d(QSIZE,
            Real_Convolution_1d(QSIZE,Real_FFT_1d(QSIZE,fata),Real_FFT_1d(QSIZE,rata)));

        Print_Real_1d(QSIZE,fata,"Convolution");
      }
  }

  { double *rata = (double *) block;
    double *fata = rata + QSIZE;
    int     i, k;

    printf("\n1-DIMENSIONAL REAL CORRELATION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < QSIZE; i++)
          { rata[i] = i%2+1;
            fata[i] = i%2+1;
          }

        Print_Real_1d(QSIZE,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_1d(QSIZE,
            Real_Correlation_1d(QSIZE,Real_FFT_1d(QSIZE,fata),Real_FFT_1d(QSIZE,rata)));

        Print_Real_1d(QSIZE,fata,"Correlation");
      }
  }

  { double *data = (double *) block;
    double *rata = data + 2*GSIZE;
    int     i, k;

    printf("\n1-DIMENSIONAL REAL CONVOLUTION TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = rata[i] = i%2+1;

        Real_FFT_Inverse_1d(k,Real_Convolution_1d(k,Real_FFT_1d(k,data),Real_FFT_1d(k,rata)));

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { if (k == 1)
                x = 1;
              else
                x = 2*k + ((i+1)%2)*(k/2);

              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { double *data = (double *) block;
    double *rata = data + 2*GSIZE;
    int     i, k;

    printf("\n1-DIMENSIONAL REAL CORRELATION TEST:\n\n");

    for (k = 1; k <= 2*GSIZE; k *= 2)
      { for (i = 0; i < k; i++)
          data[i] = rata[i] = i%2+1;

        Real_FFT_Inverse_1d(k,Real_Correlation_1d(k,Real_FFT_1d(k,data),Real_FFT_1d(k,rata)));

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k; i++)
            { if (k == 1)
                x = 1;
              else
                x = 2*k + ((i+1)%2)*(k/2);

              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %8d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                  k,mxdiff,mxdev,avg/k);
        }
      }
  }

  { Complex_d    *data = (Complex_d *) block;
    int           dims[2];
    int           i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)

      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              data[d].real = d;
              data[d].imag = base1-d;
            }

        Print_Complex_nd(2,dims,data,"Input");

        FFT_nd(2,dims,data,0);

        Print_Complex_nd(2,dims,data,"FFT");

        FFT_nd(2,dims,data,1);

        Print_Complex_nd(2,dims,data,"FFT_1(FFT)");
      }
  }

  { Complex_d    *data = (Complex_d *) block;
    int           dims[2];
    int           h, k, i, j, d;

    printf("\nMULTI-DIMENSIONAL FFT TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= DSIZE; k *= 2)
      { dims[0] = k;
        dims[1] = h;

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*k+j;
              data[d].real = (d % 256) + 1;
              data[d].imag = abs((k-d) % 256) + 1;
            }

        FFT_nd(2,dims,data,0);
        FFT_nd(2,dims,data,1);

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < dims[1]; i++)
            for (j = 0; j < dims[0]; j++)
              { d = i*k+j;
                x = (d % 256) + 1;
                dev = fabs(data[d].real - x);
                if (dev > mxdiff) mxdiff = dev;
                dev /= x;
                if (dev > mxdev) mxdev = dev;
                avg += dev;
                x = abs((k-d) % 256) + 1;
                dev = fabs(data[d].imag - x);
                if (dev > mxdiff) mxdiff = dev;
                dev /= x;
                if (dev > mxdev) mxdev = dev;
                avg += dev;
              }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { double       *data = (double *) block;
    int           dims[2];
    Complex_d    *fft;
    int           k, i;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL *REAL* FFT EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { for (i = 0; i < dims[0]*dims[1]; i++)
          data[i] = i;

        Print_Real_nd(2,dims,data,"Input");

        fft = Real_FFT_nd(2,dims,data);

        dims[0] /= 2;
        Print_Complex_nd(2,dims,fft,"Real_FFT");
        dims[0] *= 2;

        Real_FFT_Inverse_nd(2,dims,fft);

        Print_Real_nd(2,dims,data,"F_1(F(data))");
      }
  }

  { double       *data = (double *) block;
    int           dims[2];
    Complex_d    *fft;
    int           h, k, i;

    printf("\nMULTI-DIMENSIONAL *REAL* FFT TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= DSIZE; k *= 2)
      { dims[0] = k;
        dims[1] = h;

        for (i = 0; i < k*h; i++)
          data[i] = (i % 256) + 1;

        fft = Real_FFT_nd(2,dims,data);
        Real_FFT_Inverse_nd(2,dims,fft);

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < k*h; i++)
            { x = (i % 256) + 1;
              dev = fabs(data[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { double       *rata = (double *) block;
    double       *fata = rata + HSIZE*HSIZE;
    int           dims[2];
    int           i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL REAL CONVOLUTION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              rata[d] = d % 2 + 1;
              fata[d] = d % 2 + 1;
            }

        Print_Real_nd(2,dims,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_nd(2,dims,
            Real_Convolution_nd(2,dims,Real_FFT_nd(2,dims,fata),Real_FFT_nd(2,dims,rata)));

        Print_Real_nd(2,dims,fata,"Convolution");
      }
  }

  { double       *rata = (double *) block;
    double       *fata = rata + HSIZE*HSIZE;
    int           dims[2];
    int           i, j, k, d;

    dims[0] = HSIZE;
    dims[1] = HSIZE;

    printf("\nMULTI-DIMENSIONAL REAL CORRELATION EXAMPLE:\n");

    for (k = 0; k < 1; k++)
      { int base1 = dims[0];

        for (i = 0; i < dims[1]; i++)
          for (j = 0; j < dims[0]; j++)
            { d = i*base1+j;
              rata[d] = d % 2 + 1;
              fata[d] = d % 2 + 1;
            }

        Print_Real_nd(2,dims,fata,"Inputs 1 & 2");

        Real_FFT_Inverse_nd(2,dims,
            Real_Correlation_nd(2,dims,Real_FFT_nd(2,dims,fata),Real_FFT_nd(2,dims,rata)));

        Print_Real_nd(2,dims,fata,"Correlation");
      }
  }

  { double       *rata = (double *) block;
    double       *fata = rata + 2*DSIZE*DSIZE;
    int           dims[2];
    int           h, k, i;

    printf("\nMULTI-DIMENSIONAL REAL CONVOLUTION TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= 2*DSIZE; k *= 2)
      { int size;

        dims[0] = k;
        dims[1] = h;
        size    = k*h;

        for (i = 0; i < size; i++)
          { rata[i] = (i % 2) + 1;
            fata[i] = (i % 2) + 1;
          }

        Real_FFT_Inverse_nd(2,dims,
            Real_Convolution_nd(2,dims,Real_FFT_nd(2,dims,fata),Real_FFT_nd(2,dims,rata)));

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < size; i++)
            { if (size == 1)
                x = 1;
              else
                { x = 4 + (i+1)%2;
                  x = x * k * h / 2;
                }
              dev = fabs(fata[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  { double       *rata = (double *) block;
    double       *fata = rata + 2*DSIZE*DSIZE;
    int           dims[2];
    int           h, k, i;

    printf("\nMULTI-DIMENSIONAL REAL CORRELATION TEST:\n\n");

    for (h = 1; h <= DSIZE; h *= 2)
     for (k = 1; k <= 2*DSIZE; k *= 2)
      { int size;

        dims[0] = k;
        dims[1] = h;
        size    = k*h;

        for (i = 0; i < size; i++)
          { rata[i] = (i % 2) + 1;
            fata[i] = (i % 2) + 1;
          }

        Real_FFT_Inverse_nd(2,dims,
            Real_Correlation_nd(2,dims,Real_FFT_nd(2,dims,fata),Real_FFT_nd(2,dims,rata)));

        { double avg, dev, mxdev, mxdiff;
          int    x;

          avg = 0.;
          mxdiff = 0.;
          mxdev  = 0.;
          for (i = 0; i < size; i++)
            { if (size == 1)
                x = 1;
              else
                { x = 4 + (i+1)%2;
                  x = x * k * h / 2;
                }
              dev = fabs(fata[i] - x);
              if (dev > mxdiff) mxdiff = dev;
              dev /= x;
              if (dev > mxdev) mxdev = dev;
              avg += dev;
            }
          printf("For %4d,%4d: max diff = %.3e, max dev = %.3e, avg dev = %.3e\n",
                 h,k,mxdiff,mxdev,avg/(2.*k*h));
        }
      }
  }

  exit (0);

error:
  printf("Encounted an error (%s) %s\n",FFT_Error_Source(),FFT_Error_String());

  exit (0);
}
