#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fft.F.h"
#include "fft.D.h"

static double idata[10000];
static double tdata[10000];
static double cdata[10000];

static int idims[3];
static int tdims[3];
static int cdims[3];

static int cprod[3];

int compute_size(int k, int idx)
{ int p, n, s;

  p = idx / cprod[k];
  s = 0;
  n = p+tdims[k];
  if (n > cdims[k])
    { s = n-cdims[k];
      if (s > idims[k])
        s = idims[k];
    }
  if (p < idims[k])
    { if (n > idims[k])
        n = idims[k];
      s += (n-p);
    }
  if (k > 0)
    return (s*compute_size(k-1,idx%cprod[k]));
  else
    return (s);
  return (0);
}

double compute_iave(int k, int size, int idx, double *data)
{ double ave;
  int    i, p, n, s;

  p     = idx / cprod[k];
  size /= idims[k];
  ave   = 0.;
  n = p+tdims[k];
  if (n > cdims[k])
    { s = n-cdims[k];
      if (s > idims[k])
        s = idims[k];
      for (i = 0; i < s; i++)
        if (k > 0)
          ave += compute_iave(k-1,size,idx%cprod[k],data+i*size);
        else
          ave += data[i];
    }
  if (p < idims[k])
    { if (n > idims[k])
        n = idims[k];
      for (i = p; i < n; i++)
        if (k > 0)
          ave += compute_iave(k-1,size,idx%cprod[k],data+i*size);
        else
          ave += data[i];
    }
  return (ave);
}

double compute_idev(int k, int size, int idx, double *data, double ave)
{ double dev;
  int    i, p, n, s;

  p     = idx / cprod[k];
  size /= idims[k];
  dev   = 0.;
  n = p+tdims[k];
  if (n > cdims[k])
    { s = n-cdims[k];
      if (s > idims[k])
        s = idims[k];
      for (i = 0; i < s; i++)
        if (k > 0)
          dev += compute_idev(k-1,size,idx%cprod[k],data+i*size,ave);
        else
          dev += (data[i] - ave) * (data[i] - ave);
    }
  if (p < idims[k])
    { if (n > idims[k])
        n = idims[k];
      for (i = p; i < n; i++)
        if (k > 0)
          dev += compute_idev(k-1,size,idx%cprod[k],data+i*size,ave);
        else
          dev += (data[i] - ave) * (data[i] - ave);
    }
  return (dev);
}

double compute_tave(int k, int size, int idx, double *data)
{ double ave;
  int    i, p, n, s, q;

  p     = idx / cprod[k];
  size /= tdims[k];
  ave   = 0.;
  n = p+tdims[k];
  if (n > cdims[k])
    { s = n-cdims[k];
      if (s > idims[k])
        q = tdims[k] - (s-idims[k]);
      else
        q = tdims[k];
      for (i = tdims[k]-s; i < q; i++)
        if (k > 0)
          ave += compute_tave(k-1,size,idx%cprod[k],data+i*size);
        else
          ave += data[i];
    }
  if (p < idims[k])
    { if (n > idims[k])
        n = idims[k];
      for (i = 0; i < n-p; i++)
        if (k > 0)
          ave += compute_tave(k-1,size,idx%cprod[k],data+i*size);
        else
          ave += data[i];
    }
  return (ave);
}

double compute_tdev(int k, int size, int idx, double *data, double ave)
{ double dev;
  int    i, p, n, s, q;

  p     = idx / cprod[k];
  size /= tdims[k];
  dev   = 0.;
  n = p+tdims[k];
  if (n > cdims[k])
    { s = n-cdims[k];
      if (s > idims[k])
        q = tdims[k] - (s-idims[k]);
      else
        q = tdims[k];
      for (i = tdims[k]-s; i < q; i++)
        if (k > 0)
          dev += compute_tdev(k-1,size,idx%cprod[k],data+i*size,ave);
        else
          dev += (data[i] - ave) * (data[i] - ave);
    }
  if (p < idims[k])
    { if (n > idims[k])
        n = idims[k];
      for (i = 0; i < n-p; i++)
        if (k > 0)
          dev += compute_tdev(k-1,size,idx%cprod[k],data+i*size,ave);
        else
          dev += (data[i] - ave) * (data[i] - ave);
    }
  return (dev);
}

void test_result(int ndim)
{ int    k, den;
  int    isize, tsize, csize;
  double iave, idev;
  double tave, tdev;
  double rez;

  isize = 1;
  tsize = 1;
  csize = 1;
  for (k = 0; k < ndim; k++)
    { isize *= idims[k];
      tsize *= tdims[k];
      cprod[k] = csize;
      csize *= cdims[k];
    }

  for (k = 0; k < csize; k++)
    cdata[k] = 0;

  Normalize_nd(ndim,idims,idata,tdims,tdata,cdims,cdata);

  for (k = 0; k < csize; k++)

    { den = compute_size(ndim-1,k);
    
      iave = compute_iave(ndim-1,isize,k,idata)/den;
    
      idev = compute_idev(ndim-1,isize,k,idata,iave)/den;
    
      tave = compute_tave(ndim-1,tsize,k,tdata)/den;
    
      tdev = compute_tdev(ndim-1,tsize,k,tdata,tave)/den;
    
      if (tdev == 0. || idev == 0.)
        rez = 0.;
      else
        rez = -(den*iave*tave) / ((den-1)*sqrt(idev*tdev));
    
      if (fabs(rez-cdata[k]) > 1e-10)
        printf("** DIFF %d: [%g %g %g %g] %g %g (%g)\n",
               k,iave,idev,tave,tdev,cdata[k],rez,fabs(rez-cdata[k]));
    }
}

int main(int argc, char *argv[])
{
  (void) argc;
  (void) argv;

  { int k;

    for (k = 0; k < 10000; k++)
      idata[k] = tdata[k] = (double) ((k%5)+1);

    idims[0] = 3;
    idims[1] = 3;
    idims[2] = 3;

    tdims[0] = 6;
    tdims[1] = 6;
    tdims[2] = 6;

    cdims[0] = 9;
    cdims[1] = 9;
    cdims[2] = 9;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 10;
    cdims[1] = 10;
    cdims[2] = 10;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 8;
    cdims[1] = 8;
    cdims[2] = 8;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 6;
    cdims[1] = 6;
    cdims[2] = 6;

    test_result(1);
    test_result(2);
    test_result(3);

    tdims[0] = 3;
    tdims[1] = 3;
    tdims[2] = 3;

    idims[0] = 6;
    idims[1] = 6;
    idims[2] = 6;

    cdims[0] = 9;
    cdims[1] = 9;
    cdims[2] = 9;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 10;
    cdims[1] = 10;
    cdims[2] = 10;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 8;
    cdims[1] = 8;
    cdims[2] = 8;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 6;
    cdims[1] = 6;
    cdims[2] = 6;

    test_result(1);
    test_result(2);
    test_result(3);

    idims[0] = 5;
    idims[1] = 5;
    idims[2] = 5;

    tdims[0] = 5;
    tdims[1] = 5;
    tdims[2] = 5;

    cdims[0] = 10;
    cdims[1] = 10;
    cdims[2] = 10;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 11;
    cdims[1] = 11;
    cdims[2] = 11;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 6;
    cdims[1] = 6;
    cdims[2] = 6;

    test_result(1);
    test_result(2);
    test_result(3);

    cdims[0] = 5;
    cdims[1] = 5;
    cdims[2] = 5;

    test_result(1);
    test_result(2);
    test_result(3);

    idims[0] = 1;
    idims[1] = 1;
    idims[2] = 1;

    tdims[0] = 1;
    tdims[1] = 1;
    tdims[2] = 1;

    cdims[0] = 2;
    cdims[1] = 2;
    cdims[2] = 2;

    test_result(1);
    test_result(2);
    test_result(3);
  }

  exit (0);
}
