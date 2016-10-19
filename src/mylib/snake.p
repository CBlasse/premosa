#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utilities.h"
#include "array.h"
#include "linear.algebra.h"

typedef enum { FIXED_SNAKE, VARIABLE_SNAKE } Snake_Type;

typedef struct
  { Snake_Type    type;
    int           len;
    Float_Matrix *curve;
    Float_Matrix *params;
  } Snake;

double *get_dist_vector(int n, char *routine)
{ static double *Distance = NULL;
  static int     Dist_Max = 0;

  if (n > Dist_Max)
    { Dist_Max = (int) (1.2*n + 100);
      Distance = (double *) Guarded_Realloc(Distance,sizeof(double)*((size_t) Dist_Max),routine);
    }
  else if (n <= 0)
    { free(Distance);
      Distance = NULL;
      Dist_Max = 0;
    }
  return (Distance);
}

static int Vector_Top = 0;

MANAGER -IO Snake curve*Array params*Array

void resample(Double_Matrix *old, Double_Matrix *new, double *dist, double grid)
{ double *xn = AFLOAT64(new);
  double *yn = xn + new->dims[0];
  double *xo = AFLOAT64(old);
  double *yo = xo + old->dims[0];
  int     np = new->dims[0]-1;
  double  cx, cy;
  double  nx, ny;
  double  p, g;
  int     i, k;

  k = 0;
  xn[0] = xn[np] = cx = xo[-1];
  yn[0] = yn[np] = cy = yo[-1];
  nx = xo[0];
  ny = yo[0];
  for (i = 1, p = grid; i < np; i++, p += grid)
    { while (dist[k] < p)
        { p -= dist[k];
          k += 1;
          cx = nx;
          cy = ny;
          nx = xo[k];
          ny = yo[k];
        }
      g = p / dist[k];
      xn[i] = cx + g * (nx - cx);
      yn[i] = cy + g * (ny - cy);
    }
}

Snake *resample_snake(Snake *snake, double grid)
{ static Dimn_Type dims[2];

  int     n = snake->len;
  double *x = AFLOAT64(snake->curve);
  double *y = x + n;

  double *dist, len;
  int     np;
  Snake  *resamp;

  dist = get_dist_vector(n,"resample_snake");

  { double  lx, ly;
    double  cx, cy;
    double  dx, dy;
    double *dm;
    int     i;

    dm  = dist-1;
    lx  = x[0];
    ly  = y[0];
    len = 0.;
    for (i = 1; i <= n; i++)
      { cx = x[i];
        cy = y[i];
        dx = cx - lx;
        dy = cy - ly;
        lx = cx;
        ly = cy;
        len += dm[i] = sqrt(dx*dx + dy*dy);
      }
  }

  np   = (int) (len / grid);
  grid = len / np;

  if (np+1 > Vector_Top)
    Vector_Top = (int) (1.1*(np+1) + 100);

  resamp = new_snake("resample_snake");
  dims[0] = np+1;
  dims[1] = 2;
  resamp->curve = Make_Array(PLAIN_KIND,FLOAT64_TYPE,2,dims);

  if (snake->type == VARIABLE_SNAKE)
    { resamp->params = Make_Array(PLAIN_KIND,FLOAT64_TYPE,2,dims);
      resamp->type   = VARIABLE_SNAKE;
    }
  else
    { resamp->params = Inc_Array(snake->params);
      resamp->type   = FIXED_SNAKE;
    }

  resample(snake->curve,resamp->curve,dist,grid);
  if (snake->type == VARIABLE_SNAKE)
    resample(snake->params,resamp->params,dist,grid);

  return (resamp);
}

void deform_snake(Snake *snake, double gamma, double kappa, Double_Matrix *force)
{ static Double_Matrix *penta = NULL;
  static int            pmax = 0;
  static Dimn_Type      dims[2];

  int width = force->dims[0];
  int area  = force->dims[1] * width;

  double *fx = AFLOAT64(force);
  double *fy = fx + area;
  int      n = snake->len;

  double *a, *b, *c, *d, *e;

  if (n > pmax)
    { pmax = (int) (1.2*n + 300);
      if (penta != NULL)
        Free_Array(penta);
      dims[0] = pmax;
      dims[1] = 5;
      penta = Make_Array(PLAIN_KIND,FLOAT64_TYPE,2,dims);
    }

  penta->dims[0] = n;
  a = AFLOAT64(penta);
  b = a + n;
  c = b + n;
  d = c + n;
  e = d + n;

  if (snake->type == VARIABLE_SNAKE)
    { double *alpha = AFLOAT64(snake->params);
      double *beta  = alpha + snake->len;
      int     i;

      a[0] = alpha[n-1];
      b[0] = - (alpha[0] + 2*beta[n-1] + 2*beta[0]);
      c[0] = alpha[0] + alpha[1] + beta[n-1] + 4*beta[0] + beta[1];
      d[0] = - (alpha[1] + 2*beta[0] + 2*beta[1]);
      e[0] = beta[1];
      for (i = 1; i < n; i++)
        { a[i] = beta[i-1];
          b[i] = - (alpha[i] + 2*beta[i-1] + 2*beta[i]);
          c[i] = alpha[i] + alpha[i+1] + beta[i-1] + 4*beta[i] + beta[i+1];
          d[i] = - (alpha[i+1] + 2*beta[i] + 2*beta[i+1]);
          e[i] = beta[i+1];
        }
    }
  else
    { double alpha = AFLOAT64(snake->params)[0];
      double beta  = AFLOAT64(snake->params)[1];
      int   i;

      for (i = 0; i < n; i++)
        { a[i] = e[i] = beta;
          b[i] = d[i] = - (alpha + 4*beta);
          c[i] = 2*alpha + 6*beta;
        }
    }

  { Band_Factor  *lu = Pentaband_Decompose(penta);
    double       *x  = AFLOAT64(snake->curve);
    double       *y  = x + snake->len;
    Array_Bundle  abundle;
    int           i;

    for (i = 0; i < n; i++)
      { double u, v;
        double xi, yi;
        double a, b;
        int    r, s, p;

        r = (int) (xi = x[i]);
        s = (int) (yi = y[i]);

        p = s*width + r;

        a = (xi - r);
        b = (yi - s);

        u = (fx[p] * a + fx[p+1] * (1.-a)) * b
          + (fx[p+width] * a + fx[p+width+1] * (1.-a)) * (1.-b);
        v = (fy[p] * a + fy[p+1] * (1.-a)) * b
          + (fy[p+width] * a + fy[p+width+1] * (1.-a)) * (1.-b);

        x[i] = gamma * xi + kappa * u;
        y[i] = gamma * yi + kappa * v;
      }

    abundle = *(snake->curve);
    Pentaband_Solve(Get_Array_Plane(&abundle,0),lu);
    abundle = *(snake->curve);
    Pentaband_Solve(Get_Array_Plane(&abundle,1),lu);

    Free_Band_Factor(lu);
  }
}
