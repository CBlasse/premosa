/*****************************************************************************************\
*                                                                                         *
*  Histogram Data Abstraction and Array Statistics Routines                               *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  December 2008                                                                 *
*                                                                                         *
*  (c) December 20, '09, Dr. Gene Myers and Howard Hughes Medical Institute               *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>

#include "utilities.h"
#include "array.h"
#include "histogram.h"

#LISTDEF @BLUE = 1 3 42 -3 6

#LISTDEF @TYPES  =  UINT8 UINT16 UINT32 UINT64  INT8 INT16 INT32 INT64 FLOAT32 FLOAT64
#LISTDEF @UNION  =   uval   uval   uval   uval  ival  ival  ival  ival    fval    fval

#define SIZEOF(x) ((int) sizeof(x))


/****************************************************************************************
 *                                                                                      *
 *  HISTOGRAM SPACE MANAGEMENT ROUTINES AND PRIMARY GENERATOR                           *
 *                                                                                      *
 ****************************************************************************************/

// Awk-generated (manager.awk) Array memory management

typedef struct
  { Value_Kind kind;    //  Exactly the same as a Histogram save the extra field for mean.
    Value      binsize;
    Value      offset;
    int        nbins;
    Size_Type  total;
    Size_Type *counts;
    double     mean;    //  histogram mean
  } Histofull;

static Value_Kind type2kind[] = { UVAL, UVAL, UVAL, UVAL, IVAL, IVAL, IVAL, IVAL, FVAL, FVAL };

static inline int histofull_nsize(Histofull *h)
{ return (SIZEOF(Size_Type)*h->nbins); }

MANAGER -IO Histogram(Histofull) counts:nsize

Histogram *G(Make_Histogram)(Value_Kind kind, int nbins, Value binsize, Value offset)
{ Histofull *h     = new_histofull(nbins*SIZEOF(Size_Type),"Histogram_Array");
  Size_Type *count = h->counts;
  Indx_Type  p;

  h->kind    = kind;
  h->binsize = binsize;
  h->offset  = offset;
  h->nbins   = nbins;

  h->mean    = DBL_MIN;
  h->total   = 0;
  for (p = 0; p < nbins; p++)
    count[p] = 0;

  return ((Histogram *) h);
}

Histogram *Empty_Histogram(Histogram *h)
{ int        nbins = h->nbins;
  Size_Type *count = h->counts;
  Indx_Type  p;

  for (p = 0; p < nbins; p++)
    count[p] = 0;
  h->total = 0;
  ((Histofull *) h)->mean = DBL_MIN;
  return (h);
}

  /*  Generate a histogram of array a with nbins of width binsize where the smallest bin's
      lower boundary is offset (see data descriptor comments in .h file).  The type of values   
      given for binsize and offset should be congruent with the type of a.  When nbins or 
      binsize or both are zero then the histogram bins are set up as follows based on the 
      range [min,max] of values in a: 

         nbins = 0 & binsize = 0:
            bins are of length *1* and cover the range of values
              in a starting at *floor(min)*

         nbins = 0 & binsize != 0:
            bins are of length *binsize* and cover the range of values
              in a starting at *floor(min)*

         nbins != 0 & binsize = 0:
            The bin size is the smallest number of the form [1,2,5]*10^a for which nbins
            of this size cover the range [min,max].  The bins start at the first multiple
            of the bin size <= min and nbins is adjusted downwards so that the binning
            just covers the required range [min,max].

         nbins != 0 & binsize != 0
            The implied bining is used as specified and any values not in the implied range
            are not added to the histogram, i.e. the total count of the histogram can be less
            then the size of a.
  */

#define FILL_BINS_F( EXPR, ADVANCE )	\
  if (clip)				\
    for (p = 0; p < size; p++)		\
      { int i = (int) ((EXPR-o)/b);	\
        if (0 <= i && i < nbins)	\
          count[i] += 1;		\
        ADVANCE				\
      }					\
  else if (b == 1)			\
    if (o < 0)				\
      for (p = 0; p < size; p++)	\
        { count[(int) (EXPR - o)] += 1;	\
          ADVANCE			\
        }				\
    else				\
      { count -= (int) o;		\
        for (p = 0; p < size; p++)	\
          { count[(int) EXPR] += 1;	\
            ADVANCE			\
          }				\
      }					\
  else if (o == 0)			\
    for (p = 0; p < size; p++)		\
      { count[(int) (EXPR/b)] += 1;	\
        ADVANCE				\
      }					\
  else					\
    for (p = 0; p < size; p++)		\
      { count[(int) ((EXPR-o)/b)] += 1;	\
        ADVANCE				\
      }

#define FILL_BINS_D( EXPR, ADVANCE )	\
  if (clip)				\
    for (p = 0; p < size; p++)		\
      { int i = (int) ((EXPR-o)/b);	\
        if (0 <= i && i < nbins)	\
          count[i] += 1;		\
        ADVANCE				\
      }					\
  else if (b == 1)			\
    { count -= o;			\
      for (p = 0; p < size; p++)	\
        { count[EXPR] += 1;		\
          ADVANCE			\
        }				\
    }					\
  else if (o == 0)			\
    for (p = 0; p < size; p++)		\
      { count[EXPR/b] += 1;		\
        ADVANCE				\
      }					\
  else					\
    for (p = 0; p < size; p++)		\
      { count[(EXPR-o)/b] += 1;		\
        ADVANCE				\
      }

#define FILL_BINS_FREG				\
  if (clip)					\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          { int i = (int) ((v[p]-o)/b);		\
            if (0 <= i && i < nbins)		\
              count[i] += 1;			\
          }					\
      }						\
  else if (b == 1)				\
    if (o < 0)					\
      for (k = 0; k < len; k += 2)		\
        { vr = raster[k];			\
          wr = raster[k+1];			\
          for (p = vr; p <= wr; p++)		\
            count[(int) (v[p] - o)] += 1;	\
        }					\
    else					\
      { count -= (int) o;			\
        for (k = 0; k < len; k += 2)		\
          { vr = raster[k];			\
            wr = raster[k+1];			\
            for (p = vr; p <= wr; p++)		\
              count[(int) v[p]] += 1;		\
          }					\
      }						\
  else if (o == 0)				\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          count[(int) (v[p]/b)] += 1;		\
      }						\
  else						\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          count[(int) ((v[p]-o)/b)] += 1;	\
      }

#define FILL_BINS_DREG				\
  if (clip)					\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          { int i = (int) ((v[p]-o)/b);		\
            if (0 <= i && i < nbins)		\
              count[i] += 1;			\
          }					\
      }						\
  else if (b == 1)				\
    { count -= o;				\
      for (k = 0; k < len; k += 2)		\
        { vr = raster[k];			\
          wr = raster[k+1];			\
          for (p = vr; p <= wr; p++)		\
            count[v[p]] += 1;		 	\
        }					\
    }						\
  else if (o == 0)				\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          count[v[p]/b] += 1;			\
      }						\
  else						\
    for (k = 0; k < len; k += 2)		\
      { vr = raster[k];				\
        wr = raster[k+1];			\
        for (p = vr; p <= wr; p++)		\
          count[(v[p]-o)/b] += 1;		\
      }

#GENERATE N = Array Region

#IF N == Array
Histogram *G(Histogram_Array)(AForm *form, int nbins, Value binsize, Value offset)
#ELSE
Histogram *G(Histogram_Region)(Array *form, Region *reg, int nbins, Value binsize, Value offset)
#END
{ static uint64  MinOff_uval[]   = {     0,       0 };
  static uint64  FullSpan_uval[] = { 0x100, 0x10000, 0, 0 };
  static int64   MinOff_ival[]   = { 0, 0, 0, 0,  0xff,  0xffff };
  static int64   FullSpan_ival[] = { 0, 0, 0, 0, 0x100, 0x10000, 0, 0, };
  static float64 MinOff_fval[]   = { 0 };
  static float64 FullSpan_fval[] = { 0 };

  Array        *a     = AForm_Array(form);
  boolean       clip  = 0;
  Range_Bundle  rng;

  switch (type2kind[a->type]) {
    #GENERATE t,u,f = uval ival fval , uint64 int64 float64 , llu lld g
      case <T>:
        if (nbins == 0)
          { uint64 rti;

          #IF N == Array
            Array_Range(&rng,form);
          #ELSE
            Region_Range(&rng,form,reg);
          #END
            if (binsize.<t> == 0)
              binsize.<t> = 1;
            #IF t == fval
              offset.<t> = floor(rng.minval.<t> / binsize.<t>) * binsize.<t>;
            #ELSEIF t == ival
              if (rng.minval.<t> < 0)
                offset.<t> = ((rng.minval.<t>+1) / binsize.<t> - 1) * binsize.<t>;
              else
                offset.<t> = (rng.minval.<t> / binsize.<t>) * binsize.<t>;
            #ELSE
              offset.<t> = (rng.minval.<t> / binsize.<t>) * binsize.<t>;
            #END
            rti = (uint64) ((rng.maxval.<t> - offset.<t>) / binsize.<t>);
            if (offset.<t> + ((<u>) nbins)*binsize.<t> <= rng.maxval.<t>)
              rti += 1;
            if (rti > 0x7FFFFFFFull)
              { fprintf(stderr,
                        "Implied binning requires more than 2 billion bins (Histogram_<N>)\n");
                exit (1);
              }
            nbins = (int) rti;
          }
        else if (binsize.<t> == 0)
          { double bwide;

          #IF N == Array
            Array_Range(&rng,form);
          #ELSE
            Region_Range(&rng,form,reg);
          #END
            bwide = (rng.maxval.<t> - rng.minval.<t>) / (1.*nbins);
            if (bwide == 0.)
              binsize.<t> = 1;
            else
              #IF t == fval
              { double x = pow(10.,floor(log10(bwide)));
              #ELSE
              { <u> x = 1;
                while (10*x <= bwide)
                  x = 10*x;
              #END
                if (x < bwide)
                  { if (2*x < bwide)
                      { if (5*x < bwide)
                          x = 10*x;
                        else
                          x = 5*x;
                      }
                    else
                      x = 2*x;
                  }
                binsize.<t> = x;
              }
            #IF t == fval
              offset.<t> = floor(rng.minval.<t> / binsize.<t>) * binsize.<t>;
            #ELSEIF t == ival
              if (rng.minval.<t> < 0)
                offset.<t> = ((rng.minval.<t>+1) / binsize.<t> - 1) * binsize.<t>;
              else
            #ELSE
                offset.<t> = (rng.minval.<t> / binsize.<t>) * binsize.<t>;
            #END
            nbins = (int) ((rng.maxval.<t> - offset.<t>) / binsize.<t>);
            if (offset.<t> + ((<u>) nbins)*binsize.<t> <= rng.maxval.<t>)
              nbins += 1;
          } 
        else
          { if (FullSpan_<t>[a->type] == 0 || offset.<t> > MinOff_<t>[a->type] ||
                offset.<t> + ((<u>) nbins)*binsize.<t> < FullSpan_<t>[a->type])
            #IF N == Array
              { Array_Range(&rng,form);
            #ELSE
              { Region_Range(&rng,form,reg);
            #END
                clip = (offset.<t> > rng.minval.<t> ||
                        offset.<t> + binsize.<t> * ((<u>) nbins) <= rng.maxval.<t>);
              }
          }
        break;
    #END
  }

  { Histofull *h     = new_histofull(nbins*SIZEOF(Size_Type),"Histogram_Array");
    Size_Type *count = h->counts;
    Size_Type  size  = AForm_Size(form);
    Indx_Type  p;

    h->kind    = type2kind[a->type];
    h->binsize = binsize;
    h->offset  = offset;
    h->nbins   = nbins;

    for (p = 0; p < nbins; p++)
      count[p] = 0;

    switch (a->type) {
      #GENERATE T,U = @TYPES , @UNION
        case <T>_TYPE:
          { <t> *v = A<T>(a);
            <t>  o = (<t>) offset.<U>;
            <t>  b = (<t>) binsize.<U>;
            #IF N == Array
              switch (AForm_Class(form))
              { case FRAME_CLASS:
                  if (Frame_Within_Array(form))
                    { Offs_Type *off = Frame_Offsets(form);
                      v += Frame_Index(form);
                      #IF T >= FLOAT32
                        FILL_BINS_F( v[off[p]], )
                      #ELSE
                        FILL_BINS_D( v[off[p]], )
                      #END
                      break;
                    }
                  else
                    v = Frame_Values(form);
                case ARRAY_CLASS:
                  #IF T >= FLOAT32
                    FILL_BINS_F( v[p], )
                  #ELSE
                    FILL_BINS_D( v[p], )
                  #END
                  break;
                case SLICE_CLASS:
                  { Indx_Type e = Set_Slice_To_First(form);
                    #IF T >= FLOAT32
                      FILL_BINS_F( v[e], e = Next_Slice_Index(form); )
                    #ELSE
                      FILL_BINS_D( v[e], e = Next_Slice_Index(form); )
                    #END
                    break;
                  }
              }
            #ELSE
              Indx_Type *raster = reg->raster;
              Size_Type  len    = reg->rastlen;
              Indx_Type  vr, wr, k;
              #IF T >= FLOAT32
                FILL_BINS_FREG
              #ELSE
                FILL_BINS_DREG
              #END
            #END
            break;
	  }
      #END
    }
#IF N == Array
    if (clip)
#END
      { size = 0;
        for (p = 0; p < h->nbins; p++)
          size += count[p];
      }

    h->total = size;
    h->mean  = DBL_MIN;
    return ((Histogram *) h);
  }
}

#IF N == Array
Histogram *Histagain_Array(Histogram *h, AForm *form, boolean clip)
#ELSE
Histogram *Histagain_Region(Histogram *h, Array *form, Region *reg, boolean clip)
#END
{ Value      offset  = h->offset;
  Value      binsize = h->binsize;
  Size_Type  size    = AForm_Size(form);
  Size_Type *count   = h->counts;
  int        nbins   = h->nbins;
  Array     *a       = AForm_Array(form);
  Indx_Type  p;

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *v = A<T>(a);
          <t>  o = (<t>) offset.<U>;
          <t>  b = (<t>) binsize.<U>;
          #IF N == Array
            switch (AForm_Class(form))
            { case FRAME_CLASS:
                if (Frame_Within_Array(form))
                  { Offs_Type *off = Frame_Offsets(form);
                    v += Frame_Index(form);
                    #IF T >= FLOAT32
                      FILL_BINS_F( v[off[p]], )
                    #ELSE
                      FILL_BINS_D( v[off[p]], )
                    #END
                    break;
                  }
                else
                  v = Frame_Values(form);
              case ARRAY_CLASS:
                #IF T >= FLOAT32
                  FILL_BINS_F( v[p], )
                #ELSE
                  FILL_BINS_D( v[p], )
                #END
                break;
              case SLICE_CLASS:
                { Indx_Type e = Set_Slice_To_First(form);
                  #IF T >= FLOAT32
                    FILL_BINS_F( v[e], e = Next_Slice_Index(form); )
                  #ELSE
                    FILL_BINS_D( v[e], e = Next_Slice_Index(form); )
                  #END
                  break;
                }
            }
          #ELSE
            Indx_Type *raster = reg->raster;
            Size_Type  len    = reg->rastlen;
            Indx_Type  vr, wr, k;
            #IF T >= FLOAT32
              FILL_BINS_FREG
            #ELSE
              FILL_BINS_DREG
            #END
          #END
          break;
        }
    #END
  }

  size = 0;
  for (p = 0; p < h->nbins; p++)
    size += count[p];
  h->total = size;
  ((Histofull *) h)->mean = DBL_MIN;

  return (h);
}

#END

Histogram *G(Histogram_P_Vertex)(Array *a, Partition *w, int cb,
                                 int nbins, Value binsize, Value offset)
{ Region    *c;
  Histogram *h;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Histogram_P_Vertex)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,0,1);
  h = Histogram_Region(a,c,nbins,binsize,offset);
  Free_Region(c);
  return (h);
}

Histogram *G(Histogram_Level_Set)(Array *a, Level_Tree *t, Level_Set *r,
                                  int nbins, Value binsize, Value offset)
{ Region    *c;
  Histogram *h;
  c = Record_Level_Set(t,r,0,1);
  h = Histogram_Region(a,c,nbins,binsize,offset);
  Free_Region(c);
  return (h);
}

Histogram *Histagain_P_Vertex(Histogram *R(M(h)), Array *a, Partition *w, int cb, boolean clip)
{ Region    *c;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Histagain_P_Vertex)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,0,1);
  Histagain_Region(h,a,c,clip);
  Free_Region(c);
  return (h);
}

Histogram *Histagain_Level_Set(Histogram *R(M(h)), Array *a, Level_Tree *t, Level_Set *r,
                               boolean clip)
{ Region    *c;
  c = Record_Level_Set(t,r,0,1);
  Histagain_Region(h,a,c,clip);
  Free_Region(c);
  return (h);
}

  /*  Generate a histogram based on the histogram h consisting of the bins in the interval
        [min,max).  min and max need not be between 0 and h->nbins but it must be that min < max.
        If min < 0 or max > h->nbins then g's domain will be expanded as necessary to cover
        the implied range of [Bin2Value(min),Bin2Value(max)].
  */

Histogram *G(Histogram_Slice)(Histogram *h, int min, int max)
{ Histofull *g;
  int        nbins;
  Value      off;
  int        i;

  switch (h->kind) {
    #GENERATE t,u = uval ival fval , uint64 int64 float64
      case <T>:
        off.<t> = h->offset.<t> + h->binsize.<t> * ((<u>) min);
        break;
    #END
  }

  nbins = max-min;
  if (nbins <= 0)
    { fprintf(stderr,"Requested bin interval is empty (Histogram_Slice)");
      exit (1);
    }

  g = new_histofull(nbins*SIZEOF(Size_Type),"Histogram_Slice");

  g->kind    = h->kind;
  g->binsize = h->binsize;
  g->offset  = off;
  g->nbins   = nbins;
  for (i = min; i < max; i++)
    if (i < 0 || i >= h->nbins)
      g->counts[i-min] = 0;
    else
      g->counts[i-min] = h->counts[i];

  g->mean = DBL_MIN;

  return ((Histogram *) g);
}


/****************************************************************************************
 *                                                                                      *
 *  BIN INDICES, DOMAIN VALUES, AND PERCENTILES                                         *
 *                                                                                      *
 ****************************************************************************************/

  /*  Routines to map between bin indices, domain values, and percentiles:

        Bin2Value: offset + b*binsize.
        Value2Bin: max_i s.t. Bin2Value(i) <= v

        Bin2Percentile: sum_(j>=i) count[j] / total
        Percentile2Bin: max_i s.t. Bin2Percentile(i) >= fraction

        Value2Percentile: Bin2Percentile(j) - (v-Value2Bin(j))*count[j]/total,
                            for j = Value2Bin(v)
        Percentile2Value: max_v s.t. Value2Percentile(v) >= fraction

      The bin index parameters do not need to be between 0 and h->nbins-1 and the bin index
      returned by Value2Bin may not be either depending on the value of v.  The fraction
      parameter however must be between 0 and 1, and values and bin indices returned by
      the percentile routines are always in range.
  */

Value Bin2Value(Histogram *h, int b)
{ Value v;
  switch (h->kind) {
    #GENERATE t,u = uval ival fval , uint64 int64 float64
      case <T>:
        v.<t> = h->offset.<t> + h->binsize.<t> * ((<u>) b);
        break;
    #END
  }
  return (v);
}

int Value2Bin(Histogram *h, Value v)
{ int bck;
  switch (h->kind) {
    #GENERATE t,u = uval ival fval , uint64 int64 float64
      case <T>:
        { <u> o = h->offset.<t>;
          <u> b = h->binsize.<t>;
        #IF t == uval
          if (v.<t> < o)
            bck = - (int) ((o - v.<t>) / b);
          else
        #END
            bck = (int) ((v.<t> - o) / b);
          break;
        }
    #END
      default:
        bck = 0;
        break;
  }
  return (bck);
}

double Bin2Percentile(Histogram *h, int b)
{ Size_Type *count;
  int64      sum;
  int        i;

  count = h->counts;
  sum   = 0;
  for (i = b; i < h->nbins; i++)
    sum += count[i];
  return ((1.*sum)/h->total);
}

int Percentile2Bin(Histogram *h, double fraction)
{ Size_Type *count;
  int64      cthr, sum;
  int        i;

  cthr  = (int64) (h->total * fraction);
  count = h->counts;
  sum   = 0;
  if (cthr <= 0)
    return (h->nbins);
  for (i = h->nbins-1; i > 0; i--)
    { sum += count[i];
      if (sum >= cthr) break;
    }
  return (i);
}

double Value2Percentile(Histogram *h, Value v)
{ int b    = Value2Bin(h,v);
  switch (h->kind) {
    #GENERATE t,u,v = uval ival fval , uint64 int64 float64 , int64 int64 float64
      case <T>:
        return (Bin2Percentile(h,b) - 
                   (((<v>) (v.<t>-Bin2Value(h,b).<t>))*h->counts[b])/h->total);
    #END
  }
  return (0.);
}

Value Percentile2Value(Histogram *h, double fraction)
{ Size_Type *count;
  int64      cthr, sum;
  int        i;
  Value      v;

  cthr  = (int64) (h->total * fraction);
  count = h->counts;
  sum   = 0;
  if (cthr <= 0)
    return (Bin2Value(h,h->nbins));
  if (cthr >= h->total)
    return (Bin2Value(h,0));
  for (i = h->nbins-1; i > 0; i--)
    { sum += count[i];
      if (sum >= cthr) break;
    }
  switch (h->kind) {
    #GENERATE t,u = uval ival fval , uint64 int64 float64
      case <T>:
        v.<t> = h->offset.<t> + h->binsize.<t> * ((<u>) i) +
                (<u>) (h->binsize.<t> * (1.*sum-cthr)/count[i]);
        break;
    #END
  }
  return (v);
}


/****************************************************************************************
 *                                                                                      *
 *  HISTOGRAM STATISTICS                                                                *
 *                                                                                      *
 ****************************************************************************************/

static void set_histogram_mean(Histofull *h)
{ Size_Type *count;
  int        i;
  double     sum, b, u;

  count = h->counts;
  
  switch (h->kind) {
    #GENERATE t = uval ival fval
      case <T>:
        b = (double) h->binsize.<t>;
        #IF t == fval
          u = h->offset.<t> + .5*b;
        #ELSE
          u = h->offset.<t> + .5*b - .5;
        #END
        break;
    #END
      default:
        b = u = 0.;
        break;
  }

  sum  = 0.;
  for (i = 0; i < h->nbins; i++)
    { sum  += count[i] * u;
      u    += b;
    }

  h->mean  = sum / h->total;
}

double Histogram_Mean(Histogram *h)
{ if (((Histofull *) h)->mean == DBL_MIN)
    set_histogram_mean(((Histofull *) h));
  return (((Histofull *) h)->mean);
}

double Histogram_Variance(Histogram *h)
{ int        i;
  double     sum, b, u;
  Size_Type *count;

  if (((Histofull *) h)->mean == DBL_MIN)
    set_histogram_mean(((Histofull *) h));
  
  count = h->counts;

  switch (h->kind) {
    #GENERATE t = uval ival fval
      case <T>:
        b = (double) h->binsize.<t>;
        #IF t == fval
          u = (h->offset.<t> + .5*b) - ((Histofull *) h)->mean;
        #ELSE
          u = (h->offset.<t> + .5*b - .5) - ((Histofull *) h)->mean;
        #END
        break;
    #END
      default:
        u = b = 0.;
        break;
  }

  sum = 0.;
  for (i = 0; i < h->nbins; i++)
    { sum += count[i] * u * u;
      u += b;
    }
  return (sum/h->total);
}

double Histogram_Sigma(Histogram *h)
{ return (sqrt(Histogram_Variance(h))); }

double Histogram_Central_Moment(Histogram *h, int n)
{ int        i;
  double     sum, b, u;
  Size_Type *count;

  if (((Histofull *) h)->mean == DBL_MIN)
    set_histogram_mean(((Histofull *) h));
  
  count = h->counts;

  switch (h->kind) {
    #GENERATE t = uval ival fval
      case <T>:
        b = (double) h->binsize.<t>;
        #IF t == fval
          u = (h->offset.<t> + .5*b) - ((Histofull *) h)->mean;
        #ELSE
          u = (h->offset.<t> + .5*b - .5) - ((Histofull *) h)->mean;
        #END
        break;
    #END
      default:
        u = b = 0.;
        break;
  }

  sum = 0.;
  for (i = 0; i < h->nbins; i++)
    { sum += count[i] * pow(u,n);
      u   += b;
    }
  return (sum/h->total);
}


/****************************************************************************************
 *                                                                                      *
 *  HISTOGRAM ENTROPY                                                                   *
 *                                                                                      *
 ****************************************************************************************/

   /*  Assuming the histogram h is a discrete probaility distribution as defined by the
         choice of bin size, Histogram_Entropy returns - sum_b p(b) log_2 p(b) where
         p(b) is counts[b]/total for each bin b.  Cross_Entropy is - sum_b p(b) log_2 q(b)
         where q(b) is the distribution for g.  The histograms h and g must have the same
         bin size and while their offsets can be different, the difference must be a
         multiple of the bin size.  Relative_Entropy is sum_b p(b) log_2 p(b)/q(b);
   */

static double  etable[100001];
static double  loge5, loge10, loge15;
static boolean firstentropy = 1;

static pthread_mutex_t Entropy_Mutex = PTHREAD_MUTEX_INITIALIZER;

static inline void mylog2_table()
{ int i;

#ifdef _MSC_VER

  double logfac = 1./log(2.);

#define log2(e) (log(e)*logfac);

#endif

  for (i = 0; i <= 100000; i++)
    etable[i] = log2(i/100000.);
  loge5  = log2(1.e-5);
  loge10 = 2*loge5;
  loge15 = loge5+loge10;
}

static inline double mylog2(double p)
{ if (p >= 1e-5)
    return (etable[(int) (p*1.e5)]);
  else if (p >= 1e-10)
    return (loge5 + etable[(int) (p*1.e10)]);
  else if (p >= 1e-15)
    return (loge10 + etable[(int) (p*1.e15)]);
  else
    return (loge15 + etable[(int) (p*1.e20)]);
}

double Histogram_Entropy(Histogram *h)
{ Size_Type *count = h->counts;
  double     normal = 1./h->total;
  double     entropy;
  int        i;

  pthread_mutex_lock(&Entropy_Mutex);
  if (firstentropy)
    { firstentropy = 0;
      mylog2_table();
    }
  pthread_mutex_unlock(&Entropy_Mutex);

  entropy = 0.;
  for (i = 0; i < h->nbins; i++)
    { double p = count[i] * normal;
      if (p > 1.e-20)
        entropy -= p*mylog2(p);
    }
  return (entropy);
}

double Histogram_Cross_Entropy(Histogram *h, Histogram *g)
{ Size_Type *hcount = h->counts;
  Size_Type *gcount = g->counts;
  double     hnormal = 1./h->total;
  double     gnormal = 1./g->total;
  double     entropy;
  int64      disp;
  int        i, j;

  switch (h->kind) {
    #GENERATE t,u = uval ival fval , uint64 int64 float64
      case <T>:
        { <u> delt;
          if (h->binsize.<t> != g->binsize.<t>)
            { fprintf(stderr,"Histogram do not have same bin size (Histogram_Cross_Entropy)\n");
              exit (1);
            }
      #IF t != fval
        #IF t == uval
          if (h->offset.<t> >= g->offset.<t>)
            disp = (int64) (delt = h->offset.<t> - g->offset.<t>);
          else
            disp = - (int64) (delt = g->offset.<t> - h->offset.<t>);
        #ELSE
          disp = delt = h->offset.ival - g->offset.ival;
        #END
          if (delt % h->binsize.<t> != 0)
            { fprintf(stderr,"Histogram bin offsets not in synch (Histogram_Cross_Entropy)\n");
              exit (1);
            }
          disp /= h->binsize.<t>;
      #ELSE
          delt = (h->offset.<t> - g->offset.<t>) / h->binsize.<t>;
          disp = (int64) delt;
          if (fabs(disp - delt) > 1.e-10)
            { fprintf(stderr,"Histogram bin offsets not in synch (Histogram_Cross_Entropy)\n");
              exit (1);
            }
      #END
          break;
        }
    #END
      default:
        disp = 0;
        break;
    }

  pthread_mutex_lock(&Entropy_Mutex);
  if (firstentropy)
    { firstentropy = 0;
      mylog2_table();
    }
  pthread_mutex_unlock(&Entropy_Mutex);

  entropy = 0.;
  for (i = 0, j = (int) disp; i < h->nbins; i++, j++)
    { double p = hcount[i] * hnormal;
      double q;
      if (0 <= j && j < g->nbins)
        { q = gcount[j] * gnormal;
          if (q < 1.e-20)
            q = 1.e-20;
        }
      else
        q = 1.e-20;
      entropy -= p*mylog2(q);
    }

  return (entropy);
}

double Histogram_Relative_Entropy(Histogram *h, Histogram *g)
{ return (Histogram_Cross_Entropy(h,g) - Histogram_Entropy(h)); }


/****************************************************************************************
 *                                                                                      *
 *  HISTOGRAM THRESHOLDS                                                                *
 *                                                                                      *
 ****************************************************************************************/

  /*  Compute the Otsu threshold value for an image based on its histogram h: this value is
        only to the resolution of the bin size of the histogram, so a bin index b is
        returned, the implication being that everyting >= Bin2Value(b) is considered
        foreground.
  */

int Otsu_Threshold(Histogram *h)
{ int        i, t;
  Size_Type  c, pden, tden;
  double     psum, tsum;
  double     var, max;
  double     b, u;
  Size_Type *count;

  if (((Histofull *) h)->mean == DBL_MIN)
    set_histogram_mean(((Histofull *) h));

  count = h->counts;

  switch (h->kind) {
    #GENERATE t = uval ival fval
      case <T>:
        b = (double) h->binsize.<t>;
        #IF t == fval
          u = h->offset.<t> + .5*b;
        #ELSE
          u = h->offset.<t> + .5*b - .5;
        #END
        break;
    #END
      default:
        u = b = 0.;
        break;
  }

  tden = h->total;
  tsum = tden * ((Histofull *) h)->mean;

  pden = 0;
  psum = 0.;
  max  = 0.;
  t    = 0;
  for (i = 0; i < h->nbins-1; i++)
    { c = count[i];
      pden += c;
      psum += c * u;
      tden -= c;
      u    += b;
      var = psum/pden - (tsum-psum)/tden;
      var = (pden*var)*(tden*var);
      if (var > max)
        { max = var;
          t   = i;
        }
    }

  return (t+1);
}

  /*  Similarly, slope threshold returns the inflection point (to the nearest bin boundary)
        of the histogram relative to a line from the maximum bin to the larget non-zero bin.
  */

int Triangle_Threshold(Histogram *h)
{ int        i;
  int        low, hgh;
  Size_Type  c, max;
  Size_Type *count;
  double     slope;

  count = h->counts;

  low = hgh = 0;
  max = count[0];
  for (i = 1; i < h->nbins; i++)
    { c = count[i];
      if (c > max)
        { max = count[i];
          low = i;
        }
      if (c > 0)
        hgh = i;
    }
  hgh = hgh+1;

  slope = (1.*max) / (hgh-low);
  max   = 0;
  for (i = low+1; i < hgh; i++)
    { c = (Size_Type) ((hgh - i) * slope) - count[i];
      if (c > max)
        { max = c;
          low = i;
        }
    }

  return (low+1);
}

int Intermeans_Threshold(Histogram *h)
{ int        i, t, n;
  int        low, hgh;
  Size_Type  sum1, size1;
  Size_Type  sumt, sizet;
  Size_Type *count;
  double     mean1, mean2;

  count = h->counts;

  sumt  = 0;
  sizet = 0;
  low   = -1;
  hgh   = -1;
  for (i = 0; i < h->nbins; i++)
    { Size_Type c = count[i];
      if (c > 0)
        { if (low < 0)
            low = i;
          hgh = i;
          sizet += c;
          sumt  += c * i;
        }
    }

  t = (low+hgh)/2;
  n = -1;
  while (t != n)
    { n = t;

      sum1  = 0;
      size1 = 0;
      for (i = low; i < t; i++)
        { sum1  += count[i] * i;
          size1 += count[i];
        }
      mean1 = (double) ((1.*sum1) / size1);
      mean2 = (double) ((1. * ( sumt - sum1 )) / (sizet - size1));
      t = (int) (.5 * (mean1 + mean2));
    }
  
  return (t);
}


/****************************************************************************************
 *                                                                                      *
 *  HISTOGRAM DISPLAY                                                                   *
 *                                                                                      *
 ****************************************************************************************/

  /*  Print an ascii display of histogram h on FILE output indented by indent spaces.
      The parameter flag is the bitwise or of the defined constants below which determine
      what is displayed and how it is displayed.  If binsize is not 0 then the histogram
      will be displayed in bins of the given size, with the bin boundaries being multiples
      of binsize (the counts of spanning bins in the underlying histogram are interpolated).

        BIN_COUNT            0x01   //  show bin counts
        CUMULATIVE_COUNT     0x02   //  show cumulative counts
        CUMULATIVE_PERCENT   0x04   //  show cumulative per cent of the total
        ASCENDING_HGRAM      0x08   //  display in ascending order (descending by default)
          CLIP_HGRAM_HIGH    0x10   //  do not display any 0 count bins at the top
          CLIP_HGRAM_LOW     0x20   //  do not display any 0 count bins at the bottom
        CLIP_HGRAM           0x30   //  do not display any 0 count bins at either extreme
        BIN_MIDDLE           0x40   //  display the mid-value of a bin as opposed to its range
  */

void Print_Histogram(Histogram *h, FILE *output, int indent, int flags, Value binsize)
{ Size_Type sum, pre, *count;
  double    total;
  int       i, j, top, bot, inc;
  int       rwidth, dwidth, lwidth;
  int       bflag, cflag, pflag, mflag;
  int       exp;

  count = h->counts;

  if (flags == 0)
    flags = BIN_COUNT;

  if ((flags & CLIP_HGRAM_HIGH) != 0)
    { for (top = h->nbins-1; top >= 0; top--)
        if (count[top] != 0)
          break;
    }
  else
    top = h->nbins-1;

  if ((flags & CLIP_HGRAM_LOW) != 0)
    { for (bot = 0; bot < h->nbins; bot++)
        if (count[bot] != 0)
          break;
    }
  else
    bot = 0;

  if (top < bot)
    { fprintf(output,"%*sEmpty histogram!\n",indent,"");
      return;
    }

  bflag = ((flags & BIN_COUNT) != 0);
  cflag = ((flags & CUMULATIVE_COUNT) != 0);
  pflag = ((flags & CUMULATIVE_PERCENT) != 0);
  mflag = ((flags & BIN_MIDDLE) != 0);

  if ((flags & ASCENDING_HGRAM) == 0)
    inc = -1;
  else
    inc = 1;

  total  = (double) h->total;
  dwidth = ((int) floor(log10(total))) + 1;
  total  = 100./total;

  switch (h->kind) {
    #GENERATE t,u,f = uval ival fval , uint64 int64 float64 , llu lld g

      case <T>:
        { <u> b  = h->binsize.<t>;
          <u> o  = h->offset.<t>;
          <u> f  = h->offset.<t>;

          <u> u  = o + ((<u>) bot)*b;
          <u> v  = o + ((<u>) (top+1))*b;
          <u> eps = b*1e-8;             #WHEN t == fval

          <u> B  = b;
          if (binsize.<t> != 0)
            B = binsize.<t>;

          #IF t == fval
            bot = (int) floor((u+eps)/B);
            top = (int) ceil((v-eps)/B);
          #ELSEIF t == ival
            if (u < 0)
              bot = (int) ((u+1)/B-1);
            else
              bot = (int) (u/B);
            if (v < 0)
              top = (int) (v/B);
            else
              top = (int) ((v-1)/B+1);
          #ELSE
            bot = (int) (u/B);
            top = (int) ((v-1)/B+1);
          #END
          u = ((<u>) bot)*B;
          v = ((<u>) top)*B;

          if (v != 0)
            rwidth = ((int) floor(log10(fabs((double) v)))) + 1;
          else
            rwidth = 1;
          #IF t == ival
            if (u < 0)
              lwidth = ((int) floor(log10((double) -u))) + 2;
            else
              lwidth = ((int) floor(log10((double) u))) + 1;
            if (lwidth > rwidth)
              rwidth = lwidth;
          #ELSEIF t == fval
            if (u < 0)
              { lwidth = ((int) floor(log10((double) -u))) + 1;
                if (lwidth > rwidth)
                  rwidth = lwidth;
              }
            lwidth = (int) floor(log10(B));
            if (rwidth > 9 && lwidth > 4)
              { exp    = 1;
                rwidth = (rwidth-lwidth)+5;
                lwidth = rwidth-6;
                if (u < 0) rwidth += 1;
              }
            else if (lwidth < -9 && rwidth < -4)
              { exp = 1;
                rwidth = 5+(rwidth-lwidth);
                lwidth = rwidth-6;
                if (u < 0) rwidth += 1;
              }
            else
              { exp = 0;
                if (lwidth > 0)
                  lwidth = 0;
                else
                  { int ten, w;
                    lwidth = -lwidth;
                    ten = 1;
                    for (w = 0; w < lwidth; w++)
                      ten *= 10;
                    for (w = 0; w < 8; w++)
                      { if (ten*B - ((int) (ten*B)) < .01)
                          break;
                        ten *= 10;
                      }
                    lwidth += w;
                  }
                if (rwidth <= 0)
                  { rwidth = lwidth+2;
                    if (u < 0) rwidth += 1;
                  }
                else
                  { if (u < 0)
                      { if (ceil(log10((double) -u)) >= rwidth)
                          rwidth += 1;
                      }
                    if (lwidth > 0)
                      rwidth += lwidth + 1;
                  }
              }
          #END

          sum = 0;
          if (inc < 0)
            u = ((<u>) (top-1)) * B;
          else
            u = ((<u>) bot) * B; 
          v = u + B;
          j = 0;
          for (i = bot; i != top; i++)
            { if (mflag)
                #IF t == fval
                  if (exp)
                    fprintf(output,"%*s%*.*e:",indent,"",rwidth,lwidth,u+B/2.);
                  else if (B == 1)
                    fprintf(output,"%*s%*.*f:",indent,"",rwidth,lwidth,u);
                  else
                    fprintf(output,"%*s%*.*f:",indent,"",rwidth,lwidth,u+B/2.);
                #ELSE
                  if (B == 1)
                    fprintf(output,"%*s%*<f>:",indent,"",rwidth,u);
                  else
                    fprintf(output,"%*s%*<f>:",indent,"",rwidth,u+(B-1)/2);
                #END
              else
                #IF t == fval
                  if (exp)
                    fprintf(output,"%*s%*.*e - %*.*e:",indent,"",rwidth,lwidth,u,rwidth,lwidth,v);
                  else
                    fprintf(output,"%*s%*.*f - %*.*f:",indent,"",rwidth,lwidth,u,rwidth,lwidth,v);
                #ELSE
                  if (B == 1)
                    fprintf(output,"%*s%*<f>:",indent,"",rwidth,u);
                  else
                    fprintf(output,"%*s%*<f> - %*<f>:",indent,"",rwidth,u,rwidth,u+(B-1));
                #END

              while (f > u)                            #WHEN t != fval
              while (f > u+eps)                        #WHEN t == fval
                { j -= 1;
                  f -= b;
                }
              while (f + b <= u)                       #WHEN t != fval
              while (f + b <= u+eps)                   #WHEN t == fval
                { j += 1;
                  f += b;
                }

              pre = sum;

              if (f < u && j >= 0 && j < h->nbins)      #WHEN t != fval
              if (f < u-eps && j >= 0 && j < h->nbins)  #WHEN t == fval
                sum -= (Size_Type) (count[j] * ((u - f)/(1.*b)));
               
              while (f + b <= v)                        #WHEN t != fval
              while (f + b <= v+eps)                    #WHEN t == fval
                { if (j >= 0 && j < h->nbins)
                    sum += count[j];
                  j += 1;
                  f += b;
                }
            
              if (f < v && j >= 0 && j < h->nbins)      #WHEN t != fval
              if (f < v-eps && j >= 0 && j < h->nbins)  #WHEN t == fval
                sum += (Size_Type) (count[j] * ((v - f)/(1.*b)));

              if (bflag)
                fprintf(output,"  %*llu",dwidth,sum-pre);
              if (cflag)
                fprintf(output,"  %*llu",dwidth,sum);
              if (pflag)
                fprintf(output,"  %6.1f%%",sum*total);
              fprintf(output,"\n");
              if (inc > 0)
                { u  = v;
                  v += B;
                }
              else
                { v  = u;
                  u -= B;
                }
            }
          break;
        }
    #END
  }
}
