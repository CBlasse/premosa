/*****************************************************************************************\
*                                                                                         *
*  Region Data Abstractions                                                               *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  June 2007                                                                     *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "connectivity.h"
#include "array.h"
#include "level.set.h"
#include "water.shed.h"
#include "region.h"

#undef  DEBUG_CONTOUR
#undef  DEBUG_RECORD

#LISTDEF @TYPES  = UINT8 UINT16 UINT32 UINT64 INT8 INT16 INT32 INT64 FLOAT32 FLOAT64
#LISTDEF @UNION  =  uval   uval   uval   uval ival  ival  ival  ival    fval    fval
#LISTDEF @VALUES =     U     U       U      U    I     I     I     I       F       F

static int type_size[] = { 1, 2, 4, 8, 1, 2, 4, 8, 4, 8 };

/****************************************************************************************
 *                                                                                      *
 *  2D CONTOUR ABSTRACTION: TRACE_CONTOUR, SPACE MANAGEMENT                             *
 *                                                                                      *
 ****************************************************************************************/

//  Awk-generated (manager.awk) Contour space management

typedef struct
  { int        length;          // Length of the contour
    boolean    boundary;        // Does contour touch the image boundary?
    int        dims[2];         // Width & height of image from which contour was derived
    int        iscon4;          // Is this of a 4- or 8-connected region?
    Indx_Type  seed;            // Start pixel for contour (is leftmost)
    int        direct;          // Initial direction of tour;
    uint8     *tour;            // Edge direction defining outer contour of the region
    int        nsegs;           // # of segments (if segmented)
    Indx_Type *vert;            // start point of each segment
    int       *strt;            // start index of segment in 'tour' array
  } Ctour;

static inline int ctour_tsize(Ctour *contour)
{ return (contour->length); }

static inline int ctour_vsize(Ctour *contour)
{ return ((contour->nsegs+1)*sizeof(Indx_Type)); }

static inline int ctour_ssize(Ctour *contour)
{ return ((contour->nsegs+1)*sizeof(int)); }

MANAGER -IO Contour(Ctour)  tour:tsize vert:vsize strt:ssize

#define LEGAL_MOVE(good,p,d)		\
{ switch (d)				\
    { case 0:				\
        good = (p < top);		\
        break;				\
      case 1:				\
        good = ((p % width) < rgt);	\
        break;				\
      case 2:				\
        good = (p >= bot);		\
        break;				\
      case 3:				\
        good = ((p % width) > lft);	\
        break;				\
    }					\
}

Contour *Trace_Contour(APart *image, boolean iscon4, Indx_Type seed,
                       void *arg, boolean test(Indx_Type p, void *arg))
{
#ifdef DEBUG_CONTOUR
  static char     *direction[] = { "S", "E", "N", "W" };
#endif

  Ctour  *my_cont;
  Array  *array = AForm_Array(image);

  Indx_Type  width;
  Indx_Type  rgt, lft, bot, top;
  Indx_Type  offset[4];
  Indx_Type  p, q, r;
  int        d, e;
  boolean    bnd, good;
  int        len;
  uint8     *tour;

  if (array->ndims != 2)
    { fprintf(stderr,"Image is not a 2D array (Trace_Contour)\n");
      exit (1);
    }

  width = array->dims[0];
  if (Is_Slice(image))
    { Dimn_Type *bcrd = ADIMN(Slice_First(image));
      Dimn_Type *ecrd = ADIMN(Slice_Last(image));

      lft  = bcrd[0];
      rgt  = ecrd[0];
      bot  = width * (bcrd[1]+1);
      top  = width * ecrd[1];
    }
  else
    { lft  = 0;
      rgt  = width-1;
      bot  = width;
      top  = width * (array->dims[1]-1);
    }

  offset[0] =  width;
  offset[1] =  1;
  offset[2] = -width;
  offset[3] = -1;

  len  = 1;
  bnd  = 0;
  p    = seed;
  d    = 0;
  good = 0;
#ifdef DEBUG_CONTOUR
  printf("\nContour:\n  (%3lld,%3lld) -> %s\n",p%width,p/width,direction[d]);
#endif

  if (iscon4)   // 4-connected contour
    do
      { LEGAL_MOVE(good,p,d);
        q = p + offset[d];
        if (good && test(q,arg))
          { e = (d+3) % 4;
            LEGAL_MOVE(good,q,e);
            r = q + offset[e];
            if (good && test(r,arg))
              { p = r;
                d = e;
#ifdef DEBUG_CONTOUR
                printf("  (%3lld,%3lld) -> (%3lld,%3lld) %s\n",q%width,q/width,
                                                       r%width,r/width,direction[d]);
                fflush(stdout);
#endif
              }
            else
              { p = q;
#ifdef DEBUG_CONTOUR
                printf("  (%3lld,%3lld) %s\n",q%width,q/width,direction[d]);
                fflush(stdout);
#endif
              }
          }
        else
          { d = (d+1) % 4;
            bnd = 1;
#ifdef DEBUG_CONTOUR
            printf("  . %s\n",direction[d]);
            fflush(stdout);
#endif
          }
        len += 1;
      }
    while (p != seed || d != 0);

  else                // 8-connected contour
    do
      { LEGAL_MOVE(good,p,d);
        q = p + offset[d];
        if (good)
          { e = (d+3) % 4;
            LEGAL_MOVE(good,q,e)
            r = q + offset[e];
            if (good && test(r,arg))
              { p = r;
                d = e;
#ifdef DEBUG_CONTOUR
                printf("  (%3lld,%3lld) %s\n",r%width,r/width,direction[d]);
                fflush(stdout);
#endif
              }
            else if (test(q,arg))
              { p = q;
#ifdef DEBUG_CONTOUR
                printf("  (%3lld,%3lld) %s\n",q%width,q/width,direction[d]);
                fflush(stdout);
#endif
              }
            else
              { d = (d+1) % 4;
#ifdef DEBUG_CONTOUR
                printf("  . %s\n",direction[d]);
                fflush(stdout);
#endif
              }
          }
        else
          { d = (d+1) % 4;
#ifdef DEBUG_CONTOUR
            printf("  . %s\n",direction[d]);
            fflush(stdout);
#endif
          }
        len += 1;
      }
    while (p != seed || d != 0);

  my_cont = new_ctour(len,0,0,"Trace_Contour");

  my_cont->length   = len;
  my_cont->boundary = bnd;
  my_cont->dims[0]  = width;
  my_cont->dims[1]  = array->dims[1];
  my_cont->iscon4   = iscon4;
  my_cont->seed     = seed;
  my_cont->direct   = 0;
  my_cont->nsegs    = 0;

  tour = my_cont->tour;

  len  = 0;
  p = seed;
  d = 0;

  tour[len++] = d;
  if (iscon4)   // 4-connected contour
    do
      { LEGAL_MOVE(good,p,d);
        q = p + offset[d];
        if (good && test(q,arg))
          { e = (d+3) % 4;
            LEGAL_MOVE(good,q,e);
            r = q + offset[e];
            if (good && test(r,arg))
              { p = r;
                d = e;
              }
            else
              p = q;
          }
        else
          d = (d+1) % 4;
        tour[len++] = d;
      }
    while (p != seed || d != 0);

  else                // 8-connected contour
    do
      { LEGAL_MOVE(good,p,d);
        q = p + offset[d];
        if (good)
          { e = (d+3) % 4;
            LEGAL_MOVE(good,q,e);
            r = q + offset[e];
            if (good && test(r,arg))
              { p = r;
                d = e;
              }
            else if (test(q,arg))
              p = q;
            else
              d = (d+1) % 4;
          }
        else
          d = (d+1) % 4;
        tour[len++] = d;
      }
    while (p != seed || d != 0);

  return ((Contour *) my_cont);
}

void For_Contour(Contour *ctour, void *arg, boolean handler(Indx_Type p, Indx_Type q, void *arg))
{ Ctour  *my_cont = (Ctour *) ctour;

  Indx_Type  width;
  Indx_Type  rgt, lft, bot, top;
  Indx_Type  offset[4];
  Indx_Type  p, q;
  int        d, e, x;
  boolean    good = 0;
  int        i, len;
  uint8     *tour;

  width = my_cont->dims[0];
  len   = my_cont->length - 1;
  tour  = my_cont->tour;

  lft  = 0;
  rgt  = width-1;
  bot  = width;
  top  = width * (my_cont->dims[1]-1);

  offset[0] =  width;
  offset[1] =  1;
  offset[2] = -width;
  offset[3] = -1;

  p = my_cont->seed;
  LEGAL_MOVE(good,p,3);
  if (good)
    q = p-1;
  else
    q = -1;
  d = 0;
  if (handler(p,q,arg))
    return;
  for (i = 1; i < len; i++)
    { e = tour[i];
      x = (4+(e-d)) % 4;
      if (x == 0)
        p += offset[e];
      else if (x == 3)
        p += offset[d] + offset[e];
      d = (e+3)%4;
      LEGAL_MOVE(good,p,d)
      if (good)
        q = p+offset[d];
      else
        q = -1;
      if (handler(p,q,arg))
        break;
      d = e;
    }
}

int Segment_Contour(Contour *ctour, void *arg, boolean handler(Indx_Type p, Indx_Type q, void *))
{ Ctour  *my_cont = (Ctour *) ctour;

  Indx_Type  width;
  Indx_Type  rgt, lft, bot, top;
  Indx_Type  offset[4];
  Indx_Type  p, q;
  int        d, e, x;
  boolean    good = 0;
  int        i, len, nsegs;
  uint8     *tour;
  Indx_Type *vert;
  int       *strt;

  width = my_cont->dims[0];
  len   = my_cont->length - 1;
  tour  = my_cont->tour;

  lft  = 0;
  rgt  = width-1;
  bot  = width;
  top  = width * (my_cont->dims[1]-1);

  offset[0] =  width;
  offset[1] =  1;
  offset[2] = -width;
  offset[3] = -1;

  p = my_cont->seed;
  LEGAL_MOVE(good,p,3);
  if (good)
    q = p-1;
  else
    q = -1;
  d = 0;
  handler(p,q,arg);
  nsegs = 0;
  for (i = 1; i <= len; i++)
    { e = tour[i];
      x = (4+(e-d)) % 4;
      if (x == 0)
        p += offset[e];
      else if (x == 3)
        p += offset[d] + offset[e];
      d = (e+3)%4;
      LEGAL_MOVE(good,p,d)
      if (good)
        q = p+offset[d];
      else
        q = -1;
      if (handler(p,q,arg))
        nsegs += 1;
      d = e;
    }

//printf("nsegs = %d\n",nsegs);

  allocate_ctour_vert(my_cont,sizeof(Indx_Type)*(nsegs+1),"Segment_Contour");
  allocate_ctour_strt(my_cont,sizeof(int)*(nsegs+1),"Segment_Contour");
  my_cont->nsegs = nsegs;

  vert = my_cont->vert;
  strt = my_cont->strt;

  p = my_cont->seed;
  LEGAL_MOVE(good,p,3);
  if (good)
    q = p-1;
  else
    q = -1;
  d = 0;
  handler(p,q,arg);
  strt[0] = 0;
  vert[0] = p;
  nsegs = 0;
  for (i = 1; i <= len; i++)
    { e = tour[i];
      x = (4+(e-d)) % 4;
      if (x == 0)
        p += offset[e];
      else if (x == 3)
        p += offset[d] + offset[e];
      d = (e+3)%4;
      LEGAL_MOVE(good,p,d)
      if (good)
        q = p+offset[d];
      else
        q = -1;
      if (handler(p,q,arg))
        { nsegs += 1;
          strt[nsegs] = i;
          vert[nsegs] = p;
        }
      d = e;
    }
  if (strt[nsegs] < len)
    { strt[0] = strt[nsegs];
      vert[0] = vert[nsegs];
      //strt[nsegs] = strt[1];
      //vert[nsegs] = vert[1];
    }
/*
{ Array *image = (Array *) arg;
  for (i = 0; i < nsegs; i++)
    { printf(" %2d: %3d (",i,strt[i]);
      Print_Coord(stdout,Idx2CoordA(image,vert[i]));
      printf(")\n");
    }
}*/

  return (nsegs);
}

void For_Contour_Segment(Contour *ctour, int seg, void *arg,
                         boolean handler(Indx_Type p, Indx_Type q, void *))
{ Ctour  *my_cont = (Ctour *) ctour;

  Indx_Type  width;
  Indx_Type  rgt, lft, bot, top;
  Indx_Type  offset[4];
  Indx_Type  p, q;
  int        d, e, x;
  boolean    good = 0;
  int        i, len;
  uint8     *tour;

  if (my_cont->nsegs == 0)
    { fprintf(stderr,"Contour has not been segmented (For_Contour_Segment)\n");
      exit (1);
    }
  if (seg >= my_cont->nsegs)
    { fprintf(stderr,"Contour has less than %d segments (For_Contour_Segment)\n",seg+1);
      exit (1);
    }

  width = my_cont->dims[0];
  tour  = my_cont->tour;

  lft  = 0;
  rgt  = width-1;
  bot  = width;
  top  = width * (my_cont->dims[1]-1);
  if (seg == 0 && my_cont->strt[0] != 0)
    len = my_cont->length;
  else
    len  = my_cont->strt[seg+1];

  offset[0] =  width;
  offset[1] =  1;
  offset[2] = -width;
  offset[3] = -1;

  p = my_cont->vert[seg];
  i = my_cont->strt[seg];
  e = tour[i];
  d = (e+3)%4;
  LEGAL_MOVE(good,p,d);
  if (good)
    q = p+offset[d];
  else
    q = -1;
  handler(p,q,arg);
  d = e;
repeat:
  for (i = i+1; i < len; i++)
    { e = tour[i];
      x = (4+(e-d)) % 4;
      if (x == 0)
        p += offset[e];
      else if (x == 3)
        p += offset[d] + offset[e];
      d = (e+3)%4;
      LEGAL_MOVE(good,p,d)
      if (good)
        q = p+offset[d];
      else
        q = -1;
      if (handler(p,q,arg))
        break;
      d = e;
    }
  if (seg == 0 && my_cont->strt[0] != 0)
    { len = my_cont->strt[seg+1];
      i   = 0;
      seg = -1;
      goto repeat;
    }
}

void Rev_Contour_Segment(Contour *ctour, int seg, void *arg,
                         boolean handler(Indx_Type p, Indx_Type q, void *))
{ Ctour  *my_cont = (Ctour *) ctour;

  Indx_Type  width;
  Indx_Type  rgt, lft, bot, top;
  Indx_Type  offset[4];
  Indx_Type  p, q;
  int        d, e, x;
  boolean    good = 0;
  int        i, len;
  uint8     *tour;

  if (my_cont->nsegs == 0)
    { fprintf(stderr,"Contour has not been segmented (For_Contour_Segment)\n");
      exit (1);
    }
  if (seg >= my_cont->nsegs)
    { fprintf(stderr,"Contour has less than %d segments (For_Contour_Segment)\n",seg+1);
      exit (1);
    }

  width = my_cont->dims[0];
  tour  = my_cont->tour;

  lft  = 0;
  rgt  = width-1;
  bot  = width;
  top  = width * (my_cont->dims[1]-1);

  len  = my_cont->strt[seg];
  if (seg == 0 && my_cont->strt[0] != 0)
    len = 1;
  else
    len = my_cont->strt[seg];

  offset[0] =  width;
  offset[1] =  1;
  offset[2] = -width;
  offset[3] = -1;

  p = my_cont->vert[seg+1];
  i = my_cont->strt[seg+1];
  d = tour[i];
repeat:
  for (i--; i >= len; i--)
    { e = tour[i];
      x = (4+(e-d)) % 4;
      if (x == 0)
        p -= offset[e];
      else if (x == 1)
        p -= offset[e] + offset[d];
      d = (e+3)%4;
      LEGAL_MOVE(good,p,d)
      if (good)
        q = p+offset[d];
      else
        q = -1;
      if (handler(p,q,arg))
        break;
      d = e;
    }
  if (seg == 0 && my_cont->strt[0] != 0)
    { len = my_cont->strt[seg];
      i   = my_cont->length;
      seg = -1;
      goto repeat;
    }
}

Dimn_Type *Get_Contour_Dimensions(Contour *tour)
{ Ctour *cont = (Ctour *) tour;
  return (cont->dims);
}

boolean Get_Contour_Connectivity(Contour *tour)
{ Ctour *cont = (Ctour *) tour;
  return (cont->iscon4);
}

int Get_Contour_Segments(Contour *tour)
{ Ctour *cont = (Ctour *) tour;
  return (cont->nsegs);
}

boolean Boundary_Countour(Contour *tour)
{ Ctour *cont = (Ctour *) tour;
  return (cont->boundary);
}


/****************************************************************************************
 *                                                                                      *
 *  REGION ABSTRACTION: RECORD_REGION, SPACE MANAGEMENT                                 *
 *                                                                                      *
 ****************************************************************************************/

//  Awk-generated (manager.awk) Region space management

typedef struct
  { Size_Type  rastlen;      // Length of the raster pair vector (always even)
    Indx_Type *raster;       // [0..surflen-1] of all surface pixels, with the sublist
    Size_Type  surflen;      // Length of all surface pixels (surflen >= rastlen)
    boolean    iscon2n;      // Is this of a 2n- or (3^n-1)-connected region?
    int        ndims;        // Dimensionality
    Dimn_Type *dims;         // Dims
                             // [0..rastlen-1] of 0-dim extreme pixels in sorted order (paired)
    uint8     *ishole;       // Is the space between a pair part of a hole?
    Size_Type  area;         // Surface area of the object
  } RasterCon;

#ifdef DEBUG_REGION

static void Print_Region(RasterCon *region, FILE *output)
{ int i;

  fprintf(output,"\nRegion [0..%llu..%llu] is2n = %d ndims = %d area = %llu\n",
                 region->rastlen,region->surflen,region->iscon2n,region->ndims,region->area);
  for (i = 0; i < region->rastlen; i += 2)
    fprintf(output,"  r %llu %llu (%d)\n",region->raster[i],region->raster[i+1],
                                          region->ishole[i>>1]);
  for (i = region->rastlen; i < region->surflen; i += 1)
    fprintf(output,"  s %llu\n",region->raster[i]);
}

#endif

#define SIZEOF(x) ((int) sizeof(x))

static inline int64 rastercon_rsize(RasterCon *region)
{ return (SIZEOF(Indx_Type) * region->surflen); }

static inline int64 rastercon_hsize(RasterCon *region)
{ return (SIZEOF(uint8) * (region->rastlen >> 1)); }

static inline int rastercon_dsize(RasterCon *region)
{ return (SIZEOF(int) * region->ndims); }

MANAGER -IO Region(RasterCon) raster!rsize ishole!hsize dims:dsize

typedef struct
  { Size_Type   size;
    Size_Type   area;
    Size_Type   bot;
    Size_Type   top;
    RasterCon  *reg;
    Indx_Type  *vec;
  } RegArg;

#define RA(a) ((RegArg *) (a))

boolean allocate_surface(Size_Type size, Size_Type area, Size_Type dbls, void *a)
{ allocate_rastercon_raster(RA(a)->reg,SIZEOF(Indx_Type)*(size+dbls),"Record_Region");
  RA(a)->vec  = RA(a)->reg->raster;
  RA(a)->top  = RA(a)->size = size+dbls;
  RA(a)->bot  = 0;
  RA(a)->area = area;
  return (1);
}

void add_surface(Indx_Type p, int x, void *a)
{ if (x)
    { RA(a)->vec[RA(a)->bot++] = p;
      if (x > 1)
        RA(a)->vec[RA(a)->bot++] = p;
    }
  else
    RA(a)->vec[--RA(a)->top] = p;
}

static int PSORT(const void *x, const void *y)
{ Indx_Type l = *((Indx_Type *) x);
  Indx_Type r = *((Indx_Type *) y);
  if (l < r)
    return (-1);
  else if (l > r)
    return (1);
  else
    return (0);
}

Region *G(Record_Region)(APart *source, int share, boolean iscon2n, Indx_Type leftmost,
                         boolean with_holes, void *argt, boolean (*test)(Indx_Type p, void *argt))
{ Array     *array = AForm_Array(source);
  RasterCon *region;
  Size_Type  rastop, more;
  Indx_Type *raster;
  RegArg     arg, *argp = &arg;

  region  = new_rastercon(0,0,array->ndims*SIZEOF(int),"Record_Region");
  arg.reg = region;
  Flood_Surface(source,share,iscon2n,leftmost,
                argt,test,NULL,NULL,argp,allocate_surface,argp,add_surface);

  rastop = arg.bot;
  raster = arg.vec;

  qsort(raster,(size_t) rastop,sizeof(Indx_Type),PSORT);

#ifdef DEBUG_RECORD
  { Indx_Type i;

    fprintf(stdout,"\nRegion [0..%llu..%llu] is2n = %d ndims = %d area = %llu\n",
                   rastop,arg.size,iscon2n,array->ndims,arg.area);
    for (i = 0; i < rastop; i += 2)
      fprintf(stdout,"  r %llu %llu\n",raster[i],raster[i+1]);
    for (i = rastop; i < arg.size; i += 1)
      fprintf(stdout,"  s %llu\n",raster[i]);
  }
#endif

  more = 0;
  if (with_holes)
    { Indx_Type i, v, w, p;

      for (i = 0; i < rastop; i += 2)
        { v = raster[i];
          w = raster[i+1];
          for (p = v+1; p < w; p++)
            if ( ! test(p,argt))
              { more += 2; 
                while ( ! test(p,argt))
                  p += 1;
              }
        }
    }

  { Indx_Type  i;

    region->iscon2n = iscon2n;
    region->ndims   = array->ndims;
    region->area    = arg.area;
    for (i = 0; i < region->ndims; i++)
      region->dims[i] = array->dims[i];
    region->rastlen = rastop + more;
    region->surflen = arg.size + more;

    allocate_rastercon_ishole(region,((region->rastlen)>>1)*SIZEOF(uint8),"Record_Region");

    if (with_holes)
      { Indx_Type v, w, p, j;
        uint8    *ishole = region->ishole;

        allocate_rastercon_raster(region,region->surflen*SIZEOF(Indx_Type),"Record_Region");

        raster = region->raster;

        memmove(raster+region->rastlen,raster+rastop,sizeof(Indx_Type)*((size_t)(arg.size-rastop)));

        j = region->rastlen-1;
        for (i = rastop; i-- > 0; i -= 1)
          { w = raster[i];
            v = raster[i-1];
            raster[j--]  = w;
            for (p = w-1; p > v; p--)
              if ( ! test(p,argt))
                { ishole[j>>1] = 1;
                  raster[j--]  = p+1;
                  p -= 1;
                  while ( ! test(p,argt))
                    p -= 1;
                  raster[j--] = p;
                }
            ishole[j>>1] = 0;
            raster[j--] = v;
          }
      }

    else
      memset(region->ishole, 0, ((size_t) ((region->rastlen)>>1))*sizeof(uint8));
  }

#ifdef DEBUG_RECORD
  Print_Region(region,stdout);
#endif

  return ((Region *) region);
}

typedef struct
  { void       *value;
    uint64      level_uval;
    int64       level_ival;
    float64     level_fval;
  } TestArg;

#define TA(a) ((TestArg *) (a))

#GENERATE T,U = @TYPES , @UNION
  #GENERATE C,O = LE LT EQ NE GT GE , <= < == != > >=
    static boolean is_<c>_<t>(Indx_Type p, void *a)
    { return (((<t> *) TA(a)->value)[p] <o> TA(a)->level_<u>); }
  #END
#END

static boolean (*Comparator_Table[])(Indx_Type, void *) = {
#GENERATE T = @TYPES
  #GENERATE C = LE LT EQ NE GT GE
    is_<c>_<t>,
  #END
#END
  };

Region *G(Record_Basic)(APart *source, int share, boolean iscon2n, Indx_Type leftmost,
                              boolean with_holes, Comparator cmprsn, Value level)
{ Array     *array = AForm_Array(source);
  boolean  (*test)(Indx_Type, void *);
  TestArg    argt;

  switch (array->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        argt.value = array->data;
        argt.level_<u> = level.<u>;
        break;
    #END
  }

  test = Comparator_Table[6*array->type + cmprsn];

  return (Record_Region(source,share,iscon2n,leftmost,with_holes,&argt,test));
}

Contour *G(Basic_Contour)(APart *image, boolean iscon4, Indx_Type seed,
                          Comparator cmprsn, Value level)
{ Array     *array = AForm_Array(image);
  boolean  (*test)(Indx_Type, void *);
  TestArg    argt;

  switch (array->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        argt.value = array->data;
        argt.level_<u> = level.<u>;
        break;
    #END
  }

  test = Comparator_Table[6*array->type + cmprsn];

  return (Trace_Contour(image,iscon4,seed,&argt,test));
}

int Get_Region_Dimensionality(Region *region)
{ RasterCon *cont = (RasterCon *) region;
  return (cont->ndims);
}

Dimn_Type *Get_Region_Dimensions(Region *region)
{ RasterCon *cont = (RasterCon *) region;
  return (cont->dims);
}

boolean Get_Region_Connectivity(Region *region)
{ RasterCon *cont = (RasterCon *) region;
  return (cont->iscon2n);
}

void For_Region(Region *reg, void *arg, void (*handler)(Indx_Type p, void *arg))
{ RasterCon *trace = (RasterCon *) reg; 
  Indx_Type  len;
  Indx_Type *raster;

  raster = trace->raster;
  len    = trace->rastlen;

  { Indx_Type i, v, w, p;

    for (i = 0; i < len; i += 2)
      { v = raster[i];
        w = raster[i+1];
        for (p = v; p <= w; p++)
          handler(p,arg);
      }
  }
}

void For_Region_Outline(Region *reg, void *arg, void (*handler)(Indx_Type p, void *arg))
{ RasterCon *trace = (RasterCon *) reg; 
  Indx_Type  i, len;
  Indx_Type *raster;
  uint8     *ishole;

  raster = trace->raster;
  ishole = trace->ishole;
  len    = trace->rastlen;

  handler(raster[0],arg);
  for (i = 2; i < len; i += 2)
    if (!ishole[i>>1])
      { handler(raster[i-1],arg);
        handler(raster[i],arg);
      }

  len = trace->surflen;
  for (i = trace->rastlen-1; i < len; i++)
    handler(raster[i],arg);
}

void For_Region_Holes(Region *reg, void *arg, void (*handler)(Indx_Type p, void *arg))
{ RasterCon *trace = (RasterCon *) reg; 
  Indx_Type  len;
  Indx_Type *raster;
  uint8     *ishole;

  raster = trace->raster;
  ishole = trace->ishole;
  len    = trace->rastlen;

  { Indx_Type i, v, w, p;

    for (i = 2; i < len; i += 2)
      if (ishole[i>>1])
        { v = raster[i-1];
          w = raster[i];
          for (p = v+1; p < w; p++)
            handler(p,arg);
        }
  }
}

void For_Region_Exterior(Region *reg, void *arg, void (*handler)(Indx_Type p, void *arg))
{ RasterCon *trace = (RasterCon *) reg; 
  Indx_Type  len;
  Indx_Type *raster;
  uint8     *ishole;
  Size_Type  size;

  raster = trace->raster;
  ishole = trace->ishole;
  len    = trace->rastlen;

  { Indx_Type i, v, w, p;

    size = 1;
    for (i = 0; i < trace->ndims; i++)
      size *= trace->dims[i];

    for (p = 0; p < raster[0]; p++)
      handler(p,arg);
    for (i = 2; i < len; i += 2)
      if (!ishole[i>>1])
        { v = raster[i-1];
          w = raster[i];
          for (p = v+1; p < w; p++)
            handler(p,arg);
        }
    for (p = raster[len]+1; p < size; p++)
      handler(p,arg);
  }
}

Region *G(Record_Level_Set)(Level_Tree *t, Level_Set *r, boolean with_holes, int share)
{ APart  *image   = Get_Level_Tree_APart(t);
  boolean iscon2n = Get_Level_Tree_Connectivity(t);

  return (Record_Basic(image,share,iscon2n,Level_Set_Leftmost(t,r),with_holes,
                       GE_COMP,VALU(Level_Set_Level(t,r))));
}

Contour *G(Level_Set_Contour)(Level_Tree *t, Level_Set *r)
{ APart  *image   = Get_Level_Tree_APart(t);
  boolean iscon2n = Get_Level_Tree_Connectivity(t);

  return (Basic_Contour(image,iscon2n,Level_Set_Leftmost(t,r),
                        GE_COMP,VALU(Level_Set_Level(t,r))));
}

Region *G(Record_P_Vertex)(Partition *w, int cb, boolean with_holes, int share)

{ Array    *image   = Get_Partition_Labels(w);
  boolean   iscon2n = Is_Partition_2n_Connected(w);
  Indx_Type pixel   = Get_Partition_Vertex(w,cb)->seed;
  Value     val;

  if (image == NULL)
    { fprintf(stderr,"Partition does not have a label array (Record_P_Vertex)\n");
      exit (1);
    }

  switch (image->type) {
    #GENERATE T,V = @TYPES , @VALUES
      case <T>_TYPE:
        val = VAL<V>(A<T>(image)[pixel]);
        break;
    #END
  }

  return (Record_Basic(image,share,iscon2n,pixel,with_holes,EQ_COMP,val));
}

Contour *G(P_Vertex_Contour)(Partition *w, int cb)

{ Array    *image   = Get_Partition_Labels(w);
  boolean   iscon2n = Is_Partition_2n_Connected(w);
  Indx_Type pixel   = Get_Partition_Vertex(w,cb)->seed;
  Value     val;

  if (image == NULL)
    { fprintf(stderr,"Partition does not have a label array (Record_P_Vertex)\n");
      exit (1);
    }

  switch (image->type) {
    #GENERATE T,V = @TYPES , @VALUES
      case <T>_TYPE:
        val = VAL<V>(A<T>(image)[pixel]);
        break;
    #END
  }

  return (Basic_Contour(image,iscon2n,pixel,EQ_COMP,val));
}

Region *Fill_Region_Holes(Region *R(M(cont)))
{ RasterCon *region = (RasterCon *) cont;

  Indx_Type *raster;
  Size_Type  rlen, slen;
  uint8     *ishole;

  rlen   = region->rastlen;
  slen   = region->surflen;
  raster = region->raster;
  ishole = region->ishole;

  { Indx_Type i, j;

    j = 1;
    for (i = 2; i < rlen; i += 2)
      if ( ! ishole[i>>1]) 
        { raster[j++] = raster[i-1];
          raster[j++] = raster[i];
        }
    raster[j++] = raster[i-1];

    region->rastlen = j;

    for (i = rlen; i < slen; i++)
      raster[j++] = raster[i];

    region->surflen = j;

    rlen = (region->rastlen >> 1);
    for (i = 0; i < rlen; i++)
      ishole[i] = 0;
  }

#ifdef DEBUG_REGION
  Print_Region(region,stdout);
#endif

  return (cont);
}

//  Compute min and max values in 'a' of type 'type' with 'length' elements

Range_Bundle *Region_Range(Range_Bundle *R(O(rng)), APart *o, Region *reg)
{ RasterCon *cont = (RasterCon *) reg;
  Array     *a    = AForm_Array(o);
  Indx_Type *raster = cont->raster;
  Size_Type  len    = cont->rastlen;
  Indx_Type  k, v, w, p;

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          <t>  x, min, max;

          min = max = d[raster[0]];
          for (k = 0; k < len; k += 2)
            { v = raster[k];
              w = raster[k+1];
              for (p = v; p <= w; p++)
                { x = d[p];
                  if (x < min)
                    min = x;
                  else if (x > max)
                    max = x;
                }
            }
          rng->maxval.<U> = max;
          rng->minval.<U> = min;
          break;
        }
    #END
  }

  return (rng);
}


/****************************************************************************************
 *                                                                                      *
 *  CONTOUR ANALYSIS: EXTENT, CENTER OF MASS                                            *
 *                                                                                      *
 ****************************************************************************************/

Extent_Bundle *Region_Extent(Extent_Bundle *R(O(bundle)), Region *region)
{ RasterCon *cont = (RasterCon *) region;

  int        ndims;
  Dimn_Type *dims;
  Size_Type  rlen;
  Indx_Type *raster;

  ndims = cont->ndims;
  dims  = cont->dims;

  raster = cont->raster;
  rlen   = cont->rastlen;

  { Array *e;

    e = bundle->min;
    if (e == NULL)
      bundle->min = Make_Array_With_Shape(PLAIN_KIND,DIMN_TYPE,Coord1(ndims));
    else if (e->size*type_size[e->type] < ndims*type_size[DIMN_TYPE])
      { Free_Array(e);
        bundle->min = Make_Array_With_Shape(PLAIN_KIND,DIMN_TYPE,Coord1(ndims));
      }
    else
      { e->kind  = PLAIN_KIND;
        e->type  = DIMN_TYPE;
        e->scale = DIMN_SCALE;
        e->ndims = 1;
        e->size  = ndims;
        e->dims[0] = ndims;
      }

    e = bundle->max;
    if (e == NULL)
      bundle->max = Make_Array_With_Shape(PLAIN_KIND,DIMN_TYPE,Coord1(ndims));
    else if (e->size*type_size[e->type] < ndims*type_size[DIMN_TYPE])
      { Free_Array(e);
        bundle->max = Make_Array_With_Shape(PLAIN_KIND,DIMN_TYPE,Coord1(ndims));
      }
    else
      { e->kind  = PLAIN_KIND;
        e->type  = DIMN_TYPE;
        e->scale = DIMN_SCALE;
        e->ndims = 1;
        e->size  = ndims;
        e->dims[0] = ndims;
      }
  }

  { Indx_Type  p, i;
    Dimn_Type  o, x;
    int        d;
    Dimn_Type *min, *max;
    Dimn_Type  min0, max0, dim0;

    min = ADIMN(bundle->min);
    max = ADIMN(bundle->max);

    for (d = 0; d < ndims; d++)
      { max[d] = 0;
        min[d] = dims[d];
      }

    min0 = min[0];
    max0 = max[0];
    dim0 = dims[0];
    for (i = 0; i < rlen; i += 2)
      { o = dim0;
        p = raster[i+1];
        x = (Dimn_Type) (p%o);
        if (x > max0)
          max0 = x;
        p = raster[i];
        x = (Dimn_Type) (p%o);
        if (x < min0)
          min0 = x;
        for (d = 1; d < ndims; d++)
          { p /= o;
            o  = dims[d];
            x  = (Dimn_Type) (p%o);
            if (x > max[d])
              max[d] = x;
            if (x < min[d])
              min[d] = x;
          }
      }
    min[0] = min0;
    max[0] = max0;
  }

  return (bundle);
}

boolean Touches_Boundary(Extent_Bundle *bundle, APart *part)
{ Array     *array = AForm_Array(part);
  int        i, ndims;
  Dimn_Type *low, *hgh;
  Dimn_Type *min, *max;

  min     = ADIMN(bundle->min);
  max     = ADIMN(bundle->max);
  ndims   = array->ndims;

  if (Is_Slice(part))
    { low = ADIMN(Slice_First(part));
      hgh = ADIMN(Slice_Last(part));
      for (i = 0; i < ndims; i++)
        if (min[i] <= low[i] || max[i] >= hgh[i])
          return (1);
    }
  else
    { hgh = array->dims;
      for (i = 0; i < ndims; i++)
        if (min[i] <= 0 || max[i] >= hgh[i]-1)
          return (1);
    }

  return (0);
}

/* Return the area covered by the outer region */

Size_Type Region_Volume(Region *reg)
{ Indx_Type *raster, i;
  Size_Type  rlen, area;

  RasterCon *cont = (RasterCon *) reg;

  raster = cont->raster;
  rlen   = cont->rastlen;

  area = (rlen >> 1);
  for (i = 0; i < rlen; i += 2)
    area += raster[i+1] - raster[i];
  return (area);
}

Size_Type Region_Area(Region *reg)
{ RasterCon *cont = (RasterCon *) reg;
  return (cont->area);
}

/* Assuming all pixels have equal weight, return the sub-pixel coordinate
   that is at the center of mass of the region defined by "region" */

Double_Vector *G(Region_COM)(Double_Vector *R(O(com)), Region *reg)
{ Indx_Type     *raster;
  Size_Type      rlen, cnt;
  int            ndims;
  Dimn_Type     *dims;
  int64          Sum[10], *sum;

  RasterCon *cont = (RasterCon *) reg;

  raster = cont->raster;
  rlen   = cont->rastlen;

  ndims = cont->ndims;
  dims  = cont->dims;

  if (com == NULL)
    com = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(ndims));
  else if (com->size < ndims)
    { Free_Array(com);
      com = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(ndims));
    }
  else
    { com->kind  = PLAIN_KIND;
      com->type  = FLOAT64_TYPE;
      com->scale = 64;
      com->ndims = 1;
      com->size  = ndims;
      com->dims[0] = ndims;
    }

  if (ndims > 10)
    sum = (int64 *) Guarded_Malloc(sizeof(uint64)*((size_t) ndims),"Region_COM");
  else
    sum = Sum;

  { Indx_Type  v, w, c, i;
    Dimn_Type  o;
    int        d;
    double    *x;

    x = AFLOAT64(com);

    cnt = 0;
    for (d = 0; d < ndims; d++)
      sum[d] = 0;

    for (i = 0; i < rlen; i += 2)
      { c = 1;
        o = dims[0];

        v = raster[i];
        w = (raster[i+1] - v) + 1;

        cnt    += w;
        sum[0] += (((((v%o) << 1) + (w-1)) * w) >> 1);
        for (d = 1; d < ndims; d++)
          { c *= o;
            o  = dims[d];
            sum[d] += ((v/c)%o) * w;
          }
      }

    if (cnt == 0)
      return (0);

    for (d = 0; d < ndims; d++)
      x[d] = (1.*sum[d]) / cnt;

    if (ndims > 10)
      free(sum);

    return (com);
  }
}

/* Assuming all pixels have equal weight, return the sub-pixel coordinate
   that is at the center of mass of the pixels in "region" that are level
   or brigher.                                                     */

Double_Vector *G(Region_Select_COM)(Double_Vector *R(O(com)), Region *reg,
                                    APart *source, Comparator cmprsn, Value level)
{ Array     *array = AForm_Array(source);
  boolean  (*test)(Indx_Type, void *);
  Indx_Type *raster;
  Size_Type  rlen, cnt;
  int        ndims;
  Dimn_Type *dims;
  int64      Sum[10], *sum;
  TestArg    arg, *argp = &arg;

  RasterCon *cont = (RasterCon *) reg;

  raster = cont->raster;
  rlen   = cont->rastlen;

  ndims = cont->ndims;
  dims  = cont->dims;

  if (com == NULL)
    com = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(ndims));
  else if (com->size < ndims)
    { Free_Array(com);
      com = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(ndims));
    }
  else
    { com->kind  = PLAIN_KIND;
      com->type  = FLOAT64_TYPE;
      com->scale = 64;
      com->ndims = 1;
      com->size  = ndims;
      com->dims[0] = ndims;
    }

  if (ndims > 10)
    sum = (int64 *) Guarded_Malloc(sizeof(uint64)*((size_t) ndims),"Region_Select_COM");
  else
    sum = Sum;

  switch (array->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        arg.value = array->data;
        arg.level_<u> = level.<u>;
        break;
    #END
  }

  test = Comparator_Table[6*array->type + cmprsn];

  { Indx_Type    v, w, c, i;
    Dimn_Type    n, o, wgt;
    int          d;
    double      *x;
    int64        sumx;

    x = AFLOAT64(com);

    cnt = 0;
    for (d = 0; d < ndims; d++)
      sum[d] = 0;

    for (i = 0; i < rlen; i += 2)
      { c = 1;
        o = dims[0];

        v = raster[i];
        w = raster[i+1];
        n = (Dimn_Type) (v % o);
        w = w + n;

        wgt  = 0;
        sumx = 0;
        for (; v <= w; v++, n++)
          if (test(v,argp))
            { sumx += n;
              wgt  += 1;
            }

        cnt    += wgt;
        sum[0] += sumx;

        for (d = 1; d < ndims; d++)
          { c *= o;
            o  = dims[d];
            sum[d] += wgt * (w/c)%o;
          }
      }

    if (cnt == 0)
      return (0);

    for (d = 0; d < ndims; d++)
      x[d] = (1.*sum[d]) / cnt;

    if (ndims > 10)
      free(sum);

    return (com);
  }
}
