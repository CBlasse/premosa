/*****************************************************************************************\
*                                                                                         *
*  Array Data Abstraction                                                                 *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  December 2008                                                                 *
*                                                                                         *
*  (c) July 27, '09, Dr. Gene Myers and Howard Hughes Medical Institute                   *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/
 
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>

#include <stdint.h>
#include <float.h>

#define UINT8_MIN  0
#define UINT16_MIN 0
#define UINT32_MIN 0
#define UINT64_MIN 0

#include "mylib.h"
#include "utilities.h"
#include "array.h"
#include "linear.algebra.h"

#define BND_ZERO    0
#define BND_REFLECT 1
#define BND_WRAP    2
#define BND_EXTEND  3
#define BND_INVERT  4

extern int Boundary_Case_8qm5;

#LISTDEF @KINDS = PLAIN_KIND RGB_KIND RGBA_KIND COMPLEX_KIND

#LISTDEF @TYPES  =  UINT8 UINT16 UINT32 UINT64  INT8 INT16 INT32 INT64 FLOAT32 FLOAT64
#LISTDEF @UNION  =   uval   uval   uval   uval  ival  ival  ival  ival    fval    fval
#LISTDEF @CASTES = uint64 uint64 uint64 uint64 int64 int64 int64 int64  double  double
#LISTDEF @SIZES  =      1      2      4      8     1     2     4     8       4       8

#LISTDEF @OPS  = SET_OP ADD_OP SUB_OP MUL_OP DIV_OP POW_OP LSH_OP RSH_OP MIN_OP MAX_OP
#LISTDEF @SYMS =      =      +      -      *      /      =     <<     >>      >      <

#define SIZEOF(x) ((int) sizeof(x))

/****************************************************************************************
 *                                                                                      *
 *  ARRAY SIZE, SHAPE & TYPE COMPARISONS                                                *
 *                                                                                      *
 ****************************************************************************************/

static Size_Type array_size(Array *a)
{ Size_Type p;
  int       i;
  p = 1;
  for (i = 0; i < a->ndims; i++)
    p *= a->dims[i];
  return (p);
}

boolean Same_Shape(AForm *a, AForm *b)
{ int         i;
  boolean     eq;
  int         andims, bndims;
  Dimn_Type  *adims, *bdims;
  Coordinate *ashape, *bshape;

  if (Is_Array(a))
    { andims = ((Array *) a)->ndims;
      adims  = ((Array *) a)->dims;
      ashape = NULL;
    }
  else
    { if (Is_Slice(a))
        ashape = AForm_Shape(a);
      else
        ashape = Frame_Shape(a);
      andims = ashape->dims[0];
      adims  = ADIMN(ashape);
    }

  if (Is_Array(b))
    { bndims = ((Array *) b)->ndims;
      bdims  = ((Array *) b)->dims;
      bshape = NULL;
    }
  else
    { if (Is_Slice(b))
        bshape = AForm_Shape(b);
      else
        bshape = Frame_Shape(b);
      bndims = bshape->dims[0];
      bdims  = ADIMN(bshape);
    }
    
  eq = 0;
  if (andims == bndims)
    { for (i = 0; i < andims; i++)
        if (adims[i] != bdims[i])
          break;
      eq = (i >= andims);
    }

  if (Is_Slice(a))
    Free_Array(ashape);
  if (Is_Slice(b))
    Free_Array(bshape);

  return (eq);
}

boolean Same_Type(AForm *a, AForm *b)
{ if (AForm_Array(a)->type != AForm_Array(b)->type)
    return (0);
  return (Same_Shape(a,b));
}


/****************************************************************************************
 *                                                                                      *
 *  ARRAY SPACE MANAGEMENT ROUTINES                                                     *
 *                                                                                      *
 ****************************************************************************************/

// Awk-generated (manager.awk) Array memory management

static int type_size[] = { 1, 2, 4, 8, 1, 2, 4, 8, 4, 8 };

static int bit_size[] = { 8, 16, 32, 64, 8, 16, 32, 64, 32, 64 };

static int kind_size[] = { 1, 3, 4, 2 };

static Value_Kind type2kind[] = { UVAL, UVAL, UVAL, UVAL, IVAL, IVAL, IVAL, IVAL, FVAL, FVAL };

static inline int array_nsize(Array *a)
{ return (SIZEOF(int)*a->ndims); }

static inline int64 array_dsize(Array *a)
{ return (a->size*type_size[a->type]); }

static inline int array_tsize(Array *a)
{ return (a->tlen+1); }

MANAGER -IO Array dims:nsize data!dsize text:tsize


/****************************************************************************************
 *                                                                                      *
 *  MAKE AN ARRAY                                                                       *
 *                                                                                      *
 ****************************************************************************************/

static Array *make_start(Array_Kind kind, Value_Type type, int ndims, string routine)
{ Array  *a;

  a = new_array(SIZEOF(int)*ndims,0,1,routine);
  a->type    = type;
  a->kind    = kind;
  a->scale   = bit_size[type];
  a->text[0] = '\0';
  a->tlen    = 0;
  a->ndims   = ndims;
  return (a);
}

static Array *make_shape(Array_Kind kind, Value_Type type, int ndims, Dimn_Type *dims,
                         string routine)
{ Array *a;
  int    i, o, v;

  o = (kind != PLAIN_KIND);
  v = (kind == COMPLEX_KIND);
  a = make_start(kind,type,ndims+o,routine);
  if (v)
    a->dims[0] = kind_size[kind];
  for (i = 0; i < ndims; i++)
    a->dims[i+v] = dims[i];
  if (o & !v)
    a->dims[ndims] = kind_size[kind];
  a->size = array_size(a);
  return (a);
}

Array *G(Make_Array)(Array_Kind kind, Value_Type type, int ndims, Dimn_Type *dims)
{ Array *a;

  a = make_shape(kind,type,ndims,dims,"Make_Array");

  allocate_array_data(a,array_dsize(a),"Make_Array");

  return (a);
}

Array *G(Make_Array_With_Shape)(Array_Kind kind, Value_Type type, Coordinate *F(shape))
{ Array *a;

  a = make_shape(kind,type,shape->dims[0],ADIMN(shape),"Make_Array_With_Shape");

  allocate_array_data(a,array_dsize(a),"Make_Array_With_Shape");

  Free_Array(shape);

  return (a);
}

Array *G(Make_Array_Of_Data)(Array_Kind kind, Value_Type type, int ndims, Dimn_Type *dims,
                             void *data)
{ Array *a;

  a = make_shape(kind,type,ndims,dims,"Make_Array_Of_Data");

  { _Array *object = (_Array *) (((char *) a) - Array_Offset);
    if (object->dsize > 0)
      free(a->data);
    object->dsize  = array_dsize(a);
    a->data = data;
  }

  return (a);
}

Array *G(Make_Array_From_Arrays)(Array_Kind kind, int n, Array **arrays)
{ Array    *a, *a0;
  int       i, ndims;
  Size_Type dsize;

  a0 = arrays[0];
  if (n > 1)
    { for (i = 1; i < n; i++)
        if ( ! Same_Type(arrays[i],a0))
          { fprintf(stderr,
               "Arrays are not all of the same type and shape (Make_Array_From_Arrays)\n");
            exit (1);
          }
    }
  if (kind == COMPLEX_KIND && a0->ndims > 0)
    { for (i = 0; i < n; i++)
        if (arrays[i]->kind != COMPLEX_KIND)
          { fprintf(stderr,"Arrays must all be COMPLEX (Make_Array_From_Arrays)\n");
            exit (1);
          }
    }
  else
    { for (i = 0; i < n; i++)
        if (arrays[i]->kind != PLAIN_KIND)
          { fprintf(stderr,"Arrays must all be PLAIN (Make_Array_From_Arrays)\n");
            exit (1);
          }
    }
  if (kind != PLAIN_KIND && (kind != COMPLEX_KIND || a0->ndims == 0) && n != kind_size[kind])
    { fprintf(stderr,"Outer dimension and kind are not consistent (Make_Array_From_Arrays)\n");
      exit (1);
    }

  dsize = array_dsize(arrays[0]);
  ndims = a0->ndims;

  a = new_array(SIZEOF(int)*(ndims+1),n*dsize,1,"Make_Array_From_Arrays");

  a->type    = a0->type;
  a->kind    = kind;
  a->scale   = a0->scale;
  a->tlen    = 0;
  a->text[0] = '\0';
  a->ndims   = ndims+1;
  for (i = 0; i < ndims; i++)
    a->dims[i] = a0->dims[i];
  a->dims[ndims] = n;
  a->size        = n * a0->size;

  for (i = 0; i < n; i++)
    memcpy(((char *) a->data) + i*dsize, arrays[i]->data, (size_t) dsize);

  return (a);
}


/****************************************************************************************
 *                                                                                      *
 *  COORDS, GET & SET                                                                   *
 *                                                                                      *
 ****************************************************************************************/

static Dimn_Type  Coord_Ndim, Coord_Max = 0;
static Dimn_Type *Coord_Dims;
static Array_Kind Coord_Kind;

static pthread_mutex_t Coord_Mutex = PTHREAD_MUTEX_INITIALIZER;

void Set_Coord_Basis(Coordinate *F(point), Array_Kind kind)
{ pthread_mutex_lock(&Coord_Mutex);
  if (point == NULL)
    { Coord_Dims = NULL;
      Coord_Ndim = 0;
      return;
    }
  if (point->ndims != 1 || point->type != DIMN_TYPE)
    { fprintf(stderr,"Coordinate is not an unsigned integer vector (Set_Coord_Basis)\n");
      exit (1);
    }
  Coord_Ndim = point->dims[0];
  Coord_Kind = kind;
  if (Coord_Ndim > Coord_Max)
    { Coord_Max  = Coord_Ndim + 10;
      Coord_Dims = Guarded_Realloc(Coord_Dims,sizeof(Dimn_Type)*((size_t) Coord_Max),
                                   "Set_Coord_Basis");
    }
  memcpy(Coord_Dims,ADIMN(point),sizeof(Dimn_Type)*((size_t) Coord_Ndim));
  pthread_mutex_unlock(&Coord_Mutex);
  Free_Array(point);
}

Coordinate *G(Get_Coord_Basis)(Array_Kind *kind)
{ Dimn_Type  dim[1];
  Array     *point;

  if (Coord_Dims == NULL)
    return (NULL);
  pthread_mutex_lock(&Coord_Mutex);
  dim[0] = Coord_Ndim;
  point  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);
  memcpy(ADIMN(point),Coord_Dims,sizeof(Dimn_Type)*((size_t) Coord_Ndim));
  *kind = Coord_Kind;
  pthread_mutex_unlock(&Coord_Mutex);
  return (point);
}

void Use_Array_Basis(Array *a)
{ pthread_mutex_lock(&Coord_Mutex);
  Coord_Ndim = a->ndims;
  Coord_Kind = a->kind;
  if (Coord_Ndim > Coord_Max)
    { Coord_Max  = Coord_Ndim + 10;
      Coord_Dims = Guarded_Realloc(Coord_Dims,sizeof(Dimn_Type)*((size_t) Coord_Max),
                                   "Use_Array_Basis");
    }
  memcpy(Coord_Dims,a->dims,sizeof(Dimn_Type)*((size_t) Coord_Ndim));
  pthread_mutex_unlock(&Coord_Mutex);
}

Coordinate *G(Coord)(string list)
{ Array     *coord;
  Dimn_Type  dim[1], *c;
  int        i, n;
  char      *b, *e;

  n = 1;
  for (i = 0; i < (int) strlen(list); i++)
    if (list[i] == ',')
      n += 1;

  dim[0] = n;
  coord  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);

  c = ADIMN(coord);
  b = list;
  while (1)
    { c[--n] = (Dimn_Type) strtol(b,&e,10);
      if (b == e)
        { fprintf(stderr,"Not a valid constant coordinate list (Coord)\n");
          exit (1);
        }
      b = e;
      if (*b == '\0')
        break;
      if (*b++ != ',')
        { fprintf(stderr,"Not a valid constant coordinate list (Coord)\n");
          exit (1);
        }
    }

  return (coord);
}

void Print_Coord(FILE *file, Coordinate *point)
{ int  i, n;

  n = point->dims[0];
  switch (point->type) {
    #GENERATE T,F = @TYPES , %u %u %u %llu %d %d %d %lld %g %g
      case <T>_TYPE:
        { <t> *d = A<T>(point);
          fprintf(file,"<F>",d[n-1]);
          for (i = n-2; i >= 0; i--)
            fprintf(file,",<F>",d[i]);
          break;
        }
    #END
  }

  Free_Array(point);
}

Coordinate *G(Coord1)(Dimn_Type d1)
{ Array     *coord;
  Dimn_Type  dim[1], *c;

  dim[0] = 1;
  coord  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);
  c      = ADIMN(coord);
  c[0]   = d1;
  return (coord);
}

Coordinate *G(Coord2)(Dimn_Type d2, Dimn_Type d1)
{ Array     *coord;
  Dimn_Type  dim[1], *c;

  dim[0] = 2;
  coord  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);
  c      = ADIMN(coord);
  c[0]   = d1;
  c[1]   = d2;
  return (coord);
}

Coordinate *G(Coord3)(Dimn_Type d3, Dimn_Type d2, Dimn_Type d1)
{ Array     *coord;
  Dimn_Type  dim[1], *c;

  dim[0] = 3;
  coord  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);
  c      = ADIMN(coord);
  c[0]   = d1;
  c[1]   = d2;
  c[2]   = d3;
  return (coord);
}

Coordinate *G(Coord4)(Dimn_Type d4, Dimn_Type d3, Dimn_Type d2, Dimn_Type d1)
{ Array     *coord;
  Dimn_Type  dim[1], *c;

  dim[0] = 4;
  coord  = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dim);
  c      = ADIMN(coord);
  c[0]   = d1;
  c[1]   = d2;
  c[2]   = d3;
  c[3]   = d4;
  return (coord);
}

Indx_Type Coord2Idx(Coordinate *F(point))
{ Dimn_Type *x, y;
  Indx_Type  p;
  int        i;

  if (point->ndims != 1 || point->type != DIMN_TYPE)
    { fprintf(stderr,"Point is not an integer vector (Coord2Idx)\n");
      exit (1);
    }
  if (Coord_Dims == NULL)
    { fprintf(stderr,"Coordinate basis is not set (Coord2Idx)\n");
      exit (1);
    }
  pthread_mutex_lock(&Coord_Mutex);
  if (point->size != Coord_Ndim && (point->size != Coord_Ndim-1 || Coord_Kind == PLAIN_KIND))
    { fprintf(stderr,"Coordinate dimensionality doesn't match that of current basis (Coord2Idx)\n");
      exit (1);
    }
  if (point->size != Coord_Ndim)
    { if (Coord_Kind == COMPLEX_KIND)
        PrependCoord(point,0);
      else
        AppendCoord(0,point);
    }

  x = ADIMN(point);
  p = 0;
  for (i = Coord_Ndim-1; i >= 0; i--)
    { y = x[i];
      if (y >= Coord_Dims[i] || y < 0)
        { fprintf(stderr,"%d'th index out of bounds (Coord2Idx)\n",i+1);
          exit (1);
        }
      p = p*Coord_Dims[i] + y; 
    }
  pthread_mutex_unlock(&Coord_Mutex);

  Free_Array(point);

  return (p);
}

Indx_Type Coord2IdxA(Array *a, Coordinate *F(point))
{ Dimn_Type *x, *d, y;
  Indx_Type  p;
  int        n, i;

  if (point->ndims != 1 || point->type != DIMN_TYPE)
    { fprintf(stderr,"Point is not an integer vector (Coord2IdxA)\n");
      exit (1);
    }

  n = a->ndims;
  d = a->dims;

  if (point->size != n && (point->size != n-1 || a->kind == PLAIN_KIND))
   { fprintf(stderr,"Coordinate dimensionality doesn't match that of current basis (Coord2IdxA)\n");
     exit (1);
   }
  if (point->size != n)
    { if (Coord_Kind == COMPLEX_KIND)
        PrependCoord(point,0);
      else
        AppendCoord(0,point);
    }

  x = ADIMN(point);
  p = 0;
  for (i = n-1; i >= 0; i--)
    { y = x[i];
      if (y >= d[i] || y < 0)
        { fprintf(stderr,"%d'th index out of bounds (Coord2IdxA)\n",i+1);
          exit (1);
        }
      p = p*d[i] + y; 
    }

  Free_Array(point);

  return (p);
}

Coordinate *AppendCoord(Dimn_Type d, Coordinate *R(M(c)))
{ Dimn_Type n;
  n = c->dims[0];
  allocate_array_data(c,(n+1)*SIZEOF(Dimn_Type),"AppendCoord");
  c->dims[0]  = n+1;
  c->size     = n+1;
  ADIMN(c)[n] = d;
  return (c);
}

Coordinate *PrependCoord(Coordinate *R(M(c)), Dimn_Type d)
{ Dimn_Type *x = ADIMN(c);
  Dimn_Type  n;
  int        i;

  n = c->dims[0];
  allocate_array_data(c,(n+1)*SIZEOF(Dimn_Type),"PrependCoord");
  c->dims[0]  = n+1;
  c->size     = n+1;
  for (i = n; i > 0; i--)
    x[i] = x[i-1];
  x[0] = d;
  return (c);
}

static inline Coordinate *G(coord4idx)(int n, Dimn_Type *d, Indx_Type idx, string routine)
{ Dimn_Type   dims[1], m, *l;
  Coordinate *lat;
  int         i;

  dims[0] = n;
  lat     = Make_Array(PLAIN_KIND,DIMN_TYPE,1,dims);

  l = ADIMN(lat);
  for (i = 0; i < n; i++)
    { m = d[i];
      l[i] = (Dimn_Type) (idx % m);
      idx /= m;
    }
  if (idx > 0)
    { fprintf(stderr,"Index is out of array boundary (%s)\n",routine);
      exit (1);
    }
  return (lat);
}

Coordinate *G(Idx2Coord)(Indx_Type idx)
{ Coordinate *p;
  if (Coord_Dims == NULL)
    { fprintf(stderr,"Coordinate basis is not set (Idx2Coord)\n");
      exit (1);
    }
  pthread_mutex_lock(&Coord_Mutex);
  p = coord4idx(Coord_Ndim,Coord_Dims,idx,"Idx2Coord");
  pthread_mutex_unlock(&Coord_Mutex);
  return (p);
}

Coordinate *G(Idx2CoordA)(Array *a, Indx_Type idx)
{ return (coord4idx(a->ndims,a->dims,idx,"Idx2CoordA")); }

Coordinate *G(Idx2Core)(Indx_Type idx)
{ Dimn_Type  *d;
  int         n;
  Coordinate *p;

  if (Coord_Dims == NULL)
    { fprintf(stderr,"Coordinate basis is not set (Idx2Coord)\n");
      exit (1);
    }

  pthread_mutex_lock(&Coord_Mutex);
  d = Coord_Dims;
  n = Coord_Ndim;
  if (Coord_Kind != PLAIN_KIND)
    { if (Coord_Kind == COMPLEX_KIND)
        d += 1;
      n -= 1;
      idx /= kind_size[Coord_Kind];
    }

  p = coord4idx(n,d,idx,"Idx2Core");
  pthread_mutex_unlock(&Coord_Mutex);
  return (p);
}

Coordinate *G(Idx2CoreA)(Array *a, Indx_Type idx)
{ Dimn_Type  *d;
  int         n;

  d = a->dims;
  n = a->ndims;
  if (a->kind != PLAIN_KIND)
    { if (a->kind == COMPLEX_KIND)
        d += 1;
      n -= 1;
      idx /= kind_size[a->kind];
    }

  return (coord4idx(n,d,idx,"Idx2CoreA"));
}

Value Get_Array_Value(Array *a, Coordinate *coord)
{ Indx_Type p = Coord2IdxA(a,coord);
  Value     v;
  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        v.<U> = A<T>(a)[p];
        break;
    #END
  }
  return (v);
}

void Set_Array_Value(Array *M(a), Coordinate *coord, Value v)
{ Indx_Type p = Coord2IdxA(a,coord);
  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        A<T>(a)[p] = (<t>) v.<U>;
        break;
    #END
  }
}

Coordinate *G(Floor_Coord)(Double_Vector *point)
{ Array     *lat;
  double    *p;
  Dimn_Type *c, k;

  if (point->ndims != 1 || point->kind != PLAIN_KIND || point->type != FLOAT64_TYPE)
    { fprintf(stderr,"Point is not a double vector (Floor_Coord)\n");
      exit (1);
    }

  lat = Make_Array(PLAIN_KIND,DIMN_TYPE,1,point->dims);

  p = AFLOAT64(point);
  c = ADIMN(lat);
  for (k = 0; k < point->dims[0]; k++)
    c[k] = (Dimn_Type) floor(p[k]);

  return (lat);
}

Coordinate *G(Ceiling_Coord)(Double_Vector *point)
{ Array     *lat;
  double    *p;
  Dimn_Type *c, k;

  if (point->ndims != 1 || point->kind != PLAIN_KIND || point->type != FLOAT64_TYPE)
    { fprintf(stderr,"Point is not a double vector (Ceiling_Coord)\n");
      exit (1);
    }

  lat = Make_Array(PLAIN_KIND,DIMN_TYPE,1,point->dims);

  p = AFLOAT64(point);
  c = ADIMN(lat);
  for (k = 0; k < point->dims[0]; k++)
    c[k] = (Dimn_Type) ceil(p[k]);

  return (lat);
}

Coordinate *G(Nearest_Coord)(Double_Vector *point)
{ Array     *lat;
  double    *p;
  Dimn_Type *c, k;

  if (point->ndims != 1 || point->kind != PLAIN_KIND || point->type != FLOAT64_TYPE)
    { fprintf(stderr,"Point is not a double vector (Floor_Coord)\n");
      exit (1);
    }

  lat = Make_Array(PLAIN_KIND,DIMN_TYPE,1,point->dims);

  p = AFLOAT64(point);
  c = ADIMN(lat);
  for (k = 0; k < point->dims[0]; k++)
    c[k] = (Dimn_Type) floor(p[k]+.5);

  return (lat);
}


/****************************************************************************************
 *                                                                                      *
 *  SLICE ITERATOR ROUTINES                                                             *
 *                                                                                      *
 ****************************************************************************************/

#define SLICE_KIND  4

typedef struct
  { Array_Kind  kind;     //  Always SLICE_KIND in order to distinguish from an Array
    Array      *trg_ref;  //  The array the slice is in.
    Coordinate *beg;      //  References to the beg and end coordinate defining slice
    Coordinate *end;

    int         ndims;    //  Dimensionality of slice
    Size_Type   size;     //  Size of slice (# of pixels)

    Coordinate *cnt;      //  Current coordinate (= index)
    Indx_Type   p;        //  Current index
    int         clip;     //  Highest dimension currently outside of slice
    Size_Type  *dnc;      //  dnc[i] = displacement to next slice element when at an i boundary

    Dimn_Type  *acnt;     //  Direct access to elements above to speed computation
    Dimn_Type  *bcrd;
    Dimn_Type  *ecrd;
  } Slicer;

static inline int slicer_dsize(Slicer *s)
{ return (SIZEOF(Size_Type)*s->ndims); }

MANAGER -p Slice(Slicer) cnt*Array beg*Array end*Array trg_ref@Array dnc:dsize

Slice *Pack_Slice(Slice *s)
{ Slicer *o = (Slice *) s;
  boolean nok = pack_slicer(o);
  o->acnt = ADIMN(o->cnt);
  o->bcrd = ADIMN(o->beg);
  o->ecrd = ADIMN(o->end);
  if (nok) return (NULL);
  return (s);
}

Slice *G(Make_Slice)(Array *I(target), Coordinate *S(beg), Coordinate *S(end))
{ Slicer    *slice;
  Dimn_Type *ecrd = ADIMN(end);
  Dimn_Type *bcrd = ADIMN(beg);
  Dimn_Type *dims;
  int        i, ndims;

  if (beg->ndims != 1 || beg->type != DIMN_TYPE)
    { fprintf(stderr,"Beg is not a coordinate vector (Make_Slice)\n");
      exit (1);
    }
  if (end->ndims != 1 || end->type != DIMN_TYPE)
    { fprintf(stderr,"End is not a coordinate vector (Make_Slice)\n");
      exit (1);
    }
  if (Array_Refcount(beg) != 1)
    { fprintf(stderr,"Beg is not subsumable, has a reference count > 1 (Make_Slice)\n");
      exit (1);
    }
  if (Array_Refcount(end) != 1)
    { fprintf(stderr,"End is not subsumable, has a reference count > 1 (Make_Slice)\n");
      exit (1);
    }

  ndims = (int) beg->size;
  if (target->ndims != ndims && (target->ndims-1 != ndims || target->kind == PLAIN_KIND))
    { fprintf(stderr,"Target and coordinate dimensionality do not match (Make_Slice)\n");
      exit (1);
    }
  if (end->size != ndims)
    { fprintf(stderr,"Begin and end coordinate dimensionality do not match (Make_Slice)\n");
      exit (1);
    }

  if (ndims == target->ndims-1 && target->kind == COMPLEX_KIND)
    dims = target->dims+1;
  else
    dims = target->dims;
  for (i = 0; i < ndims; i++)
    { if (bcrd[i] > ecrd[i])
        { fprintf(stderr,"beg is not before end in target array (Make_Slice)\n");
          exit (1);
        }
      if (ecrd[i] >= dims[i])
        { fprintf(stderr,"end is not in basis of target array (Make_Slice)\n");
          exit (1);
        }
      if (bcrd[i] < 0)
        { fprintf(stderr,"beg is not in basis of target array (Make_Slice)\n");
          exit (1);
        }
    }
  dims = target->dims;

  if (ndims == target->ndims-1)
    { ndims += 1;
      if (target->kind == COMPLEX_KIND)
        { PrependCoord(beg,0);
          PrependCoord(end,1);
        }
      else
        { AppendCoord(0,beg);
          AppendCoord(kind_size[target->kind]-1,end);
        }
      bcrd = ADIMN(beg);
      ecrd = ADIMN(end);
    }

  slice = new_slicer(ndims*SIZEOF(Size_Type),"Make_Slice");

  slice->kind    = SLICE_KIND;
  slice->trg_ref = Inc_Array(target);
  slice->beg     = beg;
  slice->end     = end;

  slice->ndims = ndims;
  slice->size  = 1;
  for (i = 0; i < ndims; i++)
    slice->size *= (ecrd[i] - bcrd[i]) + 1;

  slice->cnt = Copy_Array(beg);

  { Size_Type *dinc = slice->dnc;
    Size_Type  offset, outer;

    offset  = 0;
    outer   = 1;
    for (i = 0; i < ndims; i++)
      { dinc[i] = outer-offset;
        offset += (ecrd[i]-bcrd[i])*outer;
        outer  *= dims[i];
      }

    slice->p = Coord2IdxA(target,Inc_Array(beg));
    slice->clip = -1;
  }

  slice->acnt = ADIMN(slice->cnt);
  slice->bcrd = bcrd;
  slice->ecrd = ecrd;

  return ((Slice *) slice);
}

Indx_Type Slice_Index(Slice *slicer)
{ return (((Slicer *) slicer)->p);  }

Coordinate *Slice_Coordinate(Slice *slicer)
{ return (((Slicer *) slicer)->cnt); }

Coordinate *Slice_First(Slice *slicer)
{ return (((Slicer *) slicer)->beg); }

Coordinate *Slice_Last(Slice *slicer)
{ return (((Slicer *) slicer)->end); }

boolean Inside_Slice(Slice *slicer)
{ return (((Slicer *) slicer)->clip < 0);  }

boolean Set_Slice_To_Index(Slice *M(slicer), Size_Type idx)
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *bcrd  = slice->bcrd;
  Dimn_Type *ecrd  = slice->ecrd;
  Dimn_Type *cnt   = slice->acnt;
  Dimn_Type *dim   = slice->trg_ref->dims;
  int        ndims = slice->ndims;
  Dimn_Type  d, c;
  int        i;

  if (idx >= slice->trg_ref->size)
    { fprintf(stderr,"Index is not in target array basis (Set_Slice_To_Index)\n");
      exit (1);
    }

  slice->p    = idx;
  slice->clip = -1;
  for (i = 0; i < ndims; i++)
    { d = dim[i];
      cnt[i] = c = (Dimn_Type) (idx % d);
      if (c < bcrd[i] || c > ecrd[i])
        slice->clip = i;
      idx = idx / d;
    }
  return (slice->clip < 0);
}

Indx_Type Set_Slice_To_First(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *bcrd  = slice->bcrd;
  Dimn_Type *dim   = slice->trg_ref->dims;
  Dimn_Type *cnt   = slice->acnt;
  int        ndims = slice->ndims;
  int        i;

  slice->p = 0;
  for (i = ndims-1; i >= 0; i--)
    { cnt[i] = bcrd[i];
      slice->p = slice->p * dim[i] + bcrd[i];
    }
  slice->clip = -1;
  return (slice->p);
}

Indx_Type Set_Slice_To_Last(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *ecrd  = slice->ecrd;
  Dimn_Type *dim   = slice->trg_ref->dims;
  Dimn_Type *cnt   = slice->acnt;
  int        ndims = slice->ndims;
  int        i;

  slice->p = 0;
  for (i = ndims-1; i >= 0; i--)
    { cnt[i] = ecrd[i];
      slice->p = slice->p * dim[i] + ecrd[i];
    }
  slice->clip = -1;
  return (slice->p);
}

Indx_Type Next_Slice_Index(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *ecrd  = slice->ecrd;
  Dimn_Type *cnt   = slice->acnt;
  int        i;

  if (slice->clip >= 0)
    { fprintf(stderr,"Must be in slice to move to next position (Next_Slice_Index)\n");
      exit (1);
    }

  if (++cnt[0] <= ecrd[0])
    return (slice->p += 1);
  else
    { int        ndims = slice->ndims;
      Dimn_Type *bcrd  = slice->bcrd;

      cnt[0] = bcrd[0];
      for (i = 1; i < ndims; i++)
        if (++cnt[i] > ecrd[i])
          cnt[i] = bcrd[i];
        else
          return (slice->p += slice->dnc[i]);
    }
  slice->p = Coord2IdxA(slice->trg_ref,Inc_Array(slice->beg));
  return (slice->p);
}

Indx_Type Prev_Slice_Index(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *bcrd  = slice->bcrd;
  Dimn_Type *cnt   = slice->acnt;
  int        i;

  if (slice->clip >= 0)
    { fprintf(stderr,"Must be in slice to move to next position (Prev_Slice_Index)\n");
      exit (1);
    }

  if (cnt[0]-- > bcrd[0])
    return (slice->p -= 1);
  else
    { int        ndims = slice->ndims;
      Dimn_Type *ecrd  = slice->ecrd;

      cnt[0] = ecrd[0];
      for (i = 1; i < ndims; i++)
        if (cnt[i]-- <= bcrd[i])
          cnt[i] = ecrd[i];
        else
          return (slice->p -= slice->dnc[i]);
    }
  slice->p = Coord2IdxA(slice->trg_ref,Inc_Array(slice->end));
  return (slice->p);
}

boolean Inc_Slice_Index(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *tcrd  = slice->trg_ref->dims;
  Dimn_Type *bcrd  = slice->bcrd;
  Dimn_Type *ecrd  = slice->ecrd;
  Dimn_Type *cnt   = slice->acnt;
  int        ndims = slice->ndims;
  int        i, nclip, clip;
  
  clip  = slice->clip;
  nclip = -1;
  for (i = 0; i < ndims; i++)
    if (++cnt[i] == tcrd[i])
      { cnt[i] = 0;
        if (cnt[i] < bcrd[i])
          nclip = i;
        else if (clip == i)
          clip = -1;
      }
    else
      { if (cnt[i] > ecrd[i])
          nclip = i;
        else if (clip == i && cnt[i] >= bcrd[i])
          clip = -1;
        break;
      }
  if (nclip > clip)
    clip = nclip;
  slice->clip = clip;
  if (i >= ndims)
    slice->p = 0;
  else
    slice->p += 1;
  return (clip < 0);
}

boolean Dec_Slice_Index(Slice *M(slicer))
{ Slicer    *slice = (Slicer *) slicer;
  Dimn_Type *tcrd  = slice->trg_ref->dims;
  Dimn_Type *bcrd  = slice->bcrd;
  Dimn_Type *ecrd  = slice->ecrd;
  Dimn_Type *cnt   = slice->acnt;
  int        ndims = slice->ndims;
  int        i, nclip, clip;
  
  clip  = slice->clip;
  nclip = -1;
  for (i = 0; i < ndims; i++)
    if (cnt[i]-- == 0)
      { cnt[i] = tcrd[i]-1;
        if (cnt[i] > ecrd[i])
          nclip = i;
        else if (clip == i)
          clip = -1;
      }
    else
      { if (cnt[i] < bcrd[i])
          nclip = i;
        else if (clip == i && cnt[i] <= ecrd[i])
          clip = -1;
        break;
      }
  if (nclip > clip)
    clip = nclip;
  slice->clip = clip;
  if (i >= ndims)
    slice->p = slice->trg_ref->size-1;
  else
    slice->p -= 1;
  return (clip < 0);
}

Array *G(Make_Array_From_Slice)(Slice *slice)
{ Array *source = AForm_Array(slice);
  Array *target = Make_Array_With_Shape(PLAIN_KIND,source->type,AForm_Shape(slice));
  Size_Type n   = AForm_Size(slice);
  Indx_Type r, p;

  target->kind = AForm_Kind(slice);
  switch (source->type) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *d = A<T>(source);
          <t> *a = A<T>(target);

          for (r = 0, p = Set_Slice_To_First(slice); r < n; p = Next_Slice_Index(slice))
            a[r++] = d[p];
          break;
        }
    #END
  }

  return (target);
}


/****************************************************************************************
 *                                                                                      *
 *  FRAME ITERATOR ROUTINES                                                             *
 *                                                                                      *
 ****************************************************************************************/

#define FRAME_KIND 5

typedef struct
  { Array_Kind  kind;     //  Always FRAME_KIND in order to distinguish from an Array
    Array      *trg_ref;  //  The array the frame is in.
    Coordinate *shape;    //  References to the shape and anchor coordinates defining frame
    Array      *anchor;

    Slice      *slice;    //  Frame's in-bound slice
    int         notempty; //  Frame can be in bounds (otherwise a fake slice and ignore in-status)

    int         ndims;    //  Dimensionality of the frame
    int64       size;     //  Size of the frame
    int64       vsize;    //  Size of value buffer

    Offs_Type  *offs;     //  Offset array of dimensions 'shape' that give offsets when within
    void       *vals;     //  Value array that holds frame values when not within
    int         mode;     //  Status of value array
  } Framer;

#define NOT_SET  0   //  Possible values for mode field
#define USE_VALS 1
#define USE_DBLS 2

static inline int64 framer_osize(Framer *f)
{ return (SIZEOF(int64)*f->size + SIZEOF(double)*f->vsize); }

MANAGER -p Frame(Framer) trg_ref@Array shape*Array anchor*Array slice*Slice offs!osize

Frame *Pack_Frame(Frame *f)
{ Framer *o = (Frame *) f;
  boolean nok = pack_framer(o);
  o->vals = (void *) (o->offs + o->size);
  if (nok) return (NULL);
  return (f);
}

/****************************************************************************************
 *                                                                                      *
 *  OFFSET MAP AND BOUNDARY VALUE DETERMINATION ROUTINES                                *
 *                                                                                      *
 ****************************************************************************************/

typedef struct
  { Dimn_Type *ffdim;    //  Used by all
    Dimn_Type *fidim;    //  Used by all
    int32     *fctr;     //  Used by frame_offsets only
    int64     *fbcrd;    //  Remainder used by other routines
    int64     *fecrd;
    int64     *fdvol;
    void     **fibuf;
  } Frame_Args;

static Offs_Type *frame_offsets(int d, Offs_Type q, Offs_Type *o, Frame_Args *args)
{ Dimn_Type j;
  Dimn_Type n = args->ffdim[d];

  q = q*args->fidim[d] - args->fctr[d];
  if (d == 0)
    for (j = 0; j < n; j++)
      *o++ = q++;
  else
    for (j = 0; j < n; j++)
      o = frame_offsets(d-1,q++,o,args);
  return (o);
}

static inline int64 min_int64(int64 a, int64 b)
{ if (a < b)
    return (a);
  else
    return (b);
}

#GENERATE T,U = @TYPES , @CASTES

static <t> *frame_zerod_<t>(int d, int64 size, <t> *o, Dimn_Type *fdim)
{ while (d >= 0)
    size *= fdim[d--];
  while (size-- > 0)
    *o++ = 0;
  return (o);
}

static <t> *frame_zero_<t>(int d, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->fidim[d];

  q *= f;
  j  = b;
  if (d == 0)
    { if (b < 0)
        { for (x = min_int64(0,e); j < x; j++)
            *o++ = 0;
          if (e <= 0)
            return (o);
        }
      for (x = min_int64(e,f); j < x; j++)
        *o++ = a[q+j];
      while (j++ < e)
        *o++ = 0;
    }
  else
    { d -= 1;
      if (b < 0)
        { o = frame_zerod_<t>(d,min_int64(0,e)-b,o,args->ffdim);
          if (e <= 0)
            return (o);
          j = 0;
        }
      for (x = min_int64(e,f); j < x; j++)
        o = frame_zero_<t>(d,q+j,o,a,args);
      if (e > f)
        o = frame_zerod_<t>(d,e-j,o,args->ffdim);
    }
  return (o);
}

static <t> *frame_wrap_<t>(int d, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, m, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->fidim[d];

  q *= f;
  if (d == 0)
    { if (b < 0)
        { m  = f-1;
          q += m;
          x  = m - min_int64(e,0);
          for (j = m-b; j > x; j--)
            *o++ = a[q - (j % f)];
          if (e <= 0)
            return (o);
          q -= m;
          j  = 0;
        }
      else
        j = b;
      for (x = min_int64(e,f); j < x; j++)
        *o++ = a[q+j];
      while (j < e)
        *o++ = a[q + ((j++) % f)];
    }
  else
    { d -= 1;
      if (b < 0)
        { m  = f-1;
          q += m;
          x  = m - min_int64(e,0);
          for (j = m-b; j > x; j--)
            o = frame_wrap_<t>(d,q - (j%f),o,a,args);
          if (e <= 0)
            return (o);
          q -= m;
          j  = 0;
        }
      else
        j = b;
      for (x = min_int64(e,f); j < x; j++)
        o = frame_wrap_<t>(d,q+j,o,a,args);
      while (j < e)
        o = frame_wrap_<t>(d,q + ((j++) % f),o,a,args);
    }
  return (o);
}

static <t> *frame_reflect_<t>(int d, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, k, m, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->fidim[d];

  q *= f;
  m  = ((f-1) << 1);
  if (d == 0)
    { if (b < 0)
        { x = -min_int64(e,0);
          for (j = -b; j > x; j--)
            { k = j % m;
              if (k >= f)
                k = m-k;
              *o++ = a[q+k];
            }
          if (e <= 0)
            return (o);
        }
      else
        j = b;
      for (x = min_int64(e,f); j < x; j++)
        *o++ = a[q+j];
      while (j < e)
        { k = (j++) % m;
          if (k >= f)
            k = m-k;
          *o++ = a[q+k];
        }
    }
  else
    { d -= 1;
      if (b < 0)
        { x = -min_int64(e,0);
          for (j = -b; j > x; j--)
            { k = j % m;
              if (k >= f)
                k = m-k;
              o = frame_reflect_<t>(d,q+k,o,a,args);
            }
          if (e <= 0)
            return (o);
        }
      else
        j = b;
      for (x = min_int64(e,f); j < x; j++)
        o = frame_reflect_<t>(d,q+j,o,a,args);
      while (j < e)
        { k = (j++) % m;
          if (k >= f)
            k = m-k;
          o = frame_reflect_<t>(d,q+k,o,a,args);
        }
    }
  return (o);
}

static <t> *frame_extend_<t>(int d, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->fidim[d];

  q *= f;
  j  = b;
  if (d == 0)
    { if (b < 0)
        { for (x = min_int64(0,e); j < x; j++)
            *o++ = a[q];
          if (e <= 0)
            return (o);
        }
      for (x = min_int64(e,f); j < x; j++)
        *o++ = a[q+j];
      q += f-1;
      while (j++ < e)
        *o++ = a[q];
    }
  else
    { d -= 1;
      if (b < 0)
        { for (x = min_int64(0,e); j < x; j++)
            o = frame_extend_<t>(d,q,o,a,args);
          if (e <= 0)
            return (o);
        }
      for (x = min_int64(e,f); j < x; j++)
        o = frame_extend_<t>(d,q+j,o,a,args);
      q += f-1;
      while (j++ < e)
        o = frame_extend_<t>(d,q,o,a,args);
    }
  return (o);
}

static void frame_invert_<t>(int d, int64 q, int64 t, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, k, m, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->fidim[d];

  q  = q*f;
  t  = t*args->ffdim[d];
  m  = ((f-1) << 1);

  if (d == 0)
    { <u>   w, ab, ae;

      a += q;
      o += t;
      ab = 2*a[0];
      ae = 2*a[f-1];

      j = b;
      if (b < 0)
        { j  = -b;
          x  = (j-1)/(f-1)+1;
          if (x % 2)
            w = ((<u>) (x/2))*(ab - ae) + ab;
          else
            w = ((<u>) (x/2))*(ab - ae);
          x = -min_int64(e,0);
          for (k = (j-1) % m + 1; j > x; j--)
            if (k >= f)
              { *o++ = (<t>) (w + a[m-k]);
                if (k-- == f)
                  w += ae;
              }
            else
              { *o++ = (<t>) (w - a[k]);
                if (k-- == 1)
                  { w -= ab;
                    k  = m;
                  }
              }
          if (e <= 0)
            return;
        }

      for (x = min_int64(e,f); j < x; j++)
        *o++ = (<t>) a[j];

      if (j < e)
        { x = (j-f)/(f-1) + 1;
          if (x % 2)
            w = ((<u>) (x/2))*(ae-ab) + ae;
          else
            w = ((<u>) (x/2))*(ae-ab);
          for (k = m - (j-1) % m; j < e; j++)
            if (k-- >= f)
              { *o++ = (<t>) (w + a[m-k]);
                if (k < f)
                  w += ae;
              }
            else
              { *o++ = (<t>) (w - a[k]);
                if (k == 0)
                  { w -= ab;
                    k  = m;
                  }
              }
        }
    }

  else
    { int64 lx, rx, y, S;
      <t>   *L, *R, *F;

      d -= 1;
      t -= b;

      S = args->fdvol[d];
      F = o;
      L = (<t> *) args->fibuf[d];
      R = L + S;

      frame_invert_<t>(d,q,0,L,a,args);
      frame_invert_<t>(d,q+(f-1),0,R,a,args);

      j = b;
      if (b < 0)
        { j  = -b;
          x  = (j-1)/(f-1)+1;
          if (x % 2)
            { lx = x+1;
              rx = -(x-1);
            }
          else
            { lx = x;
              rx = -x;
            }
          x = -min_int64(e,0);
          F = o + (t+b)*S;
          for (k = (j-1) % m + 1; j > x; j--)
            if (k >= f)
              { frame_invert_<t>(d,q + (m-k),t - j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) + ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] + F[y]);
#END
                F += S;

                if (k-- == f)
                  rx += 2;
              }
            else
              { frame_invert_<t>(d,q + k,t - j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) - ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] - F[y]);
#END
                F += S;

                if (k-- == 1)
                  { lx -= 2;
                    k   = m;
                  }
              }
          if (e <= 0)
            return;
        }

      for (x = min_int64(e,f); j < x; j++)
        frame_invert_<t>(d,q + j,t + j,o,a,args);

      if (j < e)
        { x = (j-f)/(f-1) + 1;
          if (x % 2)
            { rx = x+1;
              lx = -(x-1);
            }
          else
            { rx = x;
              lx = -x;
            }
          F = o + (t+j)*S;
          for (k = m - (j-1) % m; j < e; j++)
            if (k-- >= f)
              { frame_invert_<t>(d,q + (m-k),t + j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) + ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] + F[y]);
#END
                F += S;

                if (k < f)
                  rx += 2;
              }
            else
              { frame_invert_<t>(d,q + k,t + j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) - ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] - F[y]);
#END
                F += S;

                if (k == 0)
                  { lx -= 2;
                    k   = m;
                  }
              }
        }
    }
}

/****************************************************************************************
 *                                                                                      *
 *  BOUNDARY VALUE FILLING FOR EXPANSION (related to frame routines but subroutines     *
 *    of Pad_Array and Pad_Array_Inplace)                                               *
 *                                                                                      *
 ****************************************************************************************/

static <t> *expand_wrap_<t>(int d, int out, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, m;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->ffdim[d];

  q = q*args->fidim[d] - b;
  if (d == 0)
    { if (b < 0)
        { m  = f-1;
          q += m;
          for (j = m-b; j > m; j--)
            *o++ = a[q - (j % f)];
          q -= m;
        }
      if (out)
        for (j = 0; j < f; j++)
          *o++ = a[q+j];
      else
        { o += f;
          j = f;
        }
      while (j < e)
        *o++ = a[q + ((j++) % f)];
    }
  else
    { d -= 1;
      if (b < 0)
        { m  = f-1;
          q += m;
          for (j = m-b; j > m; j--)
            o = expand_wrap_<t>(d,1,q - (j%f),o,a,args);
          q -= m;
        }
      for (j = 0; j < f; j++)
        o = expand_wrap_<t>(d,out,q+j,o,a,args);
      while (j < e)
        o = expand_wrap_<t>(d,1,q + ((j++) % f),o,a,args);
    }
  return (o);
}

static <t> *expand_reflect_<t>(int d, int out, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, m, k;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->ffdim[d];

  q = q*args->fidim[d] - b;
  m  = ((f-1) << 1);
  if (d == 0)
    { for (j = -b; j > 0; j--)
        { k = j % m;
          if (k >= f)
            k = m-k;
          *o++ = a[q+k];
        }
      if (out)
        while (j < f)
          *o++ = a[q+(j++)];
      else
        { o += f;
          j = f;
        }
      while (j < e)
        { k = (j++) % m;
          if (k >= f)
            k = m-k;
          *o++ = a[q+k];
        }
    }
  else
    { d -= 1;
      for (j = -b; j > 0; j--)
        { k = j % m;
          if (k >= f)
            k = m-k;
          o = expand_reflect_<t>(d,1,q+k,o,a,args);
        }
      while (j < f)
        o = expand_reflect_<t>(d,out,q+(j++),o,a,args);
      while (j < e)
        { k = (j++) % m;
          if (k >= f)
            k = m-k;
          o = expand_reflect_<t>(d,1,q+k,o,a,args);
        }
    }
  return (o);
}

static <t> *expand_extend_<t>(int d, int out, int64 q, <t> *o, <t> *a, Frame_Args *args)
{ int64 j;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->ffdim[d];

  q = q*args->fidim[d] - b;
  if (d == 0)
    { for (j = b; j < 0; j++)
        *o++ = a[q];
      if (out)
        while (j < f)
          *o++ = a[q+(j++)];
      else
        { o += f;
          j = f;
        }
      q += f-1;
      while (j++ < e)
        *o++ = a[q];
    }
  else
    { d -= 1;
      for (j = b; j < 0; j++)
        o = expand_extend_<t>(d,1,q,o,a,args);
      while (j < f)
        o = expand_extend_<t>(d,out,q+(j++),o,a,args);
      q += f-1;
      while (j++ < e)
        o = expand_extend_<t>(d,1,q,o,a,args);
    }
  return (o);
}

static void expand_invert_<t>(int d, int out, int64 q, int64 t, <t> *o, <t> *a, Frame_Args *args)
{ int64 j, k, m, x;
  int64 e = args->fecrd[d];
  int64 b = args->fbcrd[d];
  int64 f = args->ffdim[d];

  q  = q*args->fidim[d] - b;
  t  = t*args->fidim[d];
  m  = ((f-1) << 1);

  if (d == 0)
    { <u>   w, ab, ae;

      a += q;
      o += t;
      ab = 2*a[0];
      ae = 2*a[f-1];

      j = b;
      if (b < 0)
        { j  = -b;
          x  = (j-1)/(f-1)+1;
          if (x % 2)
            w = ((<u>) (x/2))*(ab - ae) + ab;
          else
            w = ((<u>) (x/2))*(ab - ae);
          x = -min_int64(e,0);
          for (k = (j-1) % m + 1; j > x; j--)
            if (k >= f)
              { *o++ = (<t>) (w + a[m-k]);
                if (k-- == f)
                  w += ae;
              }
            else
              { *o++ = (<t>) (w - a[k]);
                if (k-- == 1)
                  { w -= ab;
                    k  = m;
                  }
              }
          if (e <= 0)
            return;
        }

      if (out)
        for (x = min_int64(e,f); j < x; j++)
          *o++ = (<t>) a[j];
      else
        { o += f;
          j  = f;
        }

      w = ae;
      for (k = f-1; j < e; j++)
        if (k-- >= f)
          { *o++ = (<t>) (w + a[m-k]);
            if (k < f)
              w += ae;
          }
        else
          { *o++ = (<t>) (w - a[k]);
            if (k == 0)
              { w -= ab;
                k  = m;
              }
          }
    }

  else
    { int64 lx, rx, y, S;
      <t>  *L, *R, *F;

      d -= 1;
      t -= b;

      expand_invert_<t>(d,out,q,t,o,a,args);
      expand_invert_<t>(d,out,q+(f-1),t+(f-1),o,a,args);

      S = args->fdvol[d];
      L = o + t*S;
      R = o + (t+(f-1))*S;

      j = b;
      if (b < 0)
        { j  = -b;
          x  = (j-1)/(f-1)+1;
          if (x % 2)
            { lx = x+1;
              rx = -(x-1);
            }
          else
            { lx = x;
              rx = -x;
            }
          x = -min_int64(e,0);
          F = L + b*S;
          for (k = (j-1) % m + 1; j > x; j--)
            if (k >= f)
              { expand_invert_<t>(d,1,q + (m-k),t - j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) + ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] + F[y]);
#END
                F += S;

                if (k-- == f)
                  rx += 2;
              }
            else
              { expand_invert_<t>(d,1,q + k,t - j,o,a,args);

                for (y = 0; y < S; y++)
#IF T == UINT64
                  F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) - ((int64) F[y]));
#ELSE
                  F[y] = (<t>) (lx*L[y] + rx*R[y] - F[y]);
#END
                F += S;

                if (k-- == 1)
                  { lx -= 2;
                    k   = m;
                  }
              }
          if (e <= 0)
            return;
        }

      for (j = 1; j < f-1; j++)
        expand_invert_<t>(d,out,q + j,t + j,o,a,args);

      rx = 2;
      lx = 0;
      F  = R + S;
      k  = f-1;
      for (j = f; j < e; j++)
        if (k-- >= f)
          { expand_invert_<t>(d,1,q + (m-k),t + j,o,a,args);

            for (y = 0; y < S; y++)
#IF T == UINT64
              F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) + ((int64) F[y]));
#ELSE
              F[y] = (<t>) (lx*L[y] + rx*R[y] + F[y]);
#END
            F += S;

            if (k < f)
              rx += 2;
          }
        else
          { expand_invert_<t>(d,1,q + k,t + j,o,a,args);

            for (y = 0; y < S; y++)
#IF T == UINT64
              F[y] = (<t>) (lx*((int64) L[y]) + rx*((int64) R[y]) - ((int64) F[y]));
#ELSE
              F[y] = (<t>) (lx*L[y] + rx*R[y] - F[y]);
#END
            F += S;

            if (k == 0)
              { lx -= 2;
                k   = m;
              }
          }
    }
}

#END


/****************************************************************************************
 *                                                                                      *
 *  FRAME CORE ROUTINES                                                                 *
 *                                                                                      *
 ****************************************************************************************/

Frame *G(Make_Frame)(Array *I(target), Coordinate *S(shape), Coordinate *S(anchor))
{ Framer    *frame;
  Dimn_Type *scrd = ADIMN(shape);
  int32     *ccrd = AINT32(anchor);
  int        i, ndims;
  Size_Type  size, vsize;

  if (shape->ndims != 1 || shape->type != DIMN_TYPE)
    { fprintf(stderr,"Shape is not a coordinate vector (Make_Frame)\n");
      exit (1);
    }
  if (anchor->ndims != 1 || anchor->type != DIMN_TYPE)
    { fprintf(stderr,"Center is not an integer vector (Make_Frame)\n");
      exit (1);
    }
  if (Array_Refcount(shape) != 1)
    { fprintf(stderr,"Shape is not subsumable, has a reference count > 1 (Make_Frame)\n");
      exit (1);
    }
  if (Array_Refcount(anchor) != 1)
    { fprintf(stderr,"Center is not subsumable, has a reference count > 1 (Make_Frame)\n");
      exit (1);
    }

  ndims = (int) shape->size;
  if (target->ndims != ndims && (target->ndims-1 != ndims || target->kind == PLAIN_KIND))
    { fprintf(stderr,"Target and coordinate dimensionality do not match (Make_Frame)\n");
      exit (1);
    }
  if (anchor->size != ndims)
    { fprintf(stderr,"Shape and anchor coordinate dimensionality do not match (Make_Frame)\n");
      exit (1);
    }

  if (ndims == target->ndims-1)
    { ndims += 1;
      if (target->kind == COMPLEX_KIND)
        { PrependCoord(shape,2);
          PrependCoord(anchor,0);
        }
      else
        { AppendCoord(kind_size[target->kind],shape);
          AppendCoord(0,anchor);
        }
      scrd = ADIMN(shape);
      ccrd = AINT32(anchor);
    }

  size  = 1;
  vsize = 0;
  for (i = 0; i < ndims; i++)
    { vsize += 2*size;
      size  *= scrd[i];
    }
  vsize += size-2;

  frame = new_framer(size*SIZEOF(Offs_Type)+vsize*SIZEOF(double),"Make_Frame");

  frame->kind    = FRAME_KIND;
  frame->trg_ref = Inc_Array(target);
  frame->shape   = shape;
  frame->anchor  = anchor;

  frame->ndims = ndims;
  frame->size  = size;
  frame->vsize = vsize;

  { Coordinate *low = Copy_Array(shape);
    Coordinate *hgh = Copy_Array(shape);
    Dimn_Type *hcrd = ADIMN(hgh);
    Dimn_Type *lcrd = ADIMN(low);

    frame->notempty = 1;
    for (i = 0; i < ndims; i++)
      { if (ccrd[i] < 0)
          { if (scrd[i] + ccrd[i] > target->dims[i])
              frame->notempty = 0;
          }
        else if (ccrd[i] < scrd[i])
          { if (scrd[i] > target->dims[i])
              frame->notempty = 0;
          }
        else // ccrd[i] >= scrd[i]
          { if (ccrd[i] >= target->dims[i])
              frame->notempty = 0;
          }
        if (frame->notempty)
          { if (ccrd[i] < 0)
              lcrd[i] = 0;
            else
              lcrd[i] = ccrd[i];
            if (ccrd[i] >= scrd[i])
              hcrd[i] = target->dims[i]-1;
            else
              hcrd[i] = target->dims[i] + ccrd[i] - scrd[i];
          }
        else
          lcrd[i] = hcrd[i] = 0;
        if ( ! frame->notempty)
          lcrd[i] = hcrd[i] = 0;
      }

    frame->slice = Make_Slice(target,low,hgh);
  }

  frame->vals = (void *) (frame->offs + size);

  { Frame_Args args;

    args.fctr  = ccrd;
    args.ffdim = scrd;
    args.fidim = target->dims;
    frame_offsets(ndims-1,0,frame->offs,&args);
  }

  frame->mode = NOT_SET;
  Set_Slice_To_Index(frame->slice,0);

  return ((Frame *) frame);
}

#define FRAME(f)            ((Framer *) f)
#define FRAME_SLICE(f)      (FRAME(f)->slice)
#define FRAME_INDEX(f)      (((Slicer *) FRAME_SLICE(f))->p)
#define FRAME_COORDINATE(f) (((Slicer *) FRAME_SLICE(f))->cnt)
#define FRAME_WITHIN(f)     (((Slicer *) FRAME_SLICE(f))->clip < 0 && FRAME(f)->notempty)

Coordinate *Frame_Shape(Frame *frame)
{ return (FRAME(frame)->shape); }

Coordinate *Frame_Anchor(Frame *frame)
{ return (FRAME(frame)->anchor); }

Indx_Type Frame_Index(Frame *frame)
{ return (FRAME_INDEX(frame));  }

Coordinate *Frame_Coordinate(Frame *frame)
{ return (FRAME_COORDINATE(frame));  }

boolean Place_Frame(Frame *M(frame), Indx_Type p)
{ FRAME(frame)->mode = NOT_SET;
  return (Set_Slice_To_Index(FRAME_SLICE(frame),p) && FRAME(frame)->notempty);
}

boolean Move_Frame_Forward(Frame *M(frame))
{ FRAME(frame)->mode = NOT_SET;
  return (Inc_Slice_Index(FRAME_SLICE(frame)) && FRAME(frame)->notempty);
}

boolean Move_Frame_Backward(Frame *M(frame))
{ FRAME(frame)->mode = NOT_SET;
  return (Dec_Slice_Index(FRAME_SLICE(frame)) && FRAME(frame)->notempty);
}

boolean Frame_Within_Array(Frame *frame)
{ return (FRAME_WITHIN(frame)); }

static void setup_value_compute(Framer *frame)
{ int ndims = frame->ndims;

  if (FRAME_WITHIN(frame))
    { Size_Type  n = frame->size;
      Offs_Type *o = frame->offs;
      Indx_Type j;

      switch (frame->trg_ref->type) {
        #GENERATE T = @TYPES
          case <T>_TYPE:
            { <t> *v = A<T>(frame->trg_ref) + FRAME_INDEX(frame);
              <t> *d = (<t> *) frame->vals;
              for (j = 0; j < n; j++)
                d[j] = v[o[j]];
              break;
            }
        #END
      }
    }

  else
    { Dimn_Type *pos = ADIMN(FRAME_COORDINATE(frame));
      int32     *low = AINT32(frame->anchor);
      Dimn_Type *hgh = ADIMN(frame->shape);

      int64 Fbcrd[10], *bcrd;
      int64 Fecrd[10], *ecrd;
      int64 Fdvol[10], *dvol;
      void *Fibuf[10], **ibuf;

      Frame_Args args;

      if (ndims > 10)
        { bcrd = (int64 *) Guarded_Malloc(sizeof(int64)*4*((size_t) ndims),"Frame_Offsets");
          ecrd = bcrd + ndims;
          dvol = ecrd + ndims;
          ibuf = (void **) (dvol + ndims);
        }
      else
        { bcrd = Fbcrd;
          ecrd = Fecrd;
          dvol = Fdvol;
          ibuf = Fibuf;
        }

      { int64 i, b;

        for (i = 0; i < ndims; i++)
          { bcrd[i] = b = ((int64) pos[i]) - low[i];
            ecrd[i] = b + hgh[i];
          }
      }

      args.fbcrd = bcrd;
      args.fecrd = ecrd;
      args.fidim = frame->trg_ref->dims;
      args.ffdim = hgh;

      switch (frame->trg_ref->type) {
        #GENERATE T = @TYPES
          case <T>_TYPE:
            { <t> *a = A<T>(frame->trg_ref);
              <t> *o = (<t> *) frame->vals;
              switch (Boundary_Case_8qm5)
              { case BND_ZERO:
                  frame_zero_<t>(ndims-1,0,o,a,&args);
                  break;
                case BND_REFLECT:
                  frame_reflect_<t>(ndims-1,0,o,a,&args);
                  break;
                case BND_WRAP:
                  frame_wrap_<t>(ndims-1,0,o,a,&args);
                  break;
                case BND_EXTEND:
                  frame_extend_<t>(ndims-1,0,o,a,&args);
                  break;
                case BND_INVERT:
                  { int64 b, v;
                    int   i;

                    b = 1;
                    v = frame->size;
                    for (i = 0; i < ndims; i++)
                      { ibuf[i] = o + v;
                        b *= hgh[i];
                        v += 2*b;
                        dvol[i] = b;
                      }
                    args.fibuf = ibuf;
                    args.fdvol = dvol;
                    frame_invert_<t>(ndims-1,0,0,o,a,&args);
                    break;
                  }
              }
              break;
            }
        #END
      }

      if (ndims > 10)
        free(bcrd);
    }

  frame->mode = USE_VALS;
}

void *Frame_Values(Frame *framer)
{ Framer *frame = (Framer *) framer;
  if (frame->mode != USE_VALS)
    setup_value_compute(frame);
  return ((void *) (frame->vals));
}

Offs_Type *Frame_Offsets(Frame *f)
{ return (FRAME(f)->offs); }

Array_Bundle *Frame_Array(Array_Bundle *R(O(a)), Frame *framer)
{ static char text[1] = { '\0' };
  Framer *frame = (Framer *) framer;
  a->type   = frame->trg_ref->type;
  a->ndims  = frame->trg_ref->ndims;
  a->dims   = ADIMN(frame->shape);
  a->size   = frame->size;
  a->tlen   = 0;
  a->text   = text;
  a->data   = Frame_Values(framer);
  a->scale  = frame->trg_ref->scale;
  a->kind   = AForm_Kind(frame);
  return (a);
}

Array *G(Make_Array_From_Frame)(Frame *framer)
{ Framer *frame = (Framer *) framer;
  Array  *a     = frame->trg_ref;
  
  Array  *target = Make_Array(PLAIN_KIND,a->type,a->ndims,ADIMN(frame->shape));

  target->kind = AForm_Kind(frame);

  { Size_Type n = target->size;
    Indx_Type j;

    if (FRAME_WITHIN(frame))
      { Offs_Type *o = frame->offs;
        switch (a->type) {
          #GENERATE T = @TYPES
            case <T>_TYPE:
              { <t> *d = A<T>(target);
                <t> *v = A<T>(frame->trg_ref) + FRAME_INDEX(frame);
                for (j = 0; j < n; j++)
                  d[j] = v[o[j]];
                break;
              }
          #END
        }
      }
    else
      { if (frame->mode == NOT_SET)
          setup_value_compute(frame);
        switch (a->type) {
          #GENERATE T,U = @TYPES , @UNION
            case <T>_TYPE:
              { <t> *d = A<T>(target);
                <t> *v = (<t> *) frame->vals;
                for (j = 0; j < n; j++)
                  d[j] = v[j];
                break;
              }
          #END
        }
      }
  }

  return (target);
}


/****************************************************************************************
 *                                                                                      *
 *  ARRAY FORM BASICS                                                                   *
 *                                                                                      *
 ****************************************************************************************/

static Form_Class _AForm_Class[] = { ARRAY_CLASS, ARRAY_CLASS, ARRAY_CLASS,
                                     ARRAY_CLASS, SLICE_CLASS, FRAME_CLASS };

Form_Class AForm_Class(AForm *form)
{ return (_AForm_Class[((Array *) form)->kind]); }

boolean Is_Slice(AForm *o)
{ return (((Array *) o)->kind == SLICE_KIND); }

boolean Is_Frame(AForm *o)
{ return (((Array *) o)->kind == FRAME_KIND); }

boolean Is_Array(AForm *o)
{ return (((Array *) o)->kind <= COMPLEX_KIND); }

Array *AForm_Array(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return (((Slicer *) form)->trg_ref);
    case FRAME_CLASS:
      return (((Framer *) form)->trg_ref);
    case ARRAY_CLASS:
      return ((Array *) form);
  }
  return (NULL);
}

Size_Type AForm_Size(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return (((Slicer *) form)->size);
    case FRAME_CLASS:
      return (((Framer *) form)->size);
    case ARRAY_CLASS:
    default:
      return (((Array *) form)->size);
  }
}

Array_Kind AForm_Kind(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      { Slicer    *slice = (Slicer *) form;
        int        ndims = slice->ndims;
        Array_Kind kind  = slice->trg_ref->kind;
        if (kind == RGB_KIND || kind == RGBA_KIND)
          { if (slice->bcrd[ndims-1] != 0 || slice->ecrd[ndims-1] != kind_size[kind]-1)
              kind = PLAIN_KIND;
          }
        else if (kind == COMPLEX_KIND)
          { if (slice->bcrd[0] != 0 || slice->ecrd[0] != 1)
              kind = PLAIN_KIND;
          }
        return (kind);
      }
    case FRAME_CLASS:
      { Framer    *frame = (Framer *) form;
        int        ndims = frame->ndims;
        Array_Kind kind  = frame->trg_ref->kind;
        if (kind == RGB_KIND || kind == RGBA_KIND)
          { if (ADIMN(frame->shape)[ndims-1] != kind_size[kind])
              kind = PLAIN_KIND;
          }
        else if (kind == COMPLEX_KIND)
          { if (ADIMN(frame->shape)[0] != 2)
              kind = PLAIN_KIND;
          }
        return (kind);
      }
    case ARRAY_CLASS:
      return (((Array *) form)->kind);
  }
  return (PLAIN_KIND);
}

Coordinate *G(AForm_Shape)(AForm *form)
{ Coordinate *coord;

  switch (AForm_Class(form))
  { case SLICE_CLASS:
      { Slicer     *slice = (Slicer *) form;
        Dimn_Type  *bcrd, *bval;
        int         i;

        coord = Copy_Array(slice->end);
        bval  = ADIMN(coord);
        bcrd  = slice->bcrd;
        for (i = 0; i < slice->ndims; i++)
          bval[i] = (bval[i] - bcrd[i]) + 1;
        return (coord);
      }
    case FRAME_CLASS:
      return (Copy_Array(Frame_Shape((Frame *) form)));
    case ARRAY_CLASS:
      { Array *a = (Array *) form;
        coord = Make_Array_With_Shape(PLAIN_KIND,DIMN_TYPE,Coord1(a->ndims));
        memcpy(ADIMN(coord),a->dims,sizeof(Dimn_Type)*((size_t) a->ndims));
        return (coord);
      }
  }
  return (NULL);
}

Array *G(Make_Array_From_AForm)(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return (Make_Array_From_Slice((Slice *) form));
    case FRAME_CLASS:
      return (Make_Array_From_Frame((Frame *) form));
    case ARRAY_CLASS:
      return (Copy_Array((Array *) form));
  }
  return (NULL);
}

AForm *G(Copy_AForm)(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return ((AForm *) Copy_Slice((Slice *) form));
    case FRAME_CLASS:
      return ((AForm *) Copy_Frame((Frame *) form));
    case ARRAY_CLASS:
      return ((AForm *) Copy_Array((Array *) form));
  }
  return (NULL);
}

AForm *Pack_AForm(AForm *R(M(form)))
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return ((AForm *) Pack_Slice((Slice *) form));
    case FRAME_CLASS:
      return ((AForm *) Pack_Frame((Frame *) form));
    case ARRAY_CLASS:
      return ((AForm *) Pack_Array((Array *) form));
  }
  return (NULL);
}

AForm *Inc_AForm(AForm *R(I(form)))
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return ((AForm *) Inc_Slice((Slice *) form));
    case FRAME_CLASS:
      return ((AForm *) Inc_Frame((Frame *) form));
    case ARRAY_CLASS:
      return ((AForm *) Inc_Array((Array *) form));
  }
  return (NULL);
}

void Free_AForm(AForm *F(form))
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      Free_Slice((Slice *) form);
      return;
    case FRAME_CLASS:
      Free_Frame((Frame *) form);
      return;
    case ARRAY_CLASS:
      Free_Array((Array *) form);
      return;
  }
}

void Kill_AForm(AForm *K(form))
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      Kill_Slice((Slice *) form);
      return;
    case FRAME_CLASS:
      Kill_Frame((Frame *) form);
      return;
    case ARRAY_CLASS:
      Kill_Array((Array *) form);
      return;
  }
}

void Reset_AForm()
{ Reset_Array();
  Reset_Slice();
  Reset_Frame();
}

int AForm_Usage()
{ return (Array_Usage()+Slice_Usage()+Frame_Usage()); }

void AForm_List(void (*handler)(AForm *))
{ void (*ahandler)(Array *) = (void (*)(Array *)) handler;
  Array_List(ahandler);
  Slice_List(handler);
  Frame_List(handler);
}

int AForm_Refcount(AForm *form)
{ switch (AForm_Class(form))
  { case SLICE_CLASS:
      return (Slice_Refcount((Slice *) form));
    case FRAME_CLASS:
      return (Frame_Refcount((Frame *) form));
    case ARRAY_CLASS:
      return (Array_Refcount((Array *) form));
  }
  return (0);
}

/****************************************************************************************
 *                                                                                      *
 *  DISPLAY AN ARRAY                                                                    *
 *                                                                                      *
 ****************************************************************************************/

static string type_name[] = { "uint8", "uint16", "uint32", "uint64", "int8", "int16", "int32",
                              "int64", "float32", "float64" };

static int   Uindent;
static FILE *Uoutput;

static void Uhandler(Array *a)
{ int       i;
  Size_Type bytes;

  fprintf(Uoutput,"%*s%2d : %u",Uindent,"",Array_Refcount(a),a->dims[a->ndims-1]);
  for (i = a->ndims-2; i >= 0; i--)
    fprintf(Uoutput," x %u",a->dims[i]);
  fprintf(Uoutput," %s",type_name[a->type]);
  bytes = a->size * type_size[a->type];
  if (bytes < 1.e3)
    fprintf(Uoutput," = %llu",bytes);
  else if (bytes < 1.e6)
    fprintf(Uoutput," = %.1fKb",bytes/1.e3);
  else if (bytes < 1.e9)
    fprintf(Uoutput," = %.1fMb",bytes/1.e6);
  else 
    fprintf(Uoutput," = %.2fGb",bytes/1.e9);
  if (a->text != NULL && a->text[0] != '\0')
    fprintf(Uoutput," : '%.*s'",50,a->text);
  fprintf(Uoutput,"\n");
}

static void Shandler(Slice *s)
{ int        i;
  Dimn_Type *scrd;
  Array     *a = AForm_Array(s);

  scrd = ADIMN(Slice_First(s));
  fprintf(Uoutput,"%*s%2d : (%u",Uindent,"",Slice_Refcount(a),scrd[a->ndims-1]);
  for (i = a->ndims-2; i >= 0; i--)
    fprintf(Uoutput,",%u",scrd[i]);
  scrd = ADIMN(Slice_Last(s));
  fprintf(Uoutput,") - (%u",scrd[a->ndims-1]);
  for (i = a->ndims-2; i >= 0; i--)
    fprintf(Uoutput,",%u",scrd[i]);
  fprintf(Uoutput,") %s",type_name[a->type]);
  if (a->text != NULL && a->text[0] != '\0')
    fprintf(Uoutput," @ '%.*s'",50,a->text);
  fprintf(Uoutput,"\n");
}

static void Fhandler(Frame *f)
{ int        i;
  Size_Type  bytes;
  Dimn_Type *scrd;
  Array     *a = AForm_Array(f);

  scrd = ADIMN(Frame_Shape(f));
  fprintf(Uoutput,"%*s%2d : %u",Uindent,"",Frame_Refcount(f),scrd[a->ndims-1]);
  for (i = a->ndims-2; i >= 0; i--)
    fprintf(Uoutput," x %u",scrd[i]);
  fprintf(Uoutput," %s",type_name[a->type]);
  bytes = AForm_Size(f) * 24;
  if (bytes < 1.e3)
    fprintf(Uoutput," = %llu",bytes);
  else if (bytes < 1.e6)
    fprintf(Uoutput," = %.1fKb",bytes/1.e3);
  else if (bytes < 1.e9)
    fprintf(Uoutput," = %.1fMb",bytes/1.e6);
  else 
    fprintf(Uoutput," = %.2fGb",bytes/1.e9);
  if (a->text != NULL && a->text[0] != '\0')
    fprintf(Uoutput," @ '%.*s'",50,a->text);
  fprintf(Uoutput,"\n");
}

void Print_Inuse_List(FILE *output, int indent)
{ Uoutput = output;
  Uindent = indent + 2;
  fprintf(output,"%*sArrays:\n",indent,"");
  Array_List(Uhandler);
  fprintf(output,"%*sSlices:\n",indent,"");
  Slice_List(Shandler);
  fprintf(output,"%*sFrames:\n",indent,"");
  Frame_List(Fhandler);
}

void Print_Array(AForm *o, FILE *output, int indent, string format)
{ Array      *a     = AForm_Array(o);
  int         ndims = a->ndims;
  int         kind  = AForm_Kind(o);
  Size_Type   e     = AForm_Size(o);
  Coordinate *base  = AForm_Shape(o);
  Dimn_Type  *s     = ADIMN(base);
  Array      *low   = Copy_Array(base);
  int32      *c     = AINT32(low);

  int       d0, d1, d2, od;
  Size_Type x0, x1, n, pw;
  Indx_Type i, k, p;
  int       j;
  boolean   slice;

  if (kind != PLAIN_KIND)
    ndims -= 1;
  d0 = (kind == COMPLEX_KIND);
  d1 = (ndims >= 2) + d0;
  d2 = d1+1;
  od = (ndims-1) + d0;
  e /= kind_size[kind];

  x0 = x1 = s[d0];
  if (ndims > 1)
    x1 *= s[d1];

  slice = 0;
  if (Is_Frame(o))
    { Dimn_Type *p = ADIMN(Frame_Coordinate(o));
      int32     *v = AINT32(Frame_Anchor(o));
      for (i = 0; i < a->ndims; i++)
        c[i] = p[i]-v[i];
      fprintf(output,"\n%*sFrame ",indent,"");
    }
  else if (Is_Slice(o))
    { Dimn_Type *b = ADIMN(Slice_First(o));
      for (i = 0; i < a->ndims; i++)
        c[i] = b[i];
      fprintf(output,"\n%*sSlice ",indent,"");
      slice = 1;
    }
  else
    { for (i = 0; i < a->ndims; i++)
        c[i] = 0;
      fprintf(output,"\n%*sArray ",indent,"");
    }
  if (ndims >= 2)
    fprintf(output,"[%d,%d] x ",c[d1],c[d1]+(s[d1]-1));
  fprintf(output,"[%d,%d]\n",c[d0],c[d0]+(s[d0]-1));

  switch (kind) {
    #GENERATE K = @KINDS
      case <K>:
        switch (a->type) {
          #GENERATE T,U = @TYPES , @UNION
            case <T>_TYPE:
              { <t> *v;
                <t> *w, *x;             #WHEN K == RGB_KIND || K == RGBA_KIND
                <t> *y;                 #WHEN K == RGBA_KIND

                if (Is_Frame(o))
                  v = (<t> *) Frame_Values(o);
                else
                  v = A<T>(a);

                p  = 0;
                pw = e;			#WHEN K == RGB_KIND || K == RGBA_KIND
                if (slice)
                  { p  = Set_Slice_To_First(o);
                    pw = a->size / kind_size[kind];	#WHEN K == RGB_KIND || K == RGBA_KIND
                  }

                w = v + pw;             #WHEN K == RGB_KIND || K == RGBA_KIND
                x = w + pw;             #WHEN K == RGB_KIND || K == RGBA_KIND
                y = x + pw;             #WHEN K == RGBA_KIND

                for (i = 0; i < e; i++)
                  { if (i % x1 == 0)
                      { if (i > 0)
                          { if (s[d1] == 1)
                              fprintf(output," }\n");
                            else
                              fprintf(output,"\n%*s}\n",indent,"");
                          }
                        if (ndims > 2)
                          { fprintf(output,"\n%*s(",indent,"");
                            n = e;
                            k = i;
                            for (j = od; j > d2; j--)
                              { n /= s[j];
                                fprintf(output,"%lld,",c[j]+(k/n));
                                k = (k % n);
                              }
                            n /= s[j];
                            fprintf(output,"%lld)",c[j]+(k/n));
                          }
                        fprintf(output,"\n%*s{ ",indent,"");
                      }
                    else if (i % x0 == 0)
                      fprintf(output,"\n%*s  ",indent,"");
                    else
                      fprintf(output,", ");
                    #IF K == PLAIN_KIND
                      fprintf(output,format,v[p]);
                    #ELSEIF K == RGB_KIND
                      fprintf(output,"[");
                      fprintf(output,format,v[p]);
                      fprintf(output,",");
                      fprintf(output,format,w[p]);
                      fprintf(output,",");
                      fprintf(output,format,x[p]);
                      fprintf(output,"]");
                    #ELSEIF K == RGBA_KIND
                      fprintf(output,"[");
                      fprintf(output,format,v[p]);
                      fprintf(output,",");
                      fprintf(output,format,w[p]);
                      fprintf(output,",");
                      fprintf(output,format,x[p]);
                      fprintf(output,",");
                      fprintf(output,format,y[p]);
                      fprintf(output,"]");
                    #ELSE
                      fprintf(output,format,v[p]);
                      fprintf(output," + ");
                      if (slice)
                        p = Next_Slice_Index(o);
                      else
                        p += 1;
                      fprintf(output,format,v[p]);
                      fprintf(output,"i");
                    #END
                    if (slice)
                      p = Next_Slice_Index(o);
                    else
                      p += 1;
                  }
                break;
              }
          #END
        }
        break;
    #END
  }

  if (s[d1] == 1 || d0 == d1)
    fprintf(output," }\n");
  else
    fprintf(output,"\n%*s}\n",indent,"");

  Free_Array(low);
  Free_Array(base);
}


/****************************************************************************************
 *                                                                                      *
 *  MODIFY TEXT DESCRIPTIONS                                                            *
 *                                                                                      *
 ****************************************************************************************/

void Set_Array_Text(Array *M(a), string text)
{ int len = (int) strlen(text);
  allocate_array_text(a,len+1,"Set_Array_Text");
  a->tlen = len;
  strcpy(a->text,text);
}

void Append_To_Array_Text(Array *M(a), string text)
{ int sen = (int) strlen(a->text);
  allocate_array_text(a,sen+a->tlen+1,"Append_To_Array_Text");
  a->tlen += sen;
  strcpy(a->text+sen,text);
}


/****************************************************************************************
 *                                                                                      *
 *  SUB-PLANE SELECTION                                                                 *
 *                                                                                      *
 ****************************************************************************************/


Array_Bundle *Get_Array_Plane(Array_Bundle *R(M(a)), Dimn_Type plane)
{ Dimn_Type nplanes = a->dims[a->ndims-1];
  Size_Type offset  = array_dsize(a)/nplanes;

  if (plane >= nplanes || a->ndims <= 1 /* || plane < 0 */)
    return (NULL);

  if (a->kind == COMPLEX_KIND && a->ndims > 0)
    a->kind   = COMPLEX_KIND;
  else
    a->kind   = PLAIN_KIND;
  a->ndims  = a->ndims-1;
  a->size   = a->size / nplanes;
  a->data   = ((char *) a->data) + plane*offset;

  return (a);
}


/****************************************************************************************
 *                                                                                      *
 *  COMPUTE RANGES AND SCALE IMAGES                                                     *
 *                                                                                      *
 ****************************************************************************************/

//  Compute min and max values in 'a' of type 'type' with 'length' elements

Range_Bundle *Array_Range(Range_Bundle *R(O(rng)), AForm *o)
{ Array    *a = AForm_Array(o);
  Size_Type n = AForm_Size(o);

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          <t>  x, min, max;

          switch (AForm_Class(o))
          { case SLICE_CLASS:
              min = max = d[Set_Slice_To_First(o)];
              while (n-- > 1)
                { x = d[Next_Slice_Index(o)];
                  if (x < min)
                    min = x;
                  else if (x > max)
                    max = x;
                }
              break;
            case FRAME_CLASS:
              if (Frame_Within_Array(o))
                { Offs_Type *off = Frame_Offsets(o);
                  d += Frame_Index(o);
                  min = max = d[off[0]];
                  while (n-- > 1)
                    { x = d[off[n]];
                      if (x < min)
                        min = x;
                      else if (x > max)
                        max = x;
                    }
                  break;
                }
              else
                d = (<t> *) Frame_Values(o);
            case ARRAY_CLASS:
              min = max = d[0];
              while (n-- > 1)
                { x = d[n];
                  if (x < min)
                    min = x;
                  else if (x > max)
                    max = x;
                }
              break;
            default:
              min = max = 0;
              break;
          }
          rng->maxval.<U> = max;
          rng->minval.<U> = min;
          break;
        }
    #END
  }

  return (rng);
}

APart *Scale_Array(APart *R(M(o)), double factor, double offset)
{ Size_Type n;
  Indx_Type e;
  Array    *a = AForm_Array(o);

  if (Is_Frame(o))
    { fprintf(stderr,"Scale_Array does not operate on frames\n");
      exit (1);
    }

  n = AForm_Size(o);
  switch (a->type) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          if (Is_Slice(o))
            for (e = Set_Slice_To_First(o); n-- > 0; e = Next_Slice_Index(o))
              d[e] = (<t>) (factor * (d[e] + offset));
          else
            while (n-- > 0)
              d[n] = (<t>) (factor * (d[n] + offset));
          break;
        }
    #END
  }

  return (o);
}

APart *Scale_Array_To_Range(APart *R(M(o)), Value min, Value max)
{ Range_Bundle crn;
  double       f;
  Array       *a = AForm_Array(o);

  if (Is_Frame(o))
    { fprintf(stderr,"Scale_Array_To_Range does not operate on frames\n");
      exit (1);
    }

  Array_Range(&crn,o);
  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        if (crn.maxval.<U> == crn.minval.<U>)
          break;
        if (min.<U> == max.<U>)
          break;
        f  = (double) (max.<U>-min.<U>);
        f /= crn.maxval.<U>-crn.minval.<U>;
        Scale_Array(o,f,min.<U>/f-crn.minval.<U>);
        break;
    #END
  }

  return (o);
}


/****************************************************************************************
 *                                                                                      *
 *  ARRAY SCALAR AND EL-BY-EL OPERATORS                                                 *
 *                                                                                      *
 ****************************************************************************************/

APart *Array_Op_Scalar(APart *R(M(o)), Operator op, Value_Type type, Value value)
{ Array    *a = AForm_Array(o);
  AForm    *m = o;
  uint64    uval = 0;
  int64     ival = 0;
  double    fval = 0.;
  Size_Type n;
  Indx_Type i;

  if (op == LSH_OP || op == RSH_OP)
    { if (a->type >= FLOAT32_TYPE)
        { fprintf(stderr,"Shifting floating point values is not permitted (Array_Op_Scalar)\n");
          exit (1);
        }
    }
  if (Is_Frame(m))
    { fprintf(stderr,"Array_Op_Scalar does not operate on frames\n");
      exit (1);
    }

  if (type2kind[type] == UVAL)
    uval = value.uval;
  else if (type2kind[type] == IVAL)
    ival = value.ival;
  else
    fval = value.fval;

  n = AForm_Size(m);
  switch (a->type) {
    #GENERATE T,C = @TYPES , @CASTES
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          <t>  setval;
          switch (type2kind[type]) {
            #GENERATE V = UVAL IVAL FVAL
              case <V>:
                { <c> rp = (<c>) <v>;
                  setval = (<t>) rp;
                  switch (op) {
                    #GENERATE OP,SYM = @OPS , @SYMS
                      case <OP>:
                      #GENERATE S = 0 1
                       #IF S == 0
                        if (Is_Slice(m))
                          for (i = Set_Slice_To_First(m); n-- > 0; i = Next_Slice_Index(m))
                       #ELSE
                        else
                          for (i = 0; i < n; i++)
                       #END
                          #IF OP == SET_OP
                            d[i] = setval;
                          #ELSEIF OP <= DIV_OP
                            d[i] = (<t>) (d[i] <SYM> rp);
                          #ELSEIF OP <= POW_OP
                            d[i] = (<t>) pow((double) d[i],(double) <v>);
                          #ELSEIF OP <= RSH_OP
                           #IF T < FLOAT32
                            d[i] = (<t>) (d[i] <SYM> ((int) <v>));
                           #ELSE
                            ;
                           #END
                          #ELSE
                            { if (d[i] <SYM> setval) d[i] = setval; }
                          #END
                      #END
                        break;
                    #END
                  }
                  break;
                }
            #END
          }
          break;
        }
    #END
  }

  return (o);
}

APart *Complex_Op_Scalar(APart *R(M(o)), Operator op,
                         Value_Type type, Value rpart, Value ipart)
{ Array    *a = AForm_Array(o);
  AForm    *m = o;
  Indx_Type n;
  Indx_Type i;
  double    mag, ang;
  double    mgr, agr;

  if (AForm_Kind(o) != COMPLEX_KIND)
    { fprintf(stderr,"Array form must be of COMPLEX kind (Complex_Op_Scalar)\n");
      exit (1);
    }
  if (op == LSH_OP || op == RSH_OP)
    { if (a->type >= FLOAT32_TYPE)
        { fprintf(stderr,"Shifting floating point values is not permitted (Array_Op_Scalar)\n");
          exit (1);
        }
    }
  if (Is_Frame(m))
    { fprintf(stderr,"Complex_Op_Scalar does not operate on frames\n");
      exit (1);
    }

  n = AForm_Size(m);
  switch (a->type) {
  #GENERATE T,C = @TYPES , @CASTES
    case <T>_TYPE:
      { <t> *d = A<T>(a);
        <t> *D = d+1;
        <t>  dr, di;
        <t>  setrp, setip;
        switch (type2kind[type]) {
        #GENERATE V,U = UVAL IVAL FVAL , uint64 int64 float64
          case <V>:
            { <c> rp = (<c>) rpart.<v>;
              <c> ip = (<c>) ipart.<v>;
              <c> em = rp*rp + ip*ip;
              setrp = (<t>) rp;
              setip = (<t>) ip;
              switch (op) {
              #GENERATE OP,SYM = @OPS , @SYMS
                case <OP>:
                #GENERATE S = 0 1
                #IF S == 0
                  if (Is_Slice(m))
                    for (i = Set_Slice_To_First(m); n-- > 0; n--, Next_Slice_Index(m),
                                                              i = Next_Slice_Index(m))
                #ELSE
                  else
                    for (i = 0; i < n; i += 2)
                #END
                #IF OP == SET_OP
                      { d[i] = setrp;
                        D[i] = setip;
                #ELSEIF OP <= SUB_OP
                      { d[i] = (<t>) (d[i] <SYM> rp);
                        D[i] = (<t>) (D[i] <SYM> ip);
                #ELSEIF OP <= POW_OP
                      { dr = d[i];
                        di = D[i];
                  #IF OP == MUL_OP
                        d[i] = (<t>) (dr*rp - di*ip);
                        D[i] = (<t>) (di*rp + dr*ip);
                  #ELSEIF OP == DIV_OP
                        d[i] = (<t>) ((dr*rp + di*ip) / em);
                        D[i] = (<t>) ((di*rp - dr*ip) / em);
                  #ELSE
                        mag  = (double) (di*di + dr*dr);
                        ang  = atan2((double) di,(double) dr);
                        mgr = pow(mag,.5*rp) * exp((-ang)*ip);
                        agr = .5*ip*log(mag) + rp*ang;
                        d[i] = (<t>) (mgr * cos(agr));
                        D[i] = (<t>) (mgr * sin(agr));
                  #END
                #ELSEIF OP <= RSH_OP
                  #IF T < FLOAT32
                      { d[i] = (<t>) (d[i] <SYM> ((int) rp));
                        D[i] = (<t>) (D[i] <SYM> ((int) ip));
                  #ELSE
                      { ;
                  #END
                #ELSE
                      { if (d[i] <SYM> setrp) d[i] = setrp;
                        if (D[i] <SYM> setip) D[i] = setip;
                #END
                      }
                #END
                  break;
              #END
                }
              break;
            }
        #END
          }
        break;
      }
  #END
    }

  return (o);
}

#define AOA_BUFLEN  0x40000

static pthread_mutex_t Buffer_Mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  Buffer_Cond  = PTHREAD_COND_INITIALIZER;

static int    AA_NoInit = 1;
static int    AA_BufAlloc[NUM_THREADS];
static double AA_Buffer[AOA_BUFLEN*NUM_THREADS];

static double *Get_AOA_Buffer()
{ int i;

  pthread_mutex_lock(&Buffer_Mutex);
  if (AA_NoInit)
    { AA_NoInit = 0;
      for (i = 0; i < NUM_THREADS; i++)
        AA_BufAlloc[i] = 0;
    }
  while (1)
    { for (i = 0; i < NUM_THREADS; i++)
        if (AA_BufAlloc[i] == 0)
          { AA_BufAlloc[i] = 1;
            break;
          }
      if (i < NUM_THREADS)
        break;
      pthread_cond_wait(&Buffer_Cond,&Buffer_Mutex);
    }
  pthread_mutex_unlock(&Buffer_Mutex);

  return (AA_Buffer + AOA_BUFLEN*i);
}

static void Release_AOA_Buffer(void *buf)
{ pthread_mutex_lock(&Buffer_Mutex);
  AA_BufAlloc[(((double *) buf)-AA_Buffer)/AOA_BUFLEN] = 0;
  pthread_cond_signal(&Buffer_Cond);
  pthread_mutex_unlock(&Buffer_Mutex);
}

#GENERATE R,D = UVAL IVAL FVAL , uint64 int64 float64

  static Indx_Type load_buffer_<r>(<D> *buffer, AForm *q,
                                   Size_Type crea, Indx_Type j, Size_Type brea)
  { Array    *b = AForm_Array(q);
    Indx_Type i;

    switch (b->type) {
    #GENERATE T,Q = @TYPES , 1 1 1 1 2 2 2 2 3 3
     #IF Q == R
      case <T>_TYPE:
        { <t> *e = A<T>(b);
          switch (AForm_Class(q))
          { case SLICE_CLASS:
              for (i = 0; i < crea; i++)
                { buffer[i] = e[j];
                  j = Next_Slice_Index(q);
                }
              break;
            case FRAME_CLASS:
              if (Frame_Within_Array(q))
                { Offs_Type *off = Frame_Offsets(q);
                  e += Frame_Index(q);
                  for (i = 0; i < crea; i++)
                    { buffer[i] = e[off[j++]];
                      if (j == brea) j = 0;
                    }
                }
              else
                { <t> *v = Frame_Values(q);
                  for (i = 0; i < crea; i++)
                    { buffer[i] = v[j++];
                      if (j == brea) j = 0;
                    }
                }
              break;
            case ARRAY_CLASS:
              for (i = 0; i < crea; i++)
                { buffer[i] = e[j++];
                  if (j == brea) j = 0;
                }
              break;
          }
          break;
        }
     #END
    #END
    default:
      break;
    }
    return (j);
  }
#END

APart *Array_Op_Array(APart *R(M(o)), Operator op, AForm *q)
{ Array *a = AForm_Array(o);
  Array *b = AForm_Array(q);
  AForm *m = o;

  boolean     alice;
  Size_Type   brea, crea;
  Indx_Type   i, j, k;
  Size_Type   volume;
  int         btype;

  if ((op == LSH_OP || op == RSH_OP) && a->type >= FLOAT32_TYPE)
    { fprintf(stderr,"Shifting floating point values is not permitted (Array_Op_Array)\n");
      exit (1);
    }
  if (Is_Frame(m))
    { fprintf(stderr,"Array_Op_Array does not accept a frame as it first argument\n");
      exit (1);
    }

  alice = Is_Slice(m);

  if (op == SET_OP && a->type == b->type && !alice && Is_Array(q))
    { Size_Type as = array_size(a);
      Size_Type bs = array_size(b);
      if (bs >= as)
        memcpy(a->data,b->data,(size_t) as);
      else
        { for (i = 0; i < as; i += bs)
            if (i+bs > as)
              memcpy(((char *) a->data) + i,b->data,(size_t) (as-i));
            else
              memcpy(((char *) a->data) + i,b->data,(size_t) bs);
        }
      return (o);
    }

  if (b->type <= UINT64_TYPE)
    btype = UVAL;
  else if (b->type <= INT64_TYPE)
    btype = IVAL;
  else
    btype = FVAL;

  if (Is_Slice(q))
    j = Set_Slice_To_First(q);
  else
    j = 0;

  if (alice)
    k = Set_Slice_To_First(m);
  else
    k = 0;

  volume = AForm_Size(m);
  brea   = AForm_Size(q);

  switch (btype) {
  #GENERATE R,D = UVAL IVAL FVAL , uint64 int64 float64
    case <R>:
      { <D> *buffer = (<D> *) Get_AOA_Buffer();
        while (volume != 0)
          { if (volume < AOA_BUFLEN)
              crea = volume;
            else
              crea = AOA_BUFLEN;

            j = load_buffer_<r>(buffer,q,crea,j,brea);

            switch (a->type) {
            #GENERATE T = @TYPES
              case <T>_TYPE:
                { <t> *d = A<T>(a);
                  switch (op) {
                  #GENERATE OP,SYM = @OPS , @SYMS
                    case <OP>:
                    #GENERATE S = 0 1
                    #IF S == 0
                      if (alice)
                        for (i = 0; i < crea; i++, k = Next_Slice_Index(m))
                    #ELSE
                      else
                        for (i = 0; i < crea; i++, k++)
                    #END
                      #IF OP == SET_OP
                          d[k] = (<t>) buffer[i];
                      #ELSEIF OP <= POW_OP
                        #IF T <= UINT32 && R == UVAL
                          { uint64 dr = (uint64) d[k];
                            uint64 er = (uint64) buffer[i];
                        #ELSEIF T >= FLOAT32 || R == FVAL
                          { float64 dr = (float64) d[k];
                            float64 er = (float64) buffer[i];
                        #ELSE
                          { int64 dr = (int64) d[k];
                            int64 er = (int64) buffer[i];
                        #END
                        #IF OP < POW_OP
                            d[k] = (<t>) (dr <SYM> er);
                          }
                        #ELSE
                            d[k] = (<t>) pow((double) dr,(double) er);
                          }
                        #END
                      #ELSEIF OP <= RSH_OP
                        #IF T < FLOAT32
                          d[k] = (<t>) (d[k] <SYM> (int) buffer[i]);
                        #ELSE
                          ;
                        #END
                      #ELSE
                          { <t> dr = (<t>) buffer[i];
                            if (d[k] <SYM> dr) d[k] = dr;
                          }
                      #END
                    #END
                      break;
                  #END
                  }
                  break;
                }
            #END
            }
            volume -= crea;
          }
        Release_AOA_Buffer(buffer);
        break;
      }
  #END
  }

  return (o);
}

APart *Complex_Op_Array(APart *R(M(o)), Operator op, AForm *q)
{ Array *a = AForm_Array(o);
  Array *b = AForm_Array(q);
  AForm *m = o;

  boolean     alice;
  Size_Type   brea, crea;
  Indx_Type   i, j, k;
  Size_Type   volume;
  int         btype;
  double      mag, ang;
  double      mgr, agr;

  if (AForm_Kind(o) != COMPLEX_KIND)
    { fprintf(stderr,"First array form must be complex (Complex_Op_Array)\n");
      exit (1);
    }
  if ((op == LSH_OP || op == RSH_OP) && a->type >= FLOAT32_TYPE)
    { fprintf(stderr,"Shifting floating point values is not permitted (Array_Op_Array)\n");
      exit (1);
    }
  if (Is_Frame(m))
    { fprintf(stderr,"Complex_Op_Array does not accept a frame as it first argument\n");
      exit (1);
    }

  if (b->type <= UINT64_TYPE)
    btype = UVAL;
  else if (b->type <= INT64_TYPE)
    btype = IVAL;
  else
    btype = FVAL;

  if (Is_Slice(q))
    j = Set_Slice_To_First(q);
  else
    j = 0;

  alice = Is_Slice(m);
  if (alice)
    k = Set_Slice_To_First(m);
  else
    k = 0;

  volume = AForm_Size(m);
  brea   = AForm_Size(q);

  switch (btype) {
  #GENERATE R,D = UVAL IVAL FVAL , uint64 int64 float64
    case <R>:
      { <D> *buffer = (<D> *) Get_AOA_Buffer();
        while (volume != 0)
          { if (volume < AOA_BUFLEN)
              crea = volume/2;
            else
              crea = AOA_BUFLEN/2;

            j = load_buffer_<r>(buffer,q,crea,j,brea);

            switch (a->type) {
            #GENERATE T = @TYPES
              case <T>_TYPE:
                { <t> *d = A<T>(a);
                  <t> *D = d+1;
                  switch (op) {
                  #GENERATE OP,SYM = @OPS , @SYMS
                    case <OP>:
                    #GENERATE S = 0 1
                    #IF S == 0
                      if (alice)
                        for (i = 0; i < crea; i++, Next_Slice_Index(m), k = Next_Slice_Index(m))
                    #ELSE
                      else
                        for (i = 0; i < crea; i++, k += 2)
                    #END
                      #IF OP == SET_OP
                            d[k] = D[k] = (<t>) buffer[i];
                      #ELSEIF OP <= POW_OP
                        #IF T <= UINT32 && R == UVAL
                          { uint64 dr = (uint64) d[k];
                            uint64 di = (uint64) D[k];
                            uint64 er = (uint64) buffer[i];
                        #ELSEIF T >= FLOAT32 || R == FVAL
                          { float64 dr = (float64) d[k];
                            float64 di = (float64) D[k];
                            float64 er = (float64) buffer[i];
                        #ELSE
                          { int64 dr = (int64) d[k];
                            int64 di = (int64) D[k];
                            int64 er = (int64) buffer[i];
                        #END
                        #IF OP <= DIV_OP
                            d[k] = (<t>) (dr <SYM> er);
                            D[k] = (<t>) (di <SYM> er);
                          }
                        #ELSE
                            mag  = (double) (di*di + dr*dr);
                            ang  = atan2((double) di,(double) dr);
                            mgr  = pow(mag,.5*er);
                            agr  = er*ang;
                            d[k] = (<t>) (mgr * cos(agr));
                            D[k] = (<t>) (mgr * sin(agr));
                          }
                        #END
                      #ELSEIF OP <= RSH_OP
                        #IF T < FLOAT32
                          { int er = (int) buffer[i];
                            d[k] = (<t>) (d[k] <SYM> er);
                            D[k] = (<t>) (D[k] <SYM> er);
                          }
                        #ELSE
                            ;
                        #END
                      #ELSE
                          { <t> dr = (<t>) buffer[i];
                            if (d[k] <SYM> dr) d[k] = dr;
                            if (D[k] <SYM> dr) D[k] = dr;
                          }
                        #END
                    #END
                      break;
                  #END
                    }
                  break;
                }
            #END
              }
            volume -= 2*crea;
          }
        Release_AOA_Buffer(buffer);
        break;
      }
  #END
    }

  return (o);
}

APart *Complex_Op_Complex(APart *R(M(o)), Operator op, AForm *q)
{ Array *a = AForm_Array(o);
  Array *b = AForm_Array(q);
  AForm *m = o;

  boolean     alice;
  Size_Type   brea, crea;
  Indx_Type   i, j, k;
  Size_Type   volume;
  int         btype;
  double      mag, ang;
  double      mgr, agr;

  if (AForm_Kind(o) != COMPLEX_KIND || AForm_Kind(q) != COMPLEX_KIND)
    { fprintf(stderr,"Array forms must be complex (Complex_Op_Complex)\n");
      exit (1);
    }
  if ((op == LSH_OP || op == RSH_OP) && a->type >= FLOAT32_TYPE)
    { fprintf(stderr,"Shifting floating point values is not permitted (Array_Op_Array)\n");
      exit (1);
    }
  if (Is_Frame(m))
    { fprintf(stderr,"Complex_Op_Complex does not accept a frame as it first argument\n");
      exit (1);
    }

  if (b->type <= UINT64_TYPE)
    btype = UVAL;
  else if (b->type <= INT64_TYPE)
    btype = IVAL;
  else
    btype = FVAL;

  if (Is_Slice(q))
    j = Set_Slice_To_First(q);
  else
    j = 0;

  alice = Is_Slice(m);
  if (alice)
    k = Set_Slice_To_First(m);
  else
    k = 0;

  volume = AForm_Size(m);
  brea   = AForm_Size(q);

  switch (btype) {
  #GENERATE R,D = UVAL IVAL FVAL , uint64 int64 float64
    case <R>:
      { <D> *buffer = (<D> *) Get_AOA_Buffer();
        <D> *c = buffer;
        <D> *C = c+1;
        while (volume != 0)
          { if (volume < AOA_BUFLEN)
              crea = volume;
            else
              crea = AOA_BUFLEN;

            j = load_buffer_<r>(buffer,q,crea,j,brea);

            switch (a->type) {
            #GENERATE T,C = @TYPES , @CASTES
              case <T>_TYPE:
                { <t> *d = A<T>(a);
                  <t> *D = d+1;
                  switch (op) {
                  #GENERATE OP,SYM = @OPS , @SYMS
                    case <OP>:
                    #GENERATE S = 0 1
                    #IF S == 0
                      if (alice)
                        for (i = 0; i < crea; i += 2, Next_Slice_Index(m), k = Next_Slice_Index(m))
                    #ELSE
                      else
                        for (i = 0; i < crea; i += 2, k += 2)
                    #END
                      #IF OP == SET_OP
                          { d[k] = (<t>) c[i];
                            D[k] = (<t>) C[i];
                      #ELSEIF OP <= POW_OP
                        #IF T <= UINT32 && R == UVAL
                          { uint64 dr = (uint64) d[k];
                            uint64 di = (uint64) D[k];
                            uint64 er = (uint64) c[i];
                            uint64 ei = (uint64) C[i];
                            uint64 em;				#WHEN OP == DIV_OP
                        #ELSEIF T >= FLOAT32 || R == FVAL
                          { float64 dr = (float64) d[k];
                            float64 di = (float64) D[k];
                            float64 er = (float64) c[i];
                            float64 ei = (float64) C[i];
                            float64 em;				#WHEN OP == DIV_OP
                        #ELSE
                          { int64 dr = (int64) d[k];
                            int64 di = (int64) D[k];
                            int64 er = (int64) c[i];
                            int64 ei = (int64) C[i];
                            int64 em;				#WHEN OP == DIV_OP
                        #END
                        #IF OP <= SUB_OP
                            d[k] = (<t>) (dr <SYM> er);
                            D[k] = (<t>) (di <SYM> ei);
                        #ELSEIF OP == MUL_OP
                            d[k] = (<t>) (dr*er - di*ei);
                            D[k] = (<t>) (di*er + dr*ei);
                        #ELSEIF OP == DIV_OP
                            em   = er*er + ei*ei;
                            d[k] = (<t>) ((dr*er + di*ei) / em);
                            D[k] = (<t>) ((di*er - dr*ei) / em);
                        #ELSE
                            mag  = (double) (di*di + dr*dr);
                            ang  = atan2((double) di,(double) dr);
                            mgr = pow(mag,.5*er) * exp((-ang)*ei);
                            agr = .5*ei*log(mag) + er*ang;
                            d[k] = (<t>) (mgr * cos(agr));
                            D[k] = (<t>) (mgr * sin(agr));
                        #END
                      #ELSEIF OP <= RSH_OP
                        #IF T < FLOAT32
                          { d[k] = (<t>) (d[k] <SYM> ((int) c[i]));
                            D[k] = (<t>) (D[k] <SYM> ((int) C[i]));
                        #ELSE
                          { ;
                        #END
                      #ELSE
                          { <t> dr = (<t>) c[i];
                            <t> di = (<t>) C[i];
                            if (d[k] <SYM> dr) d[k] = dr;
                            if (D[k] <SYM> di) D[k] = di;
                      #END
                          }
                    #END
                      break;
                  #END
                    }
                  break;
                }
            #END
              }
            volume -= crea;
          }
        Release_AOA_Buffer(buffer);
        break;
      }
  #END
    }

  return (o);
}

/****************************************************************************************
 *                                                                                      *
 *  APPLYING FUNCTIONS TO ARRAY ELEMENTS                                                *
 *                                                                                      *
 ****************************************************************************************/

APart *Array_Fct_Val(APart *R(M(o)), Value (*fct)(void *valp))
{ Indx_Type i, e;
  Slice    *s = (Slice *) o;
  Array    *a = AForm_Array(o);

  if (Is_Frame(o))
    { fprintf(stderr,"Array_Fct_Val does not operate on frames\n");
      exit (1);
    }

  switch (a->type) {
    #GENERATE T,C,U = @TYPES , @CASTES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          if (Is_Slice(s))
            { e = Set_Slice_To_Last(s); 
              for (i = Set_Slice_To_First(s); 1; i = Next_Slice_Index(s)) 
                { d[i] = (<t>) fct((void *) (d+i)).<U>;
                  if (i == e) break;
                }
            }
          else
            for (i = 0; i < a->size; i++)
              d[i] = (<t>) fct((void *) (d+i)).<U>;
          break;
        }
    #END
  }

  return (o);
}

APart *Array_Fct_Idx(APart *R(M(o)), Value (*fct)(Coordinate *coord))
{ Indx_Type p, e;
  Slice    *s = (Slice *) o;
  Array    *a = (Array *) AForm_Array(o);

  if (Is_Frame(o))
    { fprintf(stderr,"Array_Fct_Idx does not operate on frames\n");
      exit (1);
    }

  if (!Is_Slice(o))
    s = Make_Slice(a,Idx2CoordA(a,0),Idx2CoordA(a,a->size-1));

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          e = Set_Slice_To_Last(s); 
          for (p = Set_Slice_To_First(s); 1; p = Next_Slice_Index(s)) 
            { d[p] = (<t>) fct(Slice_Coordinate(s)).<U>;
              if (p == e) break;
            }
          break;
        }
    #END
  }

  if (s != o)
    Kill_Slice(s);

  return (o);
}

//  Threshold values less than cutoff to black, all others to white

APart *Threshold_Array(APart *R(M(o)), Value cutoff)
{ Indx_Type i, e;
  Slice    *s = (Slice *) o;
  Array    *a = AForm_Array(o);

  if (Is_Frame(o))
    { fprintf(stderr,"Threshold_Array does not operate on frames\n");
      exit (1);
    }

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          <t>  x, c = (<t>) cutoff.<U>;
          <t>  zero = 0, one = 1;

          if (Is_Slice(s))
            { e = Set_Slice_To_Last(s); 
              for (i = Set_Slice_To_First(s); 1; i = Next_Slice_Index(s)) 
                { x = d[i];
                  if (x <= c)
                    d[i] = zero;
                  else
                    d[i] = one;
                  if (i == e) break;
                }
            }
          else
            { for (i = 0; i < a->size; i++)
                { x = d[i];
                  if (x <= c)
                    d[i] = zero;
                  else
                    d[i] = one;
                }
            }
          break;
        }
    #END
  }

  return (o);
}


/****************************************************************************************
 *                                                                                      *
 *  CONVERT IMAGE TYPES                                                                 *
 *                                                                                      *
 ****************************************************************************************/

typedef struct
  { Value_Type type;
    int        scale;
    void      *data;
  } Domain;

typedef struct
  { int    tshift;
    double dscale;
  } Variable_Scale;

#define LIKE_2_LIKE  0   // conversion codes
#define CPLX_2_NORM  1
#define  RGB_2_NORM  2
#define RGBA_2_NORM  3
#define RGBA_2_RGB   4

#GENERATE C = CPLX_2_NORM RGB_2_NORM RGBA_2_NORM RGBA_2_RGB
  #GENERATE S = @TYPES

static void translate_<C>_<S>(void *sdata, void *tdata, Size_Type size, int ctype, double afactor)    
{ Indx_Type p;

  <s> *s0 = (<s> *) sdata;
  <s> *s1 = s0 + size;        #WHEN C != CPLX_2_NORM
  <s> *s2 = s1 + size;        #WHEN C != CPLX_2_NORM
  <s> *s3 = s2 + size;        #WHEN C >= RGBA_2_NORM
  
  switch (ctype) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *t0  = (<t> *) tdata;
          <t> *t1  = t0 + size;       #WHEN C == RGBA_2_RGB
          <t> *t2  = t1 + size;       #WHEN C == RGBA_2_RGB

          #IF C == RGBA_2_RGB
            for (p = 0; p < size; p++)
             t0[p] = (<t>) (s0[p] * (afactor * s3[p]));         #WHEN T >= FLOAT32
             t0[p] = (<t>) (s0[p] * (afactor * s3[p]) + .4999); #WHEN T <  FLOAT32
            for (p = 0; p < size; p++)
             t1[p] = (<t>) (s1[p] * (afactor * s3[p]));         #WHEN T >= FLOAT32
             t1[p] = (<t>) (s1[p] * (afactor * s3[p]) + .4999); #WHEN T <  FLOAT32
            for (p = 0; p < size; p++)
             t2[p] = (<t>) (s2[p] * (afactor * s3[p]));         #WHEN T >= FLOAT32
             t2[p] = (<t>) (s2[p] * (afactor * s3[p]) + .4999); #WHEN T <  FLOAT32
          #ELSE
              { Indx_Type q;
                for (q = p = 0; p < size; p++) {
                  #IF C == RGB_2_NORM
                    t0[p] = (<t>) ((.30*s0[p] + .59*s1[p] + .11*s2[p])
                                  * afactor);                       #WHEN T >= FLOAT32
                                  * afactor + .5);                  #WHEN T <  FLOAT32
                  #ELSEIF C == RGBA_2_NORM
                    t0[p] = (<t>) (((.30*s0[p] + .59*s1[p] + .11*s2[p])
                                  * (afactor * s3[p])));            #WHEN T >= FLOAT32
                                  * (afactor * s3[p])) + .5);       #WHEN T <  FLOAT32
                  #ELSE
                    { double x,y;
                      x = (double) s0[q++];
                      y = (double) s0[q++];
                      t0[p] = (<t>) (sqrt(x*x+y*y) * afactor);
                    }
                  #END
                }
              }
          #END
          return;
        }
    #END
  }
  return;
}

  #END
#END

#GENERATE C = CPLX_2_NORM RGB_2_NORM RGBA_2_NORM RGBA_2_RGB

static void translate_<C>(void *sdata, void *tdata, Size_Type size,
                          Value_Type stype, int ctype, double afactor)
{ switch (stype) {
    #GENERATE S = @TYPES
       case <S>_TYPE:
         translate_<C>_<S>(sdata,tdata,size,ctype,afactor);
         return;
    #END
  }
  return;
}

#END

static void translate(int conversion, Size_Type size, Domain *source,
                      Domain *target, Variable_Scale *var)
{ Value_Type stype,  ttype;
  void      *sdata, *tdata;
  double     dscale;
  int        shift;

  stype  = source->type;
  sdata  = source->data;

  ttype  = target->type;
  tdata  = target->data;

  dscale = var->dscale;
  shift  = var->tshift;

  if (conversion)
    { int       ctype;
      double    afactor;

      if (type_size[stype] < type_size[ttype])
        ctype = stype;
      else
        ctype = ttype;

      afactor = 1.;
      if (conversion >= RGBA_2_NORM && stype < FLOAT32_TYPE)
        { uint64 base = (((uint64) 1) << (source->scale-1));
          afactor /= ((base-1) + base);
        }
      if (type_size[stype] >= type_size[ttype])
        { if (shift > 0)
            afactor /= (1 << shift);
          else if (shift < 0)
            afactor *= (1 << -shift);
          else
            afactor *= dscale;
        }

      switch (conversion) {
        #GENERATE C = CPLX_2_NORM RGB_2_NORM RGBA_2_NORM RGBA_2_RGB
          case <C>:
            translate_<C>(sdata,tdata,size,stype,ctype,afactor);
            break;
        #END
      }

      sdata = tdata;
      if (type_size[stype] >= type_size[ttype])
        stype = ttype;
      if (conversion == RGBA_2_RGB)
        size *= 3;
  }

  if (ttype != stype || sdata != tdata)
    { Indx_Type p;

      if (shift > 0 || dscale != 1.)
        switch (stype) {
          #GENERATE S,SS = @TYPES , @SIZES
            case <S>_TYPE:
              { <s> *s = (<s> *) sdata;
                switch (ttype) {
                  #GENERATE T,ST = @TYPES , @SIZES
                    case <T>_TYPE:
                      { <t> *t  = (<t> *) tdata;
                        for (p = 0; p < size; p++)              #WHEN ST <= SS
                        for (p = size; p-- > 0; )               #WHEN ST >  SS
                          t[p] = (<t>) (s[p] >> shift);         #WHEN T <  FLOAT32 && S <  FLOAT32
                          t[p] = (<t>) (s[p] * dscale);         #WHEN T >= FLOAT32 || S >= FLOAT32
                        break;
                      }
                  #END
                }
                break;
              }
          #END
        }

      else if (shift < 0)
        switch (stype) {
          #GENERATE S,SS = @TYPES , @SIZES
            case <S>_TYPE:
              { <s> *s = (<s> *) sdata;
                shift  = - shift;
                switch (ttype) {
                  #GENERATE T,ST = @TYPES , @SIZES
                    case <T>_TYPE:
                      { <t> *t  = (<t> *) tdata;
                        for (p = 0; p < size; p++)              #WHEN ST <= SS
                        for (p = size; p-- > 0; )               #WHEN ST >  SS
                          t[p] = (<t>) (s[p] << shift);         #WHEN T <  FLOAT32 && S <  FLOAT32
                          t[p] = (<t>) (s[p] * dscale);         #WHEN T >= FLOAT32 || S >= FLOAT32
                        break;
                      }
                  #END
                }
                break;
              }
          #END
        }

      else
        switch (stype) {
          #GENERATE S,SS = @TYPES , @SIZES
            case <S>_TYPE:
              { <s> *s = (<s> *) sdata;
                switch (ttype) {
                  #GENERATE T,ST = @TYPES , @SIZES
                    case <T>_TYPE:
                      { <t> *t  = (<t> *) tdata;
                        for (p = 0; p < size; p++)                      #WHEN ST <= SS
                        for (p = size; p-- > 0; )                       #WHEN ST >  SS
                          t[p] = (<t>) s[p];
                        break;
                      }
                  #END
                }
                break;
              }
          #END
        }
    }
}

static void imaginary_fill(Size_Type size, Domain *target)
{ Indx_Type p, q;

  switch (target->type) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *d = A<T>(target);
          q = size;
          for (p = 2*size; p-- > 0; )
            { d[p--] = 0;
              d[p]   = d[--q];
            }
          break;
        }
    #END
  }
}

static void alpha_fill(Size_Type size, Domain *target, Value alpha)
{ Indx_Type p;

  switch (target->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(target) + 3*size;
          for (p = 0; p < size; p++)
            d[p] = (<t>) alpha.<u>;
          break;
        }
    #END
  }
}

static Array *convert_array(Array *sarray, Array_Kind tkind, Value_Type ttype, int scale,
                            Variable_Scale *var, int in_place)
{ Domain     source, target;
  Array     *tarray;

  string     estring;
  int        tdims;
  Size_Type  tels, pels;

  Array_Kind skind  = sarray->kind;
  Value_Type stype  = sarray->type;

  pels   = sarray->size / kind_size[skind];
  tdims  = (sarray->ndims + (tkind != PLAIN_KIND)) - (skind != PLAIN_KIND);
  tels   = pels * kind_size[tkind];

  if (in_place)
    { estring = "Convert_Array_Inplace";
      tarray  = sarray;
    }
  else
    { estring = "Convert_Array";
      tarray  = new_array(0,0,1,estring);
      tarray->tlen    = 0;
      tarray->text[0] = '\0';
    }

  if (stype < FLOAT32_TYPE && ttype < FLOAT32_TYPE)
    { if (var->tshift > bit_size[ttype])
        { fprintf(stderr,"Scale factor shift is larger than number of bits in type (%s)\n",estring);
          exit (1);
        }
    }
  else
    { if (fabs(var->dscale) < 1e-20)
        { fprintf(stderr,"Scale factor is nearly 0, i.e. < 1e-20 (%s)\n",estring);
          exit (1);
        }
    }

  allocate_array_data(tarray,tels*type_size[ttype],estring);
  allocate_array_dims(tarray,SIZEOF(Dimn_Type)*tdims,estring);

  source.type  = stype;
  source.data  = sarray->data;
  source.scale = sarray->scale;

  target.type  = ttype;
  target.data  = tarray->data;

  if (tkind == skind)
    translate(LIKE_2_LIKE,kind_size[tkind]*pels,&source,&target,var);
  else if (tkind == PLAIN_KIND || tkind == COMPLEX_KIND)
    { if (skind == RGB_KIND)
        translate(RGB_2_NORM,pels,&source,&target,var);
      else if (skind == RGBA_KIND)
        translate(RGBA_2_NORM,pels,&source,&target,var);
      else if (skind == COMPLEX_KIND)
        translate(CPLX_2_NORM,pels,&source,&target,var);
      else
        translate(LIKE_2_LIKE,kind_size[tkind]*pels,&source,&target,var);
      if (tkind == COMPLEX_KIND)
        imaginary_fill(pels,&target);
    }
  else //  tkind == RGB_KIND || tkind == RGBA_KIND
    { if (skind == COMPLEX_KIND)
        translate(CPLX_2_NORM,pels,&source,&target,var);
      else if (skind == RGBA_KIND)
        translate(RGBA_2_RGB,pels,&source,&target,var);
      else
        translate(LIKE_2_LIKE,sarray->size,&source,&target,var);
      if (skind == PLAIN_KIND || skind == COMPLEX_KIND)
        { uint64 block = (uint64) (type_size[ttype]*pels);
          memcpy(((char *) tarray->data) + block,tarray->data,(size_t) block);
          memcpy(((char *) tarray->data) + 2*block,tarray->data,(size_t) block);
        }
      if (tkind == RGBA_KIND)
        { if (ttype <= UINT64_TYPE)
            alpha_fill(pels,&target,VALU(0xFFFFFFFFull));
          else if (ttype <= INT64_TYPE)
            alpha_fill(pels,&target,VALI(0x7FFFFFFFll));
          else
            alpha_fill(pels,&target,VALF(1.));
        }
    }

  { int i, cs, ct, nd;

    cs = (skind == COMPLEX_KIND);
    ct = (tkind == COMPLEX_KIND);
    nd = tdims - (tkind != PLAIN_KIND);
    if (ct > cs)
      { for (i = nd-1; i >= 0; i--)
          tarray->dims[i+ct] = sarray->dims[i+cs];
      }
    else
      { for (i = 0; i < nd; i++)
          tarray->dims[i+ct] = sarray->dims[i+cs];
      }
    if (ct)
      tarray->dims[0] = 2;
    if (!ct && nd != tdims)
      tarray->dims[nd] = kind_size[tkind];
  }

  tarray->size  = tels;
  tarray->kind  = tkind;
  tarray->type  = ttype;
  tarray->ndims = tdims;
  tarray->scale = scale;

  return (tarray);
}

Array *Convert_Array_Inplace(Array *R(M(sarray)), Array_Kind tkind, Value_Type ttype,
                             int tscale, ...)
{ Variable_Scale var;
  va_list        ap;

  if (sarray->kind == tkind && sarray->type == ttype)
    return (sarray);

  var.dscale = 1.;
  var.tshift = 0;

  va_start(ap,tscale);
  if (sarray->type < FLOAT32_TYPE && ttype < FLOAT32_TYPE)
    var.tshift = va_arg(ap,int);
  else
    var.dscale = va_arg(ap,double);
  va_end(ap);

  return (convert_array(sarray,tkind,ttype,tscale,&var,1));
}

Array *G(Convert_Array)(Array *sarray, Array_Kind tkind, Value_Type ttype, int tscale, ...)
{ Variable_Scale var;
  va_list        ap;

  if (sarray->kind == tkind && sarray->type == ttype)
    return (Copy_Array(sarray));

  var.dscale = 1.;
  var.tshift = 0;

  va_start(ap,tscale);
  if (sarray->type < FLOAT32_TYPE && ttype < FLOAT32_TYPE)
    var.tshift = va_arg(ap,int);
  else
    var.dscale = va_arg(ap,double);
  va_end(ap);

  return (convert_array(sarray,tkind,ttype,tscale,&var,0));
}

boolean Image_Check(Array *array)
{ Indx_Type p;
  uint64    vmax;
  int64     smax;

  if (array->kind == COMPLEX_KIND)
    { fprintf(stderr,"Array cannot be COMPLEX (Image_Check)\n");
      return (1);
    }

  if (array->type <= UINT64_TYPE)
    { uint64 base = (((uint64) 1) << (array->scale-1));
      vmax = (base-1) + base;
      smax = 0;
    }
  else if (array->type <= INT64_TYPE)
    { vmax = 0;
      smax = (int64) ((((uint64) 1) << (array->scale-1)) - 1);
    }
  else
    { vmax = 0;
      smax = 0;
    }

  switch (array->type) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *v = A<T>(array);
          for (p = 0; p < array->size; p++)
            #IF T >= FLOAT32
              { if (v[p] > 1.)
                  return (1);
            #ELSEIF T >= INT8
              { if (v[p] > smax)
                  return (1);
            #ELSE
              { if (v[p] > vmax)
                  return (1);
            #END
            #IF T >= INT8
                if (v[p] < 0)
                  return (1);
            #END
              }
          break;
        }
    #END
  }

  return (0);
}

Array *convert_image(Array *array, Array_Kind kind, Value_Type type, int scale, int in_place)
{ Variable_Scale var;
  Array         *rez;
  int            atype;

  var.dscale = 1.;
  var.tshift = array->scale - scale;

  atype = array->type;
  if (type >= FLOAT32_TYPE)
    { if (atype <= UINT64_TYPE)
        { uint64 base = (((uint64) 1) << (array->scale-1));
          var.dscale = 1. / ((base-1) + base);
        }
      else if (atype <= INT64_TYPE)
        var.dscale = 1. / ((((uint64) 1) << (array->scale-1)) - 1);
    }
  else if (type >= INT8_TYPE)
    { if (atype >= FLOAT32_TYPE)
        var.dscale = (double) ((((uint64) 1) << (scale-1)) - 1);
      else if (atype <= UINT64_TYPE)
        var.tshift += 1;
    }
  else
    { if (atype >= FLOAT32_TYPE)
        { uint64 base = (((uint64) 1) << (scale-1));
          var.dscale = (double) ((base-1) + base);
        }
      else if (atype >= INT64_TYPE)
        var.tshift -= 1;
    }

  rez = convert_array(array,kind,type,scale,&var,in_place);

  return (rez);
}

Array *Convert_Image_Inplace(Array *R(M(array)), Array_Kind kind, Value_Type type, int scale)
{ return (convert_image(array,kind,type,scale,1)); }

Array *G(Convert_Image)(Array *array, Array_Kind kind, Value_Type type, int scale)
{ return (convert_image(array,kind,type,scale,0)); }

/****************************************************************************************
 *                                                                                      *
 *  ARRAY MULTIPLICATION                                                                *
 *                                                                                      *
 ****************************************************************************************/

Array *G(Array_Multiply)(Array *a, Array *b)
{ int       acmplx  = (a->kind == COMPLEX_KIND);
  int       bcmplx  = (b->kind == COMPLEX_KIND);

  Array    *prod;
  Size_Type b_off, a_off;

  { Dimn_Type *adims  = a->dims + acmplx;
    Dimn_Type *bdims  = b->dims + bcmplx;
    int        andims = a->ndims - acmplx;
    int        bndims = b->ndims - bcmplx;
    int        complex = (acmplx || bcmplx);

    if (a->dims[acmplx] != b->dims[b->ndims-1])
      { fprintf(stderr,"Corresponding inner product dimensions for array multiply don't match");
        fprintf(stderr," (Array Multiply)\n");
        exit (1);
      }

    if (a->type != b->type)
      { fprintf(stderr,"Arrays must be of the same type (Array_Multiply)\n");
        exit (1);
      }

    if (complex)
      { prod = make_start(COMPLEX_KIND,a->type,andims+bndims-1,"Array_Multiply");
        prod->dims[0] = 2;
      }
    else
      prod = make_start(PLAIN_KIND,a->type,andims+bndims-2,"Array_Multiply");

    { Dimn_Type *pdims;
      int        i;

      pdims = prod->dims + complex;
      for (i = 1; i < andims; i++)
        pdims[bndims+i-2] = adims[i];
      for (i = 0; i < bndims-1; i++)
        pdims[i] = bdims[i];
    }

    prod->size = array_size(prod);

    allocate_array_data(prod,array_dsize(prod),"Make_Array");

    b_off  = ((b->size / bdims[bndims-1]) >> bcmplx);
    a_off  = adims[0];
  }

  { Size_Type p_off;
    Indx_Type p, q, r, t, k;

    switch (acmplx*2 + bcmplx) {
      #GENERATE C,S = 0 1 2 3 , 1 2 2 2
        case <C>:
          p_off = b_off;          #WHEN C == 0
          p_off = (b_off << 1);   #WHEN C >  0
          b_off = p_off;          #WHEN C == 1 || C == 3
          a_off <<= 1;            #WHEN C >= 2
          switch (a->type) {
            #GENERATE T = @TYPES
              case <T>_TYPE:
                { <t> *z = A<T>(prod);
                  <t> *y = A<T>(b);
                  <t> *x, sr;
                  <t> si;		#WHEN C != 0

                  for (q = 0; q < p_off; q += <S>)
                    { x = A<T>(a);
                      r = q;
                      for (p = 0; p < a->size; p += a_off)
                        { sr = 0;
                          #IF C == 0
                            for (t = q, k = 0; k < a_off; k++, t += b_off)
                              sr += x[k] * y[t];
                          #ELSE
                            si = 0;
                            #IF C == 3
                              for (k = 0, t = q; k < a_off; k += 2, t += b_off)
                                { sr += x[k] * y[t]   - x[k+1] * y[t+1];
                                  si += x[k] * y[t+1] + x[k+1] * y[t];
                                }
                            #ELSEIF C == 2
                              for (k = 0, t = (q>>1); k < a_off; k += 2, t += b_off)
                                { sr += x[k]   * y[t];
                                  si += x[k+1] * y[t];
                                }
                            #ELSEIF C == 1
                              for (k = 0, t = q; k < a_off; k++, t += b_off)
                                { sr += x[k] * y[t];
                                  si += x[k] * y[t+1];
                                }
                            #END
                            z[r+1] = si;
                          #END
                          z[r] = sr;
                          r += p_off;
                          x += a_off;
                        }
                    }
                  break;
                }
            #END
          }
          break;
      #END
    }
  }

  return (prod);
}

/****************************************************************************************
 *                                                                                      *
 *  COLOR OR OTHER ARRAY MAPPING                                                        *
 *                                                                                      *
 ****************************************************************************************/


Array *G(Apply_Map)(Array *image, Array *map)
{ Array     *rez;
  Size_Type  inner;

  inner = map->dims[0];

  if (map->ndims < 1)
    { fprintf(stderr,"Map should be at least 1 dimensional (Apply_Map)\n");
      exit (1);
    }
  if (image->type > UINT64_TYPE || image->kind != PLAIN_KIND)
    { fprintf(stderr,"image must have unsigned values and be a PLAIN array (Apply_Map)\n");
      exit (1);
    }
  if (( 1ll << image->scale) != inner)
    { fprintf(stderr,"map length doesn't match image scale (Apply_Map)\n");
      exit (1);
    }

  { int i, j;

    rez = make_start(PLAIN_KIND,map->type,image->ndims+map->ndims-1,"Apply_Map");
    rez->kind  = map->kind;
    rez->scale = map->scale;
    for (i = 0; i < image->ndims; i++)
      rez->dims[i] = image->dims[i];
    for (j = 1, i = image->ndims; i < rez->ndims; i++, j++)
      rez->dims[i] = map->dims[j];
    rez->size = array_size(rez);
    allocate_array_data(rez,array_dsize(rez),"Apply_Map");
    Set_Array_Text(rez,image->text);
  }

  { Indx_Type p, q;

    switch (map->type) {
      #GENERATE M = @TYPES
        case <M>_TYPE:
          switch (image->type) {
            #GENERATE T = UINT8 UINT16 UINT32 UINT64
              case <T>_TYPE:
                { <t> *a = A<T>(image);
                  <m> *m = A<M>(map);
                  <m> *t = A<M>(rez);
                  for (q = 0; q < map->size; q += inner)
                    { for (p = 0; p < image->size; p++)
                        *t++ = m[a[p]];
                      m += inner;
                    }
                  break;
                }
            #END
              default:
                break;
          }
          break;
      #END
    }
  }

  return (rez);
}


/****************************************************************************************
 *                                                                                      *
 *  DOWN SAMPLE CODE                                                                    *
 *                                                                                      *
 ****************************************************************************************/

static void down_sample(AForm *source, void *target, Coordinate *F(voxel), Coordinate *F(basis))
{ Size_Type Sdinc[10], *dinc;
  Size_Type Stinc[10], *tinc;
  Dimn_Type Sdcnt[10], *dcnt;
  Dimn_Type Stcnt[10], *tcnt;

  Array     *src     = AForm_Array(source);
  int        isframe = Is_Frame(source);
  int        ndims   = src->ndims;
  Dimn_Type *tdims   = ADIMN(basis);
  Dimn_Type *vdims   = ADIMN(voxel);
  Size_Type  dsize;

  if (ndims > 10)
    { dinc = (Size_Type *) Guarded_Malloc((sizeof(Size_Type)+sizeof(Dimn_Type))*2*((size_t) ndims),
                                           "Down_Sample");
      tinc = dinc + ndims;
      dcnt = ((Dimn_Type *) tinc) + ndims;
      tcnt = dcnt + ndims;
    }
  else
    { dinc = Sdinc;
      tinc = Stinc;
      dcnt = Sdcnt;
      tcnt = Stcnt;
    }

  { Size_Type  offset, outer, partial;
    Dimn_Type *sdims;
    int        i;

    if (isframe)
      sdims = ADIMN(Frame_Shape(source));
    else
      sdims = src->dims;

    dcnt[0] = 0;
    tcnt[0] = 0;
    tinc[0] = vdims[0];
    partial = (tdims[0]-1) * vdims[0];
    dinc[0] = 1;
    offset  = vdims[0] - 1;
    outer   = sdims[0];
    dsize   = vdims[0];
    for (i = 1; i < ndims; i++)
      { dcnt[i] = 0;
        tcnt[i] = 0;
        tinc[i] = vdims[i] * outer - partial;
        partial += (tdims[i]-1) * vdims[i] * outer;
        dinc[i] = outer - offset;
        offset += (vdims[i]-1)*outer;
        outer  *= sdims[i];
        dsize  *= vdims[i];
      }
  }

  { Indx_Type  q, p, r;
    Dimn_Type  d0 = 0;
    Dimn_Type  c0 = vdims[0];
    int        offsetable;
    Offs_Type *off;
    int        i;

    if (Is_Slice(source))
      q = Coord2IdxA(src,Inc_Array(Slice_First(source)));
    else
      q = 0;
    offsetable = 0;
    if (isframe)
      { offsetable = Frame_Within_Array(source);
        if (offsetable)
          off = Frame_Offsets(source);
      }

  #GENERATE S = 0 1 2
    if ( ! isframe)        #WHEN S == 0
    else if (offsetable)   #WHEN S == 1
    else                   #WHEN S == 2
      switch (src->type) {
        #GENERATE T,C,U = @TYPES , @CASTES , @UNION
          case <T>_TYPE:
            { <t> *d = A<T>(src);                           #WHEN S == 0
            { <t> *d = A<T>(src) + Frame_Index(source);     #WHEN S == 1
            { <t> *d = (<t> *) Frame_Values(source);        #WHEN S == 2
              <t> *a = (<t> *) target;
              <C> sum;
              <C> den = (<C>) dsize;
  
              r = 0;
              i = 0;
              while (i != ndims)
                { p   = q;
                  sum = 0;
                  while (i != ndims)
                    { sum += d[p];                        #WHEN S != 1
                    { sum += d[off[p]];                   #WHEN S == 1
                      if (++d0 < c0)
                        p += 1;
                      else
                        { d0 = 0;
                          for (i = 1; i < ndims; i++)
                            if (++dcnt[i] == vdims[i])
                              dcnt[i] = 0;
                            else
                              { p += dinc[i];
                                break;
                              }
                        }
                    }
            
                  a[r++] = (<t>) (sum/den);
            
                  for (i = 0; i < ndims; i++)
                    if (++tcnt[i] == tdims[i])
                      tcnt[i] = 0;
                    else
                      { q += tinc[i];
                        break;
                      }
                }
              break;
            }
        #END
        }
  #END
  }

  if (ndims > 10)
    free(dinc);

  Free_Array(voxel);
  Free_Array(basis);
}
 
static void check_dsample_args(AForm *source, Coordinate *M(point), string routine)
{ Dimn_Type  *crd;
  int         i, ndims;
  Array      *a = AForm_Array(source);
 
  ndims = (int) point->size;
  if (a->ndims != ndims && (a->ndims != ndims+1 || AForm_Kind(source) == PLAIN_KIND))
    { fprintf(stderr,"Array form dimensionaliy and down sample vector don't match");
      fprintf(stderr," (%s)\n",routine);
      exit (1);
    }

  crd = ADIMN(point);
  for (i = 0; i < a->ndims; i++)
    if (crd[i] == 0)
      { fprintf(stderr,"Down sample vector has a 0-component (%s)\n",routine);
        exit (1);
      }

  if (ndims == a->ndims-1)
    { if (AForm_Kind(source) == COMPLEX_KIND)
        PrependCoord(point,1);
      else
        AppendCoord(1,point);
    }
}

Array *Down_Sample_Inplace(Array *R(M(source)), Coordinate *F(ipoint))
{ Coordinate *basis, *point;
  Dimn_Type  *dims, *crd;
  int         i;

  point = Copy_Array(ipoint);
  Free_Array(ipoint);

  check_dsample_args(source,point,"Down_Sample_Inplace");

  basis = AForm_Shape(source);
  dims  = ADIMN(basis);
  crd   = ADIMN(point);
  for (i = 0; i < source->ndims; i++)
    dims[i] /= crd[i];

  down_sample(source,source->data,point,Inc_Array(basis));

  source->size = 1;
  for (i = 0; i < source->ndims; i++)
    { source->dims[i] = dims[i];
      source->size   *= dims[i];
    }

  Free_Array(basis);
  
  return (source);
}
    

Array *G(Down_Sample)(AForm *source, Coordinate *F(ipoint))
{ Coordinate *basis, *point;
  Array      *target;
  Dimn_Type  *dims, *crd;
  int         i;
  Array      *a = AForm_Array(source);

  point = Copy_Array(ipoint);
  Free_Array(ipoint);
 
  check_dsample_args(source,point,"Down_Sample");

  basis = AForm_Shape(source);
  dims  = ADIMN(basis);
  crd   = ADIMN(point);
  for (i = 0; i < a->ndims; i++)
    dims[i] /= crd[i];

  target = Make_Array(PLAIN_KIND,a->type,a->ndims,dims);
  target->kind = a->kind;

  down_sample(source,target->data,point,basis);

  return (target);
}


/****************************************************************************************
 *                                                                                      *
 *  ARRAY CLIPPING CODE                                                                 *
 *                                                                                      *
 ****************************************************************************************/

void clip_array(Array *source, void *target, Coordinate *S(beg), Coordinate *S(end))
{ Slice    *slice;
  Indx_Type p, r, e;

  slice = Make_Slice(source,beg,end);

  switch (source->type) {
    #GENERATE T = @TYPES
      case <T>_TYPE:
        { <t> *d = A<T>(source);
          <t> *a = (<t> *) target;

          e = Set_Slice_To_Last(slice);
          r = 0;
          for (p = Set_Slice_To_First(slice); 1; p = Next_Slice_Index(slice))
            { a[r++] = d[p];
              if (p == e) break;
            }
          break;
        }
    #END
  }

  Kill_Slice(slice);
}

void clip_frame(Frame *frame, void *target, Coordinate *S(beg), Coordinate *S(end))
{ _Array       _bundle;
  Array_Bundle *bundle = &_bundle.array;
  Slice        *slice;
  Indx_Type     p, r, e;
  Array        *a = AForm_Array(frame);

  slice = Make_Slice(Inc_Array(Frame_Array(bundle,frame)),beg,end);

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *d = A<T>(a);
          <t> *a = (<t> *) target;

          if (Frame_Within_Array(frame))
            { Offs_Type *off = Frame_Offsets(frame);
              d += Frame_Index(frame);
              e = Set_Slice_To_Last(slice);
              r = 0;
              for (p = Set_Slice_To_First(slice); 1; p = Next_Slice_Index(slice))
                { a[r++] = d[off[p]];
                  if (p == e) break;
                }
            }
          else
            { <t> *v = (<t> *) Frame_Values(frame);
              e = Set_Slice_To_Last(slice);
              r = 0;
              for (p = Set_Slice_To_First(slice); 1; p = Next_Slice_Index(slice))
                { a[r++] = v[p];
                  if (p == e) break;
                }
            }
          break;
        }
    #END
  }

  Kill_Slice(slice);
}

void check_clip_args(AForm *source, Coordinate *beg, Coordinate *end, string routine)
{ Coordinate *basis = AForm_Shape(source);
  Dimn_Type  *dims  = ADIMN(basis);
  Dimn_Type  *bcrd  = ADIMN(beg);
  Dimn_Type  *ecrd  = ADIMN(end);
  int         ndims = (int) end->size;
  int         i;

  if (ndims != beg->size)
    { fprintf(stderr,"Coordinates are not of same dimensionality (%s)\n",routine);
      exit (1);
    }
  if (basis->size != ndims && (basis->size != ndims+1 || AForm_Kind(source) == PLAIN_KIND))
    { fprintf(stderr,"Source and coordinate dimensionality do not match (%s)\n",routine);
      exit (1);
    }
  
  if (ndims == basis->size-1 && AForm_Kind(source) == COMPLEX_KIND)
    dims += 1;
  for (i = 0; i < ndims; i++)
    { if (bcrd[i] > ecrd[i])
        { fprintf(stderr,"beg coordinate is not before end coordinate (%s)\n",routine);
          exit (1);
        }
      if (ecrd[i] >= dims[i])
        { fprintf(stderr,"end is not in basis of target array (%s)\n",routine);
          exit (1);
        }
    }

  Free_Array(basis);
}
 
Array *Clip_Array_Inplace(Array *R(M(source)), Coordinate *F(beg), Coordinate *F(end))
{ Coordinate *shape;
  Dimn_Type  *bcrd, *scrd, *dims, s;
  int         i, ndims;

  check_clip_args(source,beg,end,"Clip_Array_Inplace");

  ndims = (int) end->size;
  shape = Copy_Array(end);
  bcrd  = ADIMN(beg);
  scrd  = ADIMN(shape);
  for (i = 0; i < ndims; i++)
    scrd[i] -= (bcrd[i]-1);

  clip_array(source,source->data,Copy_Array(beg),Copy_Array(end));
  Free_Array(beg);
  Free_Array(end);

  dims = source->dims;
  if (source->ndims != ndims)
    { source->size = kind_size[source->kind];
      if (source->kind == COMPLEX_KIND)
        dims += 1;
    }
  else
    source->size = 1;
  for (i = 0; i < ndims; i++)
    { dims[i] = s = scrd[i];
      source->size *= s;
    }

  Free_Array(shape);

  return (source);
}
    

Array *G(Clip_Array)(AForm *source, Coordinate *F(beg), Coordinate *F(end))
{ Coordinate *shape;
  Dimn_Type  *bcrd, *scrd;
  Array      *target;
  int         ndims, i;
  boolean     islice;
  Array      *a = AForm_Array(source);

  check_clip_args(source,beg,end,"Clip_Array");

  ndims = (int) end->size;

  islice = Is_Slice(source);
  if (islice)
    { Coordinate *gpnt = Copy_Array(beg);
      Coordinate *hpnt = Copy_Array(end);
      Coordinate *fpnt = Slice_First(source);
      Dimn_Type  *gcrd = ADIMN(gpnt);
      Dimn_Type  *hcrd = ADIMN(hpnt);
      Dimn_Type  *fcrd = ADIMN(fpnt);

      if (ndims != fpnt->size && AForm_Kind(source) == COMPLEX_KIND)
        fcrd += 1;
      for (i = 0; i < ndims; i++)
        { gcrd[i] += fcrd[i];
          hcrd[i] += fcrd[i];
        }
    
      Free_Array(beg);
      Free_Array(end);
      beg = gpnt;
      end = hpnt;
    }

  shape = Copy_Array(end);
  bcrd  = ADIMN(beg);
  scrd  = ADIMN(shape);
  for (i = 0; i < ndims; i++)
    scrd[i] -= (bcrd[i]-1);

  if (a->ndims == ndims)
    { target = Make_Array_With_Shape(PLAIN_KIND,a->type,shape);
      target->kind = AForm_Kind(source);
    }
  else
    target = Make_Array_With_Shape(AForm_Kind(source),a->type,shape);

  if (Is_Frame(source))
    clip_frame(source,target->data,Copy_Array(beg),Copy_Array(end));
  else if (islice)
    clip_array(a,target->data,beg,end);
  else
    clip_array(a,target->data,Copy_Array(beg),Copy_Array(end));

  if (!islice)
    { Free_Array(beg);
      Free_Array(end);
    }

  return (target);
}


/****************************************************************************************
 *                                                                                      *
 *  ARRAY EXPANDING CODE                                                                *
 *                                                                                      *
 ****************************************************************************************/

void expand_array(Array *target, AForm *source, Size_Type size,
                  Coordinate *S(beg), Coordinate *S(end))
{ Slice     *slice;
  Indx_Type  p, q, x;
  int        in;
  Array     *base = AForm_Array(source);

  slice = Make_Slice(target,beg,end);
  in    = Set_Slice_To_Index(slice,target->size-1);

#GENERATE S = 0 1 2 3
 #IF S == 0
  if (Is_Slice(source))
    { x = Slice_Index(source);
      q = Set_Slice_To_Last(source);
 #ELSEIF S == 1
  else if (Is_Frame(source))
   if (Frame_Within_Array(source))
    { Offs_Type *off = Frame_Offsets(source) + size;
 #ELSEIF S == 2
   else
    {
 #ELSE
  else
    {
 #END
      switch (target->type) {
        #GENERATE T, U = @TYPES , @UNION
          case <T>_TYPE:
            { <t> *a = A<T>(target);
              <t> *d = A<T>(base);                           #WHEN S != 2
              <t> *d = (<t> *) Frame_Values(source);         #WHEN S == 2

              d += Frame_Index(source);                      #WHEN S == 1
              d += size;                                     #WHEN S >= 2
              for (p = target->size; p-- > 0; )
                { if (in)
                    { a[p] = d[q];                           #WHEN S == 0
                      q = Prev_Slice_Index(source);          #WHEN S == 0
                    }                                        #WHEN S == 0
                    a[p] = d[*--off];                        #WHEN S == 1
                    a[p] = *--d;                             #WHEN S >= 2
                  else
                    a[p] = 0;
                  in = Dec_Slice_Index(slice);
                }
              break;
            }
        #END
      }
      Set_Slice_To_Index(source,x);                          #WHEN S == 0
    }
#END

  if (Boundary_Case_8qm5 != BND_ZERO)
    { Coordinate *basis = AForm_Shape(source);
      int         ndims = target->ndims;
      Dimn_Type  *begc  = ADIMN(beg);
      Dimn_Type  *endc  = ADIMN(end);

      int64 Fbcrd[10], *bcrd;
      int64 Fecrd[10], *ecrd;
      int64 Fdvol[10], *dvol;

      Frame_Args args;

      if (ndims > 10)
        { bcrd = (int64 *) Guarded_Malloc(sizeof(int64)*3*((size_t) ndims),"Frame_Offsets");
          ecrd = bcrd + ndims;
          dvol = ecrd + ndims;
        }
      else
        { bcrd = Fbcrd;
          ecrd = Fecrd;
          dvol = Fdvol;
        }

      { int i;
        for (i = 0; i < ndims; i++)
          { bcrd[i] = - ((int64) begc[i]);
            ecrd[i] = target->dims[i] - ((int64) begc[i]);
          }
      }

      args.fbcrd = bcrd;
      args.fecrd = ecrd;
      args.fidim = target->dims;
      args.ffdim = ADIMN(basis);

      if (source == target)
        { int i;
          for (i = 0; i < ndims; i++)
            args.ffdim[i] = (endc[i] - begc[i]) + 1;
        }

      switch (target->type) {
        #GENERATE T = @TYPES
          case <T>_TYPE:
            { <t> *a = A<T>(target);
              <t> *o = a;
              switch (Boundary_Case_8qm5)
              { case BND_REFLECT:
                  expand_reflect_<t>(ndims-1,0,0,o,a,&args);
                  break;
                case BND_WRAP:
                  expand_wrap_<t>(ndims-1,0,0,o,a,&args);
                  break;
                case BND_EXTEND:
                  expand_extend_<t>(ndims-1,0,0,o,a,&args);
                  break;
                case BND_INVERT:
                  { int64 i, b;
                    b = 1;
                    for (i = 0; i < ndims; i++)
                      { b *= target->dims[i];
                        dvol[i] = b;
                      }
                    args.fdvol = dvol;
                    expand_invert_<t>(ndims-1,0,0,0,o,a,&args);
                    break;
                  }
                default:
                  break;
              }
              break;
            }
        #END
      }

      if (ndims > 10)
        free(bcrd);
      Free_Array(basis);
    }

  Kill_Slice(slice);
}

void check_expand_args(AForm *source, Coordinate *anchor, Coordinate *shape, string routine)
{ Coordinate *basis = AForm_Shape(source);
  Dimn_Type  *dims  = ADIMN(basis);
  Dimn_Type  *acrd  = ADIMN(anchor);
  Dimn_Type  *scrd  = ADIMN(shape);
  int         ndims = (int) anchor->size;
  int         i;

  if (anchor->size != shape->size || basis->size != shape->size)
    { fprintf(stderr,"Coordinates are not of same dimensionality (%s)\n",routine);
      fprintf(stderr," (Pad_Array_Inplace)\n");
      exit (1);
    }
  if (basis->size != ndims && (basis->size != ndims+1 || AForm_Kind(source) == PLAIN_KIND))
    { fprintf(stderr,"Source and coordinate dimensionality do not match (%s)\n",routine);
      exit (1);
    }

  if (ndims == basis->size-1 && AForm_Kind(source) == COMPLEX_KIND)
    dims += 1;
  for (i = 0; i < ndims; i++)
    { if (acrd[i] >= scrd[i])
        { fprintf(stderr,"anchor coordinate is not inside shape (%s)\n",routine);
          exit (1);
        }
      if (acrd[i] + dims[i] > scrd[i])
        { fprintf(stderr,"source does not fit into shape when placed at anchor (%s)\n",routine);
          exit (1);
        }
    }

  Free_Array(basis);
}
 
Array *Pad_Array_Inplace(Array *R(M(source)), Coordinate *F(anchor), Coordinate *F(shape))
{ Coordinate *end;
  Dimn_Type  *ecrd, *acrd, *dims;
  Size_Type   size;
  int         i, ndims;

  check_expand_args(source,anchor,shape,"Pad_Array_Inplace");
  ndims = (int) anchor->size;

  end  = Copy_Array(shape);
  Free_Array(shape);

  size = source->size;
  ecrd = ADIMN(end);
  acrd = ADIMN(anchor);
  dims = source->dims;
  if (source->ndims != ndims)
    { source->size = kind_size[source->kind];
      if (source->kind == COMPLEX_KIND)
        dims += 1;
    }
  else
    source->size = 1;
  for (i = 0; i < ndims; i++)
    { int m = ecrd[i];
      ecrd[i] = acrd[i] + (dims[i]-1);
      dims[i] = m;
      source->size *= m;
    }

  allocate_array_data(source,array_dsize(source),"Expland_Array_In_Place");

  expand_array(source,source,size,Copy_Array(anchor),end);
  Free_Array(anchor);
  
  return (source);
}


Array *G(Pad_Array)(AForm *source, Coordinate *F(anchor), Coordinate *F(shape))
{ Array      *target;
  Coordinate *end, *basis;
  Dimn_Type  *ecrd, *acrd, *dims;
  int         i, ndims;
  Array      *a = AForm_Array(source);
 
  check_expand_args(source,anchor,shape,"Pad_Array");
  ndims = (int) anchor->size;

  end  = Copy_Array(shape);
  Free_Array(shape);

  ecrd = ADIMN(end);
  acrd = ADIMN(anchor);

  if (a->ndims == ndims)
    { target = Make_Array(PLAIN_KIND,a->type,a->ndims,ecrd);
      target->kind = AForm_Kind(source);
    }
  else
    target = Make_Array(AForm_Kind(source),a->type,ndims,ecrd);

  basis = AForm_Shape(source);
  dims  = ADIMN(basis);
  if (a->ndims != ndims && AForm_Kind(source) == COMPLEX_KIND)
    dims += 1;
  for (i = 0; i < ndims; i++)
    ecrd[i] = acrd[i] + (dims[i]-1);
  Free_Array(basis);

  expand_array(target,source,AForm_Size(source),Copy_Array(anchor),end);
  Free_Array(anchor);

  return (target);
}


/****************************************************************************************
 *                                                                                      *
 *  CORRELATION AND REGRESSION ROUTINES                                                 *
 *                                                                                      *
 ****************************************************************************************/

static double inner_product_arrays(AForm *o1, AForm *o2)
{ Array      *a2 = AForm_Array(o2);
  Value_Type  t1 = AForm_Array(o1)->type;

  Size_Type n, v;
  Indx_Type p, q, i;
  double    s;

  boolean    in2  = 0;
  Offs_Type *off2 = NULL;
  Indx_Type  ctr2 = 0;
  double    *buf;

  if (Is_Slice(o1))
    p = Set_Slice_To_First(o1);
  else
    p = 0;

  q = 0;
  if (Is_Frame(o2))
    { in2 = Frame_Within_Array(o2);
      if (in2)
        { off2 = Frame_Offsets(o2);
          ctr2 = Frame_Index(o2);
        }
    }
  else if (Is_Slice(o2))
    q = Set_Slice_To_First(o2);

  buf = Get_AOA_Buffer();

  n = AForm_Size(o1);
  s = 0;
  for (n = AForm_Size(o1); n > 0; n -= v)
    { if (n < AOA_BUFLEN)
        v = n;
      else
        v = AOA_BUFLEN;

      if (t1 <= UINT32_TYPE)
        { uint64 *ubuf = (uint64 *) buf;
          p = load_buffer_uval(ubuf,o1,v,p,n);
          for (i = 0; i < v; i++)
            buf[i] = (double) ubuf[i];
        }
      else if (t1 <= INT32_TYPE)
        { int64 *ibuf = (int64 *) buf;
          p = load_buffer_ival(ibuf,o1,v,p,n);
          for (i = 0; i < v; i++)
            buf[i] = (double) ibuf[i];
        }
      else
        p = load_buffer_fval(buf,o1,v,p,n);

      switch (a2->type) {
        #GENERATE T,U = @TYPES , @UNION
          case <T>_TYPE:
            { <t> *d = A<T>(a2);
              switch (AForm_Class(o2))
              { case SLICE_CLASS:
                  for (i = 0; i < v; i++)
                    { s += buf[i] * d[q];
                      q = Next_Slice_Index(o2);
                    }
                  break;
                case FRAME_CLASS:
                  if (in2)
                    { d += ctr2;
                      for (i = 0; i < v; i++)
                        s += buf[i] * d[off2[q++]];
                      break;
                    }
                  else
                    d = (<t> *) Frame_Values(o2);
                case ARRAY_CLASS:
                  for (i = 0; i < v; i++)
                    s += buf[i] * d[q++];
                  break;
              }
              break;
            }
        #END
       }
    }

  Release_AOA_Buffer(buf);

  return (s);
}

static double sum_array(AForm *o)
{ Indx_Type p;
  double    t;

  Array    *a = AForm_Array(o);
  Size_Type n = AForm_Size(o);

  t = 0;
  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *v = A<T>(a);
          switch (AForm_Class(o))
          { case SLICE_CLASS:
              for (p = Set_Slice_To_First(o); n-- > 0; p = Next_Slice_Index(o))
                t += v[p];
              break;
            case FRAME_CLASS:
              if (Frame_Within_Array(o))
                { Offs_Type *off = Frame_Offsets(o);
                  v += Frame_Index(o);
                  while (n-- > 0)
                    t += v[off[n]];
                  return (t);
                }
              else
                v = Frame_Values(o);
            case ARRAY_CLASS:
              while (n-- > 0)
                t += v[n];
              break;
          }
        }
    #END
  }
  return (t);
}

static void accumulate_model(Double_Vector *M(model), double wgt, AForm *o)
{ Indx_Type p;
  double   *m = AFLOAT64(model);

  Array    *a = AForm_Array(o);
  Size_Type n = AForm_Size(o);

  switch (a->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        { <t> *v = A<T>(a);
          switch (AForm_Class(o))
          { case SLICE_CLASS:
              { Slice *s = (Slice *) o;
                for (p = Set_Slice_To_First(s); n-- > 0; p = Next_Slice_Index(s))
                  *m++ += wgt*v[p];
                return;
              }
            case FRAME_CLASS:
              if (Frame_Within_Array(o))
                { Offs_Type *off = Frame_Offsets(o);
                  v += Frame_Index(o);
                  while (n-- > 0)
                    *m++ += wgt*v[off[n]];
                  return;
                }
              else
                v = Frame_Values(o);
            case ARRAY_CLASS:
              while (n-- > 0)
                *m++ += wgt*v[n];
              return;
          }
        }
    #END
  }
}

static Double_Matrix *G(correlate)(int n, AForm *a1, va_list *ap, string routine)
{ AForm         *Sargs[20], **args;
  Double_Matrix *rez;

  if (n > 20)
    args = (AForm **) Guarded_Malloc(sizeof(AForm *)*((size_t) n),routine);
  else
    args = Sargs;

  { int i;

    if (n > 0)
      args[0] = a1;
    else
      { fprintf(stderr,"List of arrays is empty (%s)\n",routine);
        exit (1);
      }
    for (i = 1; i < n; i++)
      { args[i] = (AForm *) va_arg(*ap,AForm *);
        if (AForm_Size(args[i]) != AForm_Size(args[0]))
          { fprintf(stderr,"Arrays do not have the same number of elements (%s)\n",routine);
            exit (1);
          }
      }
  }

  rez = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord2(n,n));
  
  { double *a = AFLOAT64(rez);
    int     i, j;

    for (i = 0; i < n; i++)
      for (j = 0; j <= i; j++)
        a[j*n+i] = a[i*n+j] = inner_product_arrays(args[i],args[j]);
  }

  if (n > 20)
    free(args);
  
  return (rez);
}

static Double_Vector *G(array_sums)(int n, AForm *a1, va_list *ap)
{ Double_Vector *rez;

  rez = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(n));
  
  { double *a = AFLOAT64(rez);
    int     i;

    a[0] = sum_array(a1);
    for (i = 1; i < n; i++)
      a[i] = sum_array((AForm *) va_arg(*ap,AForm *));
  }
  
  return (rez);
}

double A_Correlation(AForm *a1, AForm *a2)
{ return (inner_product_arrays(a1,a2)); }

Double_Matrix *G(Correlations)(int n, AForm *a1, ...)
{ va_list        ap;
  Double_Matrix *rez;

  va_start(ap,a1);
  rez = correlate(n,a1,&ap,"Correlation");
  va_end(ap);
  return (rez);
}

double A_Covariance(AForm *a1, AForm *a2)
{ double prod = inner_product_arrays(a1,a2);
  double sum1  = sum_array(a1);
  double sum2  = sum_array(a2);
  return ((prod - sum1*sum2)/AForm_Size(a1));
}

Double_Matrix *G(Covariances)(int n, AForm *a1, ...)
{ va_list        ap;
  Double_Matrix *rez;
  Double_Vector *sum;

  va_start(ap,a1);
  rez = correlate(n,a1,&ap,"Covariance");
  va_end(ap);

  va_start(ap,a1);
  sum = array_sums(n,a1,&ap);
  va_end(ap);

  { double   *a = AFLOAT64(rez);
    double   *m = AFLOAT64(sum);
    double    mi, mj, s;
    int       i, j;
    Size_Type N = AForm_Size(a1);

    for (i = 0; i < n; i++)
      { mi = m[i]/N;
        for (j = 0; j <= i; j++)
          { s  = a[j*n+i];
            mj = m[j];
            a[j*n+i] = a[i*n+j] = (s - mj*mi)/N;
          }
      }
  }

  Free_Array(sum);
  
  return (rez);
}

double A_Pearson_Correlation(AForm *a1, AForm *a2)
{ double prod = inner_product_arrays(a1,a2);
  double sig1 = inner_product_arrays(a1,a1);
  double sig2 = inner_product_arrays(a2,a2);
  double sum1 = sum_array(a1);
  double sum2 = sum_array(a2);
  Size_Type N = AForm_Size(a1);

  sig1 -= sum1*(sum1/N);
  sig2 -= sum2*(sum2/N);
  return ((prod - sum1*(sum2/N))/sqrt(sig1*sig2));
}

Double_Matrix *G(Pearson_Correlations)(int n, AForm *a1, ...)
{ va_list        ap;
  Double_Matrix *rez;
  Double_Vector *sum;

  va_start(ap,a1);
  rez = correlate(n,a1,&ap,"Pearson_Correlation");
  va_end(ap);

  va_start(ap,a1);
  sum = array_sums(n,a1,&ap);
  va_end(ap);

  { double   *a = AFLOAT64(rez);
    double   *m = AFLOAT64(sum);
    double    mi, mj, s;
    int       i, j;
    Size_Type N = AForm_Size(a1);

    for (i = 0; i < n; i++)
      { mi = m[i]/N;
        for (j = 0; j <= i; j++)
          { s  = a[j*n+i];
            mj = m[j];
            a[j*n+i] = a[i*n+j] = (s - mj*mi)/N;
          }
      }

    for (i = n*n-1; i >= 0; i -= n+1)
      a[i] = sqrt(a[i]);

    for (i = 1; i < n; i++)
      { mi = a[i*n+i];
        for (j = 0; j < i; j++)
          a[j*n+i] = a[i*n+j] /= (mi * a[j*n+j]);
      }

    for (i = n*n-1; i >= 0; i -= n+1)
      a[i] = 1.;
  }

  Free_Array(sum);
  
  return (rez);
}

Double_Vector *Linear_Regression(int n, Regression_Bundle *O(stats), AForm *obs, AForm *in1, ...)
{ va_list        ap;
  Double_Matrix *cof;
  Double_Vector *tgt;
  LU_Factor     *lu;
  boolean        stable;
  
  if (AForm_Size(obs) != AForm_Size(in1))
    { fprintf(stderr,"Arrays do not have the same number of elements (Linear_Regression)\n");
      exit (1);
    }

  va_start(ap,in1);
  cof = correlate(n,in1,&ap,"Linear_Regression");
  va_end(ap);

  lu = LU_Decompose(cof,&stable);
  if ( ! stable)
    fprintf(stderr,"Warning: Coefficient array inversion was unstable (Linear_Regression)\n");

  tgt = Make_Array_With_Shape(PLAIN_KIND,FLOAT64_TYPE,Coord1(n));
  
  { double *a = AFLOAT64(tgt);
    int     i;

    va_start(ap,in1);
    a[0] = inner_product_arrays(in1,obs);
    for (i = 1; i < n; i++)
      a[i] = inner_product_arrays((AForm *) va_arg(ap,AForm *),obs);
    va_end(ap);
  }
  
  LU_Solve(tgt,lu);

  if (stats != NULL)
    { Size_Type      N  = AForm_Size(in1);
      double         o2 = inner_product_arrays(obs,obs);
      double         oa = sum_array(obs);
      float64       *t  = AFLOAT64(tgt);

      int            i, j;
      Array         *sval, *tval;
      Dimn_Type      dims[1];
      Double_Vector *model;

      dims[0] = (Dimn_Type) N;
      model = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      Array_Op_Scalar(model,SET_OP,FVAL,VALF(0.));

      va_start(ap,in1);
      accumulate_model(model,t[0],in1);
      for (i = 1; i < n; i++)
        accumulate_model(model,t[i],(AForm *) va_arg(ap,AForm *));
      va_end(ap);

      { double m2 = inner_product_arrays(model,model);
        double ma = sum_array(model);

        stats->tss = o2 - oa*(oa/N);
        stats->mss = m2 - (2*ma-oa)*(oa/N);
      }

      accumulate_model(model,-1.,obs);

      stats->rss = inner_product_arrays(model,model);

      stats->R2  = 1. - stats->rss / stats->tss;
      stats->aR2 = 1. - ((N-1)*(1.-stats->R2))/(N-n);
      stats->ser = sqrt(stats->rss / (N-n));

      Free_Array(model);

      dims[0] = n;
      if ((sval = stats->std_err) == NULL)
        sval = stats->std_err = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      else
        { if (sval->size < n)
            allocate_array_data(sval,SIZEOF(float64)*n,"Linear_Regression");
          sval->kind  = PLAIN_KIND;
          sval->type  = FLOAT64_TYPE;
          sval->scale = 64;
          sval->ndims = 1;
          sval->size  = n;
          sval->dims[0] = n;
        }
      if ((tval = stats->t_stat) == NULL)
        tval = stats->t_stat = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      else
        { if (tval->size < n)
            allocate_array_data(tval,SIZEOF(float64)*n,"Linear_Regression");
          tval->kind  = PLAIN_KIND;
          tval->type  = FLOAT64_TYPE;
          tval->scale = 64;
          tval->ndims = 1;
          tval->size  = n;
          tval->dims[0] = n;
        }

      { float64 *se = AFLOAT64(sval);
        float64 *ts = AFLOAT64(tval);

        for (i = 0; i < n; i++)
          { for (j = 0; j < n; j++)
              ts[j] = 0.;
            ts[i] = 1.;
            LU_Solve(tval,lu);
            se[i] = stats->ser*sqrt(ts[i]);
          }
        for (i = 0; i < n; i++)
          ts[i] = t[i] / se[i];
       }
    }

  Free_LU_Factor(lu);

  return (tgt);
}

Double_Vector *Simple_Regression(Array *R(O(vector)), Regression_Bundle *O(stats),
                                 AForm *obs, AForm *inp)
{ double    cof00, cof10, cof11;
  double    trg0, trg1;
  double    d, *data;
  Dimn_Type dims[1];

  if (AForm_Size(obs) != AForm_Size(inp))
    { fprintf(stderr,"Arrays do not have the same number of elements (Simplest_Regression)\n");
      exit (1);
    }

  dims[0] = 2;
  if (vector == NULL)
    vector = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
  else
    { if (vector->size < 2)
        allocate_array_data(vector,sizeof(float64)*2,"Simple_Regression");
      vector->kind  = PLAIN_KIND;
      vector->type  = FLOAT64_TYPE;
      vector->scale = 64;
      vector->ndims = 1;
      vector->size  = 2;
      vector->dims[0] = 2;
    }
  data = AFLOAT64(vector);

  cof00 = inner_product_arrays(inp,inp);
  cof10 = sum_array(inp);
  cof11 = (double) AForm_Size(inp);

  trg0  = inner_product_arrays(obs,inp);
  trg1  = sum_array(obs);

  d = (cof11*cof00 - cof10*cof10);
  data[0] = (cof00*trg1 - cof10*trg0) / d;
  data[1] = (cof11*trg0 - cof10*trg1) / d;

  if (stats != NULL)
    { Size_Type      N  = AForm_Size(inp);
      double         o2 = inner_product_arrays(obs,obs);
      double         oa = sum_array(obs);

      Array         *sval, *tval;
      Double_Vector *model;

      dims[0] = (Dimn_Type) N;
      model = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      Array_Op_Scalar(model,SET_OP,FVAL,VALF(data[0]));
      accumulate_model(model,data[1],inp);

      { double m2 = inner_product_arrays(model,model);
        double ma = sum_array(model);

        stats->tss = o2 - oa*(oa/N);
        stats->mss = m2 - (2*ma-oa)*(oa/N);
      }

      accumulate_model(model,-1.,obs);

      stats->rss = inner_product_arrays(model,model);

      stats->R2  = 1. - stats->rss / stats->tss;
      stats->aR2 = 1. - ((N-1)*(1.-stats->R2))/(N-2);
      stats->ser = sqrt(stats->rss / (N-2));

      Free_Array(model);

      dims[0] = 2;
      if ((sval = stats->std_err) == NULL)
        sval = stats->std_err = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      else
        { if (sval->size < 2)
            allocate_array_data(sval,sizeof(float64)*2,"Simple_Regression");
          sval->kind  = PLAIN_KIND;
          sval->type  = FLOAT64_TYPE;
          sval->scale = 64;
          sval->ndims = 1;
          sval->size  = 2;
          sval->dims[0] = 2;
        }
      if ((tval = stats->t_stat) == NULL)
        tval = stats->t_stat = Make_Array(PLAIN_KIND,FLOAT64_TYPE,1,dims);
      else
        { if (tval->size < 2)
            allocate_array_data(tval,sizeof(float64)*2,"Simple_Regression");
          tval->kind  = PLAIN_KIND;
          tval->type  = FLOAT64_TYPE;
          tval->scale = 64;
          tval->ndims = 1;
          tval->size  = 2;
          tval->dims[0] = 2;
        }

      { float64 *se = AFLOAT64(sval);
        float64 *ts = AFLOAT64(tval);

        se[0] =  stats->ser * sqrt(cof00/d);
        se[1] =  stats->ser * sqrt(cof11/d);
        ts[0] =  data[0] / se[0];
        ts[1] =  data[1] / se[1];
      }
    }

  return (vector);
}
