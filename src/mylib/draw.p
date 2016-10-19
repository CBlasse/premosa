/******************************************************************************************\
*                                                                                         *
*  Drawing Level Sets, Regions, and Basic Shapes                                          *
*                                                                                         *
*  Author  :  Gene Myers                                                                  *
*  Date    :  June 2007                                                                   *
*  Last Mod:  Aug 2008                                                                    *
*                                                                                         *
*  (c) June 19, '09, Dr. Gene Myers and Howard Hughes Medical Institute                   *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include <math.h>

#include "utilities.h"
#include "connectivity.h"
#include "draw.h"
#include "level.set.h"

//  A Brush_Bundle determines the setting of the following globals that are then
//    used by the various painter routines below to set pixels.

#LISTDEF @TYPES  =  UINT8 UINT16 UINT32 UINT64 INT8 INT16 INT32 INT64 FLOAT32 FLOAT64
#LISTDEF @UNION  =   uval   uval   uval   uval ival  ival  ival  ival    fval    fval

#LISTDEF @OP_NAME  = SET_PIX ADD_PIX SUB_PIX MUL_PIX DIV_PIX MIN_PIX MAX_PIX
#LISTDEF @OPERATOR =       =      +=      -=       *       /       <       >

typedef struct        //  Generic container for a "Paint bundle"
  { void   *red;
    void   *green;
    void   *blue;
    void   *alpha;
    double  redI;
    double  greenI;
    double  blueI;
    double  alphaI;
  } PAINT_VOID;

 //  Paint brush bundle for each type & operator combo

#GENERATE T,U = @TYPES , @UNION
  typedef struct
    { <t>    *red;       //  Typed paint bundle for non-multiplicative ops
      <t>    *green;
      <t>    *blue;
      <t>    *alpha;
      <t>     redI;
      <t>     greenI;
      <t>     blueI;
      <t>     alphaI;
    } PAINT_<T>;

  typedef struct
    { <t>    *red;       //  Typed paint bundle for multiplicative ops
      <t>    *green;
      <t>    *blue;
      <t>    *alpha;
      float64 redI;
      float64 greenI;
      float64 blueI;
      float64 alphaI;
    } MUPNT_<T>;


  #define CP<T>(a)  ((PAINT_<T> *) (a))
  #define CR<T>(a)  ((MUPNT_<T> *) (a))

#END

  //  Generate a painter for each image kind, array type, and paint operator combination!

#GENERATE T = @TYPES
  #GENERATE OP,SYM = @OP_NAME , @OPERATOR
    #IF OP < MUL_PIX

      static void PAINT_<T>_PLAIN_<OP>(Indx_Type p, void *a)
      { CP<T>(a)->red[p] <SYM> CP<T>(a)->redI; }

      static void PAINT_<T>_RGB_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green blue
          CP<T>(a)-><C>[p] <SYM> CP<T>(a)-><C>I;
        #END
      }

      static void PAINT_<T>_RGBA_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green blue alpha
          CP<T>(a)-><C>[p] <SYM> CP<T>(a)-><C>I;
        #END
      }

      static void PAINT_<T>_COMPLEX_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green
          CP<T>(a)-><C>[p] <SYM> CP<T>(a)-><C>I;
        #END
      }

    #ELSEIF OP < MIN_PIX

      static void PAINT_<T>_PLAIN_<OP>(Indx_Type p, void *a)
      { CR<T>(a)->red[p] = (<t>) (CR<T>(a)->red[p] <SYM> CR<T>(a)->redI); }

      static void PAINT_<T>_RGB_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green blue
          CR<T>(a)-><C>[p] = (<t>) (CR<T>(a)-><C>[p] <SYM> CR<T>(a)-><C>I);
        #END
      }

      static void PAINT_<T>_RGBA_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green blue alpha
          CR<T>(a)-><C>[p] = (<t>) (CR<T>(a)-><C>[p] <SYM> CR<T>(a)-><C>I);
        #END
      }

      static void PAINT_<T>_COMPLEX_<OP>(Indx_Type p, void *a) {
        #GENERATE C = red green
          CR<T>(a)-><C>[p] = (<t>) (CR<T>(a)-><C>[p] <SYM> CR<T>(a)-><C>I);
        #END
      }

    #ELSE

      static void PAINT_<T>_PLAIN_<OP>(Indx_Type p, void *a)
      { <t> x = CP<T>(a)->red[p];
        if (CP<T>(a)->redI <SYM> x) CP<T>(a)->red[p] = CP<T>(a)->redI;
      }

      static void PAINT_<T>_RGB_<OP>(Indx_Type p, void *a) {
        <t> x;
        #GENERATE C = red green blue
          x = CP<T>(a)-><C>[p];
          if (CP<T>(a)-><C>I <SYM> x) CP<T>(a)-><C>[p] = CP<T>(a)-><C>I;
        #END
      }

      static void PAINT_<T>_RGBA_<OP>(Indx_Type p, void *a) {
        <t> x;
        #GENERATE C = red green blue alpha
          x = CP<T>(a)-><C>[p];
          if (CP<T>(a)-><C>I <SYM> x) CP<T>(a)-><C>[p] = CP<T>(a)-><C>I;
        #END
      }

      static void PAINT_<T>_COMPLEX_<OP>(Indx_Type p, void *a) {
        <t> x;
        p <<= 1;
        #GENERATE C = red green
          x = CP<T>(a)-><C>[p];
          if (CP<T>(a)-><C>I <SYM> x) CP<T>(a)-><C>[p] = CP<T>(a)-><C>I;
        #END
      }

    #END
  #END
#END

static void (*Painter_Table[])(Indx_Type,void *) = {   //  Make a table of all these puppies
    #GENERATE K = PLAIN RGB RGBA COMPLEX
      #GENERATE T = @TYPES
        #GENERATE OP = @OP_NAME
          PAINT_<T>_<K>_<OP>,
        #END
      #END
    #END
  };
  
void *SETUP_PAINTER(Array *canvas, Brush_Bundle *generic, void *argp)
{ Color_Bundle *brush = (Color_Bundle  *) generic;
  int           kind  = canvas->kind;

  if (brush->op == MUL_PIX || brush->op == DIV_PIX)
    switch (canvas->type) {
      #GENERATE T,U = @TYPES , @UNION
        case <T>_TYPE:
          { MUPNT_<T> *ap = (MUPNT_<T> *) argp;
  
            ap->red  = (<t> *) (canvas->data);
            ap->redI = brush->red.fval;
            if (kind > PLAIN_KIND)
              { if (kind == COMPLEX_KIND)
                  { ap->green  = ap->red + 1;
                    ap->greenI = brush->green.fval;
                  }
                else
                  { Indx_Type n = canvas->size / canvas->dims[canvas->ndims-1];
                    ap->green  = ap->red + n;
                    ap->blue   = ap->green + n;
                    ap->greenI = brush->green.fval;
                    ap->blueI  = brush->blue.fval;
                    if (kind == RGBA_KIND)
                      { ap->alpha  = ap->blue + n;
                        ap->alphaI = brush->alpha.fval;
                      }
                  }
              }
            break;
          }
      #END
    }
  else
    switch (canvas->type) {
      #GENERATE T,U = @TYPES , @UNION
        case <T>_TYPE:
          { PAINT_<T> *ap = (PAINT_<T> *) argp;

            ap->red  = (<t> *) (canvas->data);
            ap->redI = (<t>) brush->red.<u>;
            if (kind > PLAIN_KIND)
              { if (kind == COMPLEX_KIND)
                  { ap->green  = ap->red + 1;
                    ap->greenI = (<t>) brush->green.<u>;
                  }
                else
                  { Indx_Type n = canvas->size / canvas->dims[canvas->ndims-1];
                    ap->green  = ap->red + n;
                    ap->blue   = ap->green + n;
                    ap->greenI = (<t>) brush->green.<u>;
                    ap->blueI  = (<t>) brush->blue.<u>;
                    if (kind == RGBA_KIND)
                      { ap->alpha  = ap->blue + n;
                        ap->alphaI = (<t>) brush->alpha.<u>;
                      }
                  }
              }
            break;
          }
      #END
    }

  return (Painter_Table[(canvas->kind*10 + canvas->type)*7 + brush->op]);
}

/****************************************************************************************
 *                                                                                      *
 *  DRAWING ROUTINES FOR REGIONS                                                        *
 *                                                                                      *
 ****************************************************************************************/

void check_shapes(Array *canvas, Region *reg, char *routine)
{ Array_Bundle rarray;

  rarray.kind  = PLAIN_KIND;
  rarray.ndims = Get_Region_Dimensionality(reg);
  rarray.dims  = Get_Region_Dimensions(reg);
  if (canvas->kind != PLAIN_KIND)
    { canvas->ndims -= 1;
      if (canvas->kind == COMPLEX_KIND)
        canvas->dims += 1;
    }
  if ( ! Same_Shape(canvas,&rarray))
    { fprintf(stderr,"Canvas and image from which region was computed do not have the same shape");
      fprintf(stderr," (%s)\n",routine);
      exit (1);
    }
  if (canvas->kind != PLAIN_KIND)
    { if (canvas->kind == COMPLEX_KIND)
        canvas->dims -= 1;
      canvas->ndims += 1;
    }
}

/* Color the pixels of region reg */

void Draw_Region_Outline(Array *M(canvas), Brush_Bundle *brush, Region *reg)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;

  check_shapes(canvas,reg,"Draw_Region_Outline");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Outline(reg,argp,painter);
}

/* Color the region defined by reg */

void Draw_Region(Array *M(canvas), Brush_Bundle *brush, Region *reg)
{ Indx_Type *raster;
  Indx_Type  len;
  void     (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;

  check_shapes(canvas,reg,"Draw_Region");

  painter = SETUP_PAINTER(canvas,brush,argp);

  raster = reg->raster;
  len    = reg->rastlen;

  { Indx_Type i, v, w, p;

    for (i = 0; i < len; i += 2)
      { v = raster[i];
        w = raster[i+1];
        for (p = v; p <= w; p++)
          painter(p,argp);
      }
  }
}

/* Color the complement of the region defined by reg */

void Draw_Region_Exterior(Array *M(canvas), Brush_Bundle *brush, Region *reg)
{ void     (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;

  check_shapes(canvas,reg,"Draw_Region_Exterior");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Exterior(reg,argp,painter);
}

void Draw_Region_Holes(Array *M(canvas), Brush_Bundle *brush, Region *reg)
{ void     (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;

  check_shapes(canvas,reg,"Draw_Region_Holes");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Holes(reg,argp,painter);
}

/****************************************************************************************
 *                                                                                      *
 *  DRAWING ROUTINES FOR LEVEL SETS                                                     *
 *                                                                                      *
 ****************************************************************************************/

void Draw_Level_Set_Outline(Array *M(canvas), Brush_Bundle *brush,
                            Level_Tree *t, Level_Set *r, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  c = Record_Level_Set(t,r,share,0);

  check_shapes(canvas,c,"Draw_Level_Set_Outline");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Outline(c,argp,painter);

  Free_Region(c);
}

void Draw_Level_Set_Exterior(Array *M(canvas), Brush_Bundle *brush,
                             Level_Tree *t, Level_Set *r, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  c = Record_Level_Set(t,r,share,0);

  check_shapes(canvas,c,"Draw_Level_Set_Exterior");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Exterior(c,argp,painter);

  Free_Region(c);
}

void Draw_Level_Set_Holes(Array *M(canvas), Brush_Bundle *brush,
                          Level_Tree *t, Level_Set *r, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  c = Record_Level_Set(t,r,share,1);

  check_shapes(canvas,c,"Draw_Level_Set_Holes");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Holes(c,argp,painter);

  Free_Region(c);
}


/****************************************************************************************
 *                                                                                      *
 *  DRAWING ROUTINES FOR WATERHSEDS                                                     *
 *                                                                                      *
 ****************************************************************************************/

void Draw_P_Vertex_Outline(Array *M(canvas), Brush_Bundle *brush, Partition *w, int cb, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Draw_P_Vertex_Outline)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,share,0);

  check_shapes(canvas,c,"Draw_P_Vertex_Outline");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Outline(c,argp,painter);

  Free_Region(c);
}

void Draw_P_Vertex_Exterior(Array *M(canvas), Brush_Bundle *brush, Partition *w, int cb, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Draw_P_Vertex_Exterior)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,share,0);

  check_shapes(canvas,c,"Draw_P_Vertex_Exterior");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Exterior(c,argp,painter);

  Free_Region(c);
}

void Draw_P_Vertex_Holes(Array *M(canvas), Brush_Bundle *brush, Partition *w, int cb, int share)
{ void      (*painter)(Indx_Type,void *);
  PAINT_VOID arg, *argp = &arg;
  Region    *c;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Draw_P_Vertex_Holes)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,share,1);

  check_shapes(canvas,c,"Draw_P_Vertex_Holes");

  painter = SETUP_PAINTER(canvas,brush,argp);

  For_Region_Holes(c,argp,painter);

  Free_Region(c);
}


/****************************************************************************************
 *                                                                                      *
 *  IMAGE DRAWING ROUTINES FOR REGION, LEVEL SETS, & WATERHSEDS                         *
 *                                                                                      *
 ****************************************************************************************/

void draw_region_image(Array *M(canvas), Array *image, Region *reg, string routine)
{ Indx_Type *raster;
  Indx_Type  len;

  check_shapes(canvas,reg,routine);
  if ( ! Same_Type(canvas,image))
    { fprintf(stderr,"Canvas and image must have the same shape and type (routine)\n");
      exit (1);
    }
 
  raster = reg->raster;
  len    = reg->rastlen;

  { Indx_Type i, v, w, p;

    switch (canvas->kind)
    { case PLAIN_KIND:
        switch (canvas->type) {
         #GENERATE T = @TYPES
          case <T>_TYPE:
            { <t> *val  = A<T>(canvas);
              <t> *valI = A<T>(image);
              for (i = 0; i < len; i += 2)
                { v = raster[i];
                  w = raster[i+1];
                  for (p = v; p <= w; p++)
                    val[p] = valI[p]; 
                }
              break;
            }
         #END
        }
        break;
      case RGB_KIND:
        switch (canvas->type) {
         #GENERATE T = @TYPES
          case <T>_TYPE:
            { Size_Type N = canvas->size / 3;
              <t> *red    = A<T>(canvas);
              <t> *redI   = A<T>(image);
              <t> *green  = red + N;
              <t> *greenI = redI + N;
              <t> *blue   = green + N;
              <t> *blueI  = greenI + N;
              for (i = 0; i < len; i += 2)
                { v = raster[i];
                  w = raster[i+1];
                  for (p = v; p <= w; p++)
                    { red[p]   = redI[p]; 
                      green[p] = greenI[p]; 
                      blue[p]  = blueI[p]; 
                    }
                }
              break;
            }
         #END
        }
        break;
      case RGBA_KIND:
        switch (canvas->type) {
         #GENERATE T = @TYPES
          case <T>_TYPE:
            { Size_Type N = canvas->size / 4;
              <t> *red    = A<T>(canvas);
              <t> *redI   = A<T>(image);
              <t> *green  = red + N;
              <t> *greenI = redI + N;
              <t> *blue   = green + N;
              <t> *blueI  = greenI + N;
              <t> *alpha  = blue + N;
              <t> *alphaI = blueI + N;
              for (i = 0; i < len; i += 2)
                { v = raster[i];
                  w = raster[i+1];
                  for (p = v; p <= w; p++)
                    { red[p]   = redI[p]; 
                      green[p] = greenI[p]; 
                      blue[p]  = blueI[p]; 
                      alpha[p] = alphaI[p]; 
                    }
                }
              break;
            }
         #END
        }
        break;
      case COMPLEX_KIND:
        switch (canvas->type) {
         #GENERATE T = @TYPES
          case <T>_TYPE:
            { <t> *real  = A<T>(canvas);
              <t> *realI = A<T>(image);
              <t> *imag  = real + 1;
              <t> *imagI = realI + 1;
              for (i = 0; i < len; i += 2)
                { v = raster[i];
                  w = raster[i+1];
                  for (p = v; p <= w; p++)
                    { real[p] = realI[p]; 
                      imag[p] = imagI[p]; 
                    }
                }
              break;
            }
         #END
        }
        break;
    }
  }
}

void Draw_Region_Image(Array *M(canvas), Array *image, Region *reg)
{ draw_region_image(canvas,image,reg,"Draw_Region_Image"); }

void Draw_Level_Set_Image(Array *M(canvas), Array *image, Level_Tree *t, Level_Set *r, int share)
{ Region *c;
  c = Record_Level_Set(t,r,share,1);
  draw_region_image(canvas,image,c,"Draw_Level_Set_Image");
  Free_Region(c);
}

void Draw_P_Vertex_Image(Array *M(canvas), Array *image, Partition *w, int cb, int share)
{ Region *c;

  if (Get_Partition_Labels(w) == NULL)
    { fprintf(stderr,"Partition does not have a label array (Draw_P_Vertex_Image)\n");
      exit (1);
    }

  c = Record_P_Vertex(w,cb,share,1);
  draw_region_image(canvas,image,c,"Draw_P_Vertex_Image");
  Free_Region(c);
}

/****************************************************************************************
 *                                                                                      *
 *  GENERAL FLOOD-FILL DRAWING ROUTINE                                                  *
 *                                                                                      *
 ****************************************************************************************/

typedef struct
  { void    *value;
    uint64   level_uval;
    int64    level_ival;
    float64  level_fval;
  } DrawArg;

#define DA(a) ((DrawArg *) (a))

#GENERATE T,U = @TYPES , @UNION
  #GENERATE C,O = LE LT EQ NE GT GE , <= < == != > >=
    static boolean is_<c>_<t>(Indx_Type p, void *a)
    { return (((<t> *) DA(a)->value)[p] <o> DA(a)->level_<u>); }
  #END
#END

static boolean (*Comparator_Table[])(Indx_Type,void *) = {
#GENERATE T = @TYPES
  #GENERATE C = LE LT EQ NE GT GE
    is_<c>_<t>,
  #END
#END
  };

static void check_drawing_compatibility(Array *canvas, Array *source, char *routine)
{ if (source->kind != PLAIN_KIND)
    { fprintf(stderr,"Source must be a plain array (%s)\n",routine);
      exit (1);
    }
  if (canvas->kind != PLAIN_KIND)
    { canvas->ndims -= 1;
      if (canvas->kind == COMPLEX_KIND)
        canvas->dims += 1;
    }
  if ( ! Same_Shape(canvas,source))
    { fprintf(stderr,"Canvas and source do not have the same shape! (%s)\n",routine);
      exit (1);
    }
  if (canvas->kind != PLAIN_KIND)
    { if (canvas->kind == COMPLEX_KIND)
        canvas->dims -= 1;
      canvas->ndims += 1;
    }
}

void Draw_Floodfill(Array *M(canvas), Brush_Bundle *brush,
                    APart *source, int share, boolean iscon2n,
                    Indx_Type seed, void *arg, boolean (*test)(Indx_Type p, void *arg))
{ Array      *array = AForm_Array(source);
  void      (*painter)(Indx_Type,void *);
  PAINT_VOID parg;

  check_drawing_compatibility(canvas,array,"Draw_Floodfill");

  painter = SETUP_PAINTER(canvas,brush,&parg);

  Flood_Object(source,share,iscon2n,seed,arg,test,NULL,NULL,NULL,NULL,&parg,painter);
}

void Draw_Basic(Array *M(canvas), Brush_Bundle *brush, APart *source, int share,
                boolean iscon2n, Indx_Type seed, Comparator cmprsn, Value level)
{ Array      *array = AForm_Array(source);
  boolean   (*value)(Indx_Type,void *);
  void      (*painter)(Indx_Type,void *);
  DrawArg    targ;
  PAINT_VOID parg;

  check_drawing_compatibility(canvas,array,"Draw_Basic");

  value   = Comparator_Table[6*array->type + cmprsn];
  painter = SETUP_PAINTER(canvas,brush,&parg);

  switch (array->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        targ.value = array->data;
        targ.level_<u> = level.<u>;
        break;
    #END
  }

  Flood_Object(source,share,iscon2n,seed,&targ,value,NULL,NULL,NULL,NULL,&parg,painter);
}

void Draw_P_Vertex(Array *M(canvas), Brush_Bundle *brush, Partition *shed, int cb, int share)
{ Array        *source  = Get_Partition_Labels(shed);
  boolean       iscon2n = Is_Partition_2n_Connected(shed);
  Indx_Type     seed    = Get_Partition_Vertex(shed,cb)->seed;
  boolean      (*value)(Indx_Type,void *);
  void         (*painter)(Indx_Type,void *);
  DrawArg       targ;
  PAINT_VOID    parg;

  if (source == NULL)
    { fprintf(stderr,"Partition does not have a label array (Draw_P_Vertex_Image)\n");
      exit (1);
    }

  check_drawing_compatibility(canvas,source,"Draw_P_Vertex");

  value   = Comparator_Table[6*source->type + EQ_COMP];
  painter = SETUP_PAINTER(canvas,brush,&parg);

  switch (source->type) {
    #GENERATE T,U = @TYPES , @UNION
      case <T>_TYPE:
        targ.value = source->data;
        targ.level_<u> = A<T>(source)[seed];
        break;
    #END
  }

  Flood_Object(source,share,iscon2n,seed,&targ,value,NULL,NULL,NULL,NULL,&parg,painter);
}

void Draw_Level_Set(Array *M(canvas), Brush_Bundle *brush, Level_Tree *t, Level_Set *r, int share)
{ APart        *source  = Get_Level_Tree_APart(t);
  Array        *array   = AForm_Array(source);
  boolean       iscon2n = Get_Level_Tree_Connectivity(t);
  Indx_Type     seed    = Level_Set_Leftmost(t,r);
  void         (*painter)(Indx_Type,void *);
  boolean      (*value)(Indx_Type,void *);
  DrawArg       targ;
  PAINT_VOID    parg;

  check_drawing_compatibility(canvas,array,"Draw_Level_Set");

  value   = Comparator_Table[6*array->type + 5];
  painter = SETUP_PAINTER(canvas,brush,&parg);

  targ.value      = array->data;
  targ.level_uval = (uint32) Level_Set_Level(t,r);

  Flood_Object(source,share,iscon2n,seed,&targ,value,NULL,NULL,NULL,NULL,&parg,painter);
}

/****************************************************************************************
 *                                                                                      *
 *  DRAWING ROUTINES FOR BASIC SHAPES                                                   *
 *                                                                                      *
 ****************************************************************************************/

static Array_Bundle *base_shape(Array *a, Array_Bundle *base)
{ base->ndims = a->ndims;
  base->dims  = a->dims;
  if (a->kind != PLAIN_KIND)
    { if (a->kind == COMPLEX_KIND)
        base->dims  += 1;
      base->ndims -= 1;
    }
  return (base);
}

void check_ivector(Coordinate *coord, Array *base, char *routine)
{ if (coord->kind != PLAIN_KIND || coord->type != DIMN_TYPE || coord->ndims != 1)
    { fprintf(stderr,"Coordinate is not a Dimn_Type vector (%s)\n",routine);
      exit (1);
    }
  if (coord->dims[0] != (Dimn_Type) base->ndims)
    { fprintf(stderr,"Coordinate is not of the correct dimensionality (%s)\n",routine);
      exit (1);
    }
}

/* Draw a size[0] x size[1] x ... x size[n-1] rectangle with lower left corner
   (corner[0],corner[1],...,corner[n-1]) where n = canvas->ndims-1.              */

typedef struct
  { Dimn_Type  *point1;
    Dimn_Type  *point2;
    Dimn_Type  *dims;
    void      (*painter)(Indx_Type,void *);
    PAINT_VOID *arg;
  } DrawRect; 

static void rectangle(int k, Indx_Type base, DrawRect *argp)
{ Dimn_Type beg = argp->point1[k];
  Dimn_Type end = argp->point2[k];
  Dimn_Type p;

  base *= argp->dims[k];
  if (k == 0)
    { for (p = beg; p <= end; p++)
        argp->painter(base+p,argp->arg);
    }
  else
    { for (p = beg; p <= end; p++)
        rectangle(k-1,base+p,argp);
    }
}

void Draw_Rectangle(Array *M(canvas), Brush_Bundle *brush, Coordinate *F(corner1),
                                                           Coordinate *F(corner2))
{ Array_Bundle base;
  int          i,   n;
  Dimn_Type   *c1, *c2, *d;
  PAINT_VOID   darg;
  DrawRect     rarg;
  Dimn_Type    Point1[10], *point1;
  Dimn_Type    Point2[10], *point2;

  base_shape(canvas,&base);

  check_ivector(corner1,&base,"Draw_Rectangle");
  check_ivector(corner2,&base,"Draw_Rectangle");

  n = base.ndims;
  d = base.dims;

  if (n > 10)
    { point1 = (Dimn_Type *) Guarded_Malloc(sizeof(Dimn_Type)*2*((size_t) n),"Draw_Rectangle");
      point2 = point1 + n;
    }
  else
    { point1 = Point1;
      point2 = Point2;
    }

  rarg.point1  = point1;
  rarg.point2  = point2;
  rarg.dims    = d;
  rarg.arg     = &darg;
  rarg.painter = SETUP_PAINTER(canvas,brush,&darg);

  c1 = ADIMN(corner1);
  c2 = ADIMN(corner2);
  for (i = 0; i < n; i++)
    { Dimn_Type beg = c1[i];
      Dimn_Type end = c2[i];
      if (beg < 0)
        beg  = 0;
      if (end >= d[i])
        end = d[i]-1;
      point1[i] = beg;
      point2[i] = end;
    }

  rectangle(n-1,0,&rarg);

  if (n > 10)
    free(point1);
  Free_Array(corner1);
  Free_Array(corner2);
}

/* Reset an entire image */

void Draw_Image(Array *M(canvas), Brush_Bundle *brush)
{ Array_Bundle base;
  int          i, n;
  Dimn_Type   *d;
  PAINT_VOID   darg;
  DrawRect     rarg;
  Dimn_Type    Point1[10], *point1;
  Dimn_Type    Point2[10], *point2;

  base_shape(canvas,&base);

  n = base.ndims;
  d = base.dims;

  if (n > 10)
    { point1 = (Dimn_Type *) Guarded_Malloc(sizeof(Dimn_Type)*2*((size_t) n),"Draw_Image");
      point2 = point1 + n;
    }
  else
    { point1 = Point1;
      point2 = Point2;
    }

  rarg.point1  = point1;
  rarg.point2  = point2;
  rarg.dims    = d;
  rarg.arg     = &darg;
  rarg.painter = SETUP_PAINTER(canvas,brush,&darg);

  for (i = 0; i < n; i++)
    { point1[i] = 0;
      point2[i] = d[i]-1;
    }

  rectangle(n-1,0,&rarg);

  if (n > 10)
    free(point1);
}

/* Draw a point centered a pixel point */

void Draw_Point(Array *M(canvas), Brush_Bundle *brush, Coordinate *F(point))
{ Array_Bundle base;
  Dimn_Type   *p;
  int          i;
  void        (*painter)(Indx_Type,void *);
  PAINT_VOID   arg, *argp = &arg;

  base_shape(canvas,&base);

  check_ivector(point,&base,"Draw_Cross");

  p = ADIMN(point);
  for (i = 0; i < base.ndims; i++)
    if (p[i] >= base.dims[i])
      return;

  painter = SETUP_PAINTER(canvas,brush,argp);
  painter(Coord2IdxA(&base,point),argp);
}

/* Draw a cross centered at pixel center with each arm being radius pixels long */

void Draw_Cross(Array *M(canvas), Brush_Bundle *brush, Coordinate *F(center), int radius)
{ Array_Bundle base;
  int          i, n;
  Dimn_Type   *d, *c;
  Indx_Type    p, q;
  void         (*painter)(Indx_Type,void *);
  PAINT_VOID   arg, *argp = &arg;

  base_shape(canvas,&base);

  check_ivector(center,&base,"Draw_Cross");

  n = base.ndims;
  d = base.dims;

  painter = SETUP_PAINTER(canvas,brush,argp);

  q = 1;
  c = ADIMN(center);
  p = Coord2IdxA(&base,center);
  for (i = 0; i < n; i++)
    { int64 x = c[i];
      int64 b, e, k;

      if (x < radius)
        b = -x;
      else
        b = -radius;
      if (x+radius >= d[i])
        e = (d[i] - x) - 1;
      else
        e = radius;
      for (k = b; k <= e; k++)
        if (k != 0)
          painter(p+q*k,argp);
      q *= d[i];
    }
  painter(p,argp);
}

/* Draw an n-dimensional circle centered at pixel center with radius radius,
     where n = canvas->ndims                                                  */

typedef struct
  { Dimn_Type *center;
    Dimn_Type *dims;
    void      *arg;
    void     (*painter)(Indx_Type,void *);
  } DrawCirc;

static void circle(int k, Indx_Type base, int64 rem, int64 rad, DrawCirc *argp)
{ int64 i, beg, end, rng;

  while (rad*rad > rem)
    rad -= 1;

  beg = argp->center[k] - rad;
  if (beg < 0)
    beg = 0;
  end = argp->center[k] + rad;
  if (end > argp->dims[k])
    end = argp->dims[k]-1;

  base *= argp->dims[k];
  if (k == 0)
    { for (i = beg; i <= end; i++)
        argp->painter(base+i,argp->arg);
    }
  else
    { rng = beg-argp->center[k];
      for (i = beg; i <= end; i++)
        { circle(k-1,base+i,rem-rng*rng,rad,argp);
          rng += 1;
        }
    }
}

void Draw_Circle(Array *M(canvas), Brush_Bundle *brush, Coordinate *F(center), int radius)
{ Array_Bundle base;
  PAINT_VOID   darg;
  DrawCirc     carg;

  base_shape(canvas,&base);

  check_ivector(center,&base,"Draw_Circle");

  carg.dims    = base.dims;
  carg.center  = ADIMN(center);
  carg.arg     = &darg;
  carg.painter = SETUP_PAINTER(canvas,brush,&darg);

  circle(base.ndims-1,0,radius*radius,(Dimn_Type) radius,&carg);

  Free_Array(center);
}

/*  Draw an n-dimensional line from begp to endp, clipping to the canvas as necessary */

void Draw_Line(Array *M(canvas), Brush_Bundle *brush, Coordinate *F(begp), Coordinate *F(endp))
{ int            ndims;
  Dimn_Type     *beg,  *end, *dims;
  void          (*painter)(Indx_Type,void *);
  int            kmax;
  int64          bgk, enk;
  PAINT_VOID     arg, *argp = &arg;
  double         Val[10], *val;
  double         Inc[10], *inc;

  { Array_Bundle base;

    base_shape(canvas,&base);

    check_ivector(begp,&base,"Draw_Line");
    check_ivector(endp,&base,"Draw_Line");

    painter = SETUP_PAINTER(canvas,brush,argp);

    beg  = ADIMN(begp);
    end  = ADIMN(endp);

    ndims = base.ndims;
    dims  = base.dims;

    if (ndims > 10)
      { val = (double *) Guarded_Malloc(sizeof(double)*2*((size_t) ndims),"Draw_Line");
        inc = val + ndims;
      }
    else
      { val = Val;
        inc = Inc;
      }
  }

  { Dimn_Type d, maxd = 0;           //  kmax = dimension with largest delta
    int       i;

    kmax = 0;
    for (i = 0; i < ndims; i++)
      { if (end[i] > beg[i])
          val[i] = d = end[i] - beg[i];
        else
          { d = beg[i] - end[i];
            val[i] = -1.*d;
          }
        if (d > maxd)
          { maxd = d;
            kmax = i;
          }
      }
    if (maxd == 0) goto exit_line;
  }

  { double ab, ae;
    double dkm, bkm;
    double den, dim;
    double abg, aen;
    int    i;

    ab = 0.;                   //  figure the clip interval [ab,ae] <= [0,1] for the line
    ae = 1.;
    for (i = 0; i < ndims; i++)
      { den = val[i];
        dim = dims[i] - .5;
        if (den == 0.)
          { if (beg[i] < 0 || beg[i] > dim)
              ab = 2.;
          }
        else
          { abg = - (beg[i] + .5) / den;
            aen = (dim - beg[i]) / den;
            if (abg > aen)
              { den = abg; abg = aen; aen = den; }
            if (abg > ab)
              ab = abg;
            if (aen < ae)
              ae = aen;
          }
      }

    bkm = beg[kmax];        //  then further refine to have integral start and end poinnts,
    dkm = val[kmax];        //    bgk & end, for dimension kmax
    if (dkm > 0.)
      { enk = (int64) (bkm + ae * dkm);
        bgk = (int64) ceil(bkm + ab * dkm);
      }
    else
      { bgk = (int64) (bkm + ab * dkm);
        enk = (int64) ceil(bkm + ae * dkm);
      }
    ab = (bgk - bkm) / dkm;
    ae = (enk - bkm) / dkm;

    if (ab > ae)
      goto exit_line;

    dkm = fabs(dkm);
    for (i = 0; i < ndims; i++)   // compute clipped start points and increments for every dimension
      { den = val[i];
        inc[i] = den / dkm;
        if (i == kmax)
          val[i] = (double) bgk;
        else
          val[i] = beg[i] + ab * den + .5;
      }
  }

  { int       k, step;
    int64     i;
    Indx_Type p;

    if (bgk <= enk)
      step = 1;
    else
      step = -1;
    enk += step;
    for (i = bgk; i != enk; i += step)  //  for each integer along dimension kmax,
      { k = ndims-1;                    //    paint the nearest pixel
        p = (Indx_Type) val[k];
        val[k] += inc[k];
        for (k--; k >= 0; k--)
          { p = p * dims[k] + ((int) val[k]);
            val[k] += inc[k];
          }
        painter(p,argp);
      }
  }

exit_line:
  if (ndims > 10)
    free(val);
  Free_Array(begp);
  Free_Array(endp);
}
