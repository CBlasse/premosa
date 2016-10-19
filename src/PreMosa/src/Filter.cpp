//
//  Filter.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 26.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#include "Filter.h"

#include <iostream>

extern "C" {
#include "image.h"
#include "histogram.h"
}

#include "Develop.h"
#include "HelperFunctions.h"
#include "ImageClass.h"
#include "Parameter.h"
#include "ShowArray.h"
#include "Develop.h"

using namespace HelperFunctions;
using namespace ImageClass;
using namespace Parameter;


namespace Filter {
  
  
  /***********************
   * Apply a radius 'rad' median filter to level, storing the result in median
   *    NOT USED
   **********************/
  
  Array * MedianFilter(Image & img, ImgParameter & imgPara, Array *original, int depth, int rad){ 
    
    Array  *median = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
#ifdef DEVELOP
    Array  *variance = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(original->dims[1], original->dims[0]));
    uint32 *var      = AUINT32(variance);
#endif
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*rad+1,2*rad+1),Coord2(rad,rad));
    Histogram * h = Make_Histogram(UVAL,depth,ValU(1),ValU(0));
    Place_Frame(f,0);
    
    for (Indx_Type p = 0; p < median->size; p++){ 
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      med[p] = Percentile2Bin(h,.5);
#ifdef DEVELOP
      var[p] = 10*Histogram_Sigma(h);
#endif
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "median.tif", median);
    ShowArray(img, imgPara, "variance.tif", variance);
    Free_Array(variance);
#endif
    
    if (imgPara.verbose){ 
      std::cout << "Planes filtered" << std::endl; 
    }
    
    return (median);
  }
  
  
  /***********************
   * Apply a radius 'filterRad' median filter to levels that differ from the surrounding
   *    storing the result in median
   **********************/
  
  Array * LocalMedianFilter(Image & img, ImgParameter & imgPara, Array * original, int depth, int filterRad){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
    Use_Extend_Boundary();
        
    Frame * f = Make_Frame(original,Coord2(2*filterRad+1,2*filterRad+1),Coord2(filterRad, filterRad));
    Histogram * h = Make_Histogram(UVAL,depth,ValU(1),ValU(0));
    
    
    Place_Frame(f,0);
    for (Indx_Type p = 0; p < median->size; p++){ 
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      if (abs(orig[p] - Percentile2Bin(h,.5)) > 1) {
        med[p] = Percentile2Bin(h,.5);
      } else {
        med[p] = orig[p];
      }
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
    return median;
  }
  
  
  /***********************
   * Apply a radius 'filterRad' median filter to levels that differ from the surrounding
   *    storing the result in median
   **********************/
  
  Array * LocalMedianFilter_InclVar(Image & img, ImgParameter & imgPara, Array * original, int depth, int filterRad){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*filterRad+1,2*filterRad+1),Coord2(filterRad, filterRad));
    Histogram * h = Make_Histogram(UVAL,depth,ValU(1),ValU(0));
    
    
    Place_Frame(f,0);
    for (Indx_Type p = 0; p < median->size; p++){ 
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      if (abs(orig[p] - Percentile2Bin(h,.5)) > 1) {
        med[p] = Percentile2Bin(h,.5);
      } else {
        med[p] = orig[p];
      }
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
    return median;
  }

  
  /***********************
   * Iteration of the median filter using different radii
   **********************/
  
  Array * Filtering (Image & img, ImgParameter & imgPara, Array * original, int depth, int filterSizes []) {
    Array * median = Copy_Array(original);
    
    for (int i=0; i < sizeof(filterSizes)/sizeof(*filterSizes); i++) {
      Free_Array(median);
      median = LocalMedianFilter(img, imgPara, original, depth, filterSizes[i]);
      Free_Array(original);
      original = Copy_Array(median);
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "median.tif",median);
#endif
    if (imgPara.verbose){
      std::cout << "Planes filtered" << std::endl;
    }
    
    return median;
  }

  
}