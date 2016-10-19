//
//  Interpolation.h
//
//  Created by Corinna Blasse on 27.07.12.
//

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

extern "C" {
  #include "array.h"
}

#include "RasterizedImage.h"
#include "Parameter.h"

namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computes an interpolation of the height map to provide smoother transitions between grid elements
   -> Obtained array features dimension, which are one more than the height map
   
   @param img:      input image stack
   @param imgPara:  input parameter
   @param heigthMap:height map
   */
  Array * InterpolateCorners(RasterizedImage & img, ParaSet & imgPara, Array *heightMap);
  
  
  struct ResultArrays {
    Array * interpolatedHeightMap;
    Array * projection;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (8bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  ResultArrays InterpolatePlanes_8bit(RasterizedImage & img,  ParaSet & imgPara, Array *maxim, Array *corners);
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (16bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  ResultArrays InterpolatePlanes_16bit(RasterizedImage & img,  ParaSet & imgPara, Array *maxim, Array *corners);
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map. This projection does not restrict the projection to the computed z-section, but projects the maximal intensity from the computed z-section and its neighboring z-sections
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (8bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  ResultArrays InterpolatePlanes_8bit_MaxInterpolation(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners);
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map. This projection does not restrict the projection to the computed z-section, but projects the maximal intensity from the computed z-section and its neighboring z-sections
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (16bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  ResultArrays InterpolatePlanes_16bit_MaxInterpolation(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners);
}

#endif
