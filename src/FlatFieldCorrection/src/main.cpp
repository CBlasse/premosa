//
//  main.cpp
//
//  Created by Corinna Blasse on 06/15/15.
//

#include <iostream>
#include <math.h>
#include <chrono>


extern "C" {
#include "array.h"
#include "histogram.h"
#include "image.h"
#include "utilities.h"
}

// HelperFunction Library
#include "MylibValues.h"


////////////////////////////////////////////////////////////////////////////////
/*
 Flat field correction
 -> corrected image = (image / flat field) * rescaling factor
 
 @param image:            input image (8-bit or 16-bit)
 @param flatField:        flat field image  (32-bit)
 */

Array * FlatFieldCorrection (Array * image, Array * flatField) {
  
  Range_Bundle flatFieldRange;
  Array_Range(&flatFieldRange, flatField);
  double flatFieldMaximum = ceil(flatFieldRange.maxval.fval);
  
  if (image->type == UINT8_TYPE) {
    
    uint8 * imgVals = AUINT8(image);
    float32 * noiseVals = AFLOAT32(flatField);
    Array * img8 = Make_Array(PLAIN_KIND, UINT8_TYPE, 2, image->dims);
    uint8 * img8Vals = AUINT8(img8);
    
    for (int p=0; p<image->dims[0]*image->dims[1]; ++p) {
      double val = flatFieldMaximum * ((double)imgVals[p] / (double)noiseVals[p]);     // normalization
      if (val > 255) {                                                                // clipping of oversaturated values
        img8Vals[p] = 255;
      } else {
        img8Vals[p] = (uint8) val;
      }
    }
    
    return img8;
    
  } else if (image->type == UINT16_TYPE) {
    
    uint16 * imgVals = AUINT16(image);
    float32 * noiseVals = AFLOAT32(flatField);
    Array * img16 = Make_Array(PLAIN_KIND, UINT16_TYPE, 2, image->dims);
    uint16 * img16Vals = AUINT16(img16);
    
    for (int p=0; p<image->dims[0]*image->dims[1]; ++p) {
      double val = flatFieldMaximum * ((double)imgVals[p] / (double)noiseVals[p]);   // normalization
      if (val > 65535) {                                                            // clipping of oversaturated values
        img16Vals[p] = 65535;
      } else {
        img16Vals[p] = (uint16) val;
      }
    }
    
    return img16;
    
  } else {
    std::cout << "Flat field correction is not supported for that data type " << std::endl;
  }
  
  
  return image;
}

////////////////////////////////////////////////////////////////////////////////
/*
 Approximation of the flat field by filtering the image
 
 @param image:            input image (8-bit or 16-bit)
 */

Array * FlatFieldApproximation (Array * image) {
  
  Array * approximatedFlatField = Make_Array_With_Shape(PLAIN_KIND, FLOAT32_TYPE, Coord2(image->dims[1], image->dims[0]));
  Range_Bundle range;
  Array_Range(&range, image);
  
  Use_Reflective_Boundary();
  
  // frames defines the local neighborhood
  int filterRadius = image->dims[0]/5;
  if (filterRadius%2 == 0) {
    filterRadius++;
  }
  Frame * f = Make_Frame(image,Coord2(filterRadius,filterRadius),Coord2((filterRadius-1)/2,(filterRadius-1)/2));
  Histogram * h = Make_Histogram(UVAL,range.maxval.uval,MylibValues::ValU(1), MylibValues::ValU(0));
  Place_Frame(f,0);
  
  for (Indx_Type p = 0; p < image->size; p++){
    Empty_Histogram(h);
    Histagain_Array(h,f,0);
    Set_Array_Value(approximatedFlatField, Idx2CoordA(approximatedFlatField, p), MylibValues::ValF(Percentile2Bin(h,.5)));    // median is given at 50th percentile
    Move_Frame_Forward(f);
  }
  
  Kill_Histogram(h);
  Kill_Frame(f);
  
  return approximatedFlatField;
}




////////////////////////////////////////////////////////////////////////////////
/*
 ************ MAIN ************
 
 Function to read the image, perform a flat field correction and the write the image
 
 */


int main(int argc, char * argv[])
{
  auto begin = std::chrono::high_resolution_clock::now();
  
  static char *Spec[] = { "-in <file> [-ff <file>] -out <file>", NULL };
  Process_Arguments(argc,argv,Spec,1);
  
  // Input parameter
  Array * inputImage = Read_Image(Get_String_Arg("-in"), 0);
  if (inputImage->type != UINT16_TYPE and inputImage->type != UINT8_TYPE) {
    std::cout << "Flat field correction requires a 8-bit or 16-bit input image! " << std::endl;
    Free_Array(inputImage);
    return 0;
  }
  
  Array * flatfield;
  if (Is_Arg_Matched("-ff")){
    flatfield = Read_Image(Get_String_Arg("-ff"), 0);
    if (flatfield->type != UINT32_TYPE and flatfield->type != FLOAT32_TYPE) {
      std::cout << "Flat field correction requires a 32-bit flat field image (float or uint)! " << std::endl;
      Free_Array(flatfield);
      return 0;
    }
  } else {
    flatfield = FlatFieldApproximation(inputImage);     // Approximate the flat field if no flat field is given as input
  }
  
  // Flat field correction
  Array * output = FlatFieldCorrection(inputImage, flatfield);
  Write_Image(Get_String_Arg("-out"), output, DONT_PRESS);
  
  // Cleaning
  Free_Array(inputImage);
  Free_Array(flatfield);
  Free_Array(output);
  
  auto end = std::chrono::high_resolution_clock::now();
  auto dur = end - begin;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
  std::cout << "Runtime: " << ms << " ms"<< std::endl;
  
  
  return 0;
}

