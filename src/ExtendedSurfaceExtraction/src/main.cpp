//
//  main.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>


extern "C" {
#include "utilities.h"
#include "histogram.h"
#include "image.h"
#include "mylib.h"
}

#include "ExtendedProjectionAlgorithm.h"
#include "Parameter.h"
#include "VectorFunctions.h"

#include "ImageClass.h"
#include "ZProjection.h"
#include "Filter.h"
#include "ZPlaneSelection.h"
#include "Interpolation.h"
#include "ShowArray.h"
#include "Refinement.h"
#include "Analysis.h"



#define DEVELOP

using namespace ProjectionMethod;


//Array * Project_8bit_UsingHeightMap (Image img, ParaSet imgPara, Array * median) {
//  
//  Array * maxim = MaxProjection(img, imgPara);
//  Array * index = MaxIndices(img, imgPara, maxim);
//  
//  // ---------- version 1.0 -------------
//  
//  // 1. Select reference height map with a big window (radius)
//  //      Selecting the z planes using the occurrences of the brightest pixel
//  
//  Array * level = PlaneSelection(img, imgPara, maxim, index);
//  Free_Array(index);
//  
//  
//  
//  // 2. Decomposing the compuational window based on the quad tree principle
//  //    new window size: 0.5 * radius
//  //    Smoothing of the resulting height map using local median filtering
//  
//  imgPara.radius = imgPara.radius/2;
//  img.UpdateGrid(2);
//  
//  
//  // 3. Select the final height map using a subset of the image stack
//  //      defined by the reference height map and a distance parameter
//  
//  
//  // ------ method 1 ---- variance -------------
//  
//  Array * substack = Substack_UsingMask(img, imgPara, median, true);
//  Free_Array(median);
//  level = SelectPlanes_Variance(img, imgPara, substack);
//  Free_Array(substack);
//  
//  //median = Copy_Array(level);
//  median = ExtremeLocalMedianFilter(img, imgPara, level, img.GetDepth());
//  
//  if (imgPara.printHeightMap) {
//    std::string outputpath_hm = std::string(Get_String_Arg("out")) + "_heightMap.tif";
//    Write_Image((char *)outputpath_hm.data(),median,DONT_PRESS);
//  }
//  
//  
//  Array * corners = InterpolateCorners(img, imgPara, median);
//  Array * result = InterpolatePlanes_8bit(img, imgPara, maxim, corners);
//  
//  if ((img.GetWidth()-img.GetRealWidth() > 0) || (img.GetHeight() - img.GetRealHeight() > 0)) {
//    Clip_Array_Inplace(result, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
//    //int diffGridX = img.GetWidth()-img.GetRealWidth();
//    //int diffGridY = img.GetHeight() - img.GetRealHeight();
//    
//    //Clip_Array_Inplace(result, Coord2(diffGridY, diffGridX), Coord2(img.GetRealHeight()+diffGridY-1, img.GetRealWidth()+diffGridX-1));
//  }
//  
//  
//  Free_Array(level);
//  Free_Array(maxim);
//  Free_Array(corners);
//  
//  return result;
//}
//
//Array * Project_16bit_UsingHeightMap (Image img, ParaSet imgPara, Array * median) {
//  
//  Array * maxim = MaxProjection_16bit(img, imgPara);
//  
//  imgPara.radius = imgPara.radius/2;
//  img.UpdateGrid(2);
//  
//  
//  Array * substack = Substack_UsingMask_16bit(img, imgPara, median, true);
//  Array * level = SelectPlanes_Variance_16bit(img, imgPara, substack);
//  Free_Array(substack);
//  
//  median = Copy_Array(level);
//  
//  if (imgPara.printHeightMap) {
//    std::string outputpath_hm = std::string(Get_String_Arg("out")) + "_heightMap.tif";
//    Write_Image((char *)outputpath_hm.data(),median,DONT_PRESS);
//  }
//  
//  
//  Array * corners = InterpolateCorners(img, imgPara, median);
//  Array * result = InterpolatePlanes_16bit(img, imgPara, maxim, corners);
//  
//  if ((img.GetWidth()-img.GetRealWidth() > 0) || (img.GetHeight() - img.GetRealHeight() > 0)) {
//    Clip_Array_Inplace(result, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
//    //int diffGridX = img.GetWidth()-img.GetRealWidth();
//    //int diffGridY = img.GetHeight() - img.GetRealHeight();
//    
//    //Clip_Array_Inplace(result, Coord2(diffGridY, diffGridX), Coord2(img.GetRealHeight()+diffGridY-1, img.GetRealWidth()+diffGridX-1));
//  }
//  
//  
//  Free_Array(median);
//  Free_Array(level);
//  Free_Array(maxim);
//  Free_Array(corners);
//  
//  return result;
//}
//
//Array * Fast_Project_16bit_UsingHeightMap (Image img, ParaSet imgPara, Array * median) {
//  
//  Array * img_16bit = Copy_Array(img.GetImage());
//  Array * maxim_16bit = MaxProjection_16bit(img, imgPara);
//  
//  // Works on 8-bit to figure out the height map
//  img.ScaleTo8bit();
//  
//  Array * maxim = MaxProjection(img, imgPara);
//  Array * index = MaxIndices(img, imgPara, maxim);
//  
//  imgPara.radius = imgPara.radius/2;
//  img.UpdateGrid(2);
//  
//  Array * substack = Substack_UsingMask(img, imgPara, median, true);
//  Free_Array(median);
//  Array * level = SelectPlanes_Variance(img, imgPara, substack);
//  Free_Array(substack);
//  
//  median = Copy_Array(level);
//  
//  // Projections happens on 16-bit
//  img.SetImage(img_16bit);
//  Array * corners = InterpolateCorners(img, imgPara, median);
//  Array * result = InterpolatePlanes_16bit(img, imgPara, maxim, corners);
//  
//  
//  if ((img.GetWidth()-img.GetRealWidth() > 0) || (img.GetHeight() - img.GetRealHeight() > 0)) {
//    Clip_Array_Inplace(result, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
//  }
//  
//  
//  Free_Array(median);
//  Free_Array(level);
//  Free_Array(maxim);
//  Free_Array(corners);
//  
//  return result;
//  
//}
//

Array * RasterizeHeightMap (Array * heightMap, int gridSize, int imgWidth, int imgHeight, int imgDepth) {
  
  int gridX = imgWidth / gridSize;
  int gridY = imgHeight / gridSize;

  Use_Reflective_Boundary();
  Pad_Array_Inplace(heightMap, Coord2(0, 0), Coord2(imgHeight, imgWidth));
  uint32 * hmVals = AUINT32(heightMap);
  
  Array * rasterizedHeightMap = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(gridY, gridX));
  uint32 * rasVals = AUINT32(rasterizedHeightMap);
  
  Frame * frame = Make_Frame(heightMap, Coord2(gridSize, gridSize), Coord2(0, 0));
  Histogram * h = Make_Histogram(UVAL,imgDepth,MylibValues::ValU(1),MylibValues::ValU(0));
  
  for (int y=0; y<gridY; y++) {
    for (int x=0; x<gridX; x++) {
      
      std::vector<int> zLevels;
      
      for (int r=0; r<gridSize; r++) {
        for (int s=0; s<gridSize; s++) {
          Indx_Type p = (y*gridSize + r)*imgWidth + (x*gridSize + s);
          zLevels.push_back(static_cast<int>(hmVals[p]));
        }
      }
      
      rasVals[y*gridX + x] = VectorFunctions::VectorMedian(zLevels);
    }
  }
  
  Free_Frame(frame);
  Free_Histogram(h);
  
  return rasterizedHeightMap;
}

/****************************************************************************************
 *                                                                                      *
 *  MAIN LEVEL                                                                          *
 *                                                                                      *
 ****************************************************************************************/

int main(int argc, char * argv[])
{
  static char *Spec[] = { "[-r <int(20)>] [-t <int(50)>] [-d1 <int(0)>] [-d2 <int(0)>] [-hmr] [-hmd] [-v] -in <file> -hmFile <file> -out <file>", NULL };
  
//  Process_Arguments(argc,argv,Spec,1);
  
  ProjectionMethod::ExtendedParaSet parameter (argc,argv,Spec);

  ExtendedProjectionAlgorithm projection (Get_String_Arg("-in"));
  projection.SetAllParameters(parameter);
  
  Array * heightMap = Read_Image(Get_String_Arg("-hmFile"),0);
  if (heightMap->dims[0] != projection.GetImage().GetRealWidth() or heightMap->dims[1] != projection.GetImage().GetRealHeight()) {
    std::cout << "Error: The height map needs to have the same size as the image." << std::endl;
    return 0;
  }
  if (heightMap->type != UINT32_TYPE) {
    
    Array * heightMap32 = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(heightMap->dims[1], heightMap->dims[0]));
    uint32 * newVals = AUINT32(heightMap32);
    
    if (heightMap->type == UINT8_TYPE) {
      uint8 * oldVals = AUINT8(heightMap);
      for (int p=0; p<heightMap->dims[0]*heightMap->dims[1]; p++) {
        newVals[p] = static_cast<uint32>(oldVals[p]);
      }
    } else if (heightMap->type == UINT16_TYPE) {
      uint16 * oldVals = AUINT16(heightMap);
      for (int p=0; p<heightMap->dims[0]*heightMap->dims[1]; p++) {
        newVals[p] = static_cast<uint32>(oldVals[p]);
      }
    } else if (heightMap->type == INT32_TYPE) {
      int32 * oldVals = AINT32(heightMap);
      for (int p=0; p<heightMap->dims[0]*heightMap->dims[1]; p++) {
        newVals[p] = static_cast<uint32>(oldVals[p]);
      }
    } else if (heightMap->type == FLOAT32_TYPE) {
      float32 * oldVals = AFLOAT32(heightMap);
      for (int p=0; p<heightMap->dims[0]*heightMap->dims[1]; p++) {
        newVals[p] = static_cast<uint32>(oldVals[p]);
      }
    }
    
    Free_Array(heightMap);
    heightMap = Copy_Array(heightMap32);
    Free_Array(heightMap32);
  }
  
  Array * rasterizedHeightMap = RasterizeHeightMap(heightMap, parameter.radius, projection.GetImage().GetWidth(), projection.GetImage().GetHeight(), projection.GetImage().GetDepth());
  projection.ProjectSurfaceUsingHeightMap(rasterizedHeightMap);
  projection.DrawProjection(Get_String_Arg("-out"));
  projection.DrawInterpolatedHeightMap(Get_String_Arg("-out"));
  
  
//  if (heightMap->type == FLOAT32_TYPE) {
//    Array * median2 = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(median->dims[1], median->dims[0]));
//    float32 * medVals = AFLOAT32(median);
//    uint32 * med2Vals = AUINT32(median2);
//    for (int p=0; p<median->size; p++) {
//      med2Vals[p] = static_cast<uint32>(medVals[p]);
//    }
//    median = median2;
//  }
  
//   Array * result;
//  if (img.GetType() == UINT8_TYPE) {
//    result = Project_8bit_UsingHeightMap(img, imgPara, median);
//  } else if (img.GetType() == UINT16_TYPE) {
//    result = Fast_Project_16bit_UsingHeightMap(img, imgPara, median);
//    //        //result = Project_16bit_UsingHeightMap(img, imgPara, median);  // TODO:
//    //      img.ScaleTo8bit();
//    //      result = Project_8bit_UsingHeightMap(img, imgPara, median);
//    //
//    //      Array * result16Bit = Make_Array(PLAIN_KIND, UINT16_TYPE, 2, result->dims);
//    //      uint8 * resultVals = AUINT8(result);
//    //      uint16 * result16Vals = AUINT16(result16Bit);
//    //
//    //      for (Indx_Type p=0; p<result->dims[0]*result->dims[1]; p++) {
//    //        result16Vals[p] = (uint16) resultVals[p];
//    //      }
//    //
//    //      Scale_Array_To_Range(result16Bit, HelperFunctions::ValU(0), HelperFunctions::ValU(65535));
//    //      Free_Array(result);
//    //      result = result16Bit;
//    
//    
//  } else {
//    img.ScaleTo8bit();
//    result = Project_8bit_UsingHeightMap(img, imgPara, median);
//    
//  }
  
//  std::string outputpath_p = std::string(Get_String_Arg("out"));
//  if ((int)outputpath_p.find(".tif") < 0) {
//    outputpath_p += "_projection.tif";
//  }
//  Write_Image((char *)outputpath_p.data(),result,DONT_PRESS);
//  
//  
//  Free_Array(img.GetImage());
  Free_Array(heightMap);
  Free_Array(rasterizedHeightMap);
//  Free_Array(result);
//  img.~Image();
//  
//  if (imgPara.verbose) {
//    std::cout << "Array usage: " << Array_Usage() << std::endl;
//  }
  
  return 0;
}
