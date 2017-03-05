//
//  main.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <dirent.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
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

bool StringHasSuffix(const std::string &str, const std::string &suffix)
{
  return str.size() >= suffix.size() &&
  str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void CheckOutputFolder(std::string outputPath) {
  
  if (StringHasSuffix(outputPath, "/")) {
    outputPath = outputPath.substr(0,outputPath.length()-1);
  }
  
  std::string path;
  
  if (StringHasSuffix(outputPath, ".tif")) {
    std::size_t found = outputPath.find_last_of("/");
    path = outputPath.substr(0, found);
  } else {
    path = outputPath;
  }
  
  DIR * outputDir;
  outputDir  = opendir (&(path)[0]);
  
  if (outputDir == NULL){
    mkdir(&path[0], 0777);
  }
}

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
  
  std::string outputString = Get_String_Arg("-out");
  // Verification of the output path
  CheckOutputFolder(outputString);
  if (!StringHasSuffix(outputString, ".tif")) {
    if (StringHasSuffix(outputString, "/")) {
      outputString += "Output";
    } else {
      outputString += "/Output";
    }
  } else {
    outputString = outputString.substr(0, outputString.length()-4);
  }
  
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
  projection.DrawProjection(outputString);
  projection.DrawInterpolatedHeightMap(outputString);
  
  
  Free_Array(heightMap);
  Free_Array(rasterizedHeightMap);
  
  return 0;
}
