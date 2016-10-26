//
//  main.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <chrono>


extern "C" {
#include "utilities.h"
#include "image.h"
#include "mylib.h"
}

#include "Parameter.h"
#include "ProjectionAlgorithm.h"


/****************************************************************************************
 *                                                                                      *
 *  MAIN LEVEL                                                                          *
 *                                                                                      *
 ****************************************************************************************/

int main(int argc, char * argv[])
{
  auto begin = std::chrono::high_resolution_clock::now();
  
  static char *Spec[] = { "[-r <int(20)>] [-d <int(2)>] [-l <int(1)>] [-t <int(100)>] [-hmd] [-hmr] [-mi] [-v] -in <inputFile> -out <outputFile>", NULL };
  
  ProjectionMethod::ParaSet parameter (argc,argv,Spec);
  
  std::string outputString = Get_String_Arg("-out");
  
  ProjectionMethod::ProjectionAlgorithm imgProjection (Get_String_Arg("-in"));
  imgProjection.SetParameters (parameter);
  imgProjection.ProjectImageStack();
  imgProjection.DrawProjection(outputString);
  
  if (parameter.printHeightMap == true) {
    imgProjection.DrawDownsampledHeightMap(outputString);
  }
  if (parameter.printRealHeightMap == true) {
    imgProjection.DrawInterpolatedHeightMap(outputString);
  }
  
  if (parameter.verbose == true) {
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "Runtime: " << ms << " ms"<< std::endl;
  }
  
  return 0;
}
