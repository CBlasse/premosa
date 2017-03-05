//
//  main.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <dirent.h>
#include <iostream>
#include <string>
#include <chrono>
#include <sys/stat.h>


extern "C" {
#include "utilities.h"
#include "image.h"
#include "mylib.h"
}

#include "Parameter.h"
#include "ProjectionAlgorithm.h"


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
