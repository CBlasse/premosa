//
//  FileHandling.cpp
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//
#include "FileHandling.h"

#include <unordered_map>
#include <iostream>
#include <sys/stat.h>


using namespace std;

namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Function to create all output folders
   
   @param outputDirPath:  path to the parent output folder
   */
  
  OutputDirs MakeOutputDirs (std::string outputDirPath) {
    
    OutputDirs dirs;
    
    dirs.rootPath = outputDirPath;
    dirs.rootDir  = opendir (&outputDirPath[0]);
    if (dirs.rootDir == NULL) {
      mkdir(&outputDirPath[0], 0777);
      dirs.rootDir  = opendir (&outputDirPath[0]);
    }
    
    dirs.flatFieldPath = outputDirPath + "/"+"FlatFieldCorrected_2D";
    dirs.flatFieldDir  = opendir (&dirs.flatFieldPath[0]);
    if (dirs.flatFieldDir == NULL) {
      mkdir(&dirs.flatFieldPath[0], 0777);
      dirs.flatFieldDir  = opendir (&dirs.flatFieldPath[0]);
    }
    
    dirs.contrastAdjPath = outputDirPath + "/"+"ContrastAdjusted_2D";
    dirs.constrastAdjDir  = opendir (&dirs.contrastAdjPath[0]);
    if (dirs.constrastAdjDir == NULL) {
      mkdir(&dirs.contrastAdjPath[0], 0777);
      dirs.constrastAdjDir  = opendir (&dirs.contrastAdjPath[0]);
    }

    dirs.stitchedPath = outputDirPath + "/"+"Stitched_2D";
    dirs.stitchedDir  = opendir (&dirs.stitchedPath[0]);
    if (dirs.stitchedDir == NULL) {
      mkdir(&dirs.stitchedPath[0], 0777);
      dirs.stitchedDir  = opendir (&dirs.stitchedPath[0]);
    }

    return dirs;
  }
}