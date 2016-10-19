//
//  FileHandling.cpp
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//

#include "FileHandling.h"


using namespace std;

namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Function to create all output folders
   
   @param outputDirPath:  path to the parent output folder
   */
  
  OutputDirs MakeOutputDirs (QString outputDirPath) {
    
    OutputDirs dirs;
    
    dirs.rootDir  = QDir (outputDirPath);
    if (!dirs.rootDir.exists()) {
      dirs.rootDir.mkdir(outputDirPath);
    }
    
    dirs.flatFieldDir = QDir (outputDirPath + QDir::separator()+"FlatFieldCorrected_2D");
    if (!dirs.flatFieldDir.exists()) {
      dirs.flatFieldDir.mkdir(dirs.flatFieldDir.absolutePath());
    }
    
    dirs.constrastAdjDir = QDir (outputDirPath + QDir::separator()+"ContrastAdjusted_2D");
    if (!dirs.constrastAdjDir.exists()) {
      dirs.constrastAdjDir.mkdir(dirs.constrastAdjDir.absolutePath());
    }
    
    dirs.stitchedDir = QDir (outputDirPath + QDir::separator()+"Stitched_2D");
    if (!dirs.stitchedDir.exists()) {
      dirs.stitchedDir.mkdir(dirs.stitchedDir.absolutePath());
    }
    
    return dirs;
  }
}