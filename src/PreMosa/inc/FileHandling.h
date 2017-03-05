//
//  FileHandling.h
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//

#ifndef FILEHANDLING_H
#define FILEHANDLING_H

#include <dirent.h>
#include <string>


namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   struct including all output folders
   */
  
  struct OutputDirs {
    DIR * rootDir;
    DIR * flatFieldDir;
    DIR * constrastAdjDir;
    DIR * stitchedDir;
    std::string rootPath;
    std::string flatFieldPath;
    std::string contrastAdjPath;
    std::string stitchedPath;
  };
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Function to create all output folders
   
   @param outputDirPath:  path to the parent output folder
   */
  
  OutputDirs MakeOutputDirs (std::string outputDirPath);
  
}

#endif
