//
//  FileHandling.h
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//

#ifndef FILEHANDLING_H
#define FILEHANDLING_H

#include <unordered_map>
#include <iostream>

#include <QDir>
#include <QString>


namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   struct including all output folders
   */
  
  struct OutputDirs {
    QDir rootDir;
    QDir flatFieldDir;
    QDir constrastAdjDir;
    QDir stitchedDir;
  };
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Function to create all output folders
   
   @param outputDirPath:  path to the parent output folder
   */
  
  OutputDirs MakeOutputDirs (QString outputDirPath);
  
}

#endif
