//
//  ConversionTo8bit.cpp
//  CMakeProject
//
//  Created by blasse on 7/5/13.
//
//

#include "ConversionTo8bit.h"

#include <iostream>

extern "C" {
#include "image.h"
#include "utilities.h"
}


namespace Conversion {
  
  Array * ConversionTo8Bit (QFileInfo entryInfo, QFileInfoList entries){
    
    Array * timePointStack = Read_Image("/Users/blasse/Documents/PhD/Data/20111102/rawData_8bit_selection/Stacks_t015/20111102_Pos6_p000001t00000015c01.tif",0);
    
    foreach ( QFileInfo entryInfo, entries){
        
        QString path = entryInfo.absoluteFilePath();
        std::cout << qPrintable(QString("%1").arg(entryInfo.filePath()))<< std::endl;
      
        if (entryInfo.isDir()) {
          
          QDir timePointDir(path);
          
          QFileInfoList timePointEntries = timePointDir.entryInfoList(QDir::Files);
          
          foreach (QFileInfo timePointInfo, timePointEntries){
            
            std::cout << qPrintable(timePointInfo.filePath()) << std::endl;
            
            std::string name = qPrintable(timePointInfo.fileName());
            Array * input = Read_Image(const_cast<char * > (qPrintable(timePointInfo.filePath())),0);
            
            Array ** args = (Array **) Guarded_Malloc(sizeof(Array *)*((Size_Type) 1),Program_Name());
            
            args[0] = timePointStack;
            args[1] = input;
            
            Array * newA = Make_Array_From_Arrays(PLAIN_KIND, 2, args);
            
            std::cout << "Size " << newA->dims[0] << " " << newA->dims[1] << " " << (int)newA->dims[2] << std::endl;
            
            Write_Image("/Users/blasse/Desktop/Alle.tif", newA, DONT_PRESS);
          }
        }
    }
      

    return timePointStack;
    
  }
  
}