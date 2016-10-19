//
//  ConversionTo8bit.h
//  CMakeProject
//
//  Created by blasse on 7/5/13.
//
//

#ifndef CONVERSIONTO8BIT_H
#define CONVERSIONTO8BIT_H

#include <QDir>

extern "C"{
#include "array.h"
}



namespace Conversion {
  Array * ConversionTo8Bit (QFileInfo entryInfo, QFileInfoList entries);
  
}
#endif 
