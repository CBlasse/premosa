////
////  Analysis.cpp
////  CMakeProject
////
////  Created by Corinna Blasse on 30.08.12.
////  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
////
//
//#include <iostream>
//
//#include "Analysis.h"
//
//
//using namespace ImageClass;
//
//
//namespace HeightMap {
//  
//  void GetLevelRange(Array * level){
//    
//    uint32 * lev = AUINT32(level);
//    uint32 maxVal;
//    uint32 minVal = 100;
//    
//    for (int i=0; i<level->size; i++) {
//      if (lev[i] > maxVal) {
//        maxVal = lev[i];
//      } else if (lev[i] < minVal){
//        minVal = lev[i];
//      }
//    }
//    
//    std::cout << "Max = " << maxVal << "\nMin = " << minVal << std::endl;
//  }
//}
//
