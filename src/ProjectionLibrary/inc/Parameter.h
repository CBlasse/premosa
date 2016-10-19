//
//  Parameter.h
//
//  Created by Corinna Blasse on 20.08.12.
//

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>

namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parameter required for the projection algorithm
  */
  struct ParaSet{
    
    ParaSet()
    : verbose(true)
    , grid(false)
    , radius(20)
    , distance(2)
    , threshold(50)
    , accurate(false)
    , maxInterpolation(false)
    , printHeightMap(false)
    , printRealHeightMap(false)
    , layer(1)
    {
    };
    
    ParaSet (int argc, char * argv[], char * specifications []);
    ~ParaSet();
    
    bool verbose;
    bool grid;
    int radius;
    int distance;
    int threshold;
    bool accurate;
    bool maxInterpolation;
    bool printHeightMap;
    bool printRealHeightMap;
    bool version;
    bool range;
    int layer;
    
    std::string outputFolder;
  };

  
  struct ExtendedParaSet : public ParaSet {
    ExtendedParaSet()
    : distance1(0)
    , distance2(0)
    {
    };
    
    ExtendedParaSet (int argc, char * argv[], char * specifications []);
    ~ExtendedParaSet();
    
    int distance1;
    int distance2;
  };
}

#endif
