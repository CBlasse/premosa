//
//  RasterizedImage.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef RASTERIZEDIMAGE_H
#define RASTERIZEDIMAGE_H


#include <cmath>
#include <memory>

extern "C" {
  #include "array.h"
  #include "image.h"
}


namespace ProjectionMethod {

class RasterizedImage
{
  public:    
    RasterizedImage (Array * image);
//    RasterizedImage ();
    ~RasterizedImage ();
  
    Array * GetImage ();
    int GetWidth () const;
    int GetRealWidth () const;
    int GetRealHeight () const;
    int GetHeight () const;
    int GetDepth () const;
    int GetGridX () const;
    int GetGridY () const;
    int GetArea () const;
    Value_Type GetType () const; 
  
    void SetImage (Array * image);
    void SetWidth (const int width);
    void SetHeight (const int height);
    void SetDepth (const int depth);
    void SetGrid (const int radius);
    void SetGrid_Left (const int radius);
  
    void ScaleTo8bit ();
  
    void UpdateGrid(const double scale);
      
    
  private:
    std::shared_ptr<Array *>image_;
    int realWidth_;
    int realHeight_;
    int width_;
    int height_;
    int depth_;
    int gridX_;
    int gridY_;
  
    Value_Type type_;
  };
  
   
}



#endif
