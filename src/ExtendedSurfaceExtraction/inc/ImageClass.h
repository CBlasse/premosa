////
////  ImageClass.h
////  proj_temp_cpp
////
////  Created by Corinna Blasse on 25.07.12.
////  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
////
//
//#ifndef IMAGECLASS_H
//#define IMAGECLASS_H
//
//
//#include <cmath>
//
//extern "C" {
//  #include "array.h"
//  #include "image.h"
//}
//
//#include "HelperFunctions.h"
//#include "ImageClass.h"
//
//namespace ImageClass {
//
//class Image
//{
//  public:    
//    Image (Array * image);
//    ~Image ();
//    
//    Array * GetImage ();
//    int GetWidth () const;
//    int GetRealWidth () const;
//    int GetRealHeight () const;
//    int GetHeight () const;
//    int GetDepth () const;
//    int GetGridX () const;
//    int GetGridY () const;
//    int GetArea () const;
//    Value_Type GetType () const; 
//  
//    void SetImage (Array * image);
//    void SetWidth (const int width);
//    void SetHeight (const int height);
//    void SetDepth (const int depth);
//    void SetGrid (const int radius);
//    void SetGrid_Left (const int radius);
//  
//    void ScaleTo8bit ();
//  
//    void UpdateGrid(const int scale);
//      
//    
//  private:
//    Array * image_;
//    int realWidth_;
//    int realHeight_;
//    int width_;
//    int height_;
//    int depth_;
//    int gridX_;
//    int gridY_;
//  
//    Value_Type type_;
//  };
//  
//   
//}
//
//
//
//#endif
