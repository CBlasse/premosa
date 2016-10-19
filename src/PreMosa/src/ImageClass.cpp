//
//  ImageClass.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "ImageClass.h"

extern "C" {
#include "image.h"
}

#include "HelperFunctions.h"



extern bool verbose;
extern int radius;

using namespace HelperFunctions;

namespace ImageClass{

  /***********************
   * Constructor / Destructor
   **********************/
  
  Image::Image (Array * image) {
    image_ = image;
  
    realWidth_ = image->dims[0];
    realHeight_ = image->dims[1];
    width_ = image->dims[0];
    height_ = image->dims[1];
    depth_ = image->dims[2];
    
    type_ = image->type;
  }

  Image::~Image () {
    //Free_Array(image_);
  }

  /***********************
   * Get functions
   **********************/
  
  Array * Image::GetImage(){
    return image_;
  }

  int Image::GetWidth () const {
    return width_;
  }  
  
  int Image::GetRealWidth () const {
    return realWidth_;
  }  
  
  int Image::GetRealHeight () const {
    return realHeight_;
  }  
  
  int Image::GetHeight () const {
    return height_;
  }  
  
  int Image::GetDepth () const {
    return depth_;
  }

  int Image::GetGridX () const {
    return gridX_;
  }
  
  int Image::GetGridY () const {
    return gridY_;
  }
  
  int Image::GetArea () const {
    return width_*height_;
  }
  
  Value_Type Image::GetType () const {
    return type_;
  }

  /***********************
   * Set functions
   **********************/
  
  void Image::SetImage(Array * image){
    image_ = image;
  }
  
  void Image::SetWidth(const int width){
    width_ = width;
  }
  
  void Image::SetHeight(const int height){
    height_ = height;
  }
  
  void Image::SetDepth(const int depth){
    depth_ = depth;
  }
  
  void Image::SetGrid(const int radius){
    gridX_ = ceil ((double)realWidth_/(double)radius);
    gridY_ = ceil ((double)realHeight_/(double)radius);
        
    Use_Extend_Boundary();
    Pad_Array_Inplace(image_, Coord3(0, 0, 0), Coord3(depth_, gridY_*radius, gridX_*radius));
    SetWidth(image_->dims[0]);
    SetHeight(image_->dims[1]);
  }
  
  void Image::SetGrid_Left(const int radius){
    gridX_ = ceil ((double)realWidth_/(double)radius);
    gridY_ = ceil ((double)realHeight_/(double)radius);
    
    Use_Extend_Boundary();    
    
    Pad_Array_Inplace(image_, Coord3(0, 0, gridX_*radius - realWidth_), Coord3(depth_, gridY_*radius, gridX_*radius));
    
    SetWidth(image_->dims[0]);
    SetHeight(image_->dims[1]);
  }
  
  /***********************
   * Other functions
   **********************/
  
  void Image::ScaleTo8bit(){
    Scale_Array_To_Range(image_,ValU(0),ValU(255));
    Convert_Array_Inplace(image_,PLAIN_KIND,UINT8_TYPE,8,0);
  }
  
  void Image::UpdateGrid(const int scale){    
    gridX_ = scale*gridX_;
    gridY_ = scale*gridY_;
  }
  
}