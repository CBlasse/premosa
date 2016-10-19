//
//  RasterizedImage.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "RasterizedImage.h"

extern bool verbose;
extern int radius;

#include "MylibValues.h"
using namespace MylibValues;

namespace ProjectionMethod {

  /***********************
   * Constructor / Destructor
   **********************/
  
  RasterizedImage::RasterizedImage (Array * image){
    image_ =  std::make_shared<Array * >(image);
    realWidth_ = image->dims[0];
    realHeight_ = image->dims[1];
    width_ = image->dims[0];
    height_ = image->dims[1];
    depth_ = image->dims[2];
    
    type_ = image->type;
  }
  
//  RasterizedImage::RasterizedImage () {
//  }

   RasterizedImage::~RasterizedImage () {
     
//     Free_Array(*image_);
  }

  /***********************
   * Get functions
   **********************/
  
  Array * RasterizedImage::GetImage(){
    return *image_;
  }

  int RasterizedImage::GetWidth () const {
    return width_;
  }  
  
  int RasterizedImage::GetRealWidth () const {
    return realWidth_;
  }  
  
  int RasterizedImage::GetRealHeight () const {
    return realHeight_;
  }  
  
  int RasterizedImage::GetHeight () const {
    return height_;
  }  
  
  int RasterizedImage::GetDepth () const {
    return depth_;
  }

  int RasterizedImage::GetGridX () const {
    return gridX_;
  }
  
  int RasterizedImage::GetGridY () const {
    return gridY_;
  }
  
  int RasterizedImage::GetArea () const {
    return width_*height_;
  }
  
  Value_Type RasterizedImage::GetType () const {
    return type_;
  }

  /***********************
   * Set functions
   **********************/
  
  void RasterizedImage::SetImage(Array * image) {
    image_ = std::make_shared<Array *>(image);
//    image_(image);
  }
  
  void RasterizedImage::SetWidth(const int width){
    width_ = width;
  }
  
  void RasterizedImage::SetHeight(const int height){
    height_ = height;
  }
  
  void RasterizedImage::SetDepth(const int depth){
    depth_ = depth;
  }
  
  void RasterizedImage::SetGrid(const int radius){
    gridX_ = ceil ((double)realWidth_/(double)radius);
    gridY_ = ceil ((double)realHeight_/(double)radius);
    Use_Extend_Boundary();
    Pad_Array_Inplace(*image_, Coord3(0, 0, 0), Coord3(depth_, gridY_*radius, gridX_*radius));

    SetWidth((*image_)->dims[0]);
    SetHeight((*image_)->dims[1]);
  }
  
  void RasterizedImage::SetGrid_Left(const int radius){
    gridX_ = ceil ((double)realWidth_/(double)radius);
    gridY_ = ceil ((double)realHeight_/(double)radius);
    
    Use_Extend_Boundary();    
    
    Pad_Array_Inplace(*image_, Coord3(0, 0, gridX_*radius - realWidth_), Coord3(depth_, gridY_*radius, gridX_*radius));
    
    SetWidth((*image_)->dims[0]);
    SetHeight((*image_)->dims[1]);
  }
  
  /***********************
   * Other functions
   **********************/
  
  void RasterizedImage::ScaleTo8bit(){
    Scale_Array_To_Range((*image_),ValU(0),ValU(255));
    Convert_Array_Inplace((*image_),PLAIN_KIND,UINT8_TYPE,8,0);
  }
  
  void RasterizedImage::UpdateGrid(const double scale){
    gridX_ = scale*gridX_;
    gridY_ = scale*gridY_;
  }
  
}