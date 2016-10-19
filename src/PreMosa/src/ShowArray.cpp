//
//  ShowArray.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>

#include "ShowArray.h"

extern "C" {
#include "image.h"
}

#include "ImageClass.h"
#include "Parameter.h"


using namespace HelperFunctions;
using namespace ImageClass;
using namespace Parameter;

void ShowArray(Image & img, ImgParameter & imgPara, char * name, Array * arr) 
{ 
  Dimn_Type x, y, z;
  Array    * a, * b;
    
  if (imgPara.grid == false){ 
    Write_Image((char *) name, arr, DONT_PRESS);
    return;
  }
  
  if (arr->type != UINT8_TYPE){ 
    b = Copy_Array(arr);
    Scale_Array_To_Range(b,ValU(0),ValU(255));
    Convert_Array_Inplace(b,b->kind,UINT8_TYPE,8,0);
  } else {
    b = arr;
  }
  
  if (b->kind != RGB_KIND){
    if (b == arr){
      a = Convert_Array(b,RGB_KIND,b->type,b->scale,0);
    } else {
      a = Convert_Array_Inplace(b,RGB_KIND,b->type,b->scale,0);
    }
  } else {
      a = b;
  }
  
  if (a->ndims == 4)
    for (z = 0; z < img.GetDepth(); z++){ 
      for (y = 0; y < img.GetGridY(); y++) {
        if (y % 50 == 0 && y != 0){
          Draw_Line(a,&YELLOW,Coord3(z,y*imgPara.radius,0),Coord3(z,y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
        } else if (y % 10 == 0 && y != 0){
          Draw_Line(a,&PURPLE,Coord3(z,y*imgPara.radius,0),Coord3(z,y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
        } else {
          Draw_Line(a,&CYAN,Coord3(z,y*imgPara.radius,0),Coord3(z,y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
        }
        for (x = 0; x < img.GetGridX(); x++){
          if (x % 50 == 0 && x != 0){
            Draw_Line(a,&YELLOW,Coord3(z,0,x*imgPara.radius),Coord3(z,img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
          } else if (x % 10 == 0 && x != 0){
            Draw_Line(a,&PURPLE,Coord3(z,0,x*imgPara.radius),Coord3(z,img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
          } else {
            Draw_Line(a,&CYAN,Coord3(z,0,x*imgPara.radius),Coord3(z,img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
          }
        }
      }
    }
  else{ 
    for (y = 0; y < img.GetGridY(); y++){
      if (y % 50 == 0 && y != 0){
        Draw_Line(a,&YELLOW,Coord2(y*imgPara.radius,0),Coord2(y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
      } else if (y % 10 == 0 && y != 0){
        Draw_Line(a,&PURPLE,Coord2(y*imgPara.radius,0),Coord2(y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
      } else{
        Draw_Line(a,&CYAN,Coord2(y*imgPara.radius,0),Coord2(y*imgPara.radius,img.GetGridX()*imgPara.radius-1));
      }
      for (x = 0; x < img.GetGridX(); x++){
        if (x % 50 == 0 && x != 0){
          Draw_Line(a,&YELLOW,Coord2(0,x*imgPara.radius),Coord2(img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
        } else if (x % 10 == 0 && x != 0){
          Draw_Line(a,&PURPLE,Coord2(0,x*imgPara.radius),Coord2(img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
        } else{
          Draw_Line(a,&CYAN,Coord2(0,x*imgPara.radius),Coord2(img.GetGridY()*imgPara.radius-1,x*imgPara.radius));
        }
      }
    }
  }
  
  Write_Image(name, a, DONT_PRESS);
  
  if (arr != a)
    Free_Array(a);
}


