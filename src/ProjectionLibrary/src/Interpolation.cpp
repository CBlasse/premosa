
//
//  Interpolation.cpp
//
//  Created by Corinna Blasse on 27.07.12.
//
#include "Interpolation.h"

#include <iostream>

extern "C" {
#include "array.h"
#include "image.h"
#include "histogram.h"
#include "utilities.h"
}

#include "ShowArray.h"


namespace ProjectionMethod {
  
  static Color_Bundle PURPLE = { MAX_PIX, ValU(150), ValU(0), ValU(150) };

  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computes an interpolation of the height map to provide smoother transitions between grid elements
   -> Obtained array features dimension, which are one more than the height map
   
   @param img:      input image
   @param imgPara:  input parameter
   @param heigthMap:height map
   */
   
  Array *InterpolateCorners(RasterizedImage & img, ParaSet & imgPara, Array *heightMap){

    Indx_Type p, q;
    
    int plusX = img.GetGridX() + 1;
    int plusY = img.GetGridY() + 1;
    
    Array  *corners = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(plusY, plusX));
    uint32 *map     = AUINT32(corners);
    uint32 *med     = AUINT32(heightMap);
    
    for (Dimn_Type y = 1; y < img.GetGridY(); y++){ 
      p = y*plusX + 1;
      q = y*img.GetGridX() + 1;
      for (Dimn_Type x = 1; x < img.GetGridX(); x++, q++){
        map[p++] = med[q] + med[q-1] + med[q-img.GetGridX()] + med[q-plusX];
      }
    }
    
    p = img.GetGridY()*plusX+1;
    q = (img.GetGridY()-1)*img.GetGridX()+1;
    
    for (Dimn_Type x = 1; x < img.GetGridX(); x++, q++){
      map[x] = (med[x] + med[x-1]) << 1;
      map[p++] = (med[q] + med[q-1]) << 1;
    }
    
    p = plusX;
    q = img.GetGridX();
    
    for (Dimn_Type y = 1; y < img.GetGridY(); y++){
      map[p] = (med[q] + med[q-img.GetGridX()]) << 1;
      map[p+img.GetGridX()] = (med[q+img.GetGridX()-1] + med[q-1]) << 1;
      p += plusX;
      q += img.GetGridX();
    }
    
    map[0] = med[0] << 2;
    map[img.GetGridY()*plusX] = med[(img.GetGridY()-1)*img.GetGridX()] << 2;
    map[img.GetGridX()] = med[img.GetGridX()-1] << 2;
    map[img.GetGridY()*plusX+img.GetGridX()] = med[img.GetGridY()*img.GetGridX()-1] << 2;
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "corners.tif",corners);
#endif

    if (imgPara.verbose){
      std::cout << "Plane corners ready" << std::endl;
    }
    
    return (corners);
  }

  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (8bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  
  ResultArrays InterpolatePlanes_8bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners){
    
    ResultArrays results;
    results.projection = Make_Array(PLAIN_KIND,UINT8_TYPE,2, img.GetImage()->dims);
    uint8  *res    = AUINT8(results.projection);
    results.interpolatedHeightMap = Make_Array(PLAIN_KIND,UINT32_TYPE,2, img.GetImage()->dims);
    uint32 * hmVals = AUINT32(results.interpolatedHeightMap);
    
    uint32 *map    = AUINT32(corners);
    uint8  *val    = AUINT8(img.GetImage());
    uint8  *wal    = val + img.GetArea();
    
    int plusX = img.GetGridX() + 1;
    

#ifdef DEVELOP
    Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
    Array *residue = Make_Array(PLAIN_KIND,UINT8_TYPE,2,img.GetImage()->dims);
    
    uint8 *grn = AUINT8(marked) + img.GetDepth()*img.GetArea();
    uint8 *blu = grn + img.GetDepth()*img.GetArea();
    uint8 *rem = AUINT8(residue);
    uint8 *max = AUINT8(maxim);
#endif
    
    int n = imgPara.radius-1;
    
    for (Dimn_Type y = 0; y < img.GetGridY(); y++){
      
      uint32 *dep = map + y*plusX;
      
      for (Dimn_Type x = 0; x < img.GetGridX(); x++){
        double  ul = dep[0] / 4.;
        double  ur = dep[1] / 4.;
        double  ll = dep[plusX] / 4.;
        double  lr = dep[plusX+1] / 4.;
        
        for (int r = 0; r <= n; r++){       // linaer interpolation of the z-sections
          double ru = ul*(n-r)/n + ur*r/n;
          double rl = ll*(n-r)/n + lr*r/n;
          
          for (int s = 0; s <= n; s++){
            double    wt = ru*(n-s)/n + rl*s/n;
            int       wi = (int) wt;
            double    wf = wt-wi;
            Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
            Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
            
            if (wi+1 < img.GetDepth()) {
              res[p] = val[q] * (1.-wf) + wal[q] * wf;
              hmVals[p] = wi * (1.-wf) + wi * wf;
            }else {
              res[p] = val[q];
              hmVals[p] = wi;
            }

#ifdef DEVELOP
            rem[p] = max[p] - res[p];
            if (wf < 1.)
            { grn[q] *= (1.-wf)*(1.-wf);
              blu[q] *= wf*wf;
            }
            if (wf > 0.)
            { grn[q+img.GetArea()] *= wf*wf;
              blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
            }
#endif
          }
        }
        dep += 1;
      }
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "marked.tif",marked);
    ShowArray(img, imgPara, "residue.tif",residue);
    Free_Array(residue);
    Free_Array(marked);
#endif
    
    Clip_Array_Inplace(results.interpolatedHeightMap, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
    Clip_Array_Inplace(results.projection, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));


    if (imgPara.verbose){
      std::cout << "Result interpolated" << std::endl;
    }
    
    return (results);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (16bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  
  ResultArrays InterpolatePlanes_16bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners){

    ResultArrays results;
    results.projection = Make_Array(PLAIN_KIND,UINT16_TYPE,2, img.GetImage()->dims);
    uint16  *res    = AUINT16(results.projection);
    results.interpolatedHeightMap = Make_Array(PLAIN_KIND,UINT32_TYPE,2, img.GetImage()->dims);
    uint32 * hmVals = AUINT32(results.interpolatedHeightMap);
    
    uint32 *map    = AUINT32(corners);
    uint16  *val    = AUINT16(img.GetImage());
    uint16  *wal    = val + img.GetArea();
    
    int plusX = img.GetGridX() + 1;
    
#ifdef DEVELOP
    Array *residue = Make_Array(PLAIN_KIND,UINT16_TYPE,2,img.GetImage()->dims);
    uint16 *rem = AUINT16(residue);
    uint16 *max = AUINT16(maxim);
#endif
    
    int n = imgPara.radius-1;
    
    for (Dimn_Type y = 0; y < img.GetGridY(); y++){
      
      uint32 *dep = map + y*plusX;
      
      for (Dimn_Type x = 0; x < img.GetGridX(); x++){
        double  ul = dep[0] / 4.; 
        double  ur = dep[1] / 4.; 
        double  ll = dep[plusX] / 4.; 
        double  lr = dep[plusX+1] / 4.; 
        
        for (int r = 0; r <= n; r++){
          double ru = ul*(n-r)/n + ur*r/n;
          double rl = ll*(n-r)/n + lr*r/n;
          
          for (int s = 0; s <= n; s++){
            double    wt = ru*(n-s)/n + rl*s/n; 
            int       wi = (int) wt;
            double    wf = wt-wi;
            Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
            Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
            
            if (wi+1 < img.GetDepth()) {
              res[p] = val[q] * (1.-wf) + wal[q] * wf;
              hmVals[p] = wi * (1.-wf) + wi * wf;
            }else {
              res[p] = val[q];
              hmVals[p] = wi;
            }
#ifdef DEVELOP
            rem[p] = max[p] - res[p];
#endif
          }
        }
        dep += 1;
      }
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "residue.tif",residue);
    Free_Array(residue);
#endif
  
    Clip_Array_Inplace(results.interpolatedHeightMap, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
    Clip_Array_Inplace(results.projection, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
    
    if (imgPara.verbose){ 
      std::cout << "Result interpolated" << std::endl; 
    }
    
    return (results);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map. This projection does not restrict the projection to the computed z-section, but projects the maximal intensity from the computed z-section and its neighboring z-sections
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (8bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  
  ResultArrays InterpolatePlanes_8bit_MaxInterpolation(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners){
    
    ResultArrays results;
    results.projection = Make_Array(PLAIN_KIND,UINT8_TYPE,2, img.GetImage()->dims);
    uint8  *res    = AUINT8(results.projection);
    results.interpolatedHeightMap = Make_Array(PLAIN_KIND,UINT32_TYPE,2, img.GetImage()->dims);
    uint32 * hmVals = AUINT32(results.interpolatedHeightMap);
  
    uint32 *map    = AUINT32(corners);
    uint8  *val    = AUINT8(img.GetImage());
    uint8  *wal    = val + img.GetArea();
    
    int plusX = img.GetGridX() + 1;
    
    Array *maxInt  = Convert_Array(maxim, RGB_KIND, UINT8_TYPE, maxim->scale, 0);
    
#ifdef DEVELOP
    Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
    Array *residue = Make_Array(PLAIN_KIND,UINT8_TYPE,2,img.GetImage()->dims);
    
    uint8 *grn = AUINT8(marked) + img.GetDepth()*img.GetArea();
    uint8 *blu = grn + img.GetDepth()*img.GetArea();
    uint8 *rem = AUINT8(residue);
    uint8 *max = AUINT8(maxim);
#endif
    
    int n = imgPara.radius-1;
    
    for (Dimn_Type y = 0; y < img.GetGridY(); y++){
      
      uint32 *dep = map + y*plusX;
      
      for (Dimn_Type x = 0; x < img.GetGridX(); x++){
        double  ul = dep[0] / 4.;
        double  ur = dep[1] / 4.;
        double  ll = dep[plusX] / 4.;
        double  lr = dep[plusX+1] / 4.;
        
        for (int r = 0; r <= n; r++){
          double ru = ul*(n-r)/n + ur*r/n;
          double rl = ll*(n-r)/n + lr*r/n;
          
          for (int s = 0; s <= n; s++){
            double    wt = ru*(n-s)/n + rl*s/n;
            int       wi = (int) wt;
            double    wf = wt-wi;
            Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
            
            Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
            
            
            int lowerLayer = fmax(0, wi-2);
            int upperLayer = fmin(img.GetDepth()-1, wi+2);
            
            for (int currLayer = lowerLayer; currLayer <= upperLayer; currLayer++) {
              Indx_Type currQ = (Indx_Type) currLayer*img.GetArea() + p;
              if (val[currQ] > val[q]) {
                
                q = currQ;
                wi = currLayer;
                
                if (val[q] > 20) {
                  Draw_Point(maxInt, &PURPLE, Coord2(y*imgPara.radius+s, x*imgPara.radius+r ));
                }
                
              }
            }
            
            //if (wi+1 < img.GetDepth())
            if (wi+1 < img.GetDepth() && wf==wt-wi) {
              res[p] = val[q] * (1.-wf) + wal[q] * wf;
              hmVals[p] = wi * (1.-wf) + wi * wf;
            } else {
              res[p] = val[q];
              hmVals[p] = wi;
            }
            
#ifdef DEVELOP
            rem[p] = max[p] - res[p];
            if (wf < 1.)
            { grn[q] *= (1.-wf)*(1.-wf);
              blu[q] *= wf*wf;
            }
            if (wf > 0.)
            { grn[q+img.GetArea()] *= wf*wf;
              blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
            }
#endif
          }
        }
        dep += 1;
      }
    }
    Free_Array(maxInt);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "marked.tif",marked);
    ShowArray(img, imgPara, "residue.tif",residue);
    Free_Array(residue);
    Free_Array(marked);
#endif
    
    Clip_Array_Inplace(results.interpolatedHeightMap, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
    Clip_Array_Inplace(results.projection, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));

    if (imgPara.verbose){
      std::cout << "Result interpolated" << std::endl;
    }
    
    return (results);
  }
  

  ////////////////////////////////////////////////////////////////////////////////
  /*
   Actual z-Projection using the interpolated height map. This projection does not restrict the projection to the computed z-section, but projects the maximal intensity from the computed z-section and its neighboring z-sections
   -> To avoid having artifacts because of the grid structure, the z-sections are linearly interpolated between the center of each grid element
   
   @param img:      input image stack (16bit)
   @param imgPara:  input parameter
   @param maxim:    maximum projection of the stack
   @param corners:  interpolated height map
   */
  
  ResultArrays InterpolatePlanes_16bit_MaxInterpolation(RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *corners){
    
    ResultArrays results;
    results.projection = Make_Array(PLAIN_KIND,UINT16_TYPE,2, img.GetImage()->dims);
    uint8  *res    = AUINT8(results.projection);
    results.interpolatedHeightMap = Make_Array(PLAIN_KIND,UINT32_TYPE,2, img.GetImage()->dims);
    uint32 * hmVals = AUINT32(results.interpolatedHeightMap);
    
    uint32 *map    = AUINT32(corners);
    uint16  *val    = AUINT16(img.GetImage());
    uint16  *wal    = val + img.GetArea();
    
    int plusX = img.GetGridX() + 1;
        
#ifdef DEVELOP
    Array *residue = Make_Array(PLAIN_KIND,UINT16_TYPE,2,img.GetImage()->dims);
    uint16 *rem = AUINT16(residue);
    uint16 *max = AUINT16(maxim);
#endif
    
    int n = imgPara.radius-1;
    
    for (Dimn_Type y = 0; y < img.GetGridY(); y++){
      
      uint32 *dep = map + y*plusX;
      
      for (Dimn_Type x = 0; x < img.GetGridX(); x++){
        double  ul = dep[0] / 4.;
        double  ur = dep[1] / 4.;
        double  ll = dep[plusX] / 4.;
        double  lr = dep[plusX+1] / 4.;
        
        for (int r = 0; r <= n; r++){
          double ru = ul*(n-r)/n + ur*r/n;
          double rl = ll*(n-r)/n + lr*r/n;
          
          for (int s = 0; s <= n; s++){
            double    wt = ru*(n-s)/n + rl*s/n;
            int       wi = (int) wt;
            double    wf = wt-wi;
            Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
            
            Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
            
            
            int lowerLayer = fmax(0, wi-2);
            int upperLayer = fmin(img.GetDepth()-1, wi+2);
            
            for (int currLayer = lowerLayer; currLayer <= upperLayer; currLayer++) {
              Indx_Type currQ = (Indx_Type) currLayer*img.GetArea() + p;
              if (val[currQ] > val[q]) {                
                q = currQ;
                wi = currLayer;
                
              }
            }
            
            //if (wi+1 < img.GetDepth())
            if (wi+1 < img.GetDepth() && wf==wt-wi) {
              res[p] = val[q] * (1.-wf) + wal[q] * wf;
              hmVals[p] = wi * (1.-wf) + wi * wf;
            } else {
              res[p] = val[q];
              hmVals[p] = wi;

            }
            
#ifdef DEVELOP
            rem[p] = max[p] - res[p];
#endif
          }
        }
        dep += 1;
      }
    }
    
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "residue.tif",residue);
    Free_Array(residue);
#endif
    
    Clip_Array_Inplace(results.interpolatedHeightMap, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));
    Clip_Array_Inplace(results.projection, Coord2(0, 0), Coord2(img.GetRealHeight()-1, img.GetRealWidth()-1));

    if (imgPara.verbose){
      std::cout << "Result interpolated" << std::endl;
    }
    
    return (results);
  }
  

}