//
//  ZPlaneSelection.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 26.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <algorithm>
#include <iostream>

#include "ZPlaneSelection.h"

extern "C" {
#include "filters.h"
#include "histogram.h"
#include "utilities.h"
}

#include "MylibValues.h"

using namespace MylibValues;

namespace ProjectionMethod {
  
  /***********************
   * Select the z-planes with the most occurrences of bright pixel
   **********************/
  
  std::vector<Array *> PlaneSelection_SeveralLevels_8bit_ (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    std::vector<Array *> levels;
    std::vector<uint32 *> levs;

    uint8 * imgVals = AUINT8(img.GetImage());
    
    for (int i=0; i<imgPara.layer; i++) {
      Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
      uint32 *lev   = AUINT32(level);
      
      levels.push_back(level);
      levs.push_back(lev);
    }
    
    uint32 *idx   = AUINT32(index);
    uint8  *max   = AUINT8(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint8 *big     = AUINT8(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    Array_Op_Scalar(dist, SET_OP, UINT32_TYPE, ValU(0));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Array  * brightPixelDistribution = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord3(img.GetDepth(), tmpGridY,tmpGridX));
    Array_Op_Scalar(brightPixelDistribution, SET_OP, UINT32_TYPE, ValU(0));
    uint32 * brightVals = AUINT32(brightPixelDistribution);
    
    
    Frame * fimg = Make_Frame(img.GetImage(),Coord3(img.GetDepth(), imgPara.radius,imgPara.radius),Coord3(0,0,0));
    
    
    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fimg,p);

        Empty_Histogram(wnh);
        Histagain_Array(wnh,fimg,0);
        
        { int sum = 0;
          for (thr = 255; thr > 0; thr--)
          { sum += hcnt[thr];
            if (sum >= imgPara.threshold)
              break;
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        
        for (int z =0; z<img.GetDepth(); z++) {
          for (int r=0; r < imgPara.radius; r++) {
            for (int s=0; s<imgPara.radius; s++) {
              Indx_Type q = (z*img.GetHeight() + y+r ) * img.GetWidth() + (x+s);
              
              if (imgVals[q] >= thr) {
                dcnt[z] += 1;
              }
            }
          }
        }
        
        
        { std::vector<int> best;
          std::vector<int> plan;
          
          if (dcnt[0] > dcnt[1]) {
            best.push_back(dcnt[0]);
            plan.push_back(0);
          }
          
          for (int d = 1; d < img.GetDepth()-1; d++){
            if (dcnt[d] > dcnt[d-1] && dcnt[d] > dcnt[d+1]){
              best.push_back(dcnt[d]);
              plan.push_back(d);
            }
          }
          if (dcnt[img.GetDepth()-1] > dcnt[img.GetDepth()-2]) {
            best.push_back(dcnt[img.GetDepth()-1]);
            plan.push_back(img.GetDepth()-1);
          }
          
          // Assign values in level
          if (best.size() == imgPara.layer) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
              
            }
          } else if (best.size() == 0) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = 0;
              
            }
          } else if (best.size() < imgPara.layer) {
            for (int i=0; i<best.size(); i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
            }
            for (int i=best.size(); i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[best.size()-1];
            }
          } else {
            std::vector<int> sortBest (best.size());
            std::copy(best.begin(), best.end(), sortBest.begin());
            std::sort(sortBest.begin(), sortBest.end());
            
            for (int i=0; i<best.size()-imgPara.layer; i++) {
              auto itr = std::find(best.begin(), best.end(), sortBest[i]);
              *itr = -1;
            }
            int j=0;
            for (int i=0; i<best.size(); i++) {
              if (best[i] > -1) {
                lp = Coord2IdxA(levels[j], Coord2(y/imgPara.radius, x/imgPara.radius));
                levs[j][lp] = plan[i];
                j++;
              }
            }
            
          }
       
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint8)dcnt[j];
          }
#endif
        }
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fimg);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
    for (int i=0; i<levels.size(); i++) {
      std::string path = "lev"+std::to_string(i)+".tif";
      Write_Image((char *)path.data(), levels[i], DONT_PRESS);
    }
    
    return levels;
  }
  
  
  /***********************
   * Select the z-planes with the most occurrences of bright pixel
   **********************/
  
  std::vector<Array *> PlaneSelection_SeveralLevels_8bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    std::vector<Array *> levels;
    std::vector<uint32 *> levs;
    
    
    for (int i=0; i<imgPara.layer; i++) {
      Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
      uint32 *lev   = AUINT32(level);
      
      levels.push_back(level);
      levs.push_back(lev);
    }
    
    uint32 *idx   = AUINT32(index);
    uint8  *max   = AUINT8(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint8 *big     = AUINT8(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    Array_Op_Scalar(dist, SET_OP, UINT32_TYPE, ValU(0));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 255; thr > 0; thr--)
            { sum += hcnt[thr];
              if (sum >= imgPara.threshold)
                break;
            }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){
          uint8  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint8  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else {
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { std::vector<int> best;
          std::vector<int> plan;
          
          if (dcnt[0] > dcnt[1]) {
            best.push_back(dcnt[0]);
            plan.push_back(0);
          }
          
          for (int d = 1; d < img.GetDepth()-1; d++){
            if (dcnt[d] > dcnt[d-1] && dcnt[d] > dcnt[d+1]){
              best.push_back(dcnt[d]);
              plan.push_back(d);
            }
          }
          if (dcnt[img.GetDepth()-1] > dcnt[img.GetDepth()-2]) {
            best.push_back(dcnt[img.GetDepth()-1]);
            plan.push_back(img.GetDepth()-1);
          }
          
          // Assign values in level
          if (best.size() == imgPara.layer) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
              
            }
          } else if (best.size() == 0) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = 0;
              
            }
          } else if (best.size() < imgPara.layer) {
            for (int i=0; i<best.size(); i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
            }
            for (int i=best.size(); i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[best.size()-1];
            }
          } else {
            std::vector<int> sortBest (best.size());
            std::copy(best.begin(), best.end(), sortBest.begin());
            std::sort(sortBest.begin(), sortBest.end());
            
            for (int i=0; i<best.size()-imgPara.layer; i++) {
              auto itr = std::find(best.begin(), best.end(), sortBest[i]);
              *itr = -1;
            }
            int j=0;
            for (int i=0; i<best.size(); i++) {
              if (best[i] > -1) {
                lp = Coord2IdxA(levels[j], Coord2(y/imgPara.radius, x/imgPara.radius));
                levs[j][lp] = plan[i];
                j++;
              }
            }
            
          }
          
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint8)dcnt[j];
          }
#endif
        }
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
//    for (int i=0; i<levels.size(); i++) {
//      std::string path = "lev"+std::to_string(i)+".tif";
//      Write_Image((char *)path.data(), levels[i], DONT_PRESS);
//    }
    
    return levels;
  }


  
  /***********************
   * Select the z-plane with the most occurrences of bright pixel
   **********************/
  
  Array * PlaneSelection_8bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){

    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
    uint32 *lev   = AUINT32(level);
    
    uint32 *idx   = AUINT32(index);
    uint8  *max   = AUINT8(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint8 *big     = AUINT8(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    Array_Op_Scalar(dist, SET_OP, UINT32_TYPE, ValU(0));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 255; thr > 0; thr--)
          { sum += hcnt[thr];
            if (sum >= imgPara.threshold)
              break;
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){ 
          uint8  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint8  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){ 
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else { 
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { int best = 0;
          int plan;
          
          for (int d = 0; d < img.GetDepth(); d++){
            if (dcnt[d] > best){ 
              best = dcnt[d];
              plan = d;
            }
          }
          lp = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
          lev[lp] = plan;
          
#ifdef DEVELOP
          Indx_Type pDist;
                    
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint8)dcnt[j];
          }
#endif
        }      
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){ 
      std::cout << "Z-planes selected" << std::endl; 
    }
    
    return level;
  }
  
  
  /***********************
   * Select the z-planes with the most occurrences of bright pixel
   * 16-bit images
   **********************/
  
  std::vector<Array *> PlaneSelection_SeveralLevels_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    std::vector<Array *> levels;
    std::vector<uint32 *> levs;
    
    
    for (int i=0; i<imgPara.layer; i++) {
      Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
      uint32 *lev   = AUINT32(level);
      
      levels.push_back(level);
      levs.push_back(lev);
    }
    
    uint32 *idx   = AUINT32(index);
    uint16  *max   = AUINT16(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint16 *big     = AUINT16(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    Array_Op_Scalar(dist, SET_OP, UINT32_TYPE, ValU(0));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,65536,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 65535; thr > 0; thr--)
          { sum += hcnt[thr];
            if (sum >= imgPara.threshold)
              break;
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){
          uint16  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint16  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else {
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { std::vector<int> best;
          std::vector<int> plan;
          
          if (dcnt[0] > dcnt[1]) {
            best.push_back(dcnt[0]);
            plan.push_back(0);
          }
          
          for (int d = 1; d < img.GetDepth()-1; d++){
            if (dcnt[d] > dcnt[d-1] && dcnt[d] > dcnt[d+1]){
              best.push_back(dcnt[d]);
              plan.push_back(d);
            }
          }
          if (dcnt[img.GetDepth()-1] > dcnt[img.GetDepth()-2]) {
            best.push_back(dcnt[img.GetDepth()-1]);
            plan.push_back(img.GetDepth()-1);
          }
          
          // Assign values in level
          if (best.size() == imgPara.layer) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
              
            }
          } else if (best.size() == 0) {
            for (int i=0; i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = 0;
              
            }
          } else if (best.size() < imgPara.layer) {
            for (int i=0; i<best.size(); i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
            }
            for (int i=best.size(); i<imgPara.layer; i++) {
              lp = Coord2IdxA(levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[best.size()-1];
            }
          } else {
            std::vector<int> sortBest (best.size());
            std::copy(best.begin(), best.end(), sortBest.begin());
            std::sort(sortBest.begin(), sortBest.end());
            
            for (int i=0; i<best.size()-imgPara.layer; i++) {
              auto itr = std::find(best.begin(), best.end(), sortBest[i]);
              *itr = -1;
            }
            int j=0;
            for (int i=0; i<best.size(); i++) {
              if (best[i] > -1) {
                lp = Coord2IdxA(levels[j], Coord2(y/imgPara.radius, x/imgPara.radius));
                levs[j][lp] = plan[i];
                j++;
              }
            }
            
          }
          
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint16)dcnt[j];
          }
#endif
        }
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
    for (int i=0; i<levels.size(); i++) {
      std::string path = "lev"+std::to_string(i)+".tif";
      Write_Image((char *)path.data(), levels[i], DONT_PRESS);
    }
    
    return levels;
  }
  

  /***********************
   * Select the z-planes with the most occurrences of bright pixel
   * 16-bit images
   **********************/
  
  Levels PlaneSelection_2Levels_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    Levels surfaces;
    std::vector<uint32 *> levs;
    surfaces.confidence = Make_Array_With_Shape(PLAIN_KIND,UINT8_TYPE,Coord2(tmpGridY,tmpGridX));
    uint8 * confVals = AUINT8(surfaces.confidence);
    
    
    for (int i=0; i<imgPara.layer; i++) {
      Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
      uint32 *lev   = AUINT32(level);
      
      surfaces.levels.push_back(level);
      levs.push_back(lev);
    }
    
    uint32 *idx   = AUINT32(index);
    uint16  *max   = AUINT16(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint16 *big     = AUINT16(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    Array_Op_Scalar(dist, SET_OP, UINT32_TYPE, ValU(0));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,65536,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 65535; thr > 0; thr--)
          { sum += hcnt[thr];
            if (sum >= imgPara.threshold)
              break;
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){
          uint16  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint16  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else {
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { std::vector<int> best;
          std::vector<int> plan;
          
          if (dcnt[0] > dcnt[1]) {
            best.push_back(dcnt[0]);
            plan.push_back(0);
          }
          
          for (int d = 1; d < img.GetDepth()-1; d++){
            if (dcnt[d] > dcnt[d-1] && dcnt[d] > dcnt[d+1]){
              best.push_back(dcnt[d]);
              plan.push_back(d);
            }
          }
          if (dcnt[img.GetDepth()-1] > dcnt[img.GetDepth()-2]) {
            best.push_back(dcnt[img.GetDepth()-1]);
            plan.push_back(img.GetDepth()-1);
          }
          
          // Assign values in level
          //    Several cases:
          //        # peaks = 0   -> take level 0
          //        # peaks = 2   -> assign the according to the direction
          //        # peaks = 1   -> assign both to one
          //        # peaks > 2   -> determine 2 highest ones
          
          
          
          if (best.size() == 2) {
            for (int i=0; i<2; i++) {
              lp = Coord2IdxA(surfaces.levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[i];
              
            }
            confVals[p] = 2;
          } else if (best.size() == 0) {
            for (int i=0; i<2; i++) {
              lp = Coord2IdxA(surfaces.levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = 0;
              
            }
            confVals[p] = 0;
            
          } else if (best.size() == 1) {
            for (int i=0; i<2; i++) {
              lp = Coord2IdxA(surfaces.levels[i], Coord2(y/imgPara.radius, x/imgPara.radius));
              levs[i][lp] = plan[0];
            }
            confVals[p] = 1;
          } else {
            std::vector<int> sortBest (best.size());
            std::copy(best.begin(), best.end(), sortBest.begin());
            std::sort(sortBest.begin(), sortBest.end());
            
            for (int i=0; i<best.size()-imgPara.layer; i++) {
              auto itr = std::find(best.begin(), best.end(), sortBest[i]);
              *itr = -1;
            }
            int j=0;
            for (int i=0; i<best.size(); i++) {
              if (best[i] > -1) {
                lp = Coord2IdxA(surfaces.levels[j], Coord2(y/imgPara.radius, x/imgPara.radius));
                levs[j][lp] = plan[i];
                j++;
              }
            }
            confVals[p] = 1;
          }
          
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint16)dcnt[j];
          }
#endif
        }
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
    for (int i=0; i<surfaces.levels.size(); i++) {
      std::string path = "lev"+std::to_string(i)+".tif";
      Write_Image((char *)path.data(), surfaces.levels[i], DONT_PRESS);
    }
    
    return surfaces;
  }
  
  
  

  
  
  /***********************
   *
   * Select the z-plane with the most occurrences of bright pixel
   * 16-bit images
   **********************/
  
  Array * PlaneSelection_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
    uint32 *lev   = AUINT32(level);
    
    uint32 *idx   = AUINT32(index);
    uint16  *max   = AUINT16(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT16_TYPE,2,maxim->dims);
    uint16 *big     = AUINT16(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,65536,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 65535; thr > 0; thr--)
          { sum += hcnt[thr];
            if (sum >= imgPara.threshold)
              break;
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){
          uint16  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint16  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else {
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { int best = 0;
          int plan;
          
          for (int d = 0; d < img.GetDepth(); d++){
            if (dcnt[d] > best){
              best = dcnt[d];
              plan = d;
            }
          }
          lp = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
          lev[lp] = plan;
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
            dis[pDist] = (uint8)dcnt[j];
          }
#endif
        }
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
    return level;
  }

  /***********************
   * Select the z-plane with the most occurrences of bright pixel
   *    Uses a window size twice the size of the input radius
   **********************/
  
  
  Array * PlaneSelection_DoubleRadius (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index, int rad){
        
    int tmpGridX = ceil ((double)img.GetWidth()/(double)(rad/2));
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)(rad/2));
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
    uint32 *lev   = AUINT32(level);
    
    uint32 *idx   = AUINT32(index);
    uint8  *max   = AUINT8(maxim);
    
#ifdef DEVELOP
    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
    uint8 *big     = AUINT8(largest);
    
    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
    
    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(ceil ((double)img.GetHeight()/(double)(rad)), ceil ((double)img.GetWidth()/(double)(rad)), img.GetDepth()));
    uint32 *dis    = AUINT32(dist);
#endif
    
    Use_Extend_Boundary();
    
    Frame * fid = Make_Frame(index,Coord2(rad,rad),Coord2(0,0));
    Frame * fmx = Make_Frame(maxim,Coord2(rad,rad),Coord2(0,0));
    Offs_Type * off = Frame_Offsets(fid);
    Size_Type fsz = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(),Program_Name());
    
    
    Indx_Type lp = 0;
    for (Dimn_Type y = 0; y <= img.GetHeight()-rad; y += rad){
      for (Dimn_Type x = 0; x <= img.GetWidth()-rad; x += rad){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fmx,p);
        Place_Frame(fid,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fmx,0);
        
        { int sum = 0;
          for (thr = 255; thr > 0; thr--){
            sum += hcnt[thr];
            if (sum >= imgPara.threshold){
              break;
            }
          }
        }
        
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        if (Frame_Within_Array(fid)){
          uint8  *mx = max + p;
          uint32 *ix = idx + p;
#ifdef DEVELOP
          uint8  *bg = big + p;
#endif
          for (int d = 0; d < fsz; d++)
            if (mx[off[d]] >= thr){
              dcnt[ix[off[d]]] += 1;
#ifdef DEVELOP
              bg[off[d]] = mx[off[d]];
#endif
            }
        } else {
          std::cout << "Frame should never be outside array" << std::endl;
          exit (0);
        }
        
        { int plan;
          int best = 0;
          
          for (int d = 0; d < img.GetDepth(); d++){
            if (dcnt[d] > best){
              best = dcnt[d];
              plan = d;
            }
          }
          lp = Coord2IdxA(level, Coord2((2*y)/rad, (2*x)/rad));
          lev[lp] = plan;
          
          lp = Coord2IdxA(level, Coord2((2*y+rad)/rad, (2*x)/rad));
          lev[lp] = plan;
          
          lp = Coord2IdxA(level, Coord2((2*y)/rad, (2*x+rad)/rad));
          lev[lp] = plan;
          
          lp = Coord2IdxA(level, Coord2((2*y+rad)/rad, (2*x+rad)/rad));
          lev[lp] = plan;
          
          
#ifdef DEVELOP
          Indx_Type pDist;
          
          for (int j=0; j<img.GetDepth(); j++) {
            pDist = Coord2IdxA(dist, Coord3(y/rad, x/rad, j));
            dis[pDist] = (uint8)dcnt[j];
          }
#endif
        }      
        
      }
    }
    
    free(dcnt);
    Kill_Histogram(wnh);
    Kill_Frame(fid);
    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);  
    Free_Array(dist);
    
    ShowArray(img, imgPara, "level_inital.tif", level);
#endif
    
    if (imgPara.verbose){ 
      std::cout << "Z-planes selected" << std::endl; 
    }
    
    return level;
  }
  
  /***********************
   *
   * Select the z-plane with the highest entropy
   *    Uses a subset of the stack defined by a reference height map and a distance
   * 8-bit images
   **********************/
  
  
  Array * SelectPlanes_Variance_8bit(RasterizedImage & img, ParaSet & imgPara, Array * median){
        
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
    uint32 *lev   = AUINT32(level);
    uint32 * med  = AUINT32(median);
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(img.GetImage(),Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Place_Frame(f, 0);
    
    double maxVar;
    double maxmaxVar;
    int best, maxBest;
    Indx_Type p;
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        best = 0;
        maxVar = 1.5;
        maxmaxVar = 0;
        maxBest = 0;
        
        //        std::cout << x << ", " << y << "   ( " << x/imgPara.radius << ", " << y/ imgPara.radius << std::endl;
        Indx_Type medIdx = (y/imgPara.radius) * img.GetGridX()+ (x/imgPara.radius);
        
        for (Dimn_Type z=med[medIdx]-imgPara.distance ; z<=med[medIdx]+imgPara.distance ; z++) {
          if ((z < 0) or (z >= img.GetDepth())) {
            continue;
          }
          
          p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
          Place_Frame(f, p);
          
          Empty_Histogram(h);
          Histagain_Array(h,f,0);
          //          if (Histogram_Variance(h) > 0) {
          //             std::cout << "     " << z << " => " << Histogram_Variance(h) << std::endl;
          //          }
          
          if (maxVar < Histogram_Variance(h)) {
            maxVar = Histogram_Variance(h);
            best = z;
          }
          if (maxmaxVar < Histogram_Variance(h)) {
            maxmaxVar = Histogram_Variance(h);
            maxBest = z;
          }
        }
        
        double betweenVar = maxVar;
        int betweenBest = 0;
        
        if (best == med[medIdx]-imgPara.distance) {
          
          int z = best -1;
          while (z >= 0) {
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
            Place_Frame(f, p);
            
            Empty_Histogram(h);
            Histagain_Array(h,f,0);
            //            if (Histogram_Variance(h) > 1) {
            //              std::cout << "     smaller   " << z << " => " << Histogram_Variance(h) << std::endl;
            //            }
            if (maxmaxVar < Histogram_Variance(h)) {
              maxmaxVar = Histogram_Variance(h);
              maxBest = z;
            }
            
            if (betweenVar < Histogram_Variance(h)) {
              betweenVar = Histogram_Variance(h);
              betweenBest = z;
              z--;
            } else {
              z = -1;
            }
          }
          
        } else if (best == 0) {
          
          int z = med[medIdx]-imgPara.distance -1;
          while (z >= 0) {
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
            Place_Frame(f, p);
            
            Empty_Histogram(h);
            Histagain_Array(h,f,0);
            //            if (Histogram_Variance(h) > 0) {
            //              std::cout << "     smaller   " << z << " => " << Histogram_Variance(h) << std::endl;
            //            }
            if (maxmaxVar < Histogram_Variance(h)) {
              maxmaxVar = Histogram_Variance(h);
              maxBest = z;
            }
            
            if (betweenVar < Histogram_Variance(h)) {
              betweenVar = Histogram_Variance(h);
              betweenBest = z;
              z--;
            } else if (maxVar == 1.5) {
              z--;
            }else {
              z = -1;
            }
          }
          
        }
        
        
        if (best == med[medIdx]+imgPara.distance) {
          
          int z = best + 1;
          while (z < img.GetDepth()) {
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
            Place_Frame(f, p);
            
            Empty_Histogram(h);
            Histagain_Array(h,f,0);
            //            if (Histogram_Variance(h) > 0) {
            //              std::cout << "     greater   " << z << " => " << Histogram_Variance(h) << std::endl;
            //            }
            if (maxmaxVar < Histogram_Variance(h)) {
              //              maxmaxVar = Histogram_Variance(h);
              maxBest = z;
            }
            
            if (maxVar < Histogram_Variance(h)) {
              maxVar = Histogram_Variance(h);
              best = z;
              z++;
            } else {
              z = img.GetDepth();
            }
          }
          
        } else if (best == 0) {
          
          int z = med[medIdx]+imgPara.distance + 1;
          while (z < img.GetDepth()) {
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
            Place_Frame(f, p);
            
            Empty_Histogram(h);
            Histagain_Array(h,f,0);
            //            if (Histogram_Variance(h) > 1) {
            //              std::cout << "     greater   " << z << " => " << Histogram_Variance(h) << std::endl;
            //            }
            if (maxmaxVar < Histogram_Variance(h)) {
              maxmaxVar = Histogram_Variance(h);
              maxBest = z;
            }
            
            if (maxVar < Histogram_Variance(h)) {
              maxVar = Histogram_Variance(h);
              best = z;
              z++;
            } else if (maxVar == 1.5) {
              z++;
            }else {
              z = img.GetDepth();
            }
          }
          
        }
        
        if (abs(betweenBest - med[medIdx]-imgPara.distance) < abs(best-med[medIdx]+imgPara.distance)) {
          best = betweenBest;
        }
        
        if (best == 0) {
          best = maxBest;
        }
        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
        //        std::cout << "       Best " << best << "  ( " << y/imgPara.radius << " , " << x/imgPara.radius<<  ")   p=" << y/imgPara.radius * level->dims[0] + x/imgPara.radius << " = " << p << std::endl;
        lev[p] = best;
      }
    }
    
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "level_Variance.tif", level);
#endif
    if (imgPara.verbose){
      std::cout << "Optimal planes selected (showing highest variance)" << std::endl;
    }
    Free_Frame(f);
    Free_Histogram(h);
    
    return level;
  }

  
  /***********************
   *
   * Select the z-plane with the highest entropy
   *    Uses a subset of the stack defined by a reference height map and a distance
   * 16-bit images
   **********************/
  
  
  Array * SelectPlanes_Variance_16bit(RasterizedImage & img, ParaSet & imgPara, Array * substack){
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
    uint32 *lev   = AUINT32(level);
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(substack,Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
    
    /*****
     * Fast version
     *****/
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    /*****
     * More accurate but slow version
     *****/
    //Histogram * h = Make_Histogram(UVAL,2048,ValU(32),ValU(0));
    

    Place_Frame(f, 0);
    
    double maxVar;
    int best;
    Indx_Type p;
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        best = 0;
        maxVar = 0;
        
        for (Dimn_Type z=0; z<img.GetDepth(); z++) {
          p = Coord2IdxA(substack, Coord3(z, y, x));
          Place_Frame(f, p);
          
          Empty_Histogram(h);
          Histagain_Array(h,f,0);
          
          if (maxVar < Histogram_Variance(h)) {
            maxVar = Histogram_Variance(h);
            best = z;
          }
          
        }
        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
        lev[p] = best;
        
      }
    }
    
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "level_Variance.tif", level);
#endif
    if (imgPara.verbose){
      std::cout << "Optimal planes selected (showing highest variance)" << std::endl;
    }
    
    Free_Frame(f);
    Free_Histogram(h);
    return level;
  }

  
  /***********************
   * Select the z-plane with the highest entropy
   *    Uses a subset of the stack defined by a reference height map and a distance
   **********************/

  
  Array * SelectPlanes_Entropy(RasterizedImage & img, ParaSet & imgPara, Array * substack){
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
    uint32 *lev   = AUINT32(level);
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(substack,Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));  
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Place_Frame(f, 0);  
    
    double maxVar;
    int best;
    Indx_Type p;
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        best = 0;
        maxVar = 0;
        
        for (Dimn_Type z=0; z<img.GetDepth(); z++) {
          p = Coord2IdxA(substack, Coord3(z, y, x));
          Place_Frame(f, p);
          
          Empty_Histogram(h);
          Histagain_Array(h,f,0);
          
          if (maxVar < Histogram_Entropy(h)) {
            maxVar = Histogram_Entropy(h);
            best = z;
          }
          
        }
        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
        lev[p] = best;
        
      }
    }
    
    
#ifdef DEVELOP  
    ShowArray(img, imgPara, "level_Entropy.tif", level);
#endif
    if (imgPara.verbose){
      std::cout << "Optimal planes selected (showing highest entropy)" << std::endl;
    }
    
    Free_Frame(f);
    Free_Histogram(h);
    return level;
  }

  
  /***********************
   * Select the z-plane with the highest LOG score
   **********************/
  
  
  Array * SelectPlanes_LOG(RasterizedImage & img, ParaSet & imgPara, Array * substack){

    Convert_Array_Inplace(substack, PLAIN_KIND, INT32_TYPE, 32, 0);
    
    Array * logImage = LOG_Array(substack, 1.0, 3);
        
    Convert_Array_Inplace(logImage, PLAIN_KIND, UINT32_TYPE, 32, 0);

    uint32 * logVal = AUINT32(logImage);
    
        
    // Going through the image, pixel by pixel and determine the maximal absolut value for each z-column
        
    Array * indexLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
    uint32 * idxLoG = AUINT32(indexLoG);
    
    Array * maximLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
    uint32 * maxLoG = AUINT32(maximLoG);
    
    uint32 maxValue, maxIndex;
    Indx_Type q;
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){ 
      maxValue = abs(logVal[p]);
      maxIndex = 0;
      q = p + img.GetArea();
      
      for (Dimn_Type d = 1; d < img.GetDepth(); d++){ 
        if (abs(logVal[q]) >= maxValue) {
          maxValue = abs(logVal[q]);
          maxIndex = d;
        }
        q += img.GetArea();
      }
      maxLoG[p] = maxValue;
      idxLoG[p] = maxIndex;
    }
   
    RasterizedImage imgLoG = RasterizedImage(logImage);
    imgLoG.SetGrid(imgPara.radius);
    
    Array * level = PlaneSelection_8bit(imgLoG, imgPara, maximLoG, indexLoG);
        
        
#ifdef DEVELOP
     Write_Image("logImage.tif", logImage, DONT_PRESS);
    
    ShowArray(imgLoG, imgPara, "maxim_loG.tif", maximLoG);
    ShowArray(imgLoG, imgPara, "index_loG.tif", indexLoG);
    ShowArray(img, imgPara, "level_LoG.tif", level);
#endif
    
    if (imgPara.verbose) {
      std::cout << "Optimal planes selected (using LoG)" << std::endl;
    }
    
    Free_Array(maximLoG);
    Free_Array(indexLoG);
    Free_Array(logImage);
    
    return level;
  }
  
  
  /***********************
   * Select the z-plane with the highest LOG score -> basal cell layer
   **********************/
  
  
  Array * SelectPlanes_LOG_basal(RasterizedImage & img, ParaSet & imgPara, Array * substack){
    
    Convert_Array_Inplace(substack, PLAIN_KIND, INT32_TYPE, 32, 0);
    
    Array * logImage = LOG_Array(substack, 1.0, 3);
    
    Convert_Array_Inplace(logImage, PLAIN_KIND, UINT32_TYPE, 32, 0);
    
    uint32 * logVal = AUINT32(logImage);
    
    Write_Image("logImage.tif", logImage, DONT_PRESS);
    
    
    // Going through the image, pixel by pixel and determine the maximal absolut value for each z-column
    
    Array * indexLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
    uint32 * idxLoG = AUINT32(indexLoG);
    
    Array * maximLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
    uint32 * maxLoG = AUINT32(maximLoG);
    
    uint32 maxValue, maxIndex;
    Indx_Type q;
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){ 
      maxValue = abs(logVal[p]);
      maxIndex = 0;
      q = p + img.GetArea();
      
      for (Dimn_Type d = 1; d < img.GetDepth(); d++){ 
        if (abs(logVal[q]) >= maxValue) {
          maxValue = abs(logVal[q]);
          maxIndex = d;
        }
        
        q += img.GetArea();
      }
      maxLoG[p] = maxValue;
      idxLoG[p] = maxIndex;
    }
    RasterizedImage imgLoG = RasterizedImage(logImage);
    
    
#ifdef DEVELOP
    ShowArray(imgLoG, imgPara, "maxim_loG.tif", maximLoG);
    ShowArray(imgLoG, imgPara, "index_loG.tif", indexLoG);
#endif
    
        
    Free_Array(maximLoG);
    
    Free_Array(logImage);
    
    
     // Determining the level with the highest frequency of pixels
     
    Array * level = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetGridY(), img.GetGridX()));
    uint32 * lev = AUINT32(level);
     
    Use_Reflective_Boundary();
     
    Frame * f = Make_Frame(indexLoG, Coord2(imgPara.radius,imgPara.radius),Coord2(0, 0));      
    Histogram * h = Make_Histogram(UVAL, img.GetDepth(), ValU(1), ValU(0));
    Place_Frame(f, 0);  
     
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
     for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
       
       Place_Frame(f, Coord2IdxA(indexLoG, Coord2(y, x)));
       Empty_Histogram(h);
       Histagain_Array(h, f, 0);
     
       Size_Type * hValues = h->counts;
     
       int32 maxVal = 0;
       int32 maxPlane = 0;
     
       for (int i = 0; i < img.GetDepth(); i++) {
         if (hValues[i] > maxVal) {
           maxVal = hValues[i];
           maxPlane = (uint32) i;
         }
       }
     
       Indx_Type p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
       lev[p] = maxPlane;
     }
    }
     
     
#ifdef DEVELOP
  ShowArray(img, imgPara, "level_loG.tif", level);
#endif
     
  if (imgPara.verbose) {
     std::cout << "Optimal planes selected (using LoG)" << std::endl;
  }
     
     
  Free_Frame(f);
  Free_Histogram(h);
  Free_Array(indexLoG);

  return level;
  }

  
  /***********************
   * Select the z-plane with the most occurrences of bright pixel
   **********************/
  
  Array * BrightPixelDistribution (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index){
    
    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
    
    Array  *level = Make_Array_With_Shape(PLAIN_KIND,FLOAT32_TYPE,Coord3(img.GetDepth(), img.GetHeight(), img.GetWidth()));
    Array_Op_Scalar(level, SET_OP, FLOAT32_TYPE, ValF(0));
    float32 *lev   = AFLOAT32(level);
    
    Array  *level2 = Make_Array_With_Shape(PLAIN_KIND,FLOAT32_TYPE,Coord3(img.GetDepth(), tmpGridY,tmpGridX));
    Array_Op_Scalar(level2, SET_OP, FLOAT32_TYPE, ValF(0));
    float32 *lev2   = AFLOAT32(level2);
    
    Array * finalLevel = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(tmpGridY, tmpGridX));
    Array_Op_Scalar(finalLevel, SET_OP, UINT32_TYPE, ValU(0));
    uint32 * finVals = AUINT32(finalLevel);
    
    
    
    uint32 *idx   = AUINT32(index);
    uint8  *max   = AUINT8(maxim);
    uint8 *imgVals = AUINT8(img.GetImage());
    
    Use_Extend_Boundary();
    
//    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
//    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
    
    Frame * fimg = Make_Frame(img.GetImage(),Coord3(img.GetDepth(), imgPara.radius,imgPara.radius),Coord3(0,0,0));

    
//    Offs_Type * off = Frame_Offsets(fid);
//    Size_Type fsz   = AForm_Size(fid);
    
    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    Size_Type * hcnt = wnh->counts;
    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
    
    Indx_Type lp = 0;
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        int thr;
        
        Indx_Type p = y*img.GetWidth() + x;
        
        Place_Frame(fimg,p);
        
        Empty_Histogram(wnh);
        Histagain_Array(wnh,fimg,0);
        
        int sum = 0;
        int finalThreshold = 0;
        
          for (thr = 255; thr > 0; thr--)
            { sum += hcnt[thr];
              if (sum >= imgPara.threshold) {
                finalThreshold = thr;
                break;
              }
            }
        for (int d = 0; d < img.GetDepth(); d++){
          dcnt[d] = 0;
        }
        
        for (int z=0; z<img.GetDepth(); z++) {
          for (int r = 0; r<imgPara.radius; r++) {
            for (int s = 0; s<imgPara.radius; s++) {
              Indx_Type q = (z*img.GetHeight() + y+r ) * img.GetWidth() + (x+s);
              Indx_Type pq = (z*tmpGridY + y/imgPara.radius ) * tmpGridX + (x/imgPara.radius);

              if (imgVals[q] >= finalThreshold) {
                lev[q] += 1;
                lev2[pq] += 1;
              }
            }
          }
        }
        
        int maxFequency = 0;
        int bestZ = 0;
        
        for (int z=0; z<img.GetDepth(); z++) {
          Indx_Type pq = (z*tmpGridY + y/imgPara.radius ) * tmpGridX + (x/imgPara.radius);

          if ((int)lev2[pq] > maxFequency) {
            maxFequency = lev2[pq];
            bestZ = z;
          }
        }
        finVals[(y/imgPara.radius ) * tmpGridX + (x/imgPara.radius)] = bestZ;
        
//        if (Frame_Within_Array(fimg)){
//          
//          
//          uint8  *mx = max + p;
//          uint32 *ix = idx + p;
//
//          for (int d = 0; d < fsz; d++)
//            if (mx[off[d]] >= thr){
//              dcnt[ix[off[d]]] += 1;
//            }
//        } else {
//          std::cout << "Frame should never be outside array" << std::endl;
//          exit (0);
//        }
//        
//        { int best = 0;
//          int plan;
//          
//          for (int d = 0; d < img.GetDepth(); d++){
//            
//            
//            if (dcnt[d] > best){
//              best = dcnt[d];
//              plan = d;
//            }
//          }
//          lp = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//          lev[lp] = plan;
//          
//#ifdef DEVELOP
//          Indx_Type pDist;
//          
//          for (int j=0; j<img.GetDepth(); j++) {
//            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
//            dis[pDist] = (uint8)dcnt[j];
//          }
//#endif
//        }
        
      }
    }
    Write_Image("tets4.tif", level, DONT_PRESS);
    Write_Image("tets5.tif", finalLevel, DONT_PRESS);

    
    free(dcnt);
    Kill_Histogram(wnh);
//    Kill_Frame(fid);
    Kill_Frame(fimg);
//    Kill_Frame(fmx);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "largest.tif",largest);
    Free_Array(largest);
    
    Write_Image("dist.tif", dist, DONT_PRESS);
    Free_Array(dist);
#endif
    
    if (imgPara.verbose){
      std::cout << "Z-planes selected" << std::endl;
    }
    
    return level2;
  }
  
  
  Array * SelectPlanes_Variance_UsingHeightMap(RasterizedImage & img, ParaSet & imgPara, ExtendedParaSet & imgPara2, Array * heightMap) {
    
      Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
      uint32 *lev   = AUINT32(level);
      uint32 * med  = AUINT32(heightMap);
      
      Use_Reflective_Boundary();
      
      Frame * f = Make_Frame(img.GetImage(),Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
      Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      
      Place_Frame(f, 0);
      
      double maxVar;
      double maxmaxVar;
      int best, maxBest;
      Indx_Type p;
      
      for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
        for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
          best = 0;
          maxVar = 1.5;
          maxmaxVar = 0;
          maxBest = 0;
          
          //        std::cout << x << ", " << y << "   ( " << x/imgPara.radius << ", " << y/ imgPara.radius << std::endl;
          Indx_Type medIdx = (y/imgPara.radius) * img.GetGridX()+ (x/imgPara.radius);
          
          for (Dimn_Type z=med[medIdx]+imgPara2.distance1 ; z<=med[medIdx]+imgPara2.distance2 ; z++) {
            if ((z < 0) or (z >= img.GetDepth())) {
              continue;
            }
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
            Place_Frame(f, p);
            
            Empty_Histogram(h);
            Histagain_Array(h,f,0);
            //          if (Histogram_Variance(h) > 0) {
            //             std::cout << "     " << z << " => " << Histogram_Variance(h) << std::endl;
            //          }
            
            if (maxVar < Histogram_Variance(h)) {
              maxVar = Histogram_Variance(h);
              best = z;
            }
            if (maxmaxVar < Histogram_Variance(h)) {
              maxmaxVar = Histogram_Variance(h);
              maxBest = z;
            }
          }
          
          double betweenVar = maxVar;
          int betweenBest = 0;
          
          if (best == med[medIdx]-imgPara.distance) {
            
            int z = best -1;
            while (z >= 0) {
              
              p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
              Place_Frame(f, p);
              
              Empty_Histogram(h);
              Histagain_Array(h,f,0);
              //            if (Histogram_Variance(h) > 1) {
              //              std::cout << "     smaller   " << z << " => " << Histogram_Variance(h) << std::endl;
              //            }
              if (maxmaxVar < Histogram_Variance(h)) {
                maxmaxVar = Histogram_Variance(h);
                maxBest = z;
              }
              
              if (betweenVar < Histogram_Variance(h)) {
                betweenVar = Histogram_Variance(h);
                betweenBest = z;
                z--;
              } else {
                z = -1;
              }
            }
            
          } else if (best == 0) {
            
            int z = med[medIdx]-imgPara.distance -1;
            while (z >= 0) {
              
              p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
              Place_Frame(f, p);
              
              Empty_Histogram(h);
              Histagain_Array(h,f,0);
              //            if (Histogram_Variance(h) > 0) {
              //              std::cout << "     smaller   " << z << " => " << Histogram_Variance(h) << std::endl;
              //            }
              if (maxmaxVar < Histogram_Variance(h)) {
                maxmaxVar = Histogram_Variance(h);
                maxBest = z;
              }
              
              if (betweenVar < Histogram_Variance(h)) {
                betweenVar = Histogram_Variance(h);
                betweenBest = z;
                z--;
              } else if (maxVar == 1.5) {
                z--;
              }else {
                z = -1;
              }
            }
            
          }
          
          
          if (best == med[medIdx]+imgPara.distance) {
            
            int z = best + 1;
            while (z < img.GetDepth()) {
              
              p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
              Place_Frame(f, p);
              
              Empty_Histogram(h);
              Histagain_Array(h,f,0);
              //            if (Histogram_Variance(h) > 0) {
              //              std::cout << "     greater   " << z << " => " << Histogram_Variance(h) << std::endl;
              //            }
              if (maxmaxVar < Histogram_Variance(h)) {
                //              maxmaxVar = Histogram_Variance(h);
                maxBest = z;
              }
              
              if (maxVar < Histogram_Variance(h)) {
                maxVar = Histogram_Variance(h);
                best = z;
                z++;
              } else {
                z = img.GetDepth();
              }
            }
            
          } else if (best == 0) {
            
            int z = med[medIdx]+imgPara.distance + 1;
            while (z < img.GetDepth()) {
              
              p = Coord2IdxA(img.GetImage(), Coord3(z, y, x));
              Place_Frame(f, p);
              
              Empty_Histogram(h);
              Histagain_Array(h,f,0);
              //            if (Histogram_Variance(h) > 1) {
              //              std::cout << "     greater   " << z << " => " << Histogram_Variance(h) << std::endl;
              //            }
              if (maxmaxVar < Histogram_Variance(h)) {
                maxmaxVar = Histogram_Variance(h);
                maxBest = z;
              }
              
              if (maxVar < Histogram_Variance(h)) {
                maxVar = Histogram_Variance(h);
                best = z;
                z++;
              } else if (maxVar == 1.5) {
                z++;
              }else {
                z = img.GetDepth();
              }
            }
            
          }
          
          if (abs(betweenBest - med[medIdx]-imgPara.distance) < abs(best-med[medIdx]+imgPara.distance)) {
            best = betweenBest;
          }
          
          if (best == 0) {
            best = maxBest;
          }
          p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
          //        std::cout << "       Best " << best << "  ( " << y/imgPara.radius << " , " << x/imgPara.radius<<  ")   p=" << y/imgPara.radius * level->dims[0] + x/imgPara.radius << " = " << p << std::endl;
          lev[p] = best;
        }
      }
      
      
#ifdef DEVELOP
      ShowArray(img, imgPara, "level_Variance.tif", level);
#endif
      if (imgPara.verbose){
        std::cout << "Optimal planes selected (showing highest variance)" << std::endl;
      }
      Free_Frame(f);
      Free_Histogram(h);
      
      return level;

  }

}


