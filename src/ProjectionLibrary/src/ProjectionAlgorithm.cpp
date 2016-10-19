//
//  ProjectionAlgorithm.cpp
//  CMakeProject_4.0
//
//  Created by blasse on 8/26/13.
//
//

#include "ProjectionAlgorithm.h"

#include <iostream>

extern "C" {
#include "image.h"
#include "utilities.h"
}

#include "Filter.h"
#include "HeightMapAnalysis.h"
#include "Interpolation.h"
#include "Refinement.h"
#include "ShowArray.h"
#include "ZPlaneSelection.h"
#include "ZProjection.h"

using namespace ProjectionResults;

namespace ProjectionMethod {
  
  
  /***********************
   * Constructor / Destructor
   **********************/
  
   ProjectionAlgorithm::ProjectionAlgorithm (std::string imgPath) {
    image_ = std::make_shared<RasterizedImage>(Read_Image(&imgPath[0],0));
    parametersAreSet_ = false;
  }
  
  ProjectionAlgorithm::ProjectionAlgorithm (Array * image) {
    image_ = std::make_shared<RasterizedImage>(image);
//    try {
//      image_ = RasterizedImage(image);
//    } catch (const std::exception& ex) {
//      std::cout << "FELA" << std::endl;
//    }
//    
//    std::cout << "ProjectionAlgoirthm::Constructor" << std::endl;

    parametersAreSet_ = false;
  }
  
  ProjectionAlgorithm::ProjectionAlgorithm () {
    parametersAreSet_ = false;
  }
//  
  ProjectionAlgorithm::~ProjectionAlgorithm () {
//    image_.~RasterizedImage();
//    Free_Array(image_.GetImage());
    
    for (int i=0; i<imageParameter_.layer; i++) {
      Free_Array(projections_[i]);
      Free_Array(heightMaps_[i]);
      Free_Array(interpolatedHeightMaps_[i]);
    }
  }
  
  
  
  /***********************
   * Get functions
   **********************/
  
  RasterizedImage ProjectionAlgorithm::GetImage() const {
    return *image_;
  }
  
  int ProjectionAlgorithm::GetGridSize() const {
    return imageParameter_.radius;
  }

  int ProjectionAlgorithm::GetDistance() const {
    return imageParameter_.distance;
  }
  
  int ProjectionAlgorithm::GetLayerNumber() const {
    return imageParameter_.layer;
  }
  
  int ProjectionAlgorithm::GetBrightPixelCount() const {
    return imageParameter_.threshold;
  }
  
  bool ProjectionAlgorithm::IsVerbose() const {
    return imageParameter_.verbose;
  }
  
  bool ProjectionAlgorithm::UseMaximumInterpolation() const {
    return imageParameter_.maxInterpolation;
  }
  
  bool ProjectionAlgorithm::PrintHeightMap() const {
    return imageParameter_.printHeightMap;
  }
  
  bool ProjectionAlgorithm::PrintHeightMapRealSize() const {
    return imageParameter_.printRealHeightMap;
  }
  
  ParaSet ProjectionAlgorithm::GetParameterSet() const {
    return imageParameter_;
  }
  
  /***********************
   * Set functions
   **********************/
  
  void ProjectionAlgorithm::SetGridSize(int gridSize) {
    
    if (gridSize % 2 != 0) {
      gridSize++;
    }
    imageParameter_.radius = gridSize;
  }
  
  void ProjectionAlgorithm::SetDistance(int distance) {
    imageParameter_.distance = distance;
  }
  
  void ProjectionAlgorithm::SetLayerNumber(int layerNumber) {
    imageParameter_.layer = layerNumber;
  }
  
  void ProjectionAlgorithm::SetBrightPixelCount(int brightPixelCount) {
    imageParameter_.threshold = brightPixelCount;
  }
  
  void ProjectionAlgorithm::SetVerbose(bool verbose) {
    imageParameter_.verbose = verbose;
  }
  
  void ProjectionAlgorithm::SetMaximumInterpolation(bool maximumInterpolation) {
    imageParameter_.maxInterpolation = maximumInterpolation;
  }
  
  void ProjectionAlgorithm::SetPrintHeightMap(int printHeightMap) {
    imageParameter_.printHeightMap = printHeightMap;
  }
  
  void ProjectionAlgorithm::SetPrintHeightMapRealSize(int printHeightMapRealSize) {
    imageParameter_.printRealHeightMap = printHeightMapRealSize;
  }
  
  
  /***********************
   * PUBLIC: Projection Function
   **********************/
  
  void ProjectionAlgorithm::ProjectImageStack() {
    
    if (imageParameter_.verbose) {
      std::cout << "Read in " << (*image_).GetWidth() << "x" << (*image_).GetHeight() << "x" << (*image_).GetDepth() << " image" << std::endl;
      std::cout << "Image type: " << (*image_).GetType() << " (0=uint8, 1=uint16)" << std::endl;
      std::cout << "Use an " << imageParameter_.radius << "x" << imageParameter_.radius << " window (initial grid size)" << std::endl;
      std::cout << "Distance for second selection step: " << imageParameter_.distance << std::endl;
      std::cout << "Looking at the top " << imageParameter_.threshold << " pixels" << std::endl;
    }
    
    if (parametersAreSet_ == false ) {
      std::cout << "Error: Parameters are not correctly set!" << std::endl;
      return;
    }
    
    if ((*image_).GetType() == UINT8_TYPE) {
      Project_8bit();
    } else if ((*image_).GetType() == UINT16_TYPE) {
      Fast_Project_16bit();
    } else {
      std::cout << "This image type is not supported (only uint8 and unit16)!" << std::endl;
      return;
    }
  }
  
  /***********************
   * PUBLIC: Draw interpolated height map
   **********************/
  
  void ProjectionAlgorithm::DrawInterpolatedHeightMap(std::string outputPath) {
    
    for (int i=0; i<imageParameter_.layer; i++) {
      std::string outputImgPath = outputPath + "_r" + std::to_string(2*imageParameter_.radius) + "_d" + std::to_string(imageParameter_.distance) + "_t"+ std::to_string(imageParameter_.threshold) + "_layer" + std::to_string(i) +"_interpolatedHeightMap.tif";
      Write_Image((char *)outputImgPath.data(),interpolatedHeightMaps_[i],DONT_PRESS);
    }
  }
    
  void ProjectionAlgorithm::DrawDownsampledHeightMap(std::string outputPath) {
    
    for (int i=0; i<imageParameter_.layer; i++) {
      std::string outputImgPath = outputPath + "_r" + std::to_string(2*imageParameter_.radius) + "_d" + std::to_string(imageParameter_.distance) + "_t"+ std::to_string(imageParameter_.threshold) + "_layer" + std::to_string(i) +"_downsampledHeightMap.tif";
      Write_Image((char *)outputImgPath.data(),heightMaps_[i],DONT_PRESS);
    }
  }

  void ProjectionAlgorithm::DrawProjection(std::string outputPath) {
    
    for (int i=0; i<imageParameter_.layer; i++) {
      std::string outputImgPath = outputPath + "_r" + std::to_string(2*imageParameter_.radius) + "_d" + std::to_string(imageParameter_.distance) + "_t"+ std::to_string(imageParameter_.threshold) + "_layer" + std::to_string(i) +"_projection.tif";
      
      Write_Image((char *)outputImgPath.data(),projections_[i],DONT_PRESS);
    }
    
  }


  Array * ProjectionAlgorithm::GetProjection (int layerNumber) {
    return Copy_Array(projections_[layerNumber]);
  }
  
  /***********************
   * PRIVATE: Projection Function for 8bit
   **********************/
  
  void ProjectionAlgorithm::Project_8bit (){
    
    Array * maxim = MaxProjection_8bit((*image_), imageParameter_);
    Array * index = MaxIndices_8bit((*image_), imageParameter_, maxim);
    
    // ---------- version 1.0 -------------
    
    // 1. Select reference height map with a big window (radius)
    //      Selecting the z planes using the occurrences of the brightest pixel
    
    std::vector<Array *> levels = PlaneSelection_SeveralLevels_8bit((*image_), imageParameter_, maxim, index);
    Free_Array(index);
    
    imageParameter_.radius = imageParameter_.radius/2;
    (*image_).UpdateGrid(2);
    
    for (int i=0; i<imageParameter_.layer; i++) {
      
      if (imageParameter_.verbose) {
        std::cout << "---- Layer " << i << std::endl;
      }
      
      // 2. Decomposing the compuational window based on the quad tree principle
      //    new window size: 0.5 * radius
      //    Smoothing of the resulting height map using local median filtering
      
      Array * refinedLevel = DecomposeImage(levels[i]);
//      Array * brightPixels = BrightPixelDistribution(image_, imageParameter_, maxim, index);

      int filterSizes2 [] = {4,2,1};
      Array * median = Filtering((*image_), imageParameter_, refinedLevel, filterSizes2, (*image_).GetDepth()/10);
//      Array * median = Filtering_Confidence(image_, imageParameter_, refinedLevel, filterSizes2, image_.GetDepth()/10, brightPixels);
      Free_Array(refinedLevel);
      // 3. Select the final height map using a subset of the image stack
      //      defined by the reference height map and a distance parameter
      
      
      // ------ method 1 ---- variance -------------
      
//      Array * substack = Substack_UsingMask_8bit(image_, imageParameter_, median, false);

      Array * level = SelectPlanes_Variance_8bit((*image_), imageParameter_, median);

//      Array * test = LocalMedianFilter_InclVar(img, imgPara, level, img.GetDepth(), 1);
//      Write_Image("weighted.tif", test, DONT_PRESS);
      
      //median = Copy_Array(level);
      int filterSizes3 [] = {2,1};
      //median = Filtering(img, imgPara, level, img.GetDepth(), filterSizes3);
      median = LocalMedianFilter((*image_), imageParameter_, level, 1, 1);

      
      if (imageParameter_.range) {
        GetLevelRange(median);
      }                   
      
      
      Array * corners = InterpolateCorners((*image_), imageParameter_, median);
      ResultArrays result = InterpolatePlanes_8bit((*image_), imageParameter_, maxim, corners);
      
      heightMaps_.push_back(median);
      interpolatedHeightMaps_.push_back(result.interpolatedHeightMap);
      projections_.push_back(result.projection);
      
//      imageParameter_.radius = imageParameter_.radius*2;
//      image_.UpdateGrid(0.5);

      
      Free_Array(level);
      Free_Array(corners);
    }
    
//    imageParameter_.radius = imageParameter_.radius/2;
//    image_.UpdateGrid(2);
   
    for (int i=0; i<imageParameter_.layer; i++) {
      Free_Array(levels[i]);
    }
    
    Free_Array(maxim);
  }
  
  
  /*************************************************************
   *
   *    Projection Algorithm for a 16-bit image
   *
   *************************************************************/
  
  
  void ProjectionAlgorithm::Project_16bit (){
    
    Array * maxim = MaxProjection_16bit((*image_), imageParameter_);
    Array * index = MaxIndices_16bit((*image_), imageParameter_, maxim);
    
    Levels levels = PlaneSelection_2Levels_16bit((*image_), imageParameter_, maxim, index);
    //std::vector<Array *> levels = PlaneSelection_SeveralLevels_16bit((*image_), imageParameter_, maxim, index);

    Free_Array(index);
    
    for (int i=0; i<imageParameter_.layer; i++) {
      
      std::cout << "---- Layer " << i << std::endl;
    
      Array * refinedLevel = DecomposeImage(levels.levels[i]);
    
      imageParameter_.radius = imageParameter_.radius/2;
      (*image_).UpdateGrid(2);
    
      int filterSizes [] = {4, 2, 1};
      Array * median = Filtering((*image_), imageParameter_, refinedLevel, filterSizes, 1);
      Free_Array(refinedLevel);

      Array * substack = Substack_UsingMask_16bit((*image_), imageParameter_, median, true);
      Free_Array(median);
      Array * level = SelectPlanes_Variance_16bit((*image_), imageParameter_, substack);
      Free_Array(substack);
    
      median = Copy_Array(level);
      int filterSizes2 [] = {2,1};
      median = Filtering((*image_), imageParameter_, level, filterSizes2, 1);
      

      if (imageParameter_.range) {
        GetLevelRange(median);
      }

    
      Array * corners = InterpolateCorners((*image_), imageParameter_, median);
      ResultArrays result;
    
      if (imageParameter_.maxInterpolation) {
        result = InterpolatePlanes_16bit_MaxInterpolation((*image_), imageParameter_, maxim, corners);
      } else {
        result = InterpolatePlanes_16bit((*image_), imageParameter_, maxim, corners);
      }
    
      heightMaps_.push_back(median);
      interpolatedHeightMaps_.push_back(result.interpolatedHeightMap);
      projections_.push_back(result.projection);
      
      imageParameter_.radius = imageParameter_.radius*2;
      (*image_).UpdateGrid(0.5);
      
      Free_Array(level);
      Free_Array(corners);
    }
    
    imageParameter_.radius = imageParameter_.radius/2;
    (*image_).UpdateGrid(2);
    
    for (int i=0; i<imageParameter_.layer; i++) {
      Free_Array(levels.levels[i]);
    }
    Free_Array(levels.confidence);
    Free_Array(maxim);
  }
  
  /*************************************************************
   *
   *    Fast Projection Algorithm for a 16-bit image
   *
   *************************************************************/
  void ProjectionAlgorithm::Fast_Project_16bit (){
        
    //(*image_).SetGrid(imageParameter_.radius);
        
        Array * img_16bit = Copy_Array((*image_).GetImage());
        Array * maxim_16bit = MaxProjection_16bit((*image_), imageParameter_);
        
        
        // Works on 8-bit to figure out the height map
        (*image_).ScaleTo8bit();
        
        Array * maxim = MaxProjection_8bit((*image_), imageParameter_);
        Array * index = MaxIndices_8bit((*image_), imageParameter_, maxim);
      
//        Levels levels = PlaneSelection_2Levels_16bit((*image_), imageParameter_, maxim, index);     // #TODO: check 16 / 8 bit
        std::vector<Array *> levels = PlaneSelection_SeveralLevels_8bit((*image_), imageParameter_, maxim, index);
        Free_Array(index);
    
        imageParameter_.radius = imageParameter_.radius/2;
        (*image_).UpdateGrid(2);
    
        for (int i=0; i<imageParameter_.layer; i++) {
          
          //std::cout << "---- Layer " << i << std::endl;
        
          Array * refinedLevel = DecomposeImage(levels[i]);
          
         
        
          int filterSizes [] = {4, 2,1};
          Array * median = Filtering((*image_), imageParameter_, refinedLevel, filterSizes, (*image_).GetDepth()/10);
//          Array * median = Filtering_Confidence((*image_), imageParameter_, refinedLevel, filterSizes, 1, levels.confidence);
          Free_Array(refinedLevel);
        
          Array * level = SelectPlanes_Variance_8bit((*image_), imageParameter_, median);
        
          median = LocalMedianFilter((*image_), imageParameter_, level, 1, 1);
        
          if (imageParameter_.range) {
            GetLevelRange(median);
          }
        
          // Projections happens on 16-bit
          Free_Array((*image_).GetImage());
          (*image_).SetImage(img_16bit);
          img_16bit = Copy_Array((*image_).GetImage());
          
          Array * corners = InterpolateCorners((*image_), imageParameter_, median);
          ResultArrays result;

          if (imageParameter_.maxInterpolation) {
            if (imageParameter_.verbose) {
                std::cout << "Interpolation using a maximum projection" << std::endl;
            }
            result = InterpolatePlanes_16bit_MaxInterpolation((*image_), imageParameter_, maxim_16bit, corners);
          } else {
            result = InterpolatePlanes_16bit((*image_), imageParameter_, maxim_16bit, corners);
          }
          
          heightMaps_.push_back(median);
          interpolatedHeightMaps_.push_back(result.interpolatedHeightMap);
          projections_.push_back(result.projection);
          
          (*image_).ScaleTo8bit();
        
          Free_Array(level);
          Free_Array(corners);
        }
    
        for (int i=0; i<imageParameter_.layer; i++) {
          Free_Array(levels[i]);
        }
    
        Free_Array(maxim);
        Free_Array(maxim_16bit);
        Free_Array(img_16bit);
    }
  
  
  void ProjectionAlgorithm::SetParameters (ProjectionMethod::ParaSet & parameter) {
    
    imageParameter_ = parameter;
    (*image_).SetGrid(imageParameter_.radius);
    parametersAreSet_ = true;
  }


  
}