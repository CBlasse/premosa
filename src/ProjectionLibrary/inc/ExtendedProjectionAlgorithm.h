//
//  ProjectionAlgorithm.h
//  CMakeProject_4.0
//
//  Created by blasse on 8/26/13.
//
//

#ifndef EXTENDEDPROJECTIONALGORITHM_H
#define EXTENDEDPROJECTIONALGORITHM_H

#include <string>
#include <vector>

#include "ProjectionAlgorithm.h"
#include "Parameter.h"
#include "RasterizedImage.h"

namespace ProjectionMethod {
  
  class ExtendedProjectionAlgorithm : public ProjectionAlgorithm {
    
  public: 
    ExtendedProjectionAlgorithm(std::string imgPath);
    ExtendedProjectionAlgorithm();
    ~ExtendedProjectionAlgorithm();
    
    void ExtendedProject_8bit ();
    void ExtendedProject_16bit ();
    void ExtendedFast_Project_16bit ();
    
    void SetAllParameters (ProjectionMethod::ExtendedParaSet & parameter);
    void ProjectSurfaceUsingHeightMap (Array * heightMap);
    
//    RasterizedImage GetImage() const;
//    int GetGridSize() const;
//    int GetDistance() const;
//    int GetLayerNumber() const;
//    int GetBrightPixelCount() const;
//    bool IsVerbose() const;
//    bool UseMaximumInterpolation() const;
//    bool PrintHeightMap() const;
//    bool PrintHeightMapRealSize() const;
//    ParaSet GetParameterSet () const;
//
//    
//    void SetGridSize(int gridSize);
//    void SetDistance(int distance);
//    void SetLayerNumber(int layerNumber);
//    void SetBrightPixelCount(int brightPixelCount);
//    void SetVerbose(bool verbose);
//    void SetMaximumInterpolation(bool maximumInterpolation);
//    void SetPrintHeightMap(int printHeightMap);
//    void SetPrintHeightMapRealSize(int printHeightMapRealSize);
//    void SetParameters (int gridSize, int distance, int layerNumber, int brightPixelCount, bool verbose, bool maximumInterpolation, int printHeightMap, int printHeightMapRealSize);
//    
//    void ProjectImageStack();
//    void DrawInterpolatedHeightMap (std::string outputPath);
//    void DrawHeightMap (std::string outputPath);
//    void DrawDownsampledHeightMap (std::string outputPath);
//    void DrawProjection (std::string outputPath);
//
    
  private:
    
    Array * heightMap_;
    Array * rasterizedHeightMap_;
    ExtendedParaSet extendedImageParameter_;
    
    void SurfaceProjection_8bit (Array * heightMap);
    void SurfaceProjection_Fast16bit (Array * heightMap);
  };
}

#endif