//
//  ProjectionAlgorithm.h
//  CMakeProject_4.0
//
//  Created by blasse on 8/26/13.
//
//

#ifndef PROJECTIONALGORITHM_H
#define PROJECTIONALGORITHM_H

#include <memory>
#include <string>
#include <vector>

extern "C" {
#include "array.h"
}

#include "Parameter.h"
#include "RasterizedImage.h"

namespace ProjectionMethod {
  
  class ProjectionAlgorithm {
    
    
  public:
    ProjectionAlgorithm(std::string imgPath);
    ProjectionAlgorithm(Array * image);
    ProjectionAlgorithm();
    ~ProjectionAlgorithm();
    
    RasterizedImage GetImage() const;
    int GetGridSize() const;
    int GetDistance() const;
    int GetLayerNumber() const;
    int GetBrightPixelCount() const;
    bool IsVerbose() const;
    bool UseMaximumInterpolation() const;
    bool PrintHeightMap() const;
    bool PrintHeightMapRealSize() const;
    ParaSet GetParameterSet () const;

    
    void SetGridSize(int gridSize);
    void SetDistance(int distance);
    void SetLayerNumber(int layerNumber);
    void SetBrightPixelCount(int brightPixelCount);
    void SetVerbose(bool verbose);
    void SetMaximumInterpolation(bool maximumInterpolation);
    void SetPrintHeightMap(int printHeightMap);
    void SetPrintHeightMapRealSize(int printHeightMapRealSize);
    void SetParameters (ParaSet & parameter);
    
    void ProjectImageStack();
    void DrawInterpolatedHeightMap (std::string outputPath);
    void DrawDownsampledHeightMap (std::string outputPath);
    void DrawProjection (std::string outputPath);

    Array * GetProjection (int layerNumber);
    
  protected:
    
    std::shared_ptr<RasterizedImage> image_;
//    RasterizedImage image_;
    bool parametersAreSet_;
    ParaSet imageParameter_;
    std::vector<Array *> heightMaps_;
    std::vector<Array *> interpolatedHeightMaps_;
    std::vector<Array *> projections_;

    
    Array * RescaleHeightMap (Array * heightMap);
    void Project_8bit ();
    void Project_16bit ();
    void Fast_Project_16bit ();
    
  };
}

#endif 