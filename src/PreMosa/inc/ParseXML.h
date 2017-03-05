//
//  ParseXML.h
//  CMakeProject
//
//  Created by blasse on 7/25/13.
//
//

#ifndef PARSEXML_H
#define PARSEXML_H

#include <string>

#include "tinyxml.h"

namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Class to parse parameters from a XML file
   */
  
  class XmlFile
  {
  
  public:
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Constructor / Destructor
   */
  XmlFile ();
  XmlFile (std::string xmlFileName);
  ~XmlFile ();
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parse of the different parameters (required and non-required (default value available))
   */
  void ParseRequieredParameter(TiXmlDocument & xmlDoc);
  void SetDefaultParameters ();
  void VerifyRequiredParameters (bool hasStitchedTag);

  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Return of the different parameters
   */
  bool GetComplete () const;
  bool GetPrintHeightMap () const;
  bool GetIntermediateFileDeletion () const;
  std::string GetInputDir () const;
  std::string GetOutputDir () const;
  std::string GetFolderNameExpression () const;
  std::string GetFileNameExpression () const;
  std::string GetScriptLocation () const;
  std::string GetFiji () const;
  std::string GetFlatField () const;
  std::string GetMasterTileConfig () const;
  int GetTStart () const;
  int GetTEnd () const;
  int GetGridSizeX () const;
  int GetGridSizeY () const;
  int GetFirstIndex () const;
  int GetOverlap () const;
  int GetVerbose () const;
  int GetRadius () const;
  int GetDistance () const;
  int GetThreshold () const;
  int GetMinIntensity () const;
  int GetMaxIntensity () const;
  int GetNumberOfSampledPairs () const;
  double GetLambda1() const;
  double GetLambda2() const;
  
  void SetMinIntensity (int minIntensity);
  void SetMaxIntensity (int maxIntensity);
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parameter
   */
  bool complete_;
  bool printHeightMap_;
  bool intermediateFileDeletion_;
  std::string inputDir_;
  std::string outputDir_;
  std::string folderNameExpression_;
  std::string fileNameExpression_;
  std::string fiji_;
  std::string flatField_;
  std::string scriptLocation_;
  std::string masterTileConfig_;
  int tStart_;
  int tEnd_;
  int gridSizeX_;
  int gridSizeY_;
  int firstIndex_;
  int overlap_;
  int verbose_;
  int minIntensity_;
  int maxIntensity_;
  int radius_;
  int distance_;
  int threshold_;
  int numberOfSampledPairs_;
  double lamdba1_;
  double lamdba2_;
  };
  
  
}



#endif
