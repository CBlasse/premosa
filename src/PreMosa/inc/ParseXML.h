//
//  ParseXML.h
//  CMakeProject
//
//  Created by blasse on 7/25/13.
//
//

#ifndef PARSEXML_H
#define PARSEXML_H

#include <QDomDocument>
#include <QDomElement>
#include <QDomText>
#include <QString>

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
  XmlFile (QString xmlFileName);
  ~XmlFile ();
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parse of the different parameters (required and non-required (default value available))
   */
  void ParseRequieredParameter(QDomElement xmlElements);
  void ParseNonRequieredParameter(QDomElement xmlElements);
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Return of the different parameters
   */
  bool GetComplete () const;
  bool GetPrintHeightMap () const;
  bool GetIntermediateFileDeletion () const;
  QString GetInputDir () const;
  QString GetOutputDir () const;
  QString GetFolderNameExpression () const;
  QString GetFileNameExpression () const;
  QString GetScriptLocation () const;
  QString GetFiji () const;
  QString GetFlatField () const;
  QString GetMasterTileConfig () const;
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
  
  private:
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parameter
   */
  bool complete_;
  bool printHeightMap_;
  bool intermediateFileDeletion_;
  QString inputDir_;
  QString outputDir_;
  QString folderNameExpression_;
  QString fileNameExpression_;
  QString fiji_;
  QString flatField_;
  QString scriptLocation_;
  QString masterTileConfig_;
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
