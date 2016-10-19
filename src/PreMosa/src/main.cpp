//
//  main.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


#include <iostream>
#include <time.h>


extern "C" {
#include "array.h"
#include "utilities.h"
}

#include "ParseXML.h"
#include "Pipeline.h"

using namespace PreprocessingPipeline;



#define DEVELOP




/****************************************************************************************
 *                                                                                      *
 *  MAIN LEVEL                                                                          *
 *                                                                                      *
 ****************************************************************************************/


int main(int argc, char * argv[])
{
  time_t start,end;
  time (&start);

#ifdef Q_WS_WIN
  std::cout << "This preprocessing pipeline has only been developed for Linux / OsX. It might yield issues on Windows."
  return 0;
#endif
  
  static char *Spec[] = { "[!v] <xml:file>", NULL };
  Process_Arguments(argc,argv,Spec,1);

  {
  XmlFile xmlParamter = * new XmlFile(Get_String_Arg("xml"));     // Parameter parsing
  Pipeline(xmlParamter);                                          // Perform the preprocessing pipeline
  }
  
  Array_List(Free_Array);
  std::cout << "Array " << Array_Usage()  << std::endl;
  
  
  time (&end);    double dif = difftime (end,start);
  std::cout << "Time: " << dif << "s" << std::endl;
    
  
  return 0;
}
  
