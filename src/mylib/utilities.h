/*****************************************************************************************\
*                                                                                         *
*  Utilities for allocating memory, opening files, and processing command line arguments  *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  October 2005                                                                  *
*                                                                                         *
\*****************************************************************************************/

#ifndef _SR_UTILITIES

#define _SR_UTILITIES

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef _MSC_VER

#define strdup(s)  _strdup(s)
#define inline     __inline

#endif

#define ASCII 128

void *Guarded_Malloc(size_t size, char *routine);
void *Guarded_Realloc(void *array, size_t size, char *routine);
char *Guarded_Strdup(char *string, char *routine);
FILE *Guarded_Fopen(char *name, char *options, char *routine);

void Process_Arguments(int argc, char *argv[], char *spec[], int no_escapes);

char  *Program_Name();

int    Get_Repeat_Count(char *name);
int    Is_Arg_Matched(char *name, ... /* [int no] */ );

int    Get_Int_Arg   (char *name, ... /* [int no [int an]] */ );
double Get_Double_Arg(char *name, ... /* [int no [int an]] */ );
char  *Get_String_Arg(char *name, ... /* [int no [int an]] */ );

void Print_Argument_Usage(FILE *file, int no_escapes);

#ifdef __cplusplus
}
#endif

#endif
