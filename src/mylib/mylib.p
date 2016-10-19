/*****************************************************************************************\
*                                                                                         *
*  Central library declarations and initialization routine                                *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  January 2007                                                                  *
*                                                                                         *
*  (c) June 19, '09, Dr. Gene Myers and Howard Hughes Medical Institute                   *
*      Copyrighted as per the full copy in the associated 'README' file                   *
*                                                                                         *
\*****************************************************************************************/

#include <stdlib.h>

#define BND_ZERO    0
#define BND_REFLECT 1
#define BND_WRAP    2
#define BND_EXTEND  3
#define BND_INVERT  4

int Boundary_Case_8qm5 = BND_ZERO;

void Use_Zero_Boundary()
{ Boundary_Case_8qm5 = BND_ZERO; }

void Use_Reflective_Boundary()
{ Boundary_Case_8qm5 = BND_REFLECT; }

void Use_Wrap_Boundary()
{ Boundary_Case_8qm5 = BND_WRAP; }

void Use_Extend_Boundary()
{ Boundary_Case_8qm5 = BND_EXTEND; }

void Use_Inversion_Boundary()
{ Boundary_Case_8qm5 = BND_INVERT; }
