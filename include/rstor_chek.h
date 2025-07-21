#pragma once
#ifndef rstor_chek_h
#define rstor_chek_h
#include <Remov_preexisting_BEs.h>



void check_boreholes();
void addClippedElement(const Point& a, const Point& b,
    int m, int& new_numbe,int);
int  if_tip_element(int m);
void addClippedElement2(int m, int new_numbe, BoundaryElement& new_elem, int first);


#endif // !rstor_chek_h