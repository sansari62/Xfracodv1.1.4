#pragma once
#ifndef rstor_chek_h
#define rstor_chek_h

#include<CommonPara.h>

using namespace CommonPara_h::comvar;

extern bool second_clipped;
struct Point {
    float x, y;

    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }
};

struct new_element_para
{
    BoundaryElement new_el;
    BE be1;
    float ratio;
    Joint j;
    float b01, b02;
    int tip_indx;
    new_element_para() : tip_indx(-1), ratio(1),b01(0.0f),b02(0.0f) {}
};


void check_boreholes();
void addClippedElement(const Point& a, const Point& b,
    int m, int& new_numbe,int);
int  if_tip_element(int m);
void addClippedElement2(int m, int new_numbe, BoundaryElement& new_elem, int first);


#endif // !rstor_chek_h