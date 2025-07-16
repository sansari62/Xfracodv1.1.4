#pragma once
#ifndef Remov_preexisting_BEs_h
#define   Remov_preexisting_BEs_h

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
    //float b01, b02;, b01(0.0f), b02(0.0f)
    int tip_indx;
    int s4_indx;
    new_element_para() : tip_indx(-1), ratio(1), s4_indx(-1){}
    
};
void  update_s4(int m, int new_numbe);
void fix_tip_pointer1(int new_numbe, int tip_indx, int direc);

//void setting_newelement(BoundaryElement& newelement, int m, int new_numbe, int p1_stat, int p2_stat);

bool check_segmnt_valid(const Point& a, const Point& b,int m, BoundaryElement& newelem);
void one_point_intersection(Point outside_pt, Point clipped, int m, int& new_numbe, int tip_index,
    std::vector<new_element_para>& vec);
void two_intersection_case(Point p1, Point p2, Point ip1, Point ip2, int m, int& new_numbe, int tip_index, 
    std::vector<new_element_para>& vec);
#endif // !Remov_preexisting_BEs_h