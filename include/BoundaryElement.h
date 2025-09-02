#pragma once
#include<stdafx.h>

#ifndef BoundaryElement_h
#define BoundaryElement_h

/* Element coordinates, length, orientation, joint status, joint properties,
and joint surface stresses and displacement*/

class BoundaryElement
{
public:
	int mat_no;
	float xm;
	float ym;
	float sinbet;
	float cosbet;
	float a;
	int kod;
	int frac_id;

	BoundaryElement();	
	BoundaryElement(float x, float y, float am, float cos, float sin,
		int kd ,int mat);	
	void save_to_file(std::ofstream& f, int m);
	void read_from_file(ifstream& f,int m);
	void bound(int i, float& ss, float& sn, float& ustem, float& untem, float& usneg, 
		float& unneg);
};
void  MatrixB(int m);   
#endif

