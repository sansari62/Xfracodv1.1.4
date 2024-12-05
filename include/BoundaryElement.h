#pragma once
#ifndef BoundaryElement_h
#define BoundaryElement_h
#include <common.h>

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
	

	BoundaryElement();	
	BoundaryElement(float x, float y, float am, float cos, float sin,
		int kd ,int mat);
	
	void bound(int i, float& ss, float& sn, float& ustem, float& untem, float& usneg, 
		float& unneg);

};
void  MatrixB(int m);   
#endif

