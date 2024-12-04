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
	/*float sigma_s;
	float sigma_n;
	
	int ipair;   
	int jstate;  
	int jslipd;
	int jmode; 
	float force1;
	float force2;	
	float us;
	float un;
	float forces;
	float forcen;
	float aks;
	float akn;
	float phi;
	float phid;
	float coh;
	float us_neg;
	float un_neg;
	float ss_old;
	float sn_old;*/
	

	BoundaryElement();	
	BoundaryElement(float x, float y, float am, float cos, float sin,
		int kd ,int mat);
	void read_from_file(ifstream& f,int m);
	void save_to_file(ofstream& f, int m);
	void bound(int i, float& ss, float& sn, float& ustem, float& untem, float& usneg, 
		float& unneg);

};
void  MatrixB(int m);   //Sara! not sure to be part of class or not
#endif

