#pragma once
#ifndef Elliptical_opening_h
#define Elliptical_opening_h 
#include<GeologicalForm.h>
using namespace std;


//define an elleptical opening 
class Elliptical_opening:public GeologicalForm
{

public:
		float x_cent;
		float y_cent;
		float diameter1;
		float diameter2;
		//float bvs;
		//float bvn;
		//float grad_sy;
		//float grad_ny;


		Elliptical_opening();
		Elliptical_opening(float x, float y, float r1, float r2);
	//any need to init other paras? Sara!
		Elliptical_opening(int mat, int eno, int kode, int x,int y,float r1,
			float r2,float bs,float bn,float grad1,float grad2);
		void def_arch_boundary(int elem_num, float x, float y, float r1, float r2, float ang1, float ang2);
		void increment_bnd_stress_along_arch(float x, float y, float r1, float r2, float ang1, float ang2, float dss, float dnn);

		
		void save_to_file(ofstream & f);


};
#endif

