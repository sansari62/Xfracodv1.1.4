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
		

		Elliptical_opening();
		Elliptical_opening(float x, float y, float r1, float r2);		
		void save_to_file(ofstream & f);
		void read_from_file(ifstream& f);

};
#endif

