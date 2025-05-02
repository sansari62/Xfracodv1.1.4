#pragma once
#ifndef Borehole_h
#define Borehole_h
#include<GeologicalForm.h>
#include <fstream>


class Borehole:public GeologicalForm
{
//define a tunnel
public:
	
	float diameter[10];
	float x_cent[10];
	float y_cent[10];	
	float xpp[10];
	float ypp[10];
	float dd[10];
	
	Borehole();	
	
	void save_to_file(ofstream& f);
	


};


#endif

