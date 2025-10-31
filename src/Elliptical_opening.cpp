#include<stdafx.h>
#include "Elliptical_opening.h"


Elliptical_opening::Elliptical_opening():x_cent(0), y_cent(0), diameter1(0), diameter2(0) {}


Elliptical_opening::Elliptical_opening(float x, float y, float r1, float r2) :
	x_cent(x), y_cent(y), diameter1(r1), diameter2(r2) {}




void Elliptical_opening::save_to_file(ofstream& f) 
{
	f << diameter1 << " " << diameter2 << " " << x_cent
		<< " " << y_cent << std::endl;
}

void Elliptical_opening::read_from_file(ifstream& f)
{
	f >> diameter1 >> diameter2 >> x_cent >> y_cent;
}


