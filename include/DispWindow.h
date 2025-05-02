#pragma once

#ifndef DispWindow_h
#define DispWindow_h



class DispWindow
{
public:
	float xll;			// x_leftB;
	float xur;			//x_rightB;
	float yll;         //y_topB;
	float yur;         //y_bottomB;
	int numx;           // x_gridno;
	int numy;			//y_gridno;
	float xc0;         //parameters for circular display window
	float yc0;
	float radium;
	int numr;
	int numa;
	int ID_win;
	
	DispWindow();
	
	void save_to_file(ofstream& f);

};

void CheckRange();

#endif

