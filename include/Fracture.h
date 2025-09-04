#pragma once
#ifndef Fracture_h
#define Fracture_h


#include<GeologicalForm.h>

class Fracture:public GeologicalForm
{	
	
public:
	float x_beg;
	float y_beg;
	float x_end;
	float y_end;	
	int jmat;
	int id;


	Fracture() {};
	Fracture(float x1, float y1, float x2, float y2, int mat, int eno, int jmat1, int fid);
	
	void frac_reassign_values(int eno, float x1, float y1, float x2, float y2, int kode, int mat, int jmat1,int id);

	void take_xbeg(float x);
	void take_ybeg(float y);
	void take_xend(float x);
	void take_yend(float y);
	void take_xy_beg(float x,float y);
	void take_xy_end(float x, float y);

	float get_xbeg() const;
	float get_ybeg() const;
	float get_xend() const;
	float get_yend() const;	


};
#endif