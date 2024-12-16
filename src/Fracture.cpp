#include "Fracture.h"


Fracture::Fracture(float x1, float y1, float x2, float y2, int mat, int eno,int jmat1) :
	x_beg(x1), y_beg(y1), x_end(x2), y_end(y2), jmat(jmat1), GeologicalForm(mat,eno,5){}

void Fracture::frac_reassign_values(int eno, float x1, float y1, float x2, float y2,int kode, int mat, int jmat1) 
{
	x_beg = x1;
	y_beg = y1;
	x_end = x2;
	y_end = y2;
	elem_no = eno;
	bound_type = kode ;
	mat_no = mat; 
	jmat = jmat1;
}



float Fracture::get_xbeg() const { return x_beg; }
float Fracture::get_ybeg() const { return y_beg; }
float Fracture::get_xend() const { return x_end; }
float Fracture::get_yend() const { return y_end; }


void Fracture::take_xbeg(float x) { x_beg = x; }
void Fracture::take_ybeg(float y) { y_beg = y; }
void Fracture::take_xend(float x) { x_end = x; }
void Fracture::take_yend(float y) { y_end = y; }
void Fracture::take_xy_beg(float x, float y){ x_beg = x; y_beg = y;}
void Fracture::take_xy_end(float x, float y){ x_end = x;  y_end = y;
}

