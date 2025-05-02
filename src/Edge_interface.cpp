#include<stdafx.h>

#include "Edge_interface.h"


Edge_interface::Edge_interface():x_beg(0), y_beg(0), x_end(0), y_end(0), pos_mat_no(1), neg_mat_no(1), elem_no(1) {}

Edge_interface::Edge_interface(float x1, float y1, float x2, float y2, 
	int mat1, int mat2,int eno):
	x_beg(x1),y_beg(y1),x_end(x2),y_end(y2),pos_mat_no(mat1),neg_mat_no(mat2),elem_no(eno){}

int Edge_interface::get_elemno() { return elem_no; }
int Edge_interface::get_pos_matno() { return pos_mat_no; }
int Edge_interface::get_neg_matno() { return neg_mat_no; }
float Edge_interface::get_xbeg() { return x_beg; }
float Edge_interface::get_ybeg() { return y_beg; }
float Edge_interface::get_xend() {	return x_end;}
float Edge_interface::get_yend() { return y_end; }


