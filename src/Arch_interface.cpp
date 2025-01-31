#include "Arch_interface.h"

Arch_interface::Arch_interface() {}

Arch_interface::Arch_interface(float xc, float yc, float r, float ang1, float ang2,
	int eno, int mat1, int mat2):
	x_cent(xc), y_cent(yc),diameter(r),beg_angle(ang1),end_angle(ang2),
	elem_no(eno),pos_mat(mat2),neg_mat(mat1){}

int Arch_interface::get_elenum()   { return elem_no; }
float Arch_interface::get_diameter() { return diameter; }
float Arch_interface::get_xcent()   { return x_cent; }
float Arch_interface::get_ycent()   { return y_cent; }
float Arch_interface::get_beg_ang() { return beg_angle; }
float Arch_interface::get_end_ang() { return end_angle; }
int Arch_interface::get_neg_mat() { return neg_mat; }
int Arch_interface::get_pos_mat() { return pos_mat; }
