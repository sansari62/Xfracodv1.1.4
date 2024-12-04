#pragma once
#ifndef Arch_interface_h
#define Arch_interface_h

class Arch_interface
{
	private:
		int elem_no;			//no of elements for arc interface
		float diameter;
		float x_cent;
		float y_cent;
		float beg_angle;
		float end_angle;
		int neg_mat;
		int pos_mat;

	public:
		Arch_interface();
		Arch_interface(float xc, float yc,float r,float ang1, float ang2,int eno, int mat1, int mat2);
		int get_elenum();
		float get_diameter();
		float get_xcent();
		float get_ycent();
		float get_beg_ang();
		float get_end_ang(); 
		int get_neg_mat();
		int get_pos_mat();
};

#endif





