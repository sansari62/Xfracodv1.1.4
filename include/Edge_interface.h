#ifndef Edge_interface_h
#define Edge_interface_h

class Edge_interface
{
private:
	float x_beg;	//coordinations
	float y_beg;
	float x_end;
	float y_end;
	int pos_mat_no;  //materail number for positive side
	int neg_mat_no;  //materail number for negative side
	int elem_no;		//no of elements for linear interface

public:

	Edge_interface();

	Edge_interface(float x1,float y1,float x2,float y2,int mat1,int mat2,int eno);
	int get_elemno();
	int get_pos_matno();
	int get_neg_matno();
	float get_xbeg();
	float get_ybeg();
	float get_xend();
	float get_yend();

};
#endif