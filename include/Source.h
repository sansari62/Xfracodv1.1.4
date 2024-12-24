#ifndef Source_h
#define Source_h

#include <Rock.h>


int cross(float xb1, float yb1, float xe1, float ye1, float xb2, float yb2, float& xe2, float& ye2);

void Central_control(wstring);
int  check_material_id(float xp, float yp);

void Interface(int num, float xbeg, float ybeg, float xend, float yend, int mat1, int mat2);

#endif
	