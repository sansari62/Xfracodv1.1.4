#ifndef Mainb_h
#define Mainb_h


void coeff(float xi, float yi, float xj, float yj, float aj, float cosbj, float sinbj, 
	int msym, int material);

void solve(int n, int mode);
void mainb_work1_ini();
void mainb_work1(int mode);
void mainb_work0(int mode);
void mainb_ini(int mode, int, int, int, int);

#endif