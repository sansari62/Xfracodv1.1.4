#ifndef Fmax_h
#define Fmax_h


void newcoord();

void fmax1(float& f0, float& angle);

float ang_setting(float angi0, float& fi0, int mm, int mode, float& angi);

void compute_f(float& f0, float fi0, float angi, float angii, float& angle, float fii0, int m, float wi, float wii);

float call_work1_setting_fi0(float dtt, float& fi0, int mm, int mode, float& angip, float ang);

int  check_elastic_growth(int material_no);



#endif