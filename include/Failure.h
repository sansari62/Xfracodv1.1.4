#ifndef Failure_h
#define Failure_h
#include <utility>

void failureB(float xp, float yp, float r, float alpha, int im);

void failure(float xp, float yp, float r, float alpha, int im);

void setting_elem_and_tipn_failure( int m, int im);


void Sum_Failure(float xp, float yp, float r, float alpha, int im, float fos, int location);
void Choose_Failure();
int  check_material_id(float xp, float yp);
int CheckNewElement(float ac, float xc, float yc, float cosbeta, float sinbeta, int flagB);



#endif