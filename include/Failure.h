#ifndef Failure_h
#define Failure_h
#include <utility>

void failureB(float xp, float yp, float r, float alpha, int im, float& fos);

void failure(float xp, float yp, float r, float alpha, int im, float &fos);

void setting_elem_and_tipn_failure( int m, int im);

void NewFractureCentralPoint(float xp, float yp, float& sigxx, float& sigyy, float& sigxy);

void Sum_Failure(float xp, float yp, float r, float alpha, int im, float fos, int location);
void Choose_Failure();


#endif