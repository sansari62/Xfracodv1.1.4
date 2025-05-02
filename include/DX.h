#pragma once
#ifndef DX_h
#define DX_h

float dxi(float& angle, bool flagB);
int cross(float xb1, float yb1, float xe1, float ye1, float xb2, float yb2, float& xe2, float& ye2);
void NewFractureCentralPoint(float xp, float yp, float& sigxx, float& sigyy, float& sigxy);



#endif 

