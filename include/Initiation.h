#pragma once
#ifndef Initiation_h
#define Initiation_h
#include<utility>


std::pair<bool, bool> for_j_loop(int k, int m, float xt, float yt, float& x, float& y,
    int& ID, float& beta, int &jj, int numbe0,bool intern_call);

void reassigning_boundary_values(int ID, int m, int j, int k, float beta, float x, float y,
    float& st, float& bet);

bool for_i_loop(int k, float x, float y, int m, int numbe0);

void InitiationB();

void InitiationR();
void initiation();





#endif
