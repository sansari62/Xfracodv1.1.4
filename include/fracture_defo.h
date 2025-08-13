#pragma once
#ifndef fracture_defo_h
#define fracture_defo_h

#include <WinInterface.h>
#include <fstream>


using namespace WinInterface_h::winvar;


Wjoint compute_join_attr_and_add_them(int m, int& jpoint, int round);
void fracture_defo( int& jpoint);


#endif 
