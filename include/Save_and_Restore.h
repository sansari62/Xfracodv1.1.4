#pragma once

#ifndef Save_and_Restore_h
#define Save_and_Restore_h

#include<fstream>
extern std::vector<int> jmat_list; 

void save(ofstream& file10);
void restore(string fname);

#endif

