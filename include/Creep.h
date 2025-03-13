#ifndef Creep_h
#define Creep_h
#include<fstream>

class Creep
{
public:
	float creep_x[500];
	float creep_y[500];
	float creep_l[500];
	float creep_a[500];
	float growth_length[500];

	float time;
	float time0;
	float deltaT;
	float deltaT_min;
	float deltaT_max;
	float totalT;
	float v1;
	int nn1;
	float v2;
	int nn2;
	
	float vel_creep_max;
	int ID_creep;
	int ID_fast_crack;

	Creep();

	void save_to_file(std::ofstream& f);

};
#endif

