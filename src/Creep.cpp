#include "Creep.h"

Creep::Creep() : creep_a{0},creep_x{0},creep_l{0},creep_y{0}, growth_length{0}, time(0),
time0(0), deltaT_max(1), deltaT(1), deltaT_min(1), totalT(0), v1(0), nn1(0), v2(0), nn2(0), 
vel_creep_max(0), ID_creep(0), ID_fast_crack(0){}



void Creep::read_from_file(ifstream& f)
{
	for (int m = 0; m < 500; ++m) {
		f >> creep_x[m] >> creep_y[m] >> creep_l[m] >> creep_a[m] >> growth_length[m];
	}
		f >> time >> time0 >> deltaT >> deltaT_min >> deltaT_max >> totalT >> v1 >> nn1 >> v2
			>> nn2 >> vel_creep_max >> ID_creep >> ID_fast_crack;
	return;
}



void Creep::save_to_file(ofstream& f)
{
	for (int m = 0; m < 500; ++m) 
	{
		f <<  creep_x[m] << " " << creep_y[m] << " " << creep_l[m] << " " << creep_a[m] << " " <<
			growth_length[m] << std::endl;
	}

	f << time << " " << time0 << " " << deltaT << " " << deltaT_min << " " << deltaT_max << " " << totalT <<
		" " << v1 << " " << nn1 << " " << v2
		<< " " << nn2 << " " << vel_creep_max << " " << ID_creep << " " << ID_fast_crack << std::endl;

	return;
}

