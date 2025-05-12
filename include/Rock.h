
#ifndef Rock_h
#define Rock_h


class Rock
{
public:	
	float pr;						//poisson_rate
	float e;						//young_modulus
	float viscosity;
	float density;
	float perm0;                       //permeability 
	float rcoh;					//intact rock cohesion
	float rst;					    //tensile strength;
	float rphi;                    //interanl fric_angle;
	float akic;
	float akiic;
	

	
	Rock();
	void save_to_file(ofstream & f);
	void read_from_file(ifstream& f);

};

#endif


