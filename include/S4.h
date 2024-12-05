#ifndef S4_h
#define  S4_h
#include <common.h>


class S4
{
public:
	  //what is the relation between these parameter? what thet keep and show
	std::vector<float> b;
	std::vector<float> d;
	std::vector<float> b1;
	std::vector<float> d0;
	std::vector<float> df0;
	std::vector<float> df;
	std::vector<std::vector<float>> c_d;
	std::vector<std::vector<float>> c_s;
	std::vector<std::vector<float>> c;

	std::vector<float> b0;
	std::vector<float> b0_old;
	//Sara think about the length and size
	S4(int length); 
	void limit_d();     // !limit the d() value for stability
	
	void read_from_file(ifstream& f);
	void save_to_file(ofstream& f);
	
};









#endif // !S4_h

