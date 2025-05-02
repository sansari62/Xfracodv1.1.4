#ifndef S4_h
#define  S4_h


class S4
{
public:
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
	S4(int length); 
	void limit_d();     // limit the d() value for stability
	
	void save_to_file(ofstream& f);
	
};



#endif // !S4_h

