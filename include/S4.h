#ifndef S4_h
#define  S4_h


class S4
{
public:
	std::vector<double> b;
	std::vector<double> d;
	std::vector<double> b1;
	std::vector<double> d0;
	std::vector<double> df0;
	std::vector<double> df;
	std::vector<std::vector<double>> c_d;
	std::vector<std::vector<double>> c_s;
	std::vector<std::vector<double>> c;

	std::vector<double> b0;
	std::vector<double> b0_old;
	S4(int length); 
	void limit_d();     // limit the d() value for stability
	
	void save_to_file(ofstream& f);
	void read_from_file(ifstream& f);

};



#endif // !S4_h

