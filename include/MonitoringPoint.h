#ifndef MonitoringPoint_h
#define MonitoringPoint_h


#include <fstream>


class MonitoringPoint
{
public:
	float xmon;
	float ymon;

	
	MonitoringPoint() : xmon(0.0), ymon(0.0) {}	

	void save_to_file(std::ofstream& f)
	{
		f << xmon << " " << ymon << std::endl;
	}

};

#endif