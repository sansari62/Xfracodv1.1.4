#ifndef MonitoringPoint_h
#define MonitoringPoint_h


class MonitoringPoint
{
public:
	float xmon;
	float ymon;

	
	MonitoringPoint() : xmon(0.0), ymon(0.0) {}	

	void save_to_file(ofstream& f)
	{
		f << xmon << " " << ymon << std::endl;
	}

};

#endif