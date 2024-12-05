#ifndef MonitoringPoint_h
#define MonitoringPoint_h


class MonitoringPoint
{
public:
	float xmon;
	float ymon;


	MonitoringPoint();
	void read_from_file(ifstream& f);
	void save_to_file(ofstream& f);


};




inline MonitoringPoint::MonitoringPoint(): xmon(0.0), ymon(0.0){}

inline void MonitoringPoint::read_from_file(ifstream& f)
{
	f>> xmon >> ymon;
}
inline void MonitoringPoint::save_to_file(ofstream& f)
{
	f << xmon << " " << ymon << std::endl;
}
#endif