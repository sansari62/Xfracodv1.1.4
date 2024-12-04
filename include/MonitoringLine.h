#ifndef MonitoringLine_h
#define MonitoringLine_h

class MonitoringLine
{
public:
	float x1l;
	float y1l;
	float x2l;
	float y2l;
	int npl;


	MonitoringLine();


void read_from_file(ifstream& f);
void save_to_file(ofstream& f);

};



inline MonitoringLine::MonitoringLine(): x1l(0), y1l(0), x2l(0), y2l(0), npl(1){}

inline void MonitoringLine::read_from_file(ifstream& f)
{
	f >> x1l >> y1l>> x2l >> y2l >> npl;
}

inline void MonitoringLine::save_to_file(ofstream& f)
{
	f << x1l << " " << y1l << " " << x2l << " " << y2l << " " << npl << std::endl;
}
#endif