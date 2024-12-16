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


	MonitoringLine() : x1l(0), y1l(0), x2l(0), y2l(0), npl(1) {}

	void save_to_file(ofstream& f)
	{
		f << x1l << " " << y1l << " " << x2l << " " << y2l << " " << npl << std::endl;
	}
};

#endif