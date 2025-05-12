#include<stdafx.h>

#include "Borehole.h"
#include<CommonPara.h>

using namespace comvar;



Borehole::Borehole() :diameter{0.0}, x_cent{ 0.0 }, y_cent{ 0.0 }, xpp{ 0.0 }, ypp{ 0.0 }, 
dd{ 0.0 } {}





void Borehole::read_from_file(ifstream& f)
{	
	for (int m = 0; m < ntunnel; ++m)
	{
		f >> diameter[m] >> x_cent[m] >> y_cent[m]
			>> xpp[m] >> ypp[m] >> dd[m];
	}
}





void  Borehole::save_to_file(ofstream& f)
{
	f<< ntunnel << std::endl;
	for (int m = 0; m < ntunnel; ++m)
	{
		f << diameter[m] << " " << x_cent[m] << " " << y_cent[m]
			<< xpp[m] << " " << ypp[m] << " " << dd[m] << std::endl;
	}
	//Sara not sure about 10 or 20
	/*for (int m = 0; m < 10; ++m)
	{
		f << xpp[m] << " " << ypp[m] << " " << dd[m] << std::endl;
	}*/

}
