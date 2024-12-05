#include "Rock.h"
#include<Mainb.h>
#include<Source.h>
#include<CommonPara.h>
using namespace CommonPara_h::comvar;


Rock::Rock():
	e(0),pr(0),akic(0), akiic(0), rcoh(20e6), rphi(30 * pi / 180.0), rst(10e6) {}




void Rock::give_elastic_property(int no, float y1, float rate)
{
	
	e = y1;
	pr = rate;
}



void Rock::give_permeability_property(float v, float dens, float per0)
{
}





int check_elastic_growth(int m)
{

    int kk = 0; // != 0 is not checking elastic growth
    float aki = 0, akii = 0, akie = 0;   

    int material = check_material_id(elm_list[m].xm, elm_list[m].ym );
    int mm = material;

    int ms = 2 * m; //Sara because of index error in d0[2m-1 in line64 add this indexs,make sure it works correctly
    int mn = ms + 1;

    if (s4.d0[mn] > 0) 
    {
        aki = 0;
    }
    else //if (s4.d0[mn] <= 0) 
    {
        aki = abs( rock1[mm].e / (8 * pi * (1 - (rock1[mm].pr * rock1[mm].pr))) *
            sqrt(2 * pi / elm_list[m].a) * s4.d0[mn]) * 2.5; // 2.5 is a correction factor
    }

    akii = abs( rock1[mm].e / (8 * pi * (1 - (rock1[mm].pr * rock1[mm].pr))) * 
        sqrt(2 * pi / elm_list[m].a) * s4.d0[ms]) * 2.5;

    if (akii == 0)       //place for optimization Sara!9
    {
        akii = 1e-9;
    }
   

    float temp = aki / akii;
    float seta1 = atanf(0.25 * (temp + sqrt(temp * temp + 8)));
    float seta2 = atanf(0.25 * (temp - sqrt(temp * temp + 8)));

    float temp2 = seta1 / 2;  //Sara do this for sta2
    float akie1 = aki * pow(cosf(temp2), 3) - 3 * akii * pow(cosf(temp2),2) * sinf(temp2);
    float akie2 = aki * pow(cosf(seta2 / 2), 3) - 3 * akii * pow(cosf(seta2 / 2),2) *
        sinf(seta2 / 2);

    akie = max(akie1, akie2);
    
    if (akie > factors.factor_e * rock1[mm].akic || akii > factors.factor_e * rock1[mm].akiic)
    {
        kk = 1;             // 1 is checking elastic growth
    }

    return kk;
}



void Rock::save_to_file(ofstream& f)
{
    f << pi << " " << e << " " << pr << " " << akic <<
        " " << akiic << " " << irock << " " << rphi << " "
        << rcoh << " " << rst << std::endl;

    }

void Rock::read_from_file( ifstream& f)
{
    f >> e >> pr >> akic >> akiic >>
        irock >> rphi>> rcoh >> rst;


}
