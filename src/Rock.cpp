#include<stdafx.h>

#include "Rock.h"
#include<Mainb.h>
#include<Source.h>
#include<CommonPara.h>
using namespace CommonPara_h::comvar;


Rock::Rock():
	e(0),pr(0),akic(0), akiic(0), rcoh(20e6), rphi(30 * pi / 180.0), rst(10e6) {}




void Rock::save_to_file(ofstream& f)
{
    f << pi << " " << e << " " << pr << " " << akic <<
        " " << akiic << " " << irock << " " << rphi << " "
        << rcoh << " " << rst << std::endl;

    }

    void Rock::read_from_file(ifstream& f)
    {
        f >> pi >> e >>  pr >> akic >>
             akiic >> irock >> rphi >> rcoh>>rst;
    }