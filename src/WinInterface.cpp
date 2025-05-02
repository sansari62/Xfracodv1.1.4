#include<stdafx.h>

#include <WinInterface.h>
#include<CommonPara.h>

using namespace comvar;

namespace winvar
{
	WindowExchange win_exchange;
	vector<Geom> geom(ng);
	std::vector<Stress> stress(ns, Stress());
	std::vector<Wjoint> wjoint(2*ng);
	std::vector<PermeabilityS> permeability(ns);
	std::vector<AcousricE> AE(ng);
}
