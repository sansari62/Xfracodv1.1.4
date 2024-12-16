
#pragma once
#ifndef BoundaryStress_h
#define BoundaryStress_h
#include <common.h>



class BoundaryStress
{
public:
	float dsxx;
	float dsxy;
	float dsyy;
	float dss;
	float dnn;
	int incres;

	BoundaryStress();

	void save_to_file(ofstream& f);

};

inline BoundaryStress::BoundaryStress(): dsxx(0.0), dsxy(0.0), dsyy(0.0),
dss(0.0),dnn(0.0),incres(0) {}



inline void BoundaryStress::save_to_file(ofstream& f)
{
	f << incres << " " << dsxx << " " << dsyy << " " << dsxy << " " << dss << " " << dnn <<std::endl;
}
#endif