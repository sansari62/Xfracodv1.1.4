#pragma once

#ifndef Arch_h
#define Arch_h
#include<GeologicalForm.h>



class Arch :public GeologicalForm
{
private:
	float arcx;
	float arcy;
	float arcr;
	float arcbeg;
	float arcend;
	float arcsn;
	float arcss;
	float arcdn;
	float arcds;
	float arc_gradsy;
	float arc_gradny;
	

public:
	Arch() {}
	Arch(float x, float y, float r, float beg, float end, float ss, float sn, float ds,
		float dn, float gsy, float gny, int nume, int mat, int kod) :	
		arcx(x), arcy(y), arcr(r), arcbeg(beg), arcend(end), arcsn(sn), arcss(ss),
			arcdn(dn), arcds(ds), arc_gradsy(gsy), arc_gradny(gny), GeologicalForm(mat, nume, kod) {}

	float get_arcx() const;
	float get_arcy() const;
	float get_arcr() const;
	float get_arcbeg() const;
	float get_arcend() const;
	float get_arcsn() const;
	float get_arcss() const;
	float get_ds() const;
	float get_dn() const;
	float get_gsy() const;
	float get_gny() const;

	void take_arcbeg(float beg);
	void take_arcend(float end);
	void take_arcx(float x);
	void take_arcy(float y);
	void take_arcr(float r);
	void take_ss(float ss);
	void take_sn(float sn);
	void take_ds(float ds1);
	void take_dn(float dn1);
	void take_gsy(float gs);
	void take_gny(float gn);

};


#endif