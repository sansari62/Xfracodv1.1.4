#pragma once
#ifndef BE_h
#define BE_h

/* Boundary Element rest attr., length, orientation, joint status, joint properties,
and joint surface stresses and displacement*/


class BE
{
public:

	float sigma_s;
	float sigma_n;
	int ipair;
	int jstate;
	int jslipd;
	int jmode;
	float force1;
	float force2;
	float us;
	float un;
	float forces;
	float forcen;
	float aks;
	float akn;
	float phi;
	float phid;
	float coh;
	float us_neg;
	float un_neg;
	float ss_old;
	float sn_old;

	BE() {
	
		sigma_s = 0;
		sigma_n = 0;
		ipair = 0; jstate = 0; jslipd = 0; us = 0;
		un = 0; forcen = 0;  forces = 0; force1 = 0; force2 = 0; jmode = 0;
		phi = 0;
		phid = 0;
		aks = 0;
		akn = 0;
		coh = 0;
		us_neg = 0;
		un_neg = 0;
		ss_old = 0;
		sn_old = 0;
	
	}
	
	BE( float aks, float akn, float phi,float phid, float coh): phi(phi), phid(phid),
		aks(aks), akn(akn), coh(coh) {
		sigma_s = 0;
		sigma_n = 0;
		ipair = 0; jstate = 0; jslipd = 0; us = 0;
		un = 0; forcen = 0;  forces = 0; force1 = 0; force2 = 0; jmode = 0;
		us_neg = 0;
		un_neg = 0;
		ss_old = 0;
		sn_old = 0;
	}

	
};


#endif	