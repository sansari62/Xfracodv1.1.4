#pragma once
#ifndef BE_h
#define BE_h
#include <common.h>

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


	BE(){}
	BE(float ss, float sn, int ip, int js, int jslip, int jmod,float fc1, float fc2,float us, float un, float forcs, float forcn, float aks, float akn, float phi,
		float phid, float coh, float usn, float unn, float sso, float sno):sigma_s(ss), sigma_n(sn),
		ipair(ip), jstate(js), jslipd(jslip), us(us), un(un), forcen(forcn), forces(forcs), force1(fc1), force2(fc2), jmode(jmod), phi(phi), phid(phid),
		aks(aks), akn(akn), us_neg(usn), un_neg(unn), ss_old(sso), sn_old(sno), coh(coh) {} 

	
};


#endif	