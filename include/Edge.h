#pragma once
#ifndef Edge_h
#define Edge_h
#include<GeologicalForm.h>


//define the straight boundary line for the model
class Edge:public GeologicalForm
{
private:
	float x_beg;
	float y_beg;
	float x_end;
	float y_end;
	float bvn;
	float bvs;
	float ds;
	float dn;
	float grad_sy;
	float grad_ny;
public:
	Edge();
	Edge(int mat, int eno,int kode,float xr,float yr,float xl,float yl, float bs,
		float bn, float ds, float dn, float gs, float gn);
	float get_xbeg() const;
	float get_ybeg() const;
	float get_xend() const;
	float get_yend() const;
	float get_bvn() const;
	float get_bvs() const;
	float get_ds() const;
	float get_dn() const;
	float get_gsy() const;
	float get_gny() const;

	void take_xbeg(float x) ;
	void take_ybeg(float y) ;
	void take_xend(float x) ;
	void take_yend(float y) ;
	void take_bvs(float ss) ;
	void take_bvn(float sn) ;
	void take_ds(float ds1) ;
	void take_dn(float dn1) ;
	void take_gsy(float gs) ;
	void take_gny(float gn) ;


};
#endif

