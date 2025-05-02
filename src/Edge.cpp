#include<stdafx.h>

#include "Edge.h"

Edge::Edge(): GeologicalForm(), x_beg(0), y_beg(0), x_end(0), y_end(0), bvs(0), bvn(0), ds(0), dn(0), grad_ny(0), grad_sy(0) {}

Edge::Edge(int mat, int eno, int kode, float xr, float yr, float xl, float yl,float bs,float bn,
	float ds1,float dn1, float gs,float gn):
	GeologicalForm(mat,eno,kode),x_beg(xr), y_beg(yr), x_end(xl), y_end(yl), bvs(bs), bvn(bn),ds(ds1),dn(dn1), grad_ny(gn), grad_sy(gs){}

float Edge::get_xbeg() const { return x_beg; }
float Edge::get_ybeg() const { return y_beg; }
float Edge::get_xend() const {	return x_end;}
float Edge::get_yend() const { return y_end; }
float Edge::get_bvn() const {return bvn; }
float Edge::get_bvs() const { return bvn; }
float Edge::get_ds() const { return ds; }
float Edge::get_dn() const { return dn; }
float Edge::get_gsy() const { return grad_sy; }
float Edge::get_gny() const { return grad_ny; }



void Edge::take_xbeg(float x) { x_beg = x; }
void Edge::take_ybeg(float y) { y_beg = y; }
void Edge::take_xend(float x) { x_end = x; }
void Edge::take_yend(float y) { y_end = y; }
void Edge::take_bvs(float ss) { bvs = ss; }
void Edge::take_bvn(float sn) { bvn = sn; }
void Edge::take_ds(float ds1) { ds = ds1; }
void Edge::take_dn(float dn1) { dn = dn1; }
void Edge::take_gsy(float gs) { grad_sy = gs; }
void Edge::take_gny(float gn) { grad_ny = gn; }



