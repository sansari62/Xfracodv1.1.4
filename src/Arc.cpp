#include<Arc.h>



float Arch::get_arcx() const { return arcx; }
float Arch::get_arcy() const { return arcy; }
float Arch::get_arcr() const { return arcr; }
float Arch::get_arcbeg() const { return arcbeg; }
float Arch::get_arcend() const { return arcend; }
float Arch::get_arcsn() const { return arcsn; }
float Arch::get_arcss() const { return arcss; }
float Arch::get_ds() const { return arcds; }
float Arch::get_dn() const { return arcdn; }
float Arch::get_gsy() const { return arc_gradsy; }
float Arch::get_gny() const { return arc_gradny; }

void Arch::take_arcbeg(float beg) { arcbeg = beg; }
void Arch::take_arcend(float end) { arcend = end; }
void Arch::take_arcx(float x){ arcx = x; }
void Arch::take_arcy(float y) { arcy = y; }
void Arch::take_arcr(float r) { arcr = r; }
void Arch::take_ss(float ss)  { arcss = ss; }
void Arch::take_sn(float sn)  { arcsn = sn; }
void Arch::take_ds(float ds1) { arcds = ds1; }
void Arch::take_dn(float dn1) { arcdn = dn1; }
void Arch::take_gsy(float gs) { arc_gradsy = gs; }
void Arch::take_gny(float gn) { arc_gradny = gn; }