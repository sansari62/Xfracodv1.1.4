#include "Elliptical_opening.h"

Elliptical_opening::Elliptical_opening():x_cent(0), y_cent(0), diameter1(0), diameter2(0) {}

Elliptical_opening::Elliptical_opening(int mat, int eno, int kode, int x, int y,
	float r1, float r2, float bs, float bn, float grad1, float grad2):
	GeologicalForm(mat, eno,kode), x_cent(x), y_cent(y), diameter1(r1),diameter2(r2)
	//, bvs(bs), bvn(bn), grad_ny(grad1), grad_sy(grad2)
{}

Elliptical_opening::Elliptical_opening(float x, float y, float r1, float r2) :
	x_cent(x), y_cent(y), diameter1(r1), diameter2(r2) {}


void Elliptical_opening::def_arch_boundary(int elem_num, float x, float y, float r1
	, float r2, float ang1, float ang2)
{
}


void Elliptical_opening::increment_bnd_stress_along_arch(float x, float y, float r1,
	float r2, float ang1, float ang2, float dss, float dnn)
{
}


void Elliptical_opening::save_to_file(ofstream& f) 
{
	f << " " << diameter1 << " " << diameter2 << " " << x_cent
		<< " " << y_cent << std::endl;
}

