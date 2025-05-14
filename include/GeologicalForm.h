#ifndef GeologicalForm_h
#define GeologicalForm_h
#include"BoundaryElement.h"


class GeologicalForm
{
	
public:
	int mat_no;
	int elem_no;
	int bound_type;

	GeologicalForm();
	GeologicalForm(int mat, int eleno, int kode);

	int getMatno();
	int getElem_no();
	int getBoundtype();
	

void cross_current_or_save_element(float xb1, float yb1, float xe1, float ye1, float xcross, float ycross, int m,int ii, int &);

void ctl_cross_elements(int& k,int m, int numbe0,int num);
void chk_potential_crack_growth(float xbeg, float ybeg, float xend, float yend, int numbe0, int itype);
	
//create set of boundary elements for each geologicalforms
int def_boundary_elements_for_Geoform( int num, float xbeg,float ybeg, float xend, float yend,
		float bvs, float bvn, float gradsy,float gradny,int itype,int jmat);  
	
		
};
bool isNewElementUnique(const BoundaryElement& newelem);
#endif