#ifndef ExcavationCracks_h
#define ExcavationCracks_h


int CheckNewElement(float ac, float xc, float yc, float cosbeta, float sinbeta, int flagB);

void check_point_in_rock(float xp, float yp,bool flag, int& n_valid);
void check_symm_and_set_n_valid(float xp, float yp, bool flag, int& n_valid);


void stiffness_bb(float& aks_bb, float& akn_bb, float& phi_bb, float& coh_bb, float& phid_bb,
	float& ap_bb, float& apr_bb, int im, int mm);

void ExcavationCracks();


#endif


