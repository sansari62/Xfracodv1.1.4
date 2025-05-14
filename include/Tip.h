#ifndef Tip_h
#define Tip_h

/*tip elements*/

class Tip
{
public:
	float  xbe;
	float  ybe;
	float  xen;
	float  yen;
	float  angl;
	
	float  f_value;
	float  dt;
	float  dl;
	float  sintem;
	float  costem;
		
	int		mat_no;
	int		imode;
	int		ityp;
	int		ifail;
	int		mpointer;
	int     kindtip;
	int		nu;
	

	void assign_val(float x1, float y1, float x2, float y2, float dll, float cos1,
		float sin1, int itp, int mat);

	Tip();
	void save_to_file(ofstream& f);
	void read_from_file(ifstream& f);


};

void newtips(float dr);
void arrangetip();

void check_crack_growth();
void If_No_tip();

void add_crack_growth();
void input_tip_check();
void stiffness_bb(float& aks_bb, float& akn_bb, float& phi_bb, float& coh_bb, float& phid_bb,
	float& ap_bb, float& apr_bb, int im, int mm);
void specialLabel_200(int& mergtip, int i);

#endif






