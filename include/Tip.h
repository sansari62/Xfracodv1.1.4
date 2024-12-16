#ifndef Tip_h
#define Tip_h
#include "Rock.h"

/*tip elements*/

class Tip
{
	// the begin and end points of fractures
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
	void assign_val(float dll, float cos, float sin, int itp, int mat);

	Tip();
	void save_to_file(ofstream& f);


};

void newtips(float dr);
void arrangetip();

void check_crack_growth();
void If_No_tip();

void add_crack_growth();
void input_tip_check();


#endif






