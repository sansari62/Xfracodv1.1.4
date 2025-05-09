#ifndef CommonPara_h
#define CommonPara_h

#include"BoundaryElement.h"
#include<MonitoringPoint.h>
#include<MonitoringLine.h>
#include<BoundaryStress.h>
#include<Elliptical_opening.h>

#include<Tip.h>
#include<Borehole.h>
#include<Rock.h>
#include<DispWindow.h>

#include "S4.h"
#include <array>

#include<Creep.h>
#include<Arc.h>
#include <Edge.h>
#include <Fracture.h>
#include <Edge_interface.h>
#include <Arch_interface.h>
#include<map>
#include<BE.h>

#define m0 500 //2000





namespace comvar {

	/* the first 7 from following shows S9: parameters to caculate the strain energy*/
	extern int no;				//number of fracture tips much less than “numbe”
	extern int ni;              // one specific tip out of the total number
	extern int numbe;			//total noumber of boundary elements
	extern int numbe_old;		 //old no of boundary elements
	extern float delta;
	extern float w0;
	extern float w1;

	extern int na;				//no of arcs
	extern int ntunnel;			//no of tunnel
	extern int nellipse;        //no of ellipse

	extern int nb;				//no of streight boundary
	extern int nf;			// no of fracture
	extern int npli;			//no of linear interface
	extern int nq;			// no of arc interface

	extern int ng, ns;

	extern int mf;			      // numbe of identified failure points

	extern int n_it;             //no of iteration

	extern string title;
	extern int nc;
	extern int state;
	//others
	extern int nelement;

	//const para
	extern const int m1;
	extern const int m2;
	extern const int MAX_SIZE_frac;
	extern const int MAX_SIZE_arc;
	extern const float ge;
	extern bool multi_region;
	extern bool water_mod;
	extern int j_material;
	extern float pi;

	extern const float zerof;
	extern float dist_thr;

	extern std::wstring selectedFile;
	extern wstring dir;


	extern float k_num;            //stiffness for boundary element;    struct in origincode
	extern float  d_max;		     //d_max - max joint disp in each step

	extern int mcyc;            //S13 defined here ,// cycle number
	extern int mcyc0;
	

	extern int lhist;		// no of monitoring lines
	extern int ihist;		// no of monitoring points
	

	extern std::vector<std::pair<float, float>> valid;
	extern std::vector<std::pair<float, float>> valid1;


	extern bool ktipgrow;		
	extern int ID_dtip;
	extern bool StopReturn;
	extern std::string lastinput;
	extern int line;
	extern  int prenumbe;
	extern int irock;

	extern std::vector<Tip> tips;			//Array of tips 
	extern int mat_lining;

	extern std::vector<BoundaryElement> elm_list;		// list of boundary elements
	extern std::vector<BE> b_elm;								//list of boundaryelement non freq attr.

	extern std::vector<MonitoringPoint> mpoint_list;	//list of monitoring points
	extern std::vector<MonitoringLine> mline_list;		//list of monitoring lines

	extern std::vector < Elliptical_opening>ellip_list;	//list of elliptical opening 
	extern std::vector <Arch> arc_list;
	extern std::vector<Edge> bund_list;

	extern std::vector<Fracture> frac_list;
	extern std::vector<Edge_interface> lin_intrfce_list;
	extern std::vector<Arch_interface> arc_intrface_list;
	
	
	extern Borehole tunnl;                              //S11
	extern BoundaryStress insituS;					   //define boundary stress for the model ,S12

	extern Creep creep;
	extern std::vector < Rock> rock1;
	;


	extern DispWindow dispwin;         //S10
	extern S4 s4;
	extern string filename1;
	extern ofstream file7;
	extern std::ofstream file2;
	extern std::ofstream file50;
	extern std::ifstream inFile;
	extern std::fstream file9;
	extern std::ofstream logfile;

	extern wstring filepath;
	extern string filepath1;

	extern wstring stress_dir;
	extern wstring BE_dir;
	extern wstring fd_dir;
	extern wstring monit_dir;
	extern std::array<std::ofstream, 20> mon_files;
	extern std::array<std::ofstream, 10> ml_files;



	//symmetry is S3 in the original code
	extern struct symmetry
	{
		int ksym;
		float ysym;
		float xsym;
		float pxx1;
		float pxy1;
		float pyy1;

		symmetry() : ksym(0), ysym(0), xsym(0), pxx1(0),
			pxy1(0), pyy1(0) {};

		void read_from_file(ifstream& f)
		{
			f >> ksym >> ysym >> xsym >> pxx1 >> pxy1 >> pyy1;
		}

		void save_to_file(ofstream& f)
		{
			f << ksym << " " << ysym << " " << xsym << " " << pxx1 << " " << pxy1 << " " 
				<< pyy1 << std::endl;
		}

	}symm;

	//S4 is defined as class S4 


	// parameters to calculate the influence coeficcients
	extern struct s2
	{
		float sxxs;
		float sxxn;
		float syys;
		float syyn;
		float sxys;
		float sxyn;
		float uxs;
		float uxn;
		float uys;
		float uyn;

		s2() : sxxs(0.0), sxxn(0.0), syys(0.0), syyn(0.0),
			sxys(0.0), uxs(0.0), uxn(0.0), uys(0.0), uyn(0.0), sxyn(0.0) {}


		// reset all members here
		//this method considered instead of initl
		void reset()
		{
			sxxs = 0.0;
			sxxn = 0.0;
			syys = 0.0;
			syyn = 0.0;
			sxys = 0.0;			
			sxyn = 0.0;

			uxs = 0.0;
			uxn = 0.0;
			uys = 0.0;
			uyn = 0.0;

		}
		void read_from_file(ifstream & f)
		{
			f >> sxxs >> sxxn >> syys >> syyn >> sxys >> 
				sxyn >> uxs >> uxn >> uys >> uyn;		
		}
		void save_to_file(ofstream& f)
		{
			f << sxxs << " " << sxxn << " " << syys << " " << syyn << 
				" " << sxys << " " << sxyn << " " << uxs << " " << uxn << " "
				<< uys << " " << uyn << std::endl;
		}

	} s2us;			 // I''m not sure yet what this is


	/* Model display window size */
	extern struct s5
	{
		float xmax;
		float ymax;
		float xmin;
		float ymin;
		float dtol;

		s5():xmin(-100000.), xmax(100000.), ymin(-100000.), ymax(100000.), dtol(0) {}

		void read_from_file(ifstream& f)
		{
			f >> xmax >> ymax >> xmin >> ymin >> dtol;
		}

		void save_to_file(ofstream& f)
		{
			f << xmax << " " << ymax << " " << xmin << " " << 
				ymin << " " << dtol << std::endl;
		}
	}s5u;

	//S6 implemented as class BoundaryElement
	//S7 implemented as class Tip

	//fracture mechanical properties
	struct S8
	{
		float aks0;
		float akn0;
		float phi0;
		float coh0;
		float phid0;
		float apert0;
		float apert_r;

		S8(): aks0(0.0),akn0(0.0), phi0(0.0), coh0(0.0), phid0(0.0), apert0(0.0), apert_r(0.0) {}
		bool equl_zero() 
		{
			if (aks0 == 0 && akn0 == 0 && phi0 == 0 && phid0 == 0 && coh0 == 0)
				return true;
			else
				return false;	
		
		}

		void read_from_file(ifstream& f)
		{
			f >> aks0 >> akn0 >> phi0 >> phid0 >> coh0 >> phid0 >> apert0 >> apert_r ;;
		}

		void save_to_file(ofstream& f)
		{
			f << aks0 << " " << akn0 << " " << phi0 << " " << phid0 << " " << coh0 << " " 
				<< phid0 << " " << apert0 << " " << apert_r << std::endl;
		}
	};
	extern std::vector<S8> s8;


	//S9 defined above
	//S10 define as class DispWindow
	//S11 implemented as class Borhole and elliptical_opening
	//S12 implemented as BoundaryStress class and the object instiu defined from this class
	//S13 defined above as two common variables
	//S14 , all 5 variables defined above
	//monitoring points and monitoring lines are defined as classes


	//control parameters
	extern struct S15
	{
		float aa;			//maximum element SIZE
		float aaa;			//minimum element SIZE
		int i_rand;			//!i_rand=1, random fracture initiation; =0 not random
		float f_ini0;
		int  i_bound;		// 1: allow fracture initiation at boundarieds 
		int  i_intern;		//1 allow Fracture initiation in intact rock 
		float a_ini;		//Fracture initiation element size
		int l_rand;

		S15() :aa(0), aaa(10000), i_rand(0), f_ini0(0), l_rand(1), i_bound(0), 
			i_intern(0), a_ini(0){}

		void read_from_file(ifstream& f)
		{
		
			f>> aa >> aaa >> i_rand >> l_rand >> f_ini0 >> i_bound >> i_intern;		
		}

		void save_to_file(ofstream& f)
		{
			f << aa << " " << aaa << " " << i_rand << " " << l_rand << " " << f_ini0 <<
				" " <<	i_bound << " " << i_intern << std::endl;
		}
	}s15;
		


	extern struct matrix
	{
		float A_in_js;
		float A_in_jn;
		float A_is_js;
		float A_is_jn;
		float B_in_js;
		float B_in_jn;
		float B_is_js;
		float B_is_jn;

		matrix() :A_in_js(0.0), A_in_jn(0.0), A_is_js(0.0), A_is_jn(0.0),
			B_in_js(0.0), B_in_jn(0.0), B_is_js(0.0), B_is_jn(0.0) {}

	} matx;




	extern struct gravity
	{
		float skx;
		float sky;
		float y_surf;

		gravity():skx(0.0), sky(0.0), y_surf(0) {}
	} g;	 




	//water pressure boundary conditions
	//ID_range = 1 - hole, =2 - rectangular
	extern struct waterCommon
	{
		int jwater[m0];
		float pwater[m0];
		int ID_range;
		float w_xc[10];
		float w_yc[10];
		float w_d[10];
		int iwhole;					//no of circular holes
		int iwrect;					//no of rectangular holes
		float w_x1[10];
		float w_x2[10];
		float w_y1[10];
		float w_y2[10];
		float wph[10];
		float wpr[10];



		waterCommon() : jwater{0}, pwater{0}, w_d{0}, w_x1{0}, w_x2{0}, w_y1{0}, w_y2{0}, wph{0}, wpr{0}, w_xc{0},
			w_yc{0}, iwhole(0), iwrect(0),ID_range(0) {}

		void read_from_file1(ifstream& f)
		{
			for (int m = 0; m < numbe; ++m)
				f>> jwater[m] >> pwater[m];

			for (int m = 0; m < 10; ++m)
			{
				f >> ID_range >> w_xc[m] >> w_yc[m] >> w_d[m] >> iwhole >>
					iwrect >> w_x1[m] >> w_x2[m] >> w_y1[m] >> w_y2[m] >> wph[m] >> wpr[m];
			}
			return;
		}

		//optimize point for iwhole should once written and read from file Sara!
		
		void save_to_file(ofstream& f)
		{
			// Writing jwater and pwater arrays
			for (int m = 0; m < numbe; ++m)
			{
				f << jwater[m] << " " << pwater[m] << std::endl;
			}

			// Writing data for 10 elements in arrays w_xc, w_yc, w_d, etc.
			for (int m = 0; m < 10; ++m)
			{
				f << ID_range << " " << w_xc[m] << " " << w_yc[m] << " " << w_d[m] <<
					" " << iwhole << " " << iwrect << " " <<
					w_x1[m] << " " << w_x2[m] << " " << w_y1[m] << " " << w_y2[m] << 
					" " << wph[m] << " " << wpr[m] << std::endl;
			}
			return;
		}
	}watercm;



	/* Fracture initiation */
	struct  Limited_initiation_points
	{		
		float xf;
		float yf;
		float rf;
		float alphaf;
		int imf;         
		float fosf;
		int Locationf;

		Limited_initiation_points(): xf(0), yf(0), rf(0), alphaf(0), 
			imf(0), fosf(0), Locationf(0) {}

		void read_from_file(ifstream& f)
		{			
			f >> xf >> yf >> rf >> alphaf >> imf >> fosf >> Locationf;
			
		}

		void save_to_file(ofstream& f)
		{			
			f << xf << " " << yf << " " << rf << " " << alphaf << " " <<
					imf << " " << fosf << " " << Locationf << std::endl;
			
			return;

		}
	};
	extern std::vector<Limited_initiation_points>  init_point;

	//all fields for struct num_stablility defined above

	/*  joint apertures for permeability calculations */
	struct Joint
	{
		float  aperture0;
		float aperture_r;

		Joint():aperture0(10e-6) , aperture_r(10e-6){}

	};
	extern std::vector<Joint> joint;



	extern struct Permeability
	{
		float viscosity;
		float density;
		double perm0;

		Permeability():viscosity(1e-3), density(1000), perm0(1e-9) {}

	} perm; 



	/* Factor_f */
	extern struct SET
	{
		float factor_f;
		float factor_e;
		float tolerance;

		SET() : factor_e(0.5), factor_f(0.9), tolerance(1.0) {}

	}factors;




	//to consider the excavation induced cracks
	extern struct Excavation
	{
		int ID_Exca;
		float d_wall;
		float rand_e;			//"rand_e" percent of internal points that have excavation induced cracks
		Excavation() : ID_Exca(0), d_wall(0.0), rand_e(1.0) {}

	}exca;

}
#endif
