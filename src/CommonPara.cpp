#include<stdafx.h>
#include "CommonPara.h"
//#include<map>



namespace comvar {

	int no = 0;
	int ni = 0;
	int numbe = 0;
	int numbe_old = 0;
	int n_it = 1;

	int na = 0;
	int ntunnel = 0;
	int nellipse = 0;
	int nb = 0;
	int nf = 0;

	int ng = m0;    //the real value should taken from UI, ng=m0 and 
	int ns = 10000+m0;   //ns no of grid-points + m0 (100*100+3000)

	int npli = 0;
	int nq = 0;

	float k_num = 0.0;
	float d_max = 1000;

	int nc = 0;
	string title = "FRACOD application";

	float delta = 0.0;
	float w0 = 0.0;
	float w1 = 0.0;

	int mf = 0;
	std::wstring selectedFile;

	const int m1 = 1500;
	const int m2 = 3000;
	const int MAX_SIZE_frac = 2000;
	const int MAX_SIZE_arc = 50;
	float pi = 3.1415926; 
	const int mem_size = 2*m0;
	const float ge = 9.81;

	bool multi_region = false;
	bool water_mod = false;
	int j_material = 1;

	int nelement = 0;
	int mcyc = 0;
	int mcyc0 = 0;
	bool restor_flg = false;

	int ihist = 0;
	int lhist = 0;

	int ID_dtip = 0; 
	bool ktipgrow = false ; 
	bool StopReturn = false; 
	std::string lastinput = "   ";
	int line = 0;	
	int prenumbe = 0;

	int irock = 0 ;		
	int mat_lining = 10;

	float dist_thr = 1e-6;
	float ang_thr = 1;
	float stres_thr = 1;
	int state = 0;
	

	std::vector<std::pair<float, float>> valid;
	std::vector<std::pair<float, float>> valid1;	

	std::vector<BoundaryElement> elm_list(m0);	
	std::vector<BE> b_elm(m0);
	std::vector<MonitoringPoint> mpoint_list(20);
	std::vector<MonitoringLine> mline_list(10);

	std::vector < Elliptical_opening>ellip_list(10);	
	//std::vector<Borehole> tunnl_list = {};


	std::vector <Arch> arc_list(MAX_SIZE_arc);
	std::vector<Edge> bund_list(MAX_SIZE_arc);

	std::vector<Edge_interface> lin_intrfce_list(MAX_SIZE_arc);
	std::vector<Arch_interface> arc_intrface_list(MAX_SIZE_arc);

	std::vector<Fracture> frac_list(MAX_SIZE_frac);

	Borehole tunnl;
	//std::vector<Creep> creep(500);

	BoundaryStress insituS ;
	symmetry symm;	
	s2 s2us;

	gravity g;
	S4 s4(mem_size);  
	s5 s5u;
	S15 s15;
	matrix matx;

	std::vector<S8> s8(21);		
	Creep creep;

	waterCommon watercm;
	std::vector< Rock>  rock1(10, Rock());
	std::vector<Tip> tips(m0);// 500 in v2.3, Sara think about this
	DispWindow dispwin;
	bool rect_exca = false;


	std::vector<Joint> joint(m0);       //join includes aperture0 and aperture_r
	Excavation exca;
	SET factors;

	std::vector<Limited_initiation_points>  init_point(5000);    
	
	Permeability perm;	
	std::ofstream file50("Ccreep_results.dat");
	wstring filepath = L"";

	wstring dir = L"";
	wstring stress_dir;
	wstring BE_dir;
	wstring fd_dir;
	wstring monit_dir;
	std::array<std::ofstream, 20> mon_files;
	std::array<std::ofstream, 10> ml_files;
	std::ofstream file2; 
	std::ofstream file57; 
	std::ifstream inFile;

	std::fstream file9;  // ("Cbound.dat");
	std::ofstream logfile;	
	
}
 














































                             