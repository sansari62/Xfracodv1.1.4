#ifndef WinInterface_h
#define WinInterface_h
#include<common.h>



namespace winvar
{
// Define constants
//const int ID_NoMessage = 0;
//const int ID_Pause = 1;
//const int ID_Stop = 2;
//const int ID_Resize = 3;
//const int ID_LoadSavedFile = 4;
//
//const int ID_LoadInputFile = 5;
//const int ID_NoLoadFile = 6;
//const int ID_MustLoadFile = 7;
//const int ID_PLOTG = 8;
//const int ID_PLOTS = 9;
//
//const int ID_ADDFILS = 10;
//const int ID_ADDFARS = 12;
//const int ID_ADDBOUS = 13;
//const int ID_CYCLE = 14;
//const int ID_SaveRun = 15;
//
//const int ID_PlotGeomNow = 31;
//const int ID_PlotStressNow = 32;
//const int ID_PlotShearNow = 33;
//const int ID_PlotDispNow = 34;
//const int ID_SETCREEP = 35;
//const int ID_UpdatePlot = 40;
//const int ID_Status = 50;



// Define structures
extern struct WindowExchange 
{
    std::string w_title;
    float w_xmin, w_ymin, w_xmax, w_ymax;
    float w_xll, w_yll, w_xur, w_yur;
    int w_numx, w_numy;
    float w_dsxx, w_dsyy, w_dsxy, w_dss, w_dnn;
    int w_cycle, w_numbe, w_ksym;
    float w_xsym, w_ysym, w_pxx, w_pyy, w_pxy;
    int w_npoints, w_npointg;
    float w_v1, w_v2, w_totalT, w_deltaT_min, w_deltaT_max;
    int w_nn1, w_nn2;
    float w_time, w_deltaT, w_vel_creep_max;
    int w_mcyc, w_mcyc0;
    int w_jpoint;
    int w_npointp;
}win_exchange;



struct Geom
{
    float w_xbeg, w_ybeg, w_xend, w_yend;
    int w_jstat, w_kod, w_jwater, w_mat;
};

extern vector<Geom> geom;



struct Stress
{
    float w_xp, w_yp, w_sig1, w_sig2, w_bet, w_sig12, w_set, w_disp, w_zet;
    int w_mat;

    Stress(): w_xp(0), w_yp(0), w_sig1(0), w_sig2(0), w_bet(0), w_sig12(0), w_set(0), w_disp(0), w_zet(0), w_mat(1) {}    //default costructoe not getting error in defining vector
    Stress(float x, float y, float s1, float s2, float bet, float s12, float set, float disp, float z, int m) :
        w_xp(x), w_yp(y), w_sig1(s1), w_sig2(s2), w_bet(bet), w_sig12(s12), w_set(set), w_disp(disp), w_zet(z), w_mat(m) {}
};

extern std::vector<Stress> stress;



struct AcousricE
{
    float x, y, m;
};

extern std::vector<AcousricE> AE;

 struct Wjoint
{
    float w_xp, w_yp, w_ds, w_bet, w_dn, w_set, w_aperture, w_zet;


    Wjoint() : w_xp(0), w_yp(0), w_ds(0), w_bet(0), w_dn(0), w_set(0), w_aperture(0), w_zet(0) {}


    void assign_val(float xp, float yp, float ds,float bet, float dn, float wset,
        float waper, float zet)
    {    
		w_xp = xp;
		w_yp = yp;
		w_ds = ds;
		w_bet = bet;
		w_dn = dn;
		w_set = wset;
		w_aperture = waper;
		w_zet = zet;
    }
};

extern std::vector<Wjoint> wjoint;


struct PermeabilityS
{

    float w_xp, w_yp, w_permx, w_permy, w_permr, w_permseta;
};

extern std::vector<PermeabilityS> permeability;


}
#endif