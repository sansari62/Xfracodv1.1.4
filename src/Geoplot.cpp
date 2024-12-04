#include "Geoplot.h"
#include "CommonPara.h" 
#include <WinInterface.h>
#include <chrono> 
#include <thread>
#include<Internal.h>
#include<fracture_defo.h>
#include<Plot.h>
#include <Acoustic.h>


using namespace WinInterface_h::winvar;
using namespace CommonPara_h::comvar;



void geoplot() 
{
    int numbe_real = 0;             // Exclude the ghost elements

    win_exchange.w_title = title;
    win_exchange.w_xmin = dispwin.xll;
    win_exchange.w_ymin = dispwin.yll;
    win_exchange.w_xmax = dispwin.xur;
    win_exchange.w_ymax = dispwin.yur;
    win_exchange.w_xll = dispwin.xll;
    win_exchange.w_yll = dispwin.yll;
    win_exchange.w_xur = dispwin.xur;
    win_exchange.w_yur = dispwin.yur;
    win_exchange.w_numx = dispwin.numx;
    win_exchange.w_numy = dispwin.numy;
    win_exchange.w_npointg = 1000;
    win_exchange.w_npoints = 1000;
    win_exchange.w_jpoint = numbe * 2;

    win_exchange.w_time = creep.time;                      // Creep parameters
    win_exchange.w_deltaT = creep.deltaT;                 // Creep parameters
    win_exchange.w_vel_creep_max = creep.vel_creep_max;     // Creep parameters
    win_exchange.w_mcyc = mcyc;
    win_exchange.w_mcyc0 = mcyc0;

    win_exchange.w_deltaT_min = creep.deltaT_min;         // Creep parameters
    win_exchange.w_deltaT_max = creep.deltaT_max;         // Creep parameters
    win_exchange.w_totalT = creep.totalT;                // Creep parameters
    win_exchange.w_v1 = creep.v1;                        // Creep parameters
    win_exchange.w_nn1 = creep.nn1;                      // Creep parameters
    win_exchange.w_v2 = creep.v2;                        // Creep parameters
    win_exchange.w_nn2 = creep.nn2;                      // Creep parameters

    int n1 = 0;
    int n2 = n1 + numbe;
    int l = 0;

    for (int i = 0; i < numbe; ++i) 
    {
        BoundaryElement& be = elm_list[i];    //alias for element i        
        if (be.kod != 7)
        {                
            geom[l].w_xbeg = be.xm - be.a * be.cosbet;
            geom[l].w_ybeg = be.ym - be.a * be.sinbet;
            geom[l].w_xend = be.xm + be.a * be.cosbet;
            geom[l].w_yend = be.ym + be.a * be.sinbet;
            geom[l].w_jstat = be.jstate;
            geom[l].w_kod = be.kod > 10 ? be.kod - 10 : be.kod;
            geom[l].w_jwater = watercm.jwater[i];    //Sara! think jwater fpr elem pr struct
            geom[l].w_mat = be.mat_no;
            numbe_real++;
            l = numbe_real;
        }
    }

    if (symm.ksym == 1 || symm.ksym == 4)
    {
        n1 = n2 + 1;
        n2 = n1 + numbe - 1;
        for (int i = n1; i <= n2; ++i) {
            int j = i - n1 + 1;
            BoundaryElement& be2 = elm_list[j];    //alias for element j
            if (elm_list[i].kod != 7)
            {                
                geom[l].w_xbeg = 2.0 * symm.xsym - (be2.xm + be2.a * be2.cosbet);
                geom[l].w_ybeg = be2.ym + be2.a * be2.sinbet;
                geom[l].w_xend = 2.0 * symm.xsym - (be2.xm - be2.a * be2.cosbet);
                geom[l].w_yend = be2.ym - be2.a * be2.sinbet;
                geom[l].w_jstat = be2.jstate;
                geom[l].w_kod = be2.kod > 10 ? be2.kod - 10 : be2.kod;
                geom[l].w_jwater = watercm.jwater[j];
                geom[l].w_mat = be2.mat_no;
                numbe_real++;
                l = numbe_real;
                
            }
        }
    }

    if (symm.ksym == 2 || symm.ksym == 4) 
    {
        n1 = n2 + 1;
        n2 = n1 + numbe - 1;
        for (int i = n1; i <= n2; ++i) {
            int j = i - n1 + 1;
            BoundaryElement& be2 = elm_list[j];    //alias for element j
            if (elm_list[i].kod != 7)
            {               
                geom[l].w_xbeg =be2.xm + be2.a * be2.cosbet;
                geom[l].w_ybeg = 2.0 * symm.ysym - ( be2.ym + be2.a * be2.sinbet);
                geom[l].w_xend = be2.xm - be2.a * be2.cosbet;
                geom[l].w_yend = 2.0 * symm.ysym - (be2.ym - be2.a * be2.sinbet);
                geom[l].w_jstat = be2.jstate;
                geom[l].w_kod = be2.kod > 10 ? be2.kod - 10 : be2.kod;
                geom[l].w_jwater = watercm.jwater[j];  //Sara!
                geom[l].w_mat = be2.mat_no;
                numbe_real++;
                l = numbe_real;
            }
        }
    }

    if (symm.ksym == 3 || symm.ksym == 4) 
    {
        n1 = n2 + 1;
        n2 = n1 + numbe - 1;
        for (int i = n1; i <= n2; ++i)
        {
            int j = i - n1 + 1;
            BoundaryElement& be2 = elm_list[j];    //alias for element j
            if (elm_list[i].kod != 7)
            {                
                geom[l].w_xbeg = 2.0 * symm.xsym - (be2.xm- be2.a * be2.cosbet);
                geom[l].w_ybeg = 2.0 * symm.ysym - (be2.ym - be2.a * be2.sinbet);
                geom[l].w_xend = 2.0 * symm.xsym - (be2.xm + be2.a * be2.cosbet);
                geom[l].w_yend = 2.0 * symm.ysym - (be2.ym + be2.a * be2.sinbet);
                geom[l].w_jstat = be2.jstate;
                geom[l].w_kod = be2.kod > 10 ? be2.kod - 10 : be2.kod;
                geom[l].w_jwater = watercm.jwater[j];  
                geom[l].w_mat = be2.mat_no;
                numbe_real++;
                l = numbe_real;
            }
        }
    }

    // Write data to file    Sara! filename1 for file7 should be taken from the user
   // string filename1 = "file7";
    //ofstream file7("filename1");
    file7 << numbe_real<< endl;

    for (int l = 0; l < numbe_real; ++l)
        file7 << geom[l].w_xbeg << " " << geom[l].w_ybeg << " " << geom[l].w_xend << " " << geom[l].w_xend << " " <<
        geom[l].w_jstat << " " << geom[l].w_kod << " " <<
        geom[l].w_jwater << " " << geom[l].w_mat << endl;

        
    // Open output file here and write data using ofstream or other file
    //  writing mechanisms in C++

    win_exchange.w_numbe = numbe_real;
    int npoint, jpoint, npointp;
    AcousticE();
    internal(0, npoint);   
    fracture_defo(0, jpoint);
    permeability_plot(0, npointp);


    file7 << win_exchange.w_title << endl;
    file7 << win_exchange.w_pxx << "  " << win_exchange.w_pyy << "  " << win_exchange.w_pxy << endl;
    file7 << win_exchange.w_time << "  " << win_exchange.w_deltaT << "  " << win_exchange.w_vel_creep_max << endl;
    file7 << win_exchange.w_mcyc << "  " << win_exchange.w_mcyc0 << endl;
          
    // call SendWindowMessage(ID_PlotGeomNow);  //Sara! plot geom now   
    //MessageBox(nullptr, L, L"Message!", MB_OK);

    //call system_clock(nsec0, n1, n2)
       // 502 
    // Wait for 1 second to allow plot
    // //later uncomment the following
    ////auto nsec0 = std::chrono::system_clock::now().time_since_epoch() / std::chrono::seconds(1);
    ////std::this_thread::sleep_for(std::chrono::seconds(1));
    ////auto nsec1 = std::chrono::system_clock::now().time_since_epoch() / std::chrono::seconds(1);

    ////// Check if 1 second has elapsed
    ////while ((nsec1 - nsec0) < 1 * n1)
    ////    nsec1 = std::chrono::system_clock::now().time_since_epoch() / std::chrono::seconds(1);

        return;
}
