#include<stdafx.h>
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




void export_to_vtk() {
    wstring filename = BE_dir + L"/BE" + std::to_wstring(mcyc) + L".vtk";
    std::ofstream BEfile(filename);
    std::ofstream vtkfile(filename);

    if (!vtkfile.is_open()) {
        std::cerr << "Failed to open vtk file " << std::endl;
        return;
    }

    int num_lines = win_exchange.w_numbe;
    int num_points = num_lines * 2;

    vtkfile << "# vtk DataFile Version 3.0\n";
    vtkfile << "Boundary Elements with jstat coloring\n";
    vtkfile << "ASCII\n";
    vtkfile << "DATASET POLYDATA\n";
    vtkfile << "POINTS " << num_points << " float\n";

    for (int l = 0; l < num_lines; ++l)
    {
        vtkfile << geom[l].w_xbeg << " " << geom[l].w_ybeg << " 0\n";
        vtkfile << geom[l].w_xend << " " << geom[l].w_yend << " 0\n";
    }

    vtkfile << "\nLINES " << num_lines << " " << num_lines * 3 << "\n";

    // Write lines (2 point indices per line)
    for (int i = 0; i < num_lines; ++i) {
        vtkfile << "2 " << i * 2 << " " << i * 2 + 1 << "\n";
    }


    // CELL_DATA section (per line)
    vtkfile << "\nCELL_DATA " << num_lines << "\n";
    vtkfile << "SCALARS jstat int 1\n";
    vtkfile << "LOOKUP_TABLE default\n";
    for (int l = 0; l < num_lines; ++l) {
        vtkfile << geom[l].w_jstat << "\n";
    }
    vtkfile.close();
    std::cout << "VTK file written: " << std::endl;
}






void geoplot() 
{
    int numbe_real = 0;  
    wstring filename = BE_dir + L"/BE" + std::to_wstring(mcyc) + L".csv";
    std::ofstream BEfile(filename);
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
            geom[l].w_jstat = b_elm[i].jstate;
            geom[l].w_kod = be.kod > 10 ? be.kod - 10 : be.kod;
            geom[l].w_jwater = watercm.jwater[i];
            geom[l].w_mat = be.mat_no;           
            l++;
        }
    }
    numbe_real = l;
    if (symm.ksym == 3 || symm.ksym == 4)
    {
        n1 = n2;
        n2 = n1 + numbe;
        for (int i = n1; i < n2; ++i)
        {
            int j = i - n1;
            BoundaryElement& be2 = elm_list[j];    //alias for element j
            if (elm_list[j].kod != 7)
            {
                geom[l].w_xbeg = 2.0 * symm.xsym - (be2.xm - be2.a * be2.cosbet);
                geom[l].w_ybeg = 2.0 * symm.ysym - (be2.ym - be2.a * be2.sinbet);
                geom[l].w_xend = 2.0 * symm.xsym - (be2.xm + be2.a * be2.cosbet);
                geom[l].w_yend = 2.0 * symm.ysym - (be2.ym + be2.a * be2.sinbet);
                geom[l].w_jstat = b_elm[j].jstate;
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
            n1 = n2 ;
            n2 = n1 + numbe ;
            for (int i = n1; i < n2; ++i) {
                int j = i - n1;
                BoundaryElement& be2 = elm_list[j];    //alias for element j
                if (elm_list[j].kod != 7)
                {
                    geom[l].w_xbeg = be2.xm + be2.a * be2.cosbet;
                    geom[l].w_ybeg = 2.0 * symm.ysym - (be2.ym + be2.a * be2.sinbet);
                    geom[l].w_xend = be2.xm - be2.a * be2.cosbet;
                    geom[l].w_yend = 2.0 * symm.ysym - (be2.ym - be2.a * be2.sinbet);
                    geom[l].w_jstat = b_elm[j].jstate;
                    geom[l].w_kod = be2.kod > 10 ? be2.kod - 10 : be2.kod;
                    geom[l].w_jwater = watercm.jwater[j];  //Sara!
                    geom[l].w_mat = be2.mat_no;
                    numbe_real++;
                    l = numbe_real;
                }
            }
    }

    if (symm.ksym == 1 || symm.ksym == 4)
    {
        n1 = n2 ;
        n2 = n1 + numbe;
        for (int i = n1; i < n2; ++i) {
            int j = i - n1 ;
            BoundaryElement& be2 = elm_list[j];    //alias for element j
            if (elm_list[j].kod != 7)
            {                
                geom[l].w_xbeg = 2.0 * symm.xsym - (be2.xm + be2.a * be2.cosbet);
                geom[l].w_ybeg = be2.ym + be2.a * be2.sinbet;
                geom[l].w_xend = 2.0 * symm.xsym - (be2.xm - be2.a * be2.cosbet);
                geom[l].w_yend = be2.ym - be2.a * be2.sinbet;
                geom[l].w_jstat = b_elm[j].jstate;
                geom[l].w_kod = be2.kod > 10 ? be2.kod - 10 : be2.kod;
                geom[l].w_jwater = watercm.jwater[j];
                geom[l].w_mat = be2.mat_no;
                numbe_real++;
                l = numbe_real;                
            }
        }
    }  

    /*std::stringstream buffer;
    buffer << std::fixed << std::setprecision(4);

    buffer << "x1"<<","<<"y1" << "," << "x2" << "," << "y2" << "," << "jstat" << "," << "kod" << "," << "jwater" << "," << "mat" << endl;
    for (int l = 0; l < numbe_real; ++l)
        buffer << geom[l].w_xbeg << "," << geom[l].w_ybeg << "," << geom[l].w_xend << "," << geom[l].w_yend << "," <<
        geom[l].w_jstat << "," << geom[l].w_kod <<"," <<
        geom[l].w_jwater << "," << geom[l].w_mat << endl;

    BEfile << buffer.str();
    BEfile.close();*/

    // Helper to round to N decimal places
    auto roundTo = [](float value, int decimals) {
        float factor = std::pow(10.0, decimals);
        return std::round(value * factor) / factor;
        };

    std::stringstream buffer;
    buffer << std::fixed << std::setprecision(3);

    buffer << "x1" << "," << "y1" << "," << "x2" << "," << "y2" << ","
        << "jstat" << "," << "kod" << "," << "jwater" << "," << "mat" << "\n";

    for (int l = 0; l < numbe_real; ++l) {
        float x1 = roundTo(geom[l].w_xbeg, 3);
        float y1 = roundTo(geom[l].w_ybeg, 3);
        float x2 = roundTo(geom[l].w_xend, 3);
        float y2 = roundTo(geom[l].w_yend, 3);

        buffer << x1 << "," << y1 << "," << x2 << "," << y2 << ","
            << geom[l].w_jstat << "," << geom[l].w_kod << ","
            << geom[l].w_jwater << "," << geom[l].w_mat << "\n";
    }

    BEfile << buffer.str();
    BEfile.close();

    win_exchange.w_numbe = numbe_real;
    //export_to_vtk();
    int npoint, jpoint, npointp;
    AcousticE();
    internal(npoint);   
    fracture_defo(jpoint);          
    ++state;
    return;
}



