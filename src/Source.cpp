#include<stdafx.h>
#include "Source.h"
#include<ExcavationCracks.h>
#include "CommonPara.h"
#include"Work.h"
#include<Input.h>
#include<Geoplot.h>
#include<chrono>

using namespace CommonPara_h::comvar;




void finding_in_rock_iniR()
{
    int num1, num2;
    float dx, dy, ddr, dda;
    float xp = 0, yp = 0;
    int  n_valid = 0;

    if (dispwin.ID_win == 1)
    {
        dx = (dispwin.xur - dispwin.xll) / dispwin.numx;
        dy = (dispwin.yur - dispwin.yll) / dispwin.numy;
        num1 = dispwin.numx;
        num2 = dispwin.numy;
    }
    else //if (dispwin.ID_win == 2)
    {
        ddr = dispwin.radium / dispwin.numr;
        dda = 2 * pi / dispwin.numa;
        num1 = dispwin.numa;
        num2 = dispwin.numr;
    }    
   
    for (int i1 = 0; i1 <= num1; ++i1)
    {
        for (int i2 = 0; i2 <= num2; ++i2)
        {
            if (dispwin.ID_win == 1)
            {
                xp = dispwin.xll + i1 * dx;
                yp = dispwin.yll + i2 * dy;
            }
            else //if (dispwin.ID_win == 2) 
            {
                xp = dispwin.xc0 + i1 * ddr * cosf(i2 * dda);
                yp = dispwin.yc0 + i1 * ddr * sinf(i2 * dda);
            }

            if (xp > s5u.xmax || xp < s5u.xmin || yp > s5u.ymax || yp < s5u.ymin)
            {
                continue;
            }

            if ((symm.ksym == 1 || symm.ksym == 4) && xp < (symm.xsym + s5u.dtol)) {

                continue;
            }
            if ((symm.ksym == 2 || symm.ksym == 4) && yp < (symm.ysym + s5u.dtol)) {

                continue;
            }
            if (symm.ksym == 3 && yp < 0 && xp < (symm.xsym - s5u.dtol)) {
                continue;
            }
            if (symm.ksym == 3 && yp > 0 && xp < (symm.xsym + s5u.dtol)) {
                continue;
            }

            if ((symm.ksym == 1 || symm.ksym == 4) && xp == symm.xsym) {
                continue;
            }
            if ((symm.ksym == 2 || symm.ksym == 4) && yp == symm.ysym) {
                continue;
            }

            n_valid = 1;
            check_point_in_rock(xp, yp, 1, n_valid);
            if (n_valid == 1)
            {
                comvar::valid1.push_back(std::make_pair(xp, yp));
                continue;
            }

        }
    }  
}





void compute_n_vlaid_all_points() {

    int num1 = 0, num2 = 0;
    float dx = 0, dy = 0, ddr = 0, dda = 0 ,xp = 0, yp = 0;

    if (dispwin.ID_win == 1)
    {
        dx = (dispwin.xur - dispwin.xll) / dispwin.numx;
        dy = (dispwin.yur - dispwin.yll) / dispwin.numy;
        num1 = dispwin.numx;
        num2 = dispwin.numy;
    }
    else
    {
        ddr = dispwin.radium / dispwin.numr;
        dda = 2 * pi / dispwin.numa;
        num1 = dispwin.numr;
        num2 = dispwin.numa;
    }

    for (int ix = 0; ix <= num1; ++ix)
    {
        for (int iy = 0; iy <= num2; ++iy)
        {
            if (dispwin.ID_win == 1)
            {
                xp = dispwin.xll + ix * dx;
                yp = dispwin.yll + iy * dy;
            }

            else
            {
                xp = dispwin.xc0 + ix * ddr * cosf(iy * dda);
                yp = dispwin.yc0 + ix * ddr * sinf(iy * dda);
            }


            int n_valid = 1;
            check_point_in_rock(xp, yp, 0, n_valid);           

            // Add all points in rock
            if (n_valid == 1)
            {               
                comvar::valid.push_back(std::make_pair(xp, yp));
                continue;
            }           
        }
    }   
    return;
}


std::string getTimeOnly() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::tm tm;
    localtime_s(&tm, &now_c);

    char buffer[9]; // HH:MM:SS + null terminator
    std::strftime(buffer, sizeof(buffer), "%H:%M:%S", &tm);

    return std::string(buffer);
}



void main_computing_part()
{
    if (insituS.incres != 0 && ktipgrow == true)     //special treatment th final grownth element in cycle
    {
        int incres_tem = insituS.incres;
        insituS.incres = 0;
        work0(-1);              //mode = -1 is back step
        insituS.incres = incres_tem;
    }
    if (insituS.incres != 0)
        n_it = 20;    //iteration process
    work0(0);         //First calculation in every cycle - mode = 0.                      
    check_crack_growth();              //include instant and creep
    geoplot();
    if (StopReturn == true) return;
    prenumbe = numbe;
    if(mcyc != mcyc0)
        add_crack_growth();   
}




void creep_problem()
{
    do     //Creep iteration
    {
        if (creep.ID_creep == 1 && creep.ID_fast_crack == 0)
        {
            if ((creep.time + creep.deltaT) > creep.totalT)
            {
                creep.deltaT = creep.totalT - creep.time;
                creep.time = creep.totalT;
            }
            else
                creep.time += creep.deltaT;
        }
        main_computing_part();
    } while (creep.time != creep.totalT && creep.ID_creep != 0 && creep.ID_fast_crack != 1);
}






void Central_control()
{    
    auto start = std::chrono::high_resolution_clock::now();
    if (!file50 || !file2 )
    {
        cout << "Error in opening output files2 or file50!";
        return;
    }  
    mcyc = 0;
    creep.ID_creep = 0;
    n_it = 20;
    StopReturn = false; 

    input();
    CheckRange();
    logfile << "The initial number of boundary elements:"<< numbe <<", fractures:" << nf <<
        ", Archs: " <<  na << "edges: "<<nb<<endl;
    if (StopReturn == true) return;
    file2 << " Tip No.    xc       yc         Angle            W0             Wi              Wii              F\n";
    compute_n_vlaid_all_points();
    if(irock==1)
        finding_in_rock_iniR();
    geoplot();
    prenumbe = numbe;
    if (exca.ID_Exca == 1)
    {
        ExcavationCracks();    // Excavation induced cracks within a distance into wall
        geoplot();
        prenumbe = numbe;
    }
    If_No_tip();    //Check possibility of fracture initiation if no tip
    //--------------------------Creep functions---------------    
    file50 << "   time     time step    Tip no.    Creep growth length    growth angle      K/Kc      Crack velocity\n";       
    creep.time = 0;
    while (!StopReturn)
    {
        mcyc++;  
        cout << " cycle " << comvar::mcyc << " of " <<mcyc0<<"  is running"<<" starting at(" << getTimeOnly()<<").... ";
        auto starti = std::chrono::high_resolution_clock::now();
        creep.deltaT = creep.deltaT_min;         //Creep iteration
        creep.ID_creep = (creep.time < creep.totalT) ? 1 : 0;

        if (creep.ID_creep == 1)
            creep_problem();
        else
        {           
            main_computing_part();
        } 
        auto endi = std::chrono::high_resolution_clock::now();
        std::cout <<" finished! (" << 
            std::chrono::duration<double>(endi - starti).count() << " sec)\n";

        // Pause if defined cycle finishes
        if (mcyc >= mcyc0)
        {
            if (lastinput != "endf")
            {
                input();
            }
            if (lastinput == "endf")
                { 
                    auto end1 = std::chrono::high_resolution_clock::now();                                  
                    std::cout << "Run time: " << std::chrono::duration<double>(end1 - start).count() << " seconds\n";
                    logfile << "Run time=  " << std::chrono::duration<double>(end1 - start).count() << " seconds\n";
                    MessageBox(nullptr, L"End of cycle & input file! press OK to quit.", L"Message!", MB_OK);  
                    return;                             
                }
        }                    
    }     
   
return;   
}

