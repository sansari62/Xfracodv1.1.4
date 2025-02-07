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





int  check_material_id(float xp, float yp)
{

    float dist0 = 10e8;
    int mm = 0, mclosest = 0;

    int k = 0;
    float xc = 0, yc = 0, xt = 0, yt = 0, dist = 0;
    //remove prenumbe from for loop because of the issue with multi region problem

    for (int m = 0 ; m < numbe; ++m)
    {
        dist = std::sqrt(std::pow(xp - elm_list[m].xm, 2) + std::pow(yp - elm_list[m].ym, 2));

        if (dist <= dist0)
        {
            mclosest = m;
            dist0 = dist;
            k = 0;
        }

        if (symm.ksym != 0)
        {
            if (symm.ksym == 1 || symm.ksym == 4)
            {
                xc = 2.0 * symm.xsym - elm_list[m].xm;
                yc = elm_list[m].ym;
                dist = std::sqrt(std::pow(xp - xc, 2) + std::pow(yp - yc, 2));

                if (dist < dist0)
                {
                    mclosest = m;
                    k = 1;
                    dist0 = dist;
                }
            }
            if (symm.ksym == 2 || symm.ksym == 4)
            {
                xc = elm_list[m].xm;
                yc = 2.0 * symm.ysym - elm_list[m].ym;
                dist = std::sqrt(std::pow(xp - xc, 2) + std::pow(yp - yc, 2));
                if (dist < dist0)
                {
                    mclosest = m;
                    k = 2;
                    dist0 = dist;
                }
            }
            if (symm.ksym == 3 || symm.ksym == 4)
            {
                xc = 2.0 * symm.xsym - elm_list[m].xm;
                yc = 2.0 * symm.ysym - elm_list[m].ym;
                dist = std::sqrt(std::pow(xp - xc, 2) + std::pow(yp - yc, 2));

                if (dist < dist0)
                {
                    mclosest = m;
                    k = 3;
                    dist0 = dist;
                }
            }
        }
        // label_51: continue;
    }

    if (elm_list[mclosest].kod != 6)
    {
        mm = elm_list[mclosest].mat_no;       
        return mm;
    }
    else
    {
        xt = xp;
        yt = yp;

        switch (k)
        {
        case 1:
            xt = 2.0 * symm.xsym - xp;
            break;
        case 2:
            yt = 2.0 * symm.ysym - yp;
            break;
        case 3:
            xt = 2.0 * symm.xsym - xp;
            yt = 2.0 * symm.ysym - yp;

        }

        float xpprime = (xt - elm_list[mclosest].xm) * elm_list[mclosest].cosbet +
            (yt - elm_list[mclosest].ym) * elm_list[mclosest].sinbet;
        float ypprime = -(xt - elm_list[mclosest].xm) * elm_list[mclosest].sinbet +
            (yt - elm_list[mclosest].ym) * elm_list[mclosest].cosbet;

        if (ypprime <= 0)
        {
            mm = elm_list[mclosest].mat_no;
        }
        else
        {
           
            if (elm_list[mclosest - 1].xm == elm_list[mclosest].xm &&
                elm_list[mclosest - 1].ym == elm_list[mclosest].ym)
            {
                mm = elm_list[mclosest - 1].mat_no;
            }
            if (elm_list[mclosest + 1].xm == elm_list[mclosest].xm &&
                elm_list[mclosest + 1].ym == elm_list[mclosest].ym)
            {
                mm = elm_list[mclosest + 1].mat_no;
            }
           
        }
    }

    return mm;
}





void Interface(int num, float xbeg, float ybeg, float xend, float yend, int mat1, int mat2)
{
    // Call arraysize(newmem)

    float st, angd, ang0, a0, xs, ys, xd, yd, sw;
    int  m, n, mn, ms, nn, ns;
    float pxx, pyy, pxy;
    float cosb, sinb, ss, sn;


    // Define locations, size, orientations of multiregion interfaces
    st = sqrt((xend - xbeg) * (xend - xbeg) + (yend - ybeg) * (yend - ybeg));
    angd = pi / num;
    ang0 = 0.5 * angd;
    a0 = st * 0.5;

    xs = xbeg;
    ys = ybeg;

    xd = (xend - xbeg) / num;
    yd = (yend - ybeg) / num;
    sw = sqrt(xd * xd + yd * yd);

    int numbe0 = numbe;
    for (int ne = 0; ne < num; ++ne)
    {
       
        numbe += 2;
        m = numbe - 2;
        n = numbe - 1;
        BoundaryElement& be = elm_list[m];

        sw = st / num;
        xd = (xend - xbeg) * sw / st;
        yd = (yend - ybeg) * sw / st;
        xs += xd;
        ys += yd;

        be.xm = xs - 0.5 * xd;
        be.ym = ys - 0.5 * yd;
        be.a = 0.5 * sw;
        be.sinbet = yd / sw;
        be.cosbet = xd / sw;
        be.kod = 6;
        b_elm[m].ipair = 1;
        be.mat_no = mat1;

        elm_list[n].xm = be.xm;
        elm_list[n].ym = be.ym;
        elm_list[n].a = be.a;
        elm_list[n].sinbet = -be.sinbet;
        elm_list[n].cosbet = -be.cosbet;
        elm_list[n].kod = 6;
        b_elm[n].ipair = 2;
        elm_list[n].mat_no = mat2;

        //************************* Avoid double definition of elements********************************
        bool jmp_label = false;
        for (int ii = 0; ii < m-1; ++ii) 
        {
            if (be.xm == elm_list[ii].xm && be.ym == elm_list[ii].ym &&
                be.a == elm_list[ii].a && be.sinbet == elm_list[ii].sinbet && be.cosbet == elm_list[ii].cosbet) 
            {
                numbe -= 2;    // I think those two newly added elements should be removed and the size of vector elm_list controlled by vec.size               
                jmp_label = true;
                break;
            }
        }
        if (jmp_label) continue;   //jmp to label110
        ms = 2 * m;
        mn = ms + 1;
        s4.b1[ms] = 0;
        s4.b1[mn] = 0;

        ns = 2 * n;
        nn = ns + 1;
        s4.b1[ns] = 0;
        s4.b1[nn] = 0;

        //------------------ Adjust stress boundary values to account for initial stresses-------------------------
        float y0 = g.y_surf;
        pxx = symm.pxx1 + g.skx * (y0 - be.ym);    
        pyy = symm.pyy1 + g.sky * (y0 - be.ym);
        pxy = symm.pxy1;

        cosb = be.cosbet;
        sinb = be.sinbet;
        ss = (pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb);
        sn = pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb;

        // Zero at interface - it is not real boundary value; for concrete lining insitu stresses should not be included
        if (be.mat_no == mat_lining || elm_list[n].mat_no == mat_lining)
        {
            s4.b0[ms] = -ss;                //interface element
            s4.b0[mn] = -sn;
        }
        else 
        {
            s4.b0[ms] = 0;                  //interface element
            s4.b0[mn] = 0;
        }
        s4.b0[ns] = 0;                      //interface element
        s4.b0[nn] = 0;

        // Calculate forces         !Force1 etc is not function of b0, but pxx etc.
        //still work done by interface if it is between lining and rock
        if (be.mat_no == mat_lining) 
        {
            b_elm[m].force1 = 0;
            b_elm[m].force2 = 0;
        }

        //still work done by interface if it is between lining and rock
        else
        {
            b_elm[m].force1 = 2.0 * be.a * (-ss);
            b_elm[m].force2 = 2.0 * be.a * (-sn);
        }

        if (elm_list[n].mat_no == mat_lining)
        {
            b_elm[n].force1 = 0;
            b_elm[n].force2 = 0;
        }
        else
        {
            b_elm[n].force1 = 2.0 * elm_list[n].a * (-ss);
            b_elm[n].force2 = 2.0 * elm_list[n].a * (-sn);
        }
    }
    //label110
    numbe_old = numbe;
    return;
}





int cross(float xb1, float yb1, float xe1, float ye1, float xb2, float yb2, 
    float& xe2, float& ye2)
{ 
    /* Check if the tip element crosses any existing elements,
    if yes, merge tip element to existing elements
    */

    int id = 0;
    float tan1, tan2, xcross, ycross, db, de;

    if (xe1 == xb1) 
    {
        tan1 = 10e20 * (ye1 - yb1);
    }
    else
    {
        tan1 = (ye1 - yb1) / (xe1 - xb1);
    }

    if ((xb1 == xb2 && yb1 == yb2) || (xb1 == xe2 && yb1 == ye2) ||
     (xe1 == xb2 && ye1 == yb2) || (xe1 == xe2 && ye1 == ye2) ) 
        return id;
   

    if (xe2 == xb2)
    {
        tan2 = 10e20 * (ye2 - yb2);
    }
    else
    {
        tan2 = (ye2 - yb2) / (xe2 - xb2);
    }

    if (tan1 == tan2) 
    {
        return id;
    }

    xcross = (yb2 - yb1 + tan1 * xb1 - tan2 * xb2) / (tan1 - tan2);
    ycross = (tan1 * tan2 * (xb1 - xb2) + tan1 * yb2 - tan2 * yb1) / (tan1 - tan2);

    if (xcross <= min(xb1, xe1) - 0 * s5u.dtol || ycross <= min(yb1, ye1) - 0 * s5u.dtol ||
        xcross >= max(xb1, xe1) + 0 * s5u.dtol || ycross >= max(yb1, ye1) + 0 * s5u.dtol ||
        xcross <= min(xb2, xe2) - 0 * s5u.dtol || ycross <= min(yb2, ye2) - 0 * s5u.dtol ||
        xcross >= max(xb2, xe2) + 0 * s5u.dtol || ycross >= max(yb2, ye2) + 0 * s5u.dtol) 
    {
        return id;
    }

    db = sqrt(pow(xcross - xb1, 2) + pow(ycross - yb1, 2));
    de = sqrt(pow(xcross - xe1, 2) + pow(ycross - ye1, 2));

    if (db < de) 
    {
        xe2 = xb1;
        ye2 = yb1;
    }
    else
    {
        xe2 = xe1;
        ye2 = ye1;
    }

    id = 1;
    return id; 
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
    If_No_tip(selectedFile);    //Check possibility of fracture initiation if no tip

    //--------------------------Creep functions---------------    
    file50 << "   time     time step    Tip no.    Creep growth length    growth angle      K/Kc      Crack velocity\n";
       
    creep.time = 0;
    while (!StopReturn)
    {
        mcyc++;  
        cout << "\n cycle " << comvar::mcyc << " is running.... ";

        creep.deltaT = creep.deltaT_min;         //Creep iteration
        creep.ID_creep = (creep.time < creep.totalT) ? 1 : 0;
       
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

            if (insituS.incres != 0 && ktipgrow == true) //special treatment th final grownth element in cycle
            {
                int incres_tem = insituS.incres;
                insituS.incres = 0;
                work0(-1);              //mode = -1 is back step
                insituS.incres = incres_tem;
            }
            if (insituS.incres != 0)   
                n_it = 20;    //iteration process
            work0(0);     //First calculation in every cycle - mode = 0. In work0(), only - 1 or 0 is allowed

            // --------------------------------------------------------
            check_crack_growth();              //include instant and creep
            geoplot();
            if (StopReturn == true) return;

            //----------------------------------------------------------
            prenumbe = numbe;
            add_crack_growth(); 
            cout << " finished!!! ";
            if (creep.time == creep.totalT || creep.ID_creep == 0 || creep.ID_fast_crack == 1)                
                break;
            //if not a creep problem or if fast crack growth, exit
        } while (creep.time != creep.totalT && creep.ID_creep != 0 && creep.ID_fast_crack != 1); 
        //END creep time cycle

                    // Pause if defined cycle finishes
                    if (mcyc >= mcyc0)
                    {
                        if (lastinput != "endf")
                        {
                            //MessageBox(nullptr, L"Defined cycle is completed, continue from input file.", L"Message!", MB_OK);
                            input();
                        }
                        else                                       
                            {                                
                                    auto end1 = std::chrono::high_resolution_clock::now();
                                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start);
                                     cout<< "Run time=  "<< duration.count();
                                    MessageBox(nullptr, L"End of cycle & input file! press OK to quit.", L"Message!", MB_OK);  
                                    return;                             
                            }
                    }                    
    }     
   
return;   
}

