#include<stdafx.h>
#include "CommonPara.h"
#include<Mainb.h>
#include<Failure.h>

using namespace CommonPara_h::comvar;




void safety_check()
{
    //check if the values for each element are reasonable

    for (int m = 0; m < numbe; m++)
    {
         if (elm_list[m].a == 0.0)  
             MessageBox(nullptr, L"Element length is zero", L"Warning!", MB_OK);

        if (elm_list[m].kod == 5)
        {
            if (b_elm[m].aks == 0.0 || b_elm[m].akn == 0.0)
                MessageBox(nullptr,
                    L"Element stiffness is zero.\n"
                    L"Please check jmat in the input file.",
                    L"Warning!",MB_OK);
            if (joint[m].aperture0 == 0.0 || joint[m].aperture_r == 0.0)
                MessageBox(nullptr, L"Joint aperture is zero", L"Warning!", MB_OK);

        }
        for (int n = m + 1; n < numbe; n++)
        {
            if (elm_list[m].xm == elm_list[n].xm && elm_list[m].ym == elm_list[n].ym && 
                elm_list[m].cosbet == elm_list[n].cosbet && elm_list[m].sinbet == elm_list[n].sinbet)
            {  
                MessageBox(nullptr, L"Elements overlapped", L"Warning!", MB_OK);

                string messag = "element n has overlap with m " + to_string(n);
                OutputDebugStringA(messag.c_str());
                exit(EXIT_FAILURE);
            }
        }
    }
    return;  
}





void increment()
{
    /*
    This subroutine provides the up-to-date matrix s4.b0() 
    in any iteration cycle
    For growing cracks, only the old elements are reassigned.
    The newly grown elements have their s4.b0() brought from the failure process,
    and defined in newtip
   */

    int ms = 0, mn = 0;
    float dss, dnn;    
    // Far-field stress change
    if (insituS.incres == 1)
    {
        symm.pxx1 += insituS.dsxx / n_it;
        symm.pxy1 += insituS.dsxy / n_it;
        symm.pyy1 += insituS.dsyy / n_it;
    }
    // All boundaries have stress change
    else if (insituS.incres == 2)
        {
            for (int m = 0; m < numbe_old; ++m)
            {
                if (elm_list[m].kod == 5 || elm_list[m].kod == 6)
                    continue;

                ms = 2 * m;
                mn = ms + 1;

                s4.b1[ms] += insituS.dss / n_it; //negative stress is compression
                s4.b1[mn] += insituS.dnn / n_it;
            }
        }
    // Some boundaries stress change
        else if (insituS.incres == 3)
            {
                string IDD;
                file9.open(dir + L"/Cbound.dat",std::ios::in);
                
                while(file9>> IDD)
                {
                    if (IDD == "dbou")
                    {
                        float x1 = 0, x2 = 0, y1 = 0, y2 = 0;
                        file9 >> x1 >> x2 >> y1 >> y2 >> dss >> dnn;

                        for (int m = 0; m < numbe_old; ++m)
                        {
                            if (elm_list[m].kod == 5 || elm_list[m].kod == 6)
                                continue;

                            if (elm_list[m].xm < x1 || elm_list[m].xm > x2 ||
                                elm_list[m].ym < y1 || elm_list[m].ym > y2)
                                continue;
                            ms = 2 * m;
                            mn = ms + 1;
                            s4.b1[ms] += dss / n_it;
                            s4.b1[mn] += dnn / n_it;
                        }
                    }
                    else if (IDD == "darc")
                    {
                        float xcen, ycen, diam1, diam2, ang1, ang2;
                        file9 >> xcen >> ycen >> diam1 >> diam2 >> ang1 >> ang2 >> dss >> dnn;

                        for (int m = 0; m < numbe_old; ++m)
                        {
                            if (elm_list[m].kod == 5 || elm_list[m].kod == 6)
                                continue;

                            float radius = sqrt(pow(elm_list[m].xm - xcen, 2) + pow(elm_list[m].ym - ycen, 2));
                            float ang = atan2f((elm_list[m].ym - ycen) / radius,
                                (elm_list[m].xm - xcen) / radius) * 180.0 / 3.14;  

                            if (radius < diam1 / 2 || radius > diam2 / 2 || ang < ang1 || ang > ang2)
                                continue;

                            ms = 2 * m;
                            mn = ms + 1;

                            s4.b1[ms] += dss / n_it;
                            s4.b1[mn] += dnn / n_it;
                        }
                    }
                }             

                file9.close();
            }
   
    // Final treatment
    for (int m = 0; m < numbe_old; ++m)
    {
        MatrixB(m);
    }

    if (insituS.incres != 0)
        ID_dtip = 0;

    return;
}





/*----------------------------------------
       Hydraulic pressure
  ----------------------------------------*/
//only called from work0(), no water is considered for fictitious element

static void water()
{
    vector<float> jwater_old(numbe, 0.0);
    int ms = 0, mn = 0;
    float dist = 0.0;

    // Save the current state of jwater to jwater_old
    for (int m = 0; m < numbe; ++m)
    {
        jwater_old[m] = watercm.jwater[m];
    }

    // Initialize jwater and pwater
    for (int m = 0; m < numbe; ++m)
    {
        watercm.jwater[m] = 0;
        watercm.pwater[m] = 0.0;
    }

    // Check and set jwater and pwater based on conditions
    for (int i = 0; i < watercm.iwhole; ++i)
    {
        for (int m = 0; m < numbe; ++m)
        {

            dist = std::sqrt(std::pow(elm_list[m].xm - watercm.w_xc[i], 2) +
                std::pow(elm_list[m].ym - watercm.w_yc[i], 2));

            if (dist <= watercm.w_d[i] / 2.0)
            {
                watercm.jwater[m] = 1;   
                watercm.pwater[m] = watercm.wph[i];
            }
        }
    }

    for (int i = 0; i < watercm.iwrect; ++i)
    {
        for (int m = 0; m < numbe; ++m)
        {
            if ((elm_list[m].xm <= watercm.w_x2[i] && elm_list[m].xm >= watercm.w_x1[i]) &&
                (elm_list[m].ym <= watercm.w_y2[i] && elm_list[m].ym >= watercm.w_y1[i]))
            {
                watercm.jwater[m] = 1;
                watercm.pwater[m] = watercm.wpr[i];
            }
        }
    }

    int ichange = 0;
    do
    {
        ichange = 0;
        // Merge elements based on conditions

        for (int m = 0; m < numbe; ++m)
        {
            if (elm_list[m].kod != 5 || watercm.jwater[m] == 1)
                continue;

            float xb = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
            float yb = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;
            float xe = elm_list[m].xm + elm_list[m].a * elm_list[m].cosbet;
            float ye = elm_list[m].ym + elm_list[m].a * elm_list[m].sinbet;

            for (int n = 0; n < numbe; ++n)
            {
                if (watercm.jwater[n] != 1)
                    continue;

                float xb0 = elm_list[n].xm - elm_list[n].a * elm_list[n].cosbet;
                float yb0 = elm_list[n].ym - elm_list[n].a * elm_list[n].sinbet;
                float xe0 = elm_list[n].xm + elm_list[n].a * elm_list[n].cosbet;
                float ye0 = elm_list[n].ym + elm_list[n].a * elm_list[n].sinbet;

                if ((std::abs(xb - xb0) < s5u.dtol && std::abs(yb - yb0) < s5u.dtol) ||
                    (std::abs(xb - xe0) < s5u.dtol && std::abs(yb - ye0) < s5u.dtol) ||
                    (std::abs(xe - xb0) < s5u.dtol && std::abs(ye - yb0) < s5u.dtol) ||
                    (std::abs(xe - xe0) < s5u.dtol && std::abs(ye - ye0) < s5u.dtol))
                {
                    watercm.jwater[m] = 1;
                    if (watercm.pwater[m] != watercm.pwater[n])
                    {
                        watercm.pwater[m] = max(watercm.pwater[m], watercm.pwater[n]);
                        watercm.pwater[n] = max(watercm.pwater[m], watercm.pwater[n]);
                        ichange = 1;
                    }
                }
            }
        }
    } while (ichange == 1);

    // Update s4.b0 based on conditions  
    for (int m = 0; m < numbe; ++m)
    {
        ms = 2 * m;
        mn = ms + 1;

        if (watercm.jwater[m] == 1.0 && jwater_old[m] == 0)
        {
            switch (elm_list[m].kod)
            {

            case 1:
            case 3:
            case 5:
            case 11:
            case 13:
            case 15:
                s4.b0[ms] = s4.b0[ms];              //fracture elements
                s4.b0[mn] = s4.b0[mn] - watercm.pwater[m];
                break;
            default:
                break;
                //  interfcae element, nothing change from the original boundary condition
            }            
        }
    }

    return;
}






void monitoring_point(float xp, float yp, float& sigxx, float& sigyy, float& sigxy, float& ux, float& uy)
{
    /*
    * the stress and the monitoring points user defined are calculated 
    * and later saved.
    */

    int mm = j_material;
    if (multi_region)
        mm = check_material_id(xp, yp);

    s2us.reset();
    ux = 0.0;
    uy = 0.0;

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - yp);
    float pyy = symm.pyy1 + g.sky * (y0 - yp);
    float pxy = symm.pxy1;

    sigxx = pxx;
    sigyy = pyy;
    sigxy = pxy;

    if (mm == mat_lining)
    {
        sigxx = 0;
        sigyy = 0;
        sigxy = 0;
    }

    for (int j = 0; j < numbe; ++j)
    {
        BoundaryElement& b_ei = elm_list[j];
        if (mm != b_ei.mat_no) continue;
        int js = 2 * j;
        int jn = js + 1;
        s2us.reset();
        float xj = b_ei.xm;
        float yj = b_ei.ym;
        float aj = b_ei.a;

        float cosbj = b_ei.cosbet;
        float sinbj = b_ei.sinbet;

        coeff(xp, yp, xj, yj, aj, cosbj, sinbj, +1, mm);
        switch (symm.ksym + 1)
        {
        case 1: break;
        case 2:
            xj = 2.0 * symm.xsym - b_ei.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            break;

        case 3:
            yj = 2.0 * symm.ysym - b_ei.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            break;

        case 4:
            xj = 2.0 * symm.xsym - b_ei.xm;
            yj = 2.0 * symm.ysym - b_ei.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, -sinbj, +1, mm);
            break;

        case 5:
            xj = 2.0 * symm.xsym - b_ei.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            xj = b_ei.xm;
            yj = 2.0 * symm.ysym - b_ei.ym;

            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            xj = 2.0 * symm.xsym - b_ei.xm;
            coeff(xp, yp, xj, yj, aj, -cosbj, -sinbj, +1, mm);
        }

        ux += s2us.uxs * s4.d0[js] + s2us.uxn * s4.d0[jn];
        uy += s2us.uys * s4.d0[js] + s2us.uyn * s4.d0[jn];
        sigxx += s2us.sxxs * s4.d0[js] + s2us.sxxn * s4.d0[jn];
        sigyy += s2us.syys * s4.d0[js] + s2us.syyn * s4.d0[jn];
        sigxy += s2us.sxys * s4.d0[js] + s2us.sxyn * s4.d0[jn];
    }

    return;
}





void write_monitor_data(int file_type, int i, int mcyc, int it,  float xmon,
     float ymon, float sigxx, float sigyy, float sigxy, float uxneg, float uyneg)
{
    ofstream* filem = nullptr;
    if (file_type == 0)
        filem = &mon_files[i+1];
    else
        filem = &ml_files[i+1];

    *filem << "  "
        << std::setw(6) << mcyc << "(" << std::setw(2) << it << ")" << "      "
        << std::scientific << std::setprecision(3)
        << std::setw(11) << creep.time                   // e11.3
        << "  " << std::fixed << std::setprecision(3)
        << std::setw(8) << xmon                 // f8.3
        << "  " << std::setw(8) << ymon         // f8.3
        << "  " << std::scientific << std::setprecision(3)
        << std::setw(11) << sigxx
        << "  " << std::setw(11) << sigyy
        << "  " << std::setw(11) << sigxy
        << "  " << std::setw(11) << uxneg
        << "  " << std::setw(11) << uyneg
        << "  " << std::setw(11) << creep.vel_creep_max
        << std::endl;
   }





void third_correction_run(int& it)
{
    float dist;
    float ss, sn, usneg, unneg, angle0, uxneg, uyneg, sigxx=0, sigyy=0, sigxy = 0;
    for(int i=0; i < ihist; i++)
    {
        bool jump = 0;
        for (int m = 0; m < numbe; m++)
        {
            if (elm_list[m].kod != 5 && elm_list[m].kod != 6 &&
                elm_list[m].kod != 7)
            {
                dist = sqrt(pow(comvar::elm_list[m].xm - mpoint_list[i].xmon, 2) +
                    pow(elm_list[m].ym - mpoint_list[i].ymon, 2));
                if (dist <= 1.1 * elm_list[m].a)
                {
                    ss = b_elm[m].sigma_s;
                    sn = b_elm[m].sigma_n;
                    usneg = b_elm[m].us_neg;
                    unneg = b_elm[m].un_neg;
                    angle0 = atan2f(elm_list[m].sinbet, elm_list[m].cosbet) * 180 / pi;
                    uxneg = usneg * elm_list[m].cosbet - unneg * elm_list[m].sinbet;
                    uyneg = usneg * elm_list[m].sinbet + unneg * elm_list[m].cosbet;

                    float sinbt = elm_list[m].sinbet;

                    if (abs(sinbt) > 0.95)
                    {
                        sigxx = sn;
                        sigyy = 0.0;
                        sigxy = ss;
                    }
                    else if (abs(sinbt) < 0.05)
                    {
                        sigxx = 0.0;
                        sigyy = sn;
                        sigxy = ss;
                    }
                    else //if (abs(sinbt) >= 0.05 && abs(sinbt) <= 0.95)
                    {
                        sigxx = sn;
                        sigyy = sn;
                        sigxy = ss;
                    }

                    write_monitor_data(0, i, mcyc, it, mpoint_list[i].xmon, mpoint_list[i].ymon,
                        sigxx, sigyy, sigxy, uxneg, uyneg);
                    jump = 1;
                    break;
                } //end if
            }
        }
    if (jump)
        continue;
    float ux, uy = 0;
    monitoring_point(mpoint_list[i].xmon, mpoint_list[i].ymon, sigxx, sigyy, sigxy, ux, uy); 
    write_monitor_data(0,i, mcyc, it,mpoint_list[i].xmon, mpoint_list[i].ymon,
        sigxx, sigyy, sigxy, ux, uy);
    }
    return;
 }




void run_check(int mode)
{
     mainb_work0( mode);
     s4.limit_d();   
     /*redo the d0() accumulation with correct fracture state*/
    for (int k = 0; k < numbe; k++)
    {
        int tk = k * 2;
        s4.d0[tk] +=  s4.d[tk];
        s4.d0[tk+1] += s4.d[tk+1];
    }
    return;
}




void calc_bound_stress(int it)
{
    float ss = 0, sn = 0, unneg=0, usneg = 0, untem = 0, ustem = 0;
    float ssd = 0, snd = 0;
   
    for (int m = 0; m < numbe; ++m)
    {
        BE& be = b_elm[m];
        elm_list[m].bound(m, ss, sn, ustem, untem, usneg, unneg);  //not sure Sara!
        ssd = ss - be.ss_old;
        snd = sn - be.sn_old;
        be.ss_old = ss;
        be.sn_old = sn;
        be.sigma_s += ssd; //fracture stress ss, sn is the total stress at a given moment,
        be.sigma_n += snd;  //sigma_s() is the history values
        be.us_neg = usneg;  //fracture displacement is already true and accomulative value
        be.un_neg = unneg;
        be.us = ustem;
        be.un = untem;           
        //  ----------------------------------------------------------------
        //Growth element becomes old normal fracture element now, and 
        // s4.b0() and s4.b0_old()
        // will be calculated in normal way for other elements
        s4.b0_old[m * 2] = s4.b0[m * 2 ];
        s4.b0_old[m*2 + 1] = s4.b0[m * 2 + 1];
        //  !----------------------------------------------------------------
        if (it == n_it)
        {
            be.us = ustem;
            be.un = untem;
            be.forces = be.force1 + ss * 2.0 * elm_list[m].a;
            be.forcen = be.force2 + sn * 2.0 * elm_list[m].a;
            s4.df[m * 2] = s4.d0[m * 2]; // !for (AE)
            s4.df[m * 2+1] = s4.d0[m*2+1]; // !for (AE)
        }
    }
    return;
}





//mode = -1, back step; = 0 normal w0 calculation; should not be 1 or 2 there

void work0(int mode)
{
    /*-calculate the total strain enegery before fracture growth-*/

    if (water_mod)
        water();

    safety_check();   

    //first time consider it as elastic fracture
    for (int m = 0; m < numbe; m++)
        b_elm[m].jstate = 0;                

    if (insituS.incres == 0) 
        n_it = 1;
    float thre = 0;
    for (int it = 1; it <= n_it; it++)   
    {
        increment();              //only old elements are used here(numbe_old)
        //the growth element is treated differently.s4.b0() is calculated as stress in solid near tip in subroutine newtips

        for (int k = 0; k < numbe; k++)
        {
            int tk = k * 2;
            //for old grown elements, s4.b0() needs to be calculated from pxx etc.
            if ( n_it == 1 && k < numbe_old ) 
            {
                s4.b0_old[tk] = s4.b0[tk];
                s4.b0_old[tk+1] = s4.b0[tk+1];
            }
            s4.b[tk] = s4.b0[tk] - s4.b0_old[tk];   //b for old elements become 0 and for new ones I don't know, do we need else here???
            s4.b[tk+1] = s4.b0[tk+1] - s4.b0_old[tk +1];
        }
        run_check(mode);

        float ssd = 0, streng = 0, shear = 0, snd = 0;
        float Sigma_n_prime = 0, sntem = 0, beta = 0;
        float  S_error = 0, T_error = 0;
        float ustem = 0, untem = 0, ss = 0, sn = 0, usneg = 0, unneg = 0;

        for (int m = 0; m < numbe; ++m)
        {
            //for each element all values of input usually the stresses are calculated and returned here
            //using d(m), or incremental d to calculate stresses
            BE& be = b_elm[m];
            elm_list[m].bound(m, ss, sn, ustem, untem, usneg, unneg);

            /* fracture stress ss,sn is the total stress at a given moment,
            sigma_s() is the history values
            fracture displacement is already true and accomulative value*/

            ssd = ss - be.ss_old;
            snd = sn - be.sn_old;
            be.sigma_s += ssd;
            be.sigma_n += snd;
            be.us_neg = usneg;  //displacement 
            be.un_neg = unneg;
            be.us = ustem;  
            be.un = untem;

            /* open fracture
           for initially open fractures, sigma_n is nearly 0, untem is high in open (negative)
           for initially sliding fractures, sigma_n is high (negative), but untem is nearly 0 (depending on kn)
           when reverse loading, both sigma_n and untem need to be considered. Otherwise, wrong fracture openning
           state will be resulted at the initial step where fractures are assumed to be elastic.*/

            if (elm_list[m].kod == 5)
            {
                // Open fracture conditions
                Sigma_n_prime = be.sigma_n + watercm.pwater[m] * watercm.jwater[m];
                const float eps = 0;// 1e-5f;
                if (Sigma_n_prime > -1e4 + thre && untem < thre)
                {
                    be.jstate = 1;
                    be.jslipd = 0;
                    be.coh = 0;                    
                }
                else
                {
                    // shear fracture conditions
                    //water pressure/effective stress. sn - total stress. careful about + and -
                    beta = be.phi;
                    streng = max(-(be.sigma_n + watercm.pwater[m] * watercm.jwater[m]) *
                        tanf(beta), 0.0f) + be.coh;

                    if (abs(ssd) <= 1e4)   
                        ssd = 0;                        
                    //if ssd is 0, the shear direction follows the overall shear.
                    //sign() is to consider the reverse shearing case

                    shear = (ssd == 0) ? be.sigma_s * std::copysign(1.0,be.sigma_s) :
                        be.sigma_s * std::copysign(1.0, ssd);
                    
                    if (shear > streng * 0.99)
                    {
                        be.jstate = 2;
                        be.jslipd = -std::copysign(1, ssd);
                        be.coh = 0;                        
                    }
                }
            }
            //label30
            //undo the sigma_s() etc accumulation for use in the 2nd run           

            be.sigma_s -= ssd;
            be.sigma_n -= snd;
        }

        //--------------second check-run---------------------------------
        /* re-apply the same boundary values as in the 1st run*/
        for (int k = 0; k < numbe; k++)
        {
            int ms = k * 2;
            int mn = ms + 1;
            s4.b[ms ] = s4.b0[ms] - s4.b0_old[ms];
            s4.b[mn] = s4.b0[mn] - s4.b0_old[mn];

            /* undo the d0() accumulation done in the 1st run*/
            s4.d0[ms] -=  s4.d[ms];
            s4.d0[mn] -=  s4.d[mn];
        }

       /*solve the problem for the 2nd run*/
        run_check(mode);
        calc_bound_stress(it);
        //IF no loose block, disregrad 3nd checking below
        if (ID_dtip != 0)
        {
            //--------------------- possible 3nd run - correction(below) -----------------------
           // check for slippage or tension for the case of loose block

            int ncorrection = 0;  
            for (int m = 0; m < numbe; ++m)
            {
                BE& be = b_elm[m];
                if (elm_list[m].kod == 5)
                {
                    // ------------sliding joint--------------
                    if (be.jstate == 2)
                    {
                        S_error = abs(be.sigma_s) - abs(tanf(be.phi) * be.sigma_n);

                        if (abs(S_error) > 0.1 * abs(be.sigma_s + 1e4))
                        {          
                            s4.b[m * 2] = -S_error * std::copysign(1.0, be.sigma_s);
                            s4.b[m * 2 + 1] = 0;
                            ncorrection = 1;
                        }
                    }
                    // ------------open joint-----------------
                    else if (be.jstate == 1)
                    {
                        T_error = abs(be.sigma_n);
                        if (abs(T_error) > 0)
                        {
                            s4.b[m * 2] = -be.sigma_s;
                            s4.b[m * 2 + 1] = -T_error;
                            ncorrection = 1;
                        }
                    }
                }
            }
            // if no correction is needed jump to third correc_run
            if (ncorrection != 0)
            {
                run_check(mode);
                calc_bound_stress(it);
            }           
        }
        // ---------- possible 3nd run - correction(above) ---------
       //label880 in code
        third_correction_run(it);
        //for monitoring lines
        float xp, yp;
        float sigxx, sigyy, sigxy, ux, uy;
         
        for (int i = 0; i < lhist; i++)
        {
            for (int j = 0; j <=mline_list[i].npl; j++)
            {
                xp = mline_list[i].x1l + j * (mline_list[i].x2l - mline_list[i].x1l) /mline_list[i].npl;
                yp =mline_list[i].y1l + j * (mline_list[i].y2l -mline_list[i].y1l) /mline_list[i].npl;
                monitoring_point(xp, yp, sigxx, sigyy, sigxy, ux, uy);  
                write_monitor_data(1, i, mcyc, it, xp, yp,sigxx, sigyy, sigxy, ux, uy);
            }
        }
    }
    n_it = 1; //re - set after iteration finish

    w0 =  0.0;
    for (int m = 0; m < numbe; ++m)
    {
        BE& be = b_elm[m];
        switch (elm_list[m].kod)
        {
        case 1:
        case 11:
        case 5:
        case 6:
            w0 += 0.5 * be.forces * be.us + 0.5 * be.forcen * be.un;
            break;
        case 2:
        case 12:
            w0 -= 0.5 * be.forces * be.us - 0.5 * be.forcen * be.un;
            break;
        case 3:
        case 13:
            w0 -= 0.5 * be.forces * be.us + 0.5 * be.forcen * be.un;
            break;
        case 4:
        case 14:
            w0 += 0.5 * be.forces * be.us - 0.5 * be.forcen * be.un;
            break;        
       
        }
    }

    numbe_old = numbe;

    if (mode != -1)
    {
        insituS.dsxx = 0;
        insituS.dsxy = 0;
        insituS.dsyy = 0;
        insituS.incres = 0;
    }
    else
    {
        // for first run after failure, s4.b0() needs to be calculated in normal way from pxx etc.
        //if (mode == -1)        
            for (int m = 0; m < numbe; ++m)
            {
                MatrixB(m);
                s4.b0_old[m * 2] = s4.b0[m * 2];
                s4.b0_old[m * 2 +1] = s4.b0[m * 2 +1];
            }       
    }
    return;
}

       



void work1(int mode)
{
    /* calculate the total strain energy when there is a fracture growth */
   
    float streng = 0.0;
    float ss = 0.0, ustem = 0.0,usneg = 0.0,unneg = 0.0,untem = 0.0;
    float sn = 0.0;
    auto& belm = b_elm[numbe - 1];
    float thre = 0;

    belm.jstate = 0;
    float ph = 30.0 * 3.1416 / 180.0;

    for (int k = 0; k < 2 * (numbe - 1); ++k)
         s4.b[k] = 0;   

    for (int k = 2 * numbe - 2; k < 2 * numbe; ++k)
            s4.b[k] = s4.b0[k];      //s4.b0(k) has been defined in sub newcoordinate as the stresses in intact rock

    mainb_work1(mode);
    //total stress/displacement using d0(m) total, not increment d(m)   
    elm_list[numbe-1].bound(numbe-1, ss, sn, ustem, untem, usneg, unneg);

    streng = max(- (sn + watercm.pwater[numbe - 1] * watercm.jwater[numbe - 1]) * tanf(ph), 0.0f);

    if (sn + watercm.pwater[numbe - 1] * watercm.jwater[numbe - 1] > 1e4 +thre)
    {
        belm.jstate = 1;
        belm.jslipd = 0;      
    }
    else
    {
        if (abs(ss) > streng+ thre)
        {
           belm.jstate = 2;
           belm.jslipd = -copysign(1.0, ss);
        }
        //Sara set m to 0   change m to numbe-1 Sara! 9.9.2024
        float epsilon = 0;// 1e-5;
        if (sn + watercm.pwater[numbe - 1] * watercm.jwater[numbe - 1] < thre && abs(ss) < streng+ thre)    //m instead of numbe check later in test
           belm.jstate = 0;
    }
    //label30
    if (belm.jstate != 0)
        mainb_work1(mode);

    //temperatory accumulation of element displacement, undo it later 
    for (int k = 0; k < 2 * numbe; k += 2)
    {
        s4.d0[k] += s4.d[k];
        s4.d0[k + 1] += s4.d[k + 1];
    }
    for (int m = 0; m < numbe; ++m)
    {
        auto& elm = elm_list[m];
        auto& belm = b_elm[m];
        elm.bound(m, ss, sn, ustem, untem, usneg, unneg);
        belm.us = ustem;
        belm.un = untem;
        belm.forces = belm.force1 + ss * 2.0 * elm.a;
        belm.forcen = belm.force2 + sn * 2.0 * elm.a;
    }


    //undo accumulation because the fictitious element is not real element yet
    for (int k = 0; k < 2 * numbe; ++k)
        s4.d0[k] -= s4.d[k];  
    w1 = 0.0;
    for (int m = 0; m < numbe; ++m)
    {
        switch (elm_list[m].kod)
        {
        case 1:
        case 11:
            w1 += 0.5 * b_elm[m].forces * b_elm[m].us + 0.5 * b_elm[m].forcen * b_elm[m].un;
            break;

        case 2:
        case 12:
            w1 -= 0.5 * b_elm[m].forces * b_elm[m].us - 0.5 * b_elm[m].forcen * b_elm[m].un;
            break;
        case 3:
        case 13:
            w1 -= 0.5 * b_elm[m].forces * b_elm[m].us + 0.5 * b_elm[m].forcen * b_elm[m].un;
            break;
        case 4:
        case 14:
            w1 += 0.5 * b_elm[m].forces * b_elm[m].us - 0.5 * b_elm[m].forcen * b_elm[m].un;
            break;
        case 5:
        case 6:
            w1 += 0.5 * b_elm[m].forces * b_elm[m].us + 0.5 * b_elm[m].forcen * b_elm[m].un;
        }
    }

   /* float sum = 0.0;
    for (int m = 0; m < numbe; ++m)
    {
        auto& belm = b_elm[m];
        double term1 = 0.5 * belm.forces * belm.us;
        double term2 = 0.5 * belm.forcen * belm.un;

        switch (elm_list[m].kod)
        {
        case 1: case 11:
        case 5: case 6:
            sum += term1 + term2; break;

        case 2: case 12:
            sum -= term1 - term2; break;

        case 3: case 13:
            sum -= term1 + term2; break;

        case 4: case 14:
            sum += term1 - term2; break;
        }
    }
    w1 = sum;*/


    //reset the fictitous element to no water pressure   
    // watercm.jwater[numbe - 1]= 0;

    return;
}









        










