#include <Tip.h>
#include "CommonPara.h"
#include<Work.h>
#include<Source.h>
#include<Input.h>
#include<Geoplot.h>
#include<Failure.h>
#include<ExcavationCracks.h>
#include<Fmax.h>
#include<Initiation.h>

using namespace CommonPara_h::comvar;



Tip::Tip() : xbe(0), ybe(0), xen(0), yen(0), ityp(0), imode(0), ifail(0), mpointer(0),
kindtip(0), f_value(0), dt(0), dl(0), nu(0), angl(0), mat_no(0),costem(0),sintem(0) {}



void Tip::assign_val(float x1, float y1, float x2, float y2, float dll, float cos1, 
    float sin1, int itp, int mat)
{
    xbe = x1;
    ybe = y1;
    xen = x2;
    yen = y2;
    dl = dll;
    sintem = sin1;
    costem = cos1;
    ityp = itp;
    mat_no = mat;
}




void Tip::save_to_file(ofstream& f) 
{
    f << dt << " " << dl << " " << ifail << " " 
        << nu << " " << xbe << " " << ybe << " " << xen << 
        " " << yen << " " << ityp << " " << imode << " " <<
        mpointer << " " << kindtip << " " << angl << " " << mat_no
        << " " << f_value << std::endl;
}






void label400_new_coordin_for_tip(int n, int mm, int mergtip, float xt, float yt)
{
    Tip& t = tips[n];
    if (t.ityp == +1 || t.ityp == 3 || t.ityp == 4)
    {
        t.xen = xt;
        t.yen = yt;
    }
    else if (t.ityp == -1 || t.ityp == -3) 
    {
        t.xbe = xt;
        t.ybe = yt;
    }
    t.dl = sqrt(pow(t.xen - t.xbe, 2) + pow(t.yen - t.ybe, 2));
    //-------------------------Add new element--------------------------
    int m = numbe; 
    elm_list[m].kod = 5;
    elm_list[m].a = t.dl / 2.0;
    elm_list[m].xm = 0.5 * (t.xen + t.xbe);
    elm_list[m].ym = 0.5 * (t.yen + t.ybe);
    elm_list[m].mat_no = mm;
    numbe++;
    BoundaryElement& be = elm_list[m];
   
    be.cosbet = (t.xen - t.xbe) / t.dl;
    be.sinbet = (t.yen - t.ybe) / t.dl;
    float sinb = be.sinbet;
    float cosb = be.cosbet;
    float sigxx = 0, sigyy = 0, sigxy = 0;

    NewFractureCentralPoint(be.xm, be.ym, sigxx, sigyy, sigxy);
    if (be.mat_no == mat_lining)
    {
        s4.b0[m * 2] = 0;
        s4.b0[m * 2 + 1] = 0;
    }
    else
    {
        float ss = (sigyy - sigxx) * sinb * cosb + sigxy * (cosb * cosb - sinb * sinb);
        float sn = sigxx * sinb * sinb - 2.0 * sigxy * sinb * cosb + sigyy * cosb * cosb;
        s4.b0[m * 2] = -ss;
        s4.b0[m * 2 + 1] = -sn;         
        if (sn < 0 && abs(ss) > (-sn * tanf(30.0 / 180.0 * pi)))  
        {
            //additional stress applied to crack growth element
            s4.b0[m * 2] = -(abs(ss) - (-sn * tanf(0.524))) * copysign(1, ss);
            s4.b0[m * 2 + 1] = 0;
        }
    }
    s4.df[m * 2] = s4.d0[m * 2]; 
    s4.df[m * 2 + 1] = s4.d0[m * 2 + 1];

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - be.ym);
    float pyy = symm.pyy1 + g.sky * (y0 - be.ym);
    float pxy = symm.pxy1;

    b_elm[m].force1 = 2.0 * be.a * (-((pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb)));// !old b0()
    b_elm[m].force2 = 2.0 * be.a * (-(pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb));// !old b0()
    b_elm[m].jmode = t.imode;
    b_elm[m].jstate = t.imode;
    float ss = (sigyy - sigxx) * sinb * cosb + sigxy * (cosb * cosb - sinb * sinb);
    b_elm[m].jslipd = -copysign(1.0, ss);

    float aks_bb = 0, akn_bb = 0, phi_bb = 0;
    float coh_bb = 0;
    float phid_bb = 0, ap_bb = 0, apr_bb = 0;
    stiffness_bb(aks_bb, akn_bb, phi_bb, coh_bb, phid_bb, ap_bb, apr_bb, tips[ni].imode, mm);
    b_elm[m].akn = akn_bb;
    b_elm[m].aks = aks_bb;
    b_elm[m].phi = phi_bb;;
    b_elm[m].phid = phid_bb;
    b_elm[m].coh = 0;
    joint[m].aperture0 = ap_bb;
    joint[m].aperture_r = apr_bb;

    //-------------------------------------------------
        if (mergtip == 1) 
        {
            t.ityp = 0;
            t.mpointer = 0;
        }           
        else
        {
            int mt = t.ityp / abs(t.ityp);  
            t.nu = 1;
            t.xbe = be.xm + be.a * be.cosbet * (2 - mt) * mt;
            t.ybe = be.ym + be.a * be.sinbet * (2 - mt) * mt;
            t.xen = t.xbe + 2.0 * be.a * be.cosbet;
            t.yen = t.ybe + 2.0 * be.a * be.sinbet;
            t.dl = 2.0 * be.a;
            t.ityp = mt;
            t.mpointer = m;
            t.mat_no = be.mat_no;
        }       
        
    return;
}




void specialLabel_200(int & mergtip,int i)
{
    mergtip = 1;
    comvar::ID_dtip = 1; 
    //tips disappeared after link to other element
    for (int k = 0; k < no; ++k)
    {
        if (tips[k].mpointer == i)
          tips[k].ityp = 0;
    }
}




void newtips(float dr)
    {    
        int n = ni;
        int id = 0;    // return value from cross func
        float xt = 0, yt = 0, xt0 = 0, yt0 = 0, xbeg = 0, xend = 0, ybeg = 0,
            yend = 0, tol = 0, tol1 = 0, xc = 0, yc = 0, dbeg = 0, dend = 0;

        if (tips[n].ityp == 0 || tips[n].ifail == 0)
        { 
            return;
        }
        //-----------------------
        if (tips[n].ityp == +1 || tips[n].ityp == 3 || tips[n].ityp == 4)
        {
            xt = tips[n].xbe + tips[n].dl * (cosf(dr) * tips[n].costem - sinf(dr) * tips[n].sintem);
            yt = tips[n].ybe + tips[n].dl * (sinf(dr) * tips[n].costem + cosf(dr) * tips[n].sintem);
            xt0 = tips[n].xbe;
            yt0 = tips[n].ybe;
        }

        else if (tips[n].ityp == -1 || tips[n].ityp == -3)
        {
            xt = tips[n].xen - tips[n].dl * (cosf(dr) * tips[n].costem - sinf(dr) * tips[n].sintem);
            yt = tips[n].yen - tips[n].dl * (sinf(dr) * tips[n].costem + cosf(dr) * tips[n].sintem);
            xt0 = tips[n].xen;
            yt0 = tips[n].yen;
        }
        int mm = tips[ni].mat_no;

        //----------------------------------------------
        int mergtip = 0;
        int m = tips[n].mpointer;

        if (! ((abs(elm_list[m].xm - symm.xsym) <= elm_list[m].a / 1000 && abs(elm_list[m].sinbet) >
            sinf(85.0 * pi / 180.0) && (symm.ksym == 1 || symm.ksym == 4)) ||
            (abs(elm_list[m].ym - symm.ysym) <= elm_list[m].a / 1000 && abs(elm_list[m].cosbet) >
            cosf(5.0 * pi / 180.0) && (symm.ksym == 2 || symm.ksym == 4))
            || symm.ksym == 0 ) )
        {            
            if (symm.ksym == 1 || symm.ksym == 4)
            {
                //element cross the symmetry line
                if ((xt - symm.xsym) * (xt0 - symm.xsym) <= 0)  
                {
                    mergtip = 1;
                    yt = yt - (xt - symm.xsym) / (xt - xt0) * (yt - yt0);
                    xt = symm.xsym;
                    label400_new_coordin_for_tip(n, mm, mergtip, xt, yt);
                    return;
                }
            }

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                if ((yt - symm.ysym) * (yt0 - symm.ysym) <= 0)  //!element cross the symmetry line
                {
                    mergtip = 1;
                    xt = xt - (yt - symm.ysym) / (yt - yt0) * (xt - xt0);
                    yt = symm.ysym;
                    label400_new_coordin_for_tip(n, mm, mergtip, xt, yt);
                    return;
                }
            }

            if (symm.ksym == 3 || symm.ksym == 4)
            {
                if ((xt - symm.xsym) * (yt0 - symm.ysym) == -(xt0 - symm.xsym) * (yt - symm.ysym))
                {
                    mergtip = 1;
                    xt = symm.xsym;
                    yt = symm.ysym;
                    label400_new_coordin_for_tip(n, mm, mergtip, xt, yt);
                    return;
                }              
            }
        }
        float z = 0.0;
        //------------------------------------
        for (int i = 0; i < numbe; ++i)
        {
            BoundaryElement& be = elm_list[i];
            if (i == m) continue; //do not check the mother element of the new tip
            xbeg = be.xm - be.a * be.cosbet;
            ybeg = be.ym - be.a * be.sinbet;
            xend = be.xm + be.a * be.cosbet;
            yend = be.ym + be.a * be.sinbet;
            xc = be.xm;
            yc = be.ym;

            tol = factors.tolerance;
            tol1 = factors.tolerance;

            float dc = min(sqrtf(powf(xt - xc, 2) + powf(yt - yc, 2)),
                sqrtf(powf(xt0 - xc, 2) + powf(yt0 - yc, 2)));

            if (dc <= tol1 * max(be.a, z))
            {

                dbeg = sqrt(pow(xt - xbeg, 2) + pow(yt - ybeg, 2));
                dend = sqrt(pow(xt - xend, 2) + pow(yt - yend, 2));
                if (dbeg <= dend)
                {
                    xt = xbeg;
                    yt = ybeg;
                }
                else
                {
                    xt = xend;
                    yt = yend;
                }
                specialLabel_200(mergtip, i);
                continue;
            }

            dbeg = min(sqrt(pow(xt - xbeg, 2) + pow((yt - ybeg), 2)),
                sqrt(pow(xt0 - xbeg, 2) + pow(yt0 - ybeg, 2)));

            dend = min(sqrt(pow(xt - xend, 2) + pow((yt - yend), 2)),
                sqrt(pow(xt0 - xend, 2) + pow(yt0 - yend, 2)));

            if (dbeg <= dend && dbeg <= tol * max(be.a, z))
            {
                xt = xbeg;
                yt = ybeg;
                specialLabel_200(mergtip, i);
                continue;
            }

            if (dend <= dbeg && dend <= tol * max(be.a, z))
            {
                xt = xend;
                yt = yend;
                specialLabel_200(mergtip, i);
                continue;
            }
            id = cross(xbeg, ybeg, xend, yend, xt0, yt0, xt, yt);
            if (id == 1) 
            {
                specialLabel_200(mergtip, i);
                continue;
            }

            if (symm.ksym != 0)
            {
                if (symm.ksym == 1 || symm.ksym == 4)
                {
                    dc = min(sqrt(pow(xt - (2. * symm.xsym - xc), 2) + pow(yt - yc, 2)), 
                        sqrt(pow(xt0 - (2. * symm.xsym - xc), 2) + pow(yt0 - yc, 2)));

                    if (dc <= tol1 * max(be.a, z))
                    {
                        dbeg = sqrt(pow((xt - (2. * symm.xsym - xbeg)), 2) + pow((yt - ybeg), 2));
                        dend = sqrt(pow((xt - (2. * symm.xsym - xend)), 2) + pow((yt - yend), 2));

                        if (dbeg <= dend)
                        {
                            xt = 2. * symm.xsym - xbeg;
                            yt = ybeg;
                        }
                        else
                        {
                            xt = 2. * symm.xsym - xend;
                            yt = yend;
                        }

                        specialLabel_200(mergtip,i);
                        continue;
                    }

                    dbeg = min(sqrt(pow(xt - (2. * symm.xsym - xbeg), 2) + pow(yt - ybeg, 2)), 
                        sqrt(pow(xt0 - (2. * symm.xsym - xbeg), 2) + pow(yt0 - ybeg, 2)));

                    dend = min(sqrt(pow(xt - (2. * symm.xsym - xend), 2) + pow(yt - yend, 2)),
                        sqrt(pow(xt0 - (2. * symm.xsym - xend), 2) + pow(yt0 - yend, 2)));

                    if (dbeg <= dend && dbeg <= tol * max(be.a, z))
                    {
                        xt = 2. * symm.xsym - xbeg;
                        specialLabel_200(mergtip, i);
                        continue;
                    }

                    if (dend <= dbeg && dend <= tol * max(be.a, z))
                    {
                        xt = 2. * symm.xsym - xend;
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                    id = cross(2. * symm.xsym - xbeg, ybeg, 2. * symm.xsym - xend, yend, xt0, yt0, xt, yt);  
                    if (id == 1)
                    {
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                }

                if (symm.ksym == 2 || symm.ksym == 4)
                {
                    dc = min(sqrt(pow(xt - xc, 2) + pow(yt - (2. * symm.ysym - yc), 2)), 
                        sqrt(pow((xt0 - xc), 2) + pow((yt0 - (2. * symm.ysym - yc)), 2)));
                    if (dc <= tol1 * max(be.a, 0))
                    {
                        dbeg = sqrt(pow((xt - xbeg), 2) + pow((yt - (2. * symm.ysym - ybeg)), 2));
                        dend = sqrt(pow((xt - xend), 2) + pow((yt - (2. * symm.ysym - yend)), 2));
                        if (dbeg <= dend)
                        {
                            xt = xbeg;
                            yt = 2.0 * symm.ysym - ybeg;
                        }
                        else
                        {
                            xt = xend;
                            yt = 2.0 * symm.ysym - yend;
                        }
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                    dbeg = min(sqrt(pow(xt - xbeg, 2) + pow(yt - (2. * symm.ysym - ybeg), 2)),
                        sqrt(pow(xt0 - xbeg, 2) + pow(yt0 - (2. * symm.ysym - ybeg), 2)));

                    dend = min(sqrt(pow(xt - xend, 2) + pow(yt - (2. * symm.ysym - yend), 2)),
                        sqrt(pow(xt0 - xend, 2) + pow(yt0 - (2. * symm.ysym - yend), 2)));

                    if (dbeg <= dend && dbeg <= tol * max(be.a, z))
                    {
                        yt = 2.0 * symm.ysym - ybeg;
                        specialLabel_200(mergtip, i);
                        continue;
                    }

                    if (dend <= dbeg && dend <= tol * max(be.a, z))
                    {
                        yt = 2.0 * symm.ysym - yend;
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                    id = cross(xbeg, 2.0 * symm.ysym - ybeg, xend, 2. * symm.ysym - yend, xt0, yt0, xt, yt);
                    if (id == 1)
                    {
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                }

                if (symm.ksym == 3 || symm.ksym == 4)
                {
                    dc = min(sqrt(pow(xt - (2. * symm.xsym - xc), 2) + pow(yt - (2. * symm.ysym - yc), 2)), 
                        sqrt(pow(xt0 - (2. * symm.xsym - xc), 2) + pow(yt0 - (2. * symm.ysym - yc), 2)));
                    if (dc <= tol1 * max(be.a, z))
                    {
                        dbeg = sqrt(pow(xt - (2. * symm.xsym - xbeg), 2) + pow((yt - (2. * symm.ysym - ybeg)), 2));
                        dend = sqrt(pow(xt - (2. * symm.xsym - xend), 2) + pow((yt - (2. * symm.ysym - yend)), 2));
                        if (dbeg <= dend)
                        {
                            xt = 2. * symm.xsym - xbeg;
                            yt = 2. * symm.ysym - ybeg;
                        }
                        else
                        {
                            xt = 2. * symm.xsym - xend;
                            yt = 2. * symm.ysym - yend;
                        }
                        specialLabel_200(mergtip, i);
                        continue;
                    }

                    dbeg = min(sqrt(pow(xt - (2. * symm.xsym - xbeg), 2) + pow(yt - (2. * symm.ysym - ybeg), 2)),
                        sqrt(pow((xt0 - (2. * symm.xsym - xbeg)), 2) + pow(yt0 - (2. * symm.ysym - ybeg), 2)));
                    dend = min(sqrt(pow(xt - (2. * symm.xsym - xend), 2) + pow(yt - (2. * symm.ysym - yend), 2)),
                        sqrt(pow(xt0 - (2. * symm.xsym - xend), 2) + pow(yt0 - (2. * symm.ysym - yend), 2)));
                    if (dbeg <= dend && dbeg <= tol * max(be.a, z))
                    {
                        xt = 2. * symm.xsym - xbeg;
                        yt = 2. * symm.ysym - ybeg;
                        specialLabel_200(mergtip, i);
                        continue;
                    }

                    if (dend <= dbeg && dend <= tol * max(be.a, z))
                    {
                        xt = 2. * symm.xsym - xend;
                        yt = 2. * symm.ysym - yend;
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                    id = cross(2. * symm.xsym - xbeg, 2. * symm.ysym - ybeg, 2. * symm.xsym - xend,
                        2. * symm.ysym - yend, xt0, yt0, xt, yt);
                    if (id == 1)
                    {
                        specialLabel_200(mergtip, i);
                        continue;
                    }
                }
            }
        }

        label400_new_coordin_for_tip( n, mm, mergtip, xt, yt);
        return;
 }



 void arrangetip()
{
        int id =0;

        do {            
            for (int i = 0; i < no; i++) 
            {  
                id = 0;
                if (tips[i].ityp == 0) 
                {
                    for (int j = i + 1; j < no; j++)
                    {
                        tips[j - 1].nu = tips[j].nu;
                        tips[j - 1].xbe = tips[j].xbe;
                        tips[j - 1].ybe = tips[j].ybe;
                        tips[j - 1].xen = tips[j].xen;
                        tips[j - 1].yen = tips[j].yen;
                        tips[j - 1].dl = tips[j].dl;
                        tips[j - 1].ityp = tips[j].ityp;
                        tips[j - 1].mpointer = tips[j].mpointer;
                        tips[j - 1].mat_no = tips[j].mat_no;
                    }
                    id = 1;
                    no--;
                }
            }
        } while (id == 1 && no != 0);
       
        return;
}




    void check_crack_growth()
    {
        float creep_growth_max = 0;
        float fm = 0;

        float angle = 0, vel = 0;
        float vel_creep = 0;
        creep.vel_creep_max = 0;
        creep.ID_fast_crack = 0;
        file2 << "Cycle#" << mcyc << endl;
        for ( ni = 0; ni < no; ni++) 
        {
            tips[ni].ifail = 0;
            if (tips[ni].ityp == 0) continue; //ityp() = -1, 0, 1, (-3; +3 ? )
            StopReturn = false;
            fmax1(fm,angle);  

            if (fm < 0) fm = 0;
            if (StopReturn == true) return; //stop
            if (abs(angle) >= 100 * pi / 180)
            {
                //if code did not find proper max, then ignore this tip
               tips[ni].angl = 0;
               tips[ni].f_value = 0;
            }
            else
            {
               tips[ni].angl = angle;
               tips[ni].f_value = fm;
            }

            if (fm >= 1.0)    //instant growth
            {
                tips[ni].ifail = 1;
                creep.growth_length[ni] = tips[ni].dl;
                creep_growth_max = 1.0;                
                if (tips[ni].imode == 1) vel = creep.v1;
                else if (tips[ni].imode == 2) vel = creep.v2;

                file50 << std::scientific << std::setprecision(4);
                file50 << std::setw(10) << creep.time << std::setw(1) << " " << std::setw(10) <<
                    creep.deltaT << std::setw(1) << " "
                    << std::setw(4) << ni << std::setw(8) << " " << std::setw(10) << 
                    creep.growth_length[ni] << std::setw(7) << " "
                    << std::setw(8) << tips[ni].angl * 180.0 / pi <<  std::setw(6) << 
                    std::sqrt(std::abs(fm)) << std::setw(5) << " "
                    << std::setw(10) << vel << "    (Fast crack growth)" << std::endl;
                creep.ID_fast_crack = 1;                
            }

            if (creep.ID_creep == 1)
            {
                if (fm >= 0.0 && fm < 1.0)    //creep growth
                {
                    if (tips[ni].imode == 1) 
                        vel_creep = creep.v1 * pow(fm / 1.0, (float(creep.nn1 / 2.))); //mode I creep propagation
                    else if (tips[ni].imode == 2)
                        vel_creep = creep.v2 * pow(fm / 1.0, (float(creep.nn2 / 2.))); //mode II creep propagation

                    float deltaL = vel_creep * creep.deltaT;
                    creep.creep_x[ni] = creep.creep_x[ni] + deltaL * cosf(tips[ni].angl);
                    creep.creep_y[ni] = creep.creep_y[ni] + deltaL * sinf(tips[ni].angl);
                    creep.creep_l[ni] = sqrt(pow(creep.creep_x[ni], 2) + pow(creep.creep_y[ni], 2));
                    creep.creep_a[ni] = atan2f(creep.creep_y[ni], creep.creep_x[ni]);

                    if (creep.creep_l[ni] < tips[ni].dl)   
                    {
                        file50 << std::scientific << std::setprecision(4);
                        file50 << std::setw(10) << time << std::setw(1) << " " << std::setw(10) << creep.deltaT << std::setw(1) << " "
                            << std::setw(4) << ni << std::setw(8) << " " << std::setw(10) << creep.creep_l[ni] << std::setw(7) << " "
                            << std::setw(8) << creep.creep_a[ni] * 180/ 3.14 << std::setw(7) << " " << std::setw(6) <<
                            std::sqrt(std::abs(fm)) << std::setw(5) << " "
                            << std::setw(10) << vel_creep << "    (Temporary creep growth)" << std::endl;

                    }

                    else //if (creep.creep_l[ni] >= tips[ni].dl)
                    {                       
                        tips[ni].ifail = 1;
                        creep.growth_length[ni] = tips[ni].dl;
                        tips[ni].angl = creep.creep_a[ni];
                        creep.creep_l[ni] = 0;
                        creep.creep_x[ni] = 0;
                        creep.creep_y[ni] = 0;
                        
                        file50 << std::scientific << std::setprecision(4);
                        file50 << std::setw(10) << time << std::setw(1) << " " << std::setw(10) << creep.deltaT << std::setw(1) << " "
                            << std::setw(4) << ni << std::setw(8) << " " << std::setw(10) << creep.growth_length[ni] << std::setw(7) << " "
                            << std::setw(8) << tips[ni].angl * 180 / 3.14 << std::setw(7) << " " << std::setw(6) <<
                            std::sqrt(std::abs(fm)) << std::setw(5) << " "
                            << std::setw(10) << vel_creep << "    (Creep growth)" << std::endl;

                    }
                    creep_growth_max = max(creep_growth_max, deltaL / tips[ni].dl);
                }
            }

            //no creep growth - elastic fracture tips
            if (fm <= 0.0)
            {
                vel_creep = 0;
                tips[ni].ifail = 0;
            }

            (((creep.vel_creep_max) > (vel_creep)) ? (creep.vel_creep_max) : (vel_creep)); 

        }

        // - automatically adjust time step within(min, max)
        if (creep_growth_max != 0.0) 
            creep.deltaT = creep.deltaT / (creep_growth_max * 10);
        else
            creep.deltaT = creep.deltaT_max;

        if (creep.deltaT > creep.totalT) creep.deltaT = creep.totalT / 10;
        if (creep.deltaT > creep.deltaT_max) creep.deltaT = creep.deltaT_max;
        if (creep.deltaT < creep.deltaT_min) creep.deltaT = creep.deltaT_min;

        return;
    }




    void If_No_tip(wstring selectedFile)
    {
        if (no == 0)
        {
            work0(0);   
            if (irock == 1)
            {
                initiation();      //if rock strength is given
                geoplot();
            }
            if (no != 0)
            {
               // geoplot();
                //------ - Initialising the seed elements for growth---------
                for (ni = 0; ni < no; ni++)
                {
                    if (tips[ni].ityp == 0) break;
                    nelement = tips[ni].mpointer;
                  
                }
                return;
            }
            if (lastinput == "endf")
            {
                
                return;  // stop
            }
            else
                input(); //first time to call input()  
            return;

        }
            //------ - Initialising the seed elements for growth----
            for (ni = 0; ni < no; ni++)
            {
                if (tips[ni].ityp == 0) break;
                nelement = tips[ni].mpointer;
              
            }

            return;
    }
      




    void add_crack_growth()
    {       
        /* add new elements to simulate crack growth */

        ktipgrow = false;
        for ( ni = 0; ni < no; ni++)
        {
            if (tips[ni].ityp == 0) continue;
            if (tips[ni].ifail != 1) continue;            
            nelement = tips[ni].mpointer;              
            newtips(tips[ni].angl);
            creep.creep_x[ni] = 0;         //reset the creep growth to zero
            creep.creep_y[ni] = 0;          //reset the creep growth to zero
            ktipgrow = true;            
        }

        arrangetip();
        if (irock == 1)       
        {                 
            initiation();
            return;
        }
        if (creep.ID_creep == 1 && creep.time >= creep.totalT)        //creep problem,id_creep=ID_creep
        {
            geoplot();
            MessageBox(nullptr, L"Defined creep time is completed, continue cycle without creep", L"Message!", MB_OK);
            creep.ID_creep = 0;
            return;
        }

        else if (creep.ID_creep == 0 && ktipgrow == false)        //non - creep problem
        {
            geoplot();
            if (lastinput != "endf")
            {
                MessageBox(nullptr, L"No more fracture propogation, continue from input file", L"Error!", MB_OK);
                input();
            }
            else if (lastinput == "endf")
            {
               //geoplot();
                StopReturn = false;
                return;                    
               
            }
        }   

        return;
    }




void input_tip_check()
    {

        float xt = 0, yt = 0, x1 = 0, y1 = 0, x2 = 0, y2 = 0;

        for (int n = 0; n < no; ++n)
        {
            Tip& t = tips[n];       //t is alias 
            if (t.ityp == +1 || t.ityp == 3 || t.ityp == 4)
            {
                xt = t.xbe;
                yt = t.ybe;
            }
            else if (t.ityp == -1 || t.ityp == -3)
            {
                xt = t.xen;
                yt = t.yen;
            }

            for (int m = 0; m < numbe; ++m)
            {
                BoundaryElement& b = elm_list[m];
                if (m == t.mpointer) continue;      // Skip if m == mpointer[n]
                x1 = b.xm - b.a * b.cosbet;
                y1 = b.ym - b.a * b.sinbet;
                x2 = b.xm + b.a * b.cosbet;
                y2 = b.ym + b.a * b.sinbet;

                //if (std::fabs(xt - x1) <= 1e-5 && std::fabs(yt - y1) <= 1e-5)
                
                
                if (std::fabsf(static_cast<float>(xt - x2)) <= static_cast<float>(1e-5) &&
                    std::fabsf(static_cast<float>(yt - y2)) <= static_cast<float>(1e-5)) {
                    t.ityp = 0;
                }
                //if (std::fabs(xt - x2) <= 1e-5 && std::fabs(yt - y2) <= 1e-5) 
                
                if (std::fabsf(xt - x2) <= 1e-5f && std::fabsf(yt - y2) <= 1e-5f) {
                    t.ityp = 0;
                }
            }
        }
        arrangetip();
        return;
    }


