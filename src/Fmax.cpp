#include<stdafx.h>

#include <Fmax.h>
#include "CommonPara.h"
#include "Work.h"
#include <Failure.h>
#include <WinInterface.h>
#include <DX.h>


using namespace WinInterface_h::winvar;
using namespace CommonPara_h::comvar;






int check_elastic_growth(int m)
{

    int kk = 0; // = 0 is not checking elastic growth
    float aki = 0, akii = 0, akie = 0;
    int mm = elm_list[m].mat_no;
    if (multi_region)
            mm = check_material_id(elm_list[m].xm, elm_list[m].ym);
    

    int ms = 2 * m; //Sara because of index error in d0[2m-1 in line64 add this indexs,make sure it works correctly
    int mn = ms + 1;

    if (s4.d0[mn] > 0)
    {
        aki = 0;
    }
    else
    {
        aki = abs(rock1[mm].e / (8 * pi * (1 - (rock1[mm].pr * rock1[mm].pr))) *
            sqrt(2 * pi / elm_list[m].a) * s4.d0[mn]) * 2.5; // 2.5 is a correction factor
    }

    akii = abs(rock1[mm].e / (8 * pi * (1 - (rock1[mm].pr * rock1[mm].pr))) *
        sqrt(2 * pi / elm_list[m].a) * s4.d0[ms]) * 2.5;

    if (akii == 0)       
    {
        akii = 1e-9;
    }


    float temp = aki / akii;
    float seta1 = atanf(0.25 * (temp + sqrt(temp * temp + 8)));
    float seta2 = atanf(0.25 * (temp - sqrt(temp * temp + 8)));

    float temp2 = seta1 / 2;  //Sara do this for sta2
    float akie1 = aki * pow(cosf(temp2), 3) - 3 * akii * pow(cosf(temp2), 2) * sinf(temp2);
    float akie2 = aki * pow(cosf(seta2 / 2), 3) - 3 * akii * pow(cosf(seta2 / 2), 2) *
        sinf(seta2 / 2);

    akie = max(akie1, akie2);

    if (akie > factors.factor_e * rock1[mm].akic || akii > factors.factor_e * rock1[mm].akiic)
    {
        kk = 1;             // 1 is checking elastic growth
    }

    return kk;
}






void newcoord()
{    

    /*compute the coordinates of gost element*/

    float sigxx = 0, sigyy = 0, sigxy = 0;  
    float xd = 0, yd = 0, sw = 0, cosb = 0, sinb = 0, ss = 0, sn = 0,
        pxx = 0, pyy = 0, pxy = 0, y0 = 0;
    float dr = 0;

    float xb = 0.0f, yb = 0.0f, xe = 0.0f, ye = 0.0f;   
    Tip& t = tips[ni];  //alias
  
    xb = t.xbe;
    yb = t.ybe;
    xe = t.xen;
    ye = t.yen;

    if (t.ityp == -1 || t.ityp == -3)
    {
        xb = xe - t.dl * (cosf(dr) * t.costem - sinf(dr) * t.sintem);
        yb = ye - t.dl * (sinf(dr) * t.costem + cosf(dr) * t.sintem);
    }
    else
    {
        if ( t.ityp == +1 || t.ityp == +3 || t.ityp == 4 )
        {
            xe = xb + t.dl * (cosf(dr) * t.costem - sinf(dr) * t.sintem);
            ye = yb + t.dl * (sinf(dr) * t.costem + cosf(dr) * t.sintem);

        }
    }
  
    xd = xe - xb;
    yd = ye - yb;
    sw = sqrt(xd*xd + yd*yd);
    int mm = t.mat_no;

    int m = numbe -1;     // index of new element(ghost element)
    elm_list[m].xm = xb + 0.5 * xd;
    elm_list[m].ym = yb + 0.5 * yd;
    elm_list[m].a = 0.5 * sw;
    elm_list[m].sinbet = yd / sw;
    elm_list[m].cosbet = xd / sw;
    elm_list[m].kod = 5;

    elm_list[m].mat_no = mm;
    cosb = elm_list[m].cosbet;
    sinb = elm_list[m].sinbet;

    NewFractureCentralPoint(elm_list[m].xm, elm_list[m].ym, sigxx, sigyy, sigxy);

    ss = (sigyy - sigxx) * sinb * cosb + sigxy * (cosb * cosb - sinb * sinb);
    sn = sigxx * sinb * sinb - 2.0 * sigxy * sinb * cosb + sigyy * cosb * cosb;

    if (elm_list[m].mat_no == mat_lining)
    {
        s4.b0[m * 2] = 0;
        s4.b0[m * 2 + 1] = 0;
    }
    else 
    {
        s4.b0[m * 2] = -ss;
        s4.b0[m * 2 + 1] = -sn;
    }

    //increment from crack growth, the additional stress applied should be the diff ss and sn* tan(phi)
    if ((sn < 0) && (abs(ss) > (-sn * tanf(30.0 / 180.0 * pi))))
    {
        s4.b0[m * 2] = -(abs(ss) - (-sn * tanf(0.524))) * std::copysign(1, ss);
        s4.b0[m * 2 + 1] = 0;
    }

    y0 = g.y_surf;
    pxx = symm.pxx1 + g.skx * (y0 - elm_list[m].ym);
    pyy = symm.pyy1 + g.sky * (y0 - elm_list[m].ym);
    pxy = symm.pxy1;

    b_elm[m].force1 = 2.0 * elm_list[m].a * ( -(pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb) );
    b_elm[m].force2 = 2.0 * elm_list[m].a * ( -(pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb) );

    return;
}





float call_work1_setting_fi0(float dtt, float& fi0, int mm, int mode, float& angip, float ang)
{
    float  gi0 = 0.0, fi = 0.0;
    float wi = 0.0;
    if (mode == 1)
    {
        work1(1); // 0 = not first time; 1 = mode I    2 = mode II
        //cout << "tip: "<<ni<<"  mode1:" << ang << "  w1:" << w1<<endl;

        wi = (b_elm[numbe-1].jstate != 1) ? w0 : w1;       // if tip element not openp, ignore it
         if (dtt == 0)
        {
             MessageBox(nullptr,L"dtt equals 0",L"Error", MB_OK);
        }
        else
        {
            gi0 = abs(wi - w0) / dtt;
            fi = gi0 / (pow(rock1[mm].akic, 2) * (1 - pow(rock1[mm].pr, 2)) / rock1[mm].e);
        }
    }
    else
    {
        work1(2); // 0 = not first time; 1 = mode I    2 = mode II
        //cout << "tip: " << ni << "mode2: " << ang << "  w1:" << w1 << endl;

        wi = (b_elm[numbe - 1].jstate != 2) ? w0 : w1;    // if tip element not openp, ignore it
        if (dtt == 0)
        {
           MessageBox(nullptr, L"dtt equals 0", L"Error", MB_OK);
        }
        else
        {
            gi0 = abs(wi - w0) / dtt;
            fi = gi0 / (pow(rock1[mm].akiic, 2) * (1 - pow(rock1[mm].pr, 2)) / rock1[mm].e);
        }
    }
    if (fi < 0) fi = 0;
    if (fi >= fi0)
    {
        angip = ang * 180.0/ pi; 
        fi0 = fi;
    }
    
    return wi;
}






float ang_setting(float angi0, float& fi0, int mm, int mode, float& angi)
{
    float dtt = 0.0;
    
    float ang = 0.0;
    bool kpause = false;
    float wi = 0.0;
    float s_index;
    float e_index;
      
    angi = angi0;
    //0.000009 only till cycl13
    //0.00001 to cycl13  

    //0.000001works till test 11
    if(float(int(angi0)) != round(angi0))
    {
        if (abs(angi0 - 8 - round(angi0 - 8)) > 0.0000019)   //  //0.000001
            s_index = int(angi0 - 8);
        else
            s_index = round(angi0 - 8);

        if (abs(angi0 + 8 - round(angi0 + 8)) > 0.0000019)
            e_index = int(angi0 + 8);
        else
            e_index = round(angi0 + 8);
    }
    else 
    {

        s_index = int(angi0 - 8);
        e_index = int(angi0 + 8);
    }    


    for (int icyc = s_index; icyc <= e_index; icyc += 2)
    {
        ang = icyc * pi / 180.0;

        
        if (tips[ni].ityp == 4)
        {
            dtt = dxi(ang, true); // first fracture initiation from boundary
        }
        else
        {
            dtt = dxi(ang, false);
        }

        wi = call_work1_setting_fi0(dtt, fi0, mm, mode, angi, ang);
    }
  
    return wi;
}





void compute_f(float& f0, float fi0, float angi, float angii, float& angle, float fii0, int m, float wi, float wii)
{    
    float dtt = 0.0;
    if (fi0 >= fii0)
    {
        f0 = fi0;
        angle = angi;
        tips[ni].imode = 1;
    }
    else
    {
        f0 = fii0;
        angle = angii;
        tips[ni].imode = 2;
    }

    angle = angle * pi / 180.0;

    if (tips[ni].ityp == 4)
    {
        dtt = dxi(angle,true);   //dxiB 
    }
    else
    {
        dtt = dxi(angle,false);   //dxi
    }

    if ((elm_list[m].xm - symm.xsym) <= elm_list[m].a / 1000 && abs(elm_list[m].sinbet) >
        sinf(85 * pi / 180) && (symm.ksym == 1 || symm.ksym == 4)  ||
        (elm_list[m].ym - symm.ysym) <= elm_list[m].a / 1000 && abs(elm_list[m].cosbet) > 
        cosf(5 * pi / 180) && (symm.ksym == 2 || symm.ksym == 4 ))
    {
        f0 = 2.0 * f0;
    }

    float xc = 0.5 * (tips[ni].xbe + tips[ni].xen);
    float yc = 0.5 * (tips[ni].ybe + tips[ni].yen);
    
    file2 << std::setw(5) << ni+1 << std::fixed << std::setprecision(4)
        << std::setw(8) << xc << std::setw(8) << yc << std::setprecision(3)
        << std::setw(12) << angle * 180. / pi << std::scientific << std::setprecision(8)
        << std::setw(18) << w0 << std::setw(18) << wi << std::setw(18) << wii
        << std::scientific << std::setprecision(3) << std::setw(14) << f0 << std::endl;
        
    return;
}






void fmax1(float& f0, float& angle)
{
    /*determine the fmax and angel using the F - criterion */ 
    
    float fi0 = 0., fii0 = 0., angi0 = 0., angii0 = 0.;
    int message = 0, itooclose = 0;
    angle = 0;
    f0 = 0.;

    int m = tips[ni].mpointer;
    int mm = elm_list[m].mat_no;     //mat_no instead of mat
    int kk = 0;

    float sn =0;
    float dtt = 0;
    float angi = 0, angii = 0, ang = 0;     
    float wii = 0, wi = 0;

    BoundaryElement& be = elm_list[m];
    bool kpause = false;

    if (b_elm[m].jstate == 0)
    {
        kk = check_elastic_growth(m);
        if (kk == 0) return;                
    }

    numbe++;
    newcoord();    
    bool flag0 = (be.xm == symm.xsym && (symm.ksym == 1 || symm.ksym == 4) ||
        be.ym == symm.ysym && (symm.ksym == 2 || symm.ksym == 4));
    bool flag1 = b_elm[m].jstate == 2 && b_elm[m].jslipd == 1;
    bool flag2 = b_elm[m].jstate == 2 && b_elm[m].jslipd == -1;
    bool flag3 = (be.xm - symm.xsym) <= be.a / 1000 && abs(be.sinbet) > sinf(85 * pi / 180)
        && (symm.ksym == 1 || symm.ksym == 4);
    bool flag4 = (be.ym - symm.ysym) <= be.a / 1000 && abs(be.cosbet) > cosf(5 * pi / 180)
        && (symm.ksym == 2 || symm.ksym == 4);
    //float factor = pi / 180.0;   

    for (int icyc = 100; icyc >= -100; icyc -= 10)
    {
        ang = icyc * pi / 180.0;

        if (flag0)
        {
            ang = 0;
        }

        if (tips[ni].ityp == 4)
        {
            dtt = dxi(ang,true);  // first fracture initiation from boundary,dxiB flagB = true
            //cout <<" icyc: "<<icyc<< " ang: " << ang << "  dtt: " << dtt << endl;
        }
        else
        {
            dtt = dxi(ang,false);      //call dxi flagB = false
           // cout << " icyc: " << icyc << " ang: " << ang << "  dtt: " << dtt << endl;

        }

        if (!((flag1 && icyc < -20)||(flag2 && icyc > 20)))
        {
            wi = call_work1_setting_fi0(dtt, fi0, mm, 1, angi0, ang);
        }
       
        if (icyc < -50 || icyc > 50) continue;// next round of loop label400
        
        wii = call_work1_setting_fi0(dtt, fii0, mm, 2, angii0, ang);
                
        if (flag3 || flag4)
        {
            compute_f(f0, fi0, angi, angii, angle, fii0, m, wi, wii);
            numbe--;
            return;
        }
        
    } //label400
    //--- ---------------------------------------------------------

      wi = ang_setting(angi0, fi0, mm, 1, angi);

	  wii = ang_setting(angii0, fii0, mm, 2, angii);
	  compute_f(f0, fi0, angi, angii, angle, fii0, m , wi, wii);  //label600
      numbe--;
     return;
}












