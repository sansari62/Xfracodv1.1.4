#include<stdafx.h>
#include <Failure.h>
#include "CommonPara.h"
#include<DX.h>
#include<Tip.h>

using namespace CommonPara_h::comvar;




int CheckNewElement(float ac, float xc, float yc, float cosbeta, float sinbeta, int flagB)
{
    //if flagB = 1 then it works like checkNewElementB in the original code
    //one more condition added to loop, maybe need to optimized

    int legal = 1;
    float x1 = xc - ac * cosbeta;
    float y1 = yc - ac * sinbeta;
    float x2 = xc + ac * cosbeta;
    float y2 = yc + ac * sinbeta;
    int  en;
    float thr = 0.0003; //6e-4;
    //const float epsilon = 1e-6;

    if (flagB == 0)
        en = numbe - 1;
    else
        en = numbe;

    for (int i = 0; i < en; ++i)
    {
        BoundaryElement& be = elm_list[i];
        if (flagB == 1)
        {
            if (be.kod != 5) continue;
        }

        float dist1 = sqrt(powf(x1 - be.xm, 2) + powf(y1 - be.ym, 2));
        float dist2 = sqrt(powf(x2 - be.xm, 2) + powf(y2 - be.ym, 2));


        if (min(dist1, dist2) - thr <= factors.tolerance * max(ac, be.a))
        {

            return(0);
        }

        if (symm.ksym == 1 || symm.ksym == 4)
        {
            dist1 = sqrt(powf(x1 - (2.0 * symm.xsym - be.xm), 2) +
                powf(y1 - be.ym, 2));
            dist2 = sqrtf(powf(x2 - (2.0 * symm.xsym - be.xm), 2) +
                powf(y2 - be.ym, 2));
            if (min(dist1, dist2) <= factors.tolerance * max(ac, be.a))
            {

                return(0);
            }

        }

        if (symm.ksym == 2 || symm.ksym == 4)
        {
            dist1 = sqrt(powf(x1 - be.xm, 2) + powf(y1 - (2.0 * symm.ysym - be.ym), 2));
            dist2 = sqrt(pow(x2 - be.xm, 2) + powf(y2 - (2.0 * symm.ysym - be.ym), 2));
            if (min(dist1, dist2) <= factors.tolerance * max(ac, be.a))
            {

                return(0);
            }
        }

        if (symm.ksym == 3 || symm.ksym == 4)
        {
            dist1 = sqrt(powf(x1 - (2.0 * symm.xsym - be.xm), 2) + powf(y1 - (2.0 * symm.ysym - be.ym), 2));
            dist2 = sqrt(powf(x2 - (2.0 * symm.xsym - be.xm), 2) + powf(y2 - (2.0 * symm.ysym - be.ym), 2));
            if (min(dist1, dist2) <= factors.tolerance * max(ac, be.a))
            {

                return(0);
            }
        }
    }

    if (symm.ksym == 1 || symm.ksym == 4)
    {
        if ((x1 - symm.xsym) * (x2 - symm.xsym) < -s5u.dtol)
        {

            return(0);
        }
    }
    if (symm.ksym == 2 || symm.ksym == 4)
    {
        if ((y1 - symm.ysym) * (y2 - symm.ysym) < -s5u.dtol)
        {
            //legal = 0;
            return(0);
        }
    }
    if (symm.ksym == 3 || symm.ksym == 4)
    {
        if ((x1 - symm.xsym) / (y1 - symm.ysym) == -(x2 - symm.xsym) / (y2 - symm.ysym))
        {
            //legal = 0;
            return(0);
        }
    }

    return legal;
}







int  check_material_id(float xp, float yp)
{
    float dist0 = 10e8;
    int mm = 1, mclosest = 0;
    int k = 0;
    float xc = 0, yc = 0, xt = 0, yt = 0, dist = 0;
    float ypprime, xpprime;
    //remove prenumbe from for loop because of the issue with multi region problem

    for (int m = 0; m < numbe; ++m)
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






void  Sum_Failure(float xp, float yp, float r, float alpha, int im, float fos, int location)
{
    /* This function reassign value to the failure points
        mf - numbe of identified failure points
    */
    int index = mf;
    init_point[index].xf = xp;
    init_point[index].yf = yp;
    init_point[index].rf = r;
    init_point[index].alphaf = alpha;
    init_point[index].imf = im;
    init_point[index].fosf = fos;
    init_point[index].Locationf = location;
    mf += 1;
    return;
}






void setting_elem_and_tipn_failure(int m, int im)
{
    float sigxx = 0.0, sigyy = 0.0, sigxy = 0.0;
    float sinb = elm_list[m].sinbet;
    float cosb = elm_list[m].cosbet;

    if (elm_list[m].mat_no == mat_lining)
    {
        s4.b0[m * 2] = 0;
        s4.b0[m * 2 +1] = 0;
    }
    else
    {
        NewFractureCentralPoint(elm_list[m].xm, elm_list[m].ym, sigxx, sigyy, sigxy);
        float ss = (sigyy - sigxx) * sinb * cosb + sigxy * (cosb * cosb - sinb * sinb);
        float sn = sigxx * sinb * sinb - 2.0 * sigxy * sinb * cosb + sigyy * cosb * cosb;
        s4.b0[m * 2] = -ss;
        s4.b0[m * 2+1] = -sn;
        if (sn < 0 && abs(ss) > (- sn * tanf(30.0 / 180.0 * pi)))
        {
            // additional stress applied to crack growth element
            s4.b0[m * 2] = -(abs(ss) - (-sn * tanf(0.524))) * copysign(1, ss);
            s4.b0[m * 2+1] = 0;
        }
    }

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - elm_list[m].ym);
    float pyy = symm.pyy1 + g.sky * (y0 - elm_list[m].ym);
    float pxy = symm.pxy1;

    b_elm[m].force1 = 2.0 * elm_list[m].a * (-((pyy - pxx) * sinb * cosb + pxy *
        (cosb * cosb - sinb * sinb)));           // old b0()
    b_elm[m].force2 = 2.0 * elm_list[m].a * (-(pxx * sinb * sinb - 2.0 * pxy *
        sinb * cosb + pyy * cosb * cosb));        // old b0()

    int mm = 1;  // index of rock
    float aks_bb = 0, akn_bb = 0, phi_bb = 0, coh_bb = 0, phid_bb = 0, ap_bb = 0,
        apr_bb = 0;
    stiffness_bb(aks_bb, akn_bb, phi_bb, coh_bb, phid_bb, ap_bb, apr_bb, im, mm);
    b_elm[m].aks = aks_bb;
    b_elm[m].akn = akn_bb;
    b_elm[m].phi = phi_bb;
    b_elm[m].coh = coh_bb;
    b_elm[m].phid = phid_bb;

    joint[m].aperture0 = ap_bb;
    joint[m].aperture_r = apr_bb;
    b_elm[m].jstate = 0;
    b_elm[m].jmode = im;

    return;
}




void failure(float xp, float yp, float r, float alpha, int im)
{
    //IM = 1, tensile, IM = 2, shear

    int m = 0;   
    float xb, yb, xn, yn;
    int legal, material;      
    if (numbe >= m0-1)
    {
        MessageBox(nullptr, L"Maximum BE limit exceeded!", L"Message!", MB_OK);
        exit(0);        
        return;   
    }
    int n = no;
    xb = xp - 1 / 2.0 * r * cosf(alpha);     // xbe(?) etc are potential (additional) element
    yb = yp - 1 / 2.0 * r * sinf(alpha);
    xn = xp;
    yn = yp;    
    int  material = (multi_region)?check_material_id(0.5 * (xb + xn), 0.5 * (yb + yn)): j_material;
    float dl = sqrtf(std::powf(xn - xb,2) + std::powf(yn - yb,2));
    tips[n].assign_val(xb, yb, xn, yn, dl, cosf(alpha), sinf(alpha), -1, material);  
    no++;    
    //-------------------------------add new element------------
    m = numbe;  
    elm_list[m].xm= 0.5 * (xb + xn);
    elm_list[m].ym= 0.5 * (yb + yn);
    elm_list[m].a = dl / 2.0;
    elm_list[m].sinbet= (yn - yb) / dl;
    elm_list[m].cosbet= (xn - xb) / dl;
    elm_list[m].kod = 5;
    elm_list[m].mat_no = material;

    legal = CheckNewElement(elm_list[m].a, elm_list[m].xm, elm_list[m].ym,
        elm_list[m].cosbet, elm_list[m].sinbet, 0);    
    m += 1;
    if (legal == 0)
    {
        elm_list[m-1].mat_no = 0;
        tips[no-1].mat_no = 0;
        m -= 1;
        no -= 1;       
    }
    else
    {
        setting_elem_and_tipn_failure( m-1,  im);        
        tips[n].mpointer = m-1;
        tips[n].xen = xp - 1 / 2. * r * cosf(alpha);
        tips[n].yen = yp - 1 / 2. * r * sinf(alpha);
        tips[n].xbe = tips[n].xen - 1 / 2. * r * cosf(alpha);     
        tips[n].ybe = tips[n].yen - 1 / 2. * r * sinf(alpha);
        tips[n].imode = im;       
        // -----------------------------------------------
        // avoid repeating element at symmetry locations

        if (symm.ksym == 1 && abs(xp - symm.xsym) < s5u.dtol)
        {
            numbe = m;
            return;
        }
        if(symm.ksym == 2 && abs(yp - symm.ysym) < s5u.dtol)
        {
            numbe = m;
            return;
        }
        if(symm.ksym == 3 && abs(xp - symm.xsym) < s5u.dtol && abs(yp - symm.ysym) <
            s5u.dtol)
        {
            numbe = m;
            return;
        }
        //FFFF
        if(symm.ksym == 4 && (abs(xp - symm.xsym) < s5u.dtol || abs(yp - symm.ysym) <
            s5u.dtol))           
        {
            numbe = m;
            return;
        }
        //------------Add new tip----------------
        n = no;
        xb = xp;   
        yb = yp;
        xn = xp + 1 / 2.0 * r * cosf(alpha); 
        yn = yp + 1 / 2.0 * r * sinf(alpha);
        if(multi_region)
             material = check_material_id(0.5 * (xb + xn), 0.5 * (yb + yn));  
        dl = std::sqrt(std::pow(xn - xb,2) + std::pow(yn - yb, 2));
        tips[n].assign_val(xb, yb, xn, yn, dl, cosf(alpha), sinf(alpha), 1, material);
        no++;     //new tip added
        
        //-------------Add new element------------------
        float sinb = (yn - yb) / dl;
        float cosb = (xn - xb) / dl;       
        elm_list[m].xm = 0.5 * (xb + xn);
        elm_list[m].ym = 0.5 * (yb + yn);
        elm_list[m].a = dl / 2.0;
        elm_list[m].sinbet = sinb;    
        elm_list[m].cosbet = cosb;       
        elm_list[m].kod = 5;
        elm_list[m].mat_no = material;

        float a = elm_list[m].a; float xm = elm_list[m].xm; float ym = elm_list[m].ym;
        legal = CheckNewElement(a,xm, ym, cosb, sinb, 0);
        m += 1;      //new element added
        if (legal == 0)
        {
            elm_list[m-1].mat_no = 0;
            tips[no-1].mat_no = 0;
            m -= 2;
            no -= 2;
        }
        else
        {
            setting_elem_and_tipn_failure(m-1, im);
            tips[n].mpointer = m-1;
            tips[n].xbe = tips[n].xen;     
            tips[n].ybe = tips[n].yen;
            tips[n].xen = tips[n].xbe + 1 / 2. * r * cosf(alpha);
            tips[n].yen = tips[n].ybe + 1 / 2. * r * sinf(alpha);
            tips[n].imode = im;
        }
    }
    numbe = m;    
    return;
}





void failureB(float xp, float yp, float r, float alpha, int im)
{
    //IM = 1, tensile, IM = 2, shear
   
    int m = numbe;    
    float xb, yb, xe, ye;    
    if (numbe >= m0-1)
    {
        MessageBox(nullptr, L"Maximum BE limit exceeded!", L"Message!", MB_OK);
        exit(0);          
        return;  
    }
    int n = no;
    float sinb = sinf(alpha);
    float cosb = cosf(alpha);

    xb = xp + 0.5 * r * cosb;     
    yb = yp + 0.5 * r * sinb;
    xe = xp + 1.0 * r * cosb;
    ye = yp + 1.0 * r * sinb;    
    
    int material = (multi_region)? check_material_id(0.5 * (xb + xe), 0.5 * (yb + ye)):
        j_material;
    float dl = sqrt(pow(xe - xb,2) + pow(ye - yb,2));
    tips[n].assign_val(xb, yb, xe, ye, dl, cosb, sinb, 4, material);
    no++; 
    //---------------Add new element-------------------------------
    float x = 0.5 * (xb + xp);
    float y = 0.5 * (yb + yp);    
    elm_list[m].xm= 0.5 * (xb + xp);
    elm_list[m].ym = 0.5 * (yb + yp);
    elm_list[m].a = dl / 2;
    elm_list[m].sinbet = sinb;
    elm_list[m].cosbet = cosb;
    elm_list[m].kod = 5;
    elm_list[m].mat_no = material;
     
    int legal = CheckNewElement(elm_list[m].a, elm_list[m].xm, elm_list[m].ym, 
        elm_list[m].cosbet, elm_list[m].sinbet, 1);
      
    m = numbe + 1;     
    if (legal == 0)
    {
        m -= 1;
        no -= 1;
    }
    else
    {
        setting_elem_and_tipn_failure(m-1, im);
        // -----------------------------------------------       
        if (multi_region)
            material = check_material_id(0.5 * (xb + xe), 0.5 * (yb + ye));
             
        elm_list[m].xm = 0.5 * (xb + xe);
        elm_list[m].ym = 0.5 * (yb + ye);
        elm_list[m].a = dl / 2;
        elm_list[m].sinbet = sinb;
        elm_list[m].cosbet = cosb;
        elm_list[m].kod = 5;
        elm_list[m].mat_no = material;

        legal = CheckNewElement(elm_list[m].a, elm_list[m].xm, elm_list[m].ym, 
            elm_list[m].cosbet, elm_list[m].sinbet, 1);
        m = numbe + 2;
        if (legal == 0)
        {
            elm_list[m-1].mat_no = 0;
            tips[no-1].mat_no = 0;
            m -= 2;
            no -= 1;
        }
        else
        {
            setting_elem_and_tipn_failure(m-1, im);            
            tips[n].mpointer = m-1;
            ktipgrow = true;           
        }
    }
     tips[n].xbe = xe;
     tips[n].ybe = ye;
     tips[n].xen = tips[n].xbe + 0.5 * r * cosf(alpha);  
     tips[n].yen = tips[n].ybe + 0.5 * r * sinf(alpha);
     tips[n].imode = im;
     numbe = m;     
    return;
}






void Choose_Failure()
{
    /* Choose the points with FoS > 0.9*FoS_max as the failure points*/


    float fos_max = 0, fos = 0;
    int legal = 0;
    float xp, yp, r, alpha;
    int im;

    // mf - numbe of identified failure points
    for (int np = 0; np < mf; ++np)
    {        
        if (init_point[np].fosf > fos_max)
        {
            xp = init_point[np].xf;
            yp = init_point[np].yf;
            r = init_point[np].rf;
            alpha = init_point[np].alphaf;

            if (init_point[np].Locationf == 1)
            {
                legal = CheckNewElement(r, xp, yp, cosf(alpha), sinf(alpha), 1);
            }
            else //if (init_point[np].Locationf == 2)
            {
                legal = CheckNewElement(r, xp, yp, cosf(alpha), sinf(alpha), 0);
            }
            if (legal == 1)
                    fos_max = init_point[np].fosf;
        }
    }

    for (int np = 0; np < mf; ++np)
    {
        if (init_point[np].fosf> factors.factor_f * fos_max)
        {
            xp = init_point[np].xf;
            yp = init_point[np].yf;
            r = init_point[np].rf;

            alpha = init_point[np].alphaf;
            im = init_point[np].imf;
            fos = init_point[np].fosf;

            if (init_point[np].Locationf == 1)
            {
                legal = CheckNewElement( r, xp, yp, cosf(alpha), sinf(alpha),1);
                if (legal == 1)
                    failureB(xp, yp, r, alpha, im);
            }
            else
            {
                legal = CheckNewElement( r, xp, yp, cosf(alpha), sinf(alpha),0);
                if (legal == 1)
                    failure(xp, yp, r, alpha, im);
            }
        }
    }
    return;
}




