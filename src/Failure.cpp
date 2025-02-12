#include <Failure.h>
#include "CommonPara.h"
#include<ExcavationCracks.h>
#include <Mainb.h>
#include<Source.h>



using namespace CommonPara_h::comvar;




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



void NewFractureCentralPoint(float xp, float yp, float& sigxx, float& sigyy, float& sigxy)
{
    /* obtain the centre point of a new fracture and add it into the matrix */
      
    float  y0 = 0.0;
    int mm = check_material_id(xp, yp);   

    s2us.reset();
    y0 = g.y_surf;

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

    for (int j = 0; j < numbe_old; ++j)
    {
        BoundaryElement& be = elm_list[j];
        if (mm != be.mat_no)
            continue;

        int js = 2 * j;
        int jn = js + 1;
        s2us.reset();
        float xj = be.xm;
        float yj = be.ym;
        float aj = be.a;

        float cosbj = be.cosbet;
        float sinbj = be.sinbet;

        coeff(xp, yp, xj, yj, aj, cosbj, sinbj, +1, mm);

        switch (symm.ksym + 1) {
        case 1:
            break;

        case 2:
            xj = 2.0 * symm.xsym - be.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            break;

        case 3:
            yj = 2.0 * symm.ysym - be.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            break;

        case 4:
            xj = 2.0 * symm.xsym - be.xm;
            yj = 2.0 * symm.ysym - be.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, -sinbj, +1, mm);
            break;

        case 5:
            xj = 2.0 * symm.xsym - be.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            xj = be.xm;
            yj = 2.0 * symm.ysym - be.ym;

            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            xj = 2.0 * symm.xsym - be.xm;
            coeff(xp, yp, xj, yj, aj, -cosbj, -sinbj, +1, mm);
           
        }
        
        sigxx += s2us.sxxs * s4.d0[js] + s2us.sxxn * s4.d0[jn];
        sigyy += s2us.syys * s4.d0[js] + s2us.syyn * s4.d0[jn];
        sigxy += s2us.sxys * s4.d0[js] + s2us.sxyn * s4.d0[jn];
    }
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




void failure(float xp, float yp, float r, float alpha, int im, float& fos)
{
    //IM = 1, tensile, IM = 2, shear

    int m = 0;
   
    float xb, yb, xn, yn;
    int legal, material;
    if (fos > 1.0) 
    {
        fos = 1.0;
    }
    
    if (numbe >= 1498)
    {
        MessageBox(nullptr, L"Memory Overflow - Reduce In Situ Stresses.", L"Message!", MB_OK);
        exit(0);        
        return;   //instead of stop
    }
    int n = no;
    xb = xp - 1 / 2.0 * r * cosf(alpha);     // xbe(?) etc are potential (additional) element, not actual  instead of 1/2r, dl used here
    yb = yp - 1 / 2.0 * r * sinf(alpha);
    xn = xp;
    yn = yp;

    material = check_material_id(0.5 * (xb + xn), 0.5 * (yb + yn));    
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

    //legal = CheckNewElement(elm_list[m].a, elm_list[m].xm, elm_list[m].ym,
     //   elm_list[m].cosbet, elm_list[m].sinbet, 0);

    float a = elm_list[m].a; float xm = elm_list[m].xm; float ym = elm_list[m].ym;
    legal = CheckNewElement(a, xm, ym, elm_list[m].cosbet, elm_list[m].sinbet, 0);

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





void failureB(float xp, float yp, float r, float alpha, int im, float& fos)
{
    //IM = 1, tensile, IM = 2, shear
   
    int m = numbe;
    
    float xb, yb, xe, ye;

    if (fos > 1.0) 
    {
        fos = 1.0;
    }
    else
        if (fos == 0)
        {
            //WRITE(1,*) 0000          write into file 
        }

    if (numbe >= 1498)
    {
       MessageBox(nullptr, L"Memory Overflow - Reduce In Situ Stresses.", L"Message!", MB_OK);

        exit(0);          
        return;   //instead of stop
    } 

    int n = no;
    float sinb = sinf(alpha);
    float cosb = cosf(alpha);

    xb = xp + 0.5 * r * cosb;     
    yb = yp + 0.5 * r * sinb;
    xe = xp + 1.0 * r * cosb;
    ye = yp + 1.0 * r * sinb;
   
    int material = check_material_id(0.5 * (xb + xe), 0.5 * (yb + ye)); 
    float dl = sqrt((xe - xb) * (xe - xb) + (ye - yb) * (ye - yb));     
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
      
    m = numbe + 1;    //no of new elements 
    if (legal == 0)
    {
        m -= 1;
        no -= 1;
    }
    else
    {
        setting_elem_and_tipn_failure(m-1, im);
        // -----------------------------------------------       
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
                 legal = CheckNewElement(r, xp, yp, cosf(alpha), sinf(alpha), 1);

            else //if (init_point[np].Locationf == 2)
                    legal = CheckNewElement(r, xp, yp, cosf(alpha), sinf(alpha), 0);
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
                    failureB(xp, yp, r, alpha, im, fos);
            }
            else //if (init_point[np].Locationf == 2)
            {
                legal = CheckNewElement( r, xp, yp, cosf(alpha), sinf(alpha),0);
                if (legal == 1)
                    failure(xp, yp, r, alpha, im, fos);
            }
        }
    }
    return;
}




