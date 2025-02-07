#include <DX.h>
#include <CommonPara.h>
#include <Source.h>
#include <Failure.h>

using namespace comvar;





void check_all_tips(float& xt, float& yt, float xt0, float yt0,bool flagB)
{
    int n = ni;  //not sure Sara!
    int id = 0;
    float z = 0;

    for (int i = 0; i < numbe - 1; ++i)
    {
        int m = tips[ni].mpointer;
        if (i == m || (flagB == true && elm_list[i].kod != 5))   continue;     //the second condition applys to dxib

        float xbeg = elm_list[i].xm - elm_list[i].a * elm_list[i].cosbet;
        float ybeg = elm_list[i].ym - elm_list[i].a * elm_list[i].sinbet;
        float xend = elm_list[i].xm + elm_list[i].a * elm_list[i].cosbet;
        float yend = elm_list[i].ym + elm_list[i].a * elm_list[i].sinbet;
        float xc = elm_list[i].xm;
        float yc = elm_list[i].ym;

        float tol = factors.tolerance;
        float tol1 = factors.tolerance;

        float dc = min(sqrt((xt - xc) * (xt - xc) + (yt - yc) * (yt - yc)),
            sqrt((xt0 - xc) * (xt0 - xc) + (yt0 - yc) * (yt0 - yc)));
        if (dc <= tol1 * max(elm_list[i].a, 0))
        {
            float dbeg = sqrt((xt - xbeg) * (xt - xbeg) + (yt - ybeg) * (yt - ybeg));
            float dend = sqrt((xt - xend) * (xt - xend) + (yt - yend) * (yt - yend));
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
            return;
        }


        float dbeg = min(sqrt((xt - xbeg) * (xt - xbeg) + (yt - ybeg) * (yt - ybeg)),
            sqrt((xt0 - xbeg) * (xt0 - xbeg) + (yt0 - ybeg) * (yt0 - ybeg)));

        float dend = min(sqrt((xt - xend) * (xt - xend) + (yt - yend) * (yt - yend)),
            sqrt((xt0 - xend) * (xt0 - xend) + (yt0 - yend) * (yt0 - yend)));

        if (dbeg <= dend && dbeg <= tol * max(elm_list[i].a, z))
        {
            xt = xbeg;
            yt = ybeg;
            return;
        }

        if (dend <= dbeg && dend <= tol * max(elm_list[i].a, z))
        {
            xt = xend;
            yt = yend;
            return;
        }
        id = cross(xbeg, ybeg, xend, yend, xt0, yt0, xt, yt);
        if (id == 1) return;

        //Sara! need to be optimized and similar to newtips
        if (symm.ksym != 0)
        {
            if (symm.ksym == 1 || symm.ksym == 4)
            {
                dc = min(std::sqrt((xt - (2. * symm.xsym - xc)) * (xt - (2. * symm.xsym - xc)) +
                    (yt - yc) * (yt - yc)),
                    std::sqrt((xt0 - (2. * symm.xsym - xc)) * (xt0 - (2. * symm.xsym - xc)) +
                        (yt0 - yc) * (yt0 - yc)));

                    if (dc <= tol1 * max(elm_list[i].a, z))
                    {
                        dbeg = sqrt((xt - (2. * symm.xsym - xbeg)) * (xt - (2. * symm.xsym - xbeg)) +
                            (yt - ybeg) * (yt - ybeg));
                        dend = sqrt((xt - (2. * symm.xsym - xend)) * (xt - (2. * symm.xsym - xend)) +
                            (yt - yend) * (yt - yend));

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
                        return;

                    }

                dbeg = min(sqrt((xt - (2. * symm.xsym - xbeg)) * (xt - (2. * symm.xsym - xbeg))
                    + (yt - ybeg) * (yt - ybeg)),
                    sqrt((xt0 - (2. * symm.xsym - xbeg)) * (xt0 - (2. * symm.xsym - xbeg)) +
                        (yt0 - ybeg) * (yt0 - ybeg)));
                dend = min(sqrt((xt - (2. * symm.xsym - xend)) * (xt - (2. * symm.xsym - xend)) +
                    (yt - yend) * (yt - yend)),
                    sqrt((xt0 - (2. * symm.xsym - xend)) * (xt0 - (2. * symm.xsym - xend)) +
                        (yt0 - yend) * (yt0 - yend)));

                if (dbeg <= dend && dbeg <= tol * max(elm_list[i].a, z))
                {
                    xt = 2. * symm.xsym - xbeg;
                    return;
                }

                if (dend <= dbeg && dend <= tol * max(elm_list[i].a, z))
                {
                    xt = 2. * symm.xsym - xend;
                    return;
                }

                id = cross(2.0 * symm.xsym - xbeg, ybeg, 2.0 * symm.xsym - xend, yend, xt0, yt0, xt, yt);

                if (id == 1)
                    return;
            }

            /////////////////////////////////////////////////////

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                dc = min(sqrt((xt - xc) * (xt - xc) + (yt - (2. * symm.ysym - yc)) * (yt - (2. * symm.ysym - yc))),
                    sqrt((xt0 - xc) * (xt0 - xc) + (yt0 - (2. * symm.ysym - yc)) * (yt0 - (2. * symm.ysym - yc))));

                if (dc <= tol1 * max(elm_list[i].a, 0 * tips[n].dl))
                {
                    dbeg = sqrt((xt - xbeg) * (xt - xbeg) + (yt - (2. * symm.ysym - ybeg)) * (yt - (2. * symm.ysym - ybeg)));
                    dend = sqrt((xt - xend) * (xt - xend) + (yt - (2. * symm.ysym - yend)) * (yt - (2. * symm.ysym - yend)));

                    if (dbeg <= dend) 
                    {
                        xt = xbeg;
                        yt = 2. * symm.ysym - ybeg;
                    }
                    else
                    {
                        xt = xend;
                        yt = 2. * symm.ysym - yend;
                    }
                    return;
                }

                dbeg = min(sqrt((xt - xbeg) * (xt - xbeg) +
                    (yt - (2. * symm.ysym - ybeg)) * (yt - (2. * symm.ysym - ybeg))),
                    sqrt((xt0 - xbeg) * (xt0 - xbeg) +
                        (yt0 - (2. * symm.ysym - ybeg)) * (yt0 - (2. * symm.ysym - ybeg))));

                dend = min(sqrt((xt - xend) * (xt - xend) +
                    (yt - (2. * symm.ysym - yend)) * (yt - (2. * symm.ysym - yend))),
                    sqrt((xt0 - xend) * (xt0 - xend) +
                        (yt0 - (2. * symm.ysym - yend)) * (yt0 - (2. * symm.ysym - yend))));

                if (dbeg <= dend && dbeg <= tol * max(elm_list[i].a, z)) {
                    yt = 2. * symm.ysym - ybeg;
                    return;
                }

                if (dend <= dbeg && dend <= tol * max(elm_list[i].a, z)) {
                    yt = 2. * symm.ysym - yend;
                    return;
                }

                id = cross(xbeg, 2. * symm.ysym - ybeg, xend, 2. * symm.ysym - yend, xt0, yt0, xt, yt);

                if (id == 1)
                    return;
            }
            //////////////////////////////

            if (symm.ksym == 3 || symm.ksym == 4) 
            {
                dc = min(sqrt((xt - (2. * symm.xsym - xc)) * (xt - (2. * symm.xsym - xc)) +
                    (yt - (2. * symm.ysym - yc)) * (yt - (2. * symm.ysym - yc))),
                    sqrt((xt0 - (2. * symm.xsym - xc)) * (xt0 - (2. * symm.xsym - xc)) +
                        (yt0 - (2. * symm.ysym - yc)) * (yt0 - (2. * symm.ysym - yc))));

                if (dc <= tol1 * max(elm_list[i].a, 0))
                {
                    dbeg = sqrt((xt - (2. * symm.xsym - xbeg)) * (xt - (2. * symm.xsym - xbeg)) +
                        (yt - (2. * symm.ysym - ybeg)) * (yt - (2. * symm.ysym - ybeg)));
                    dend = sqrt((xt - (2. * symm.xsym - xend)) * (xt - (2. * symm.xsym - xend)) +
                        (yt - (2. * symm.ysym - yend)) * (yt - (2. * symm.ysym - yend)));

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
                    return;

                }

                dbeg = min(sqrt((xt - (2. * symm.xsym - xbeg)) * (xt - (2. * symm.xsym - xbeg)) +
                    (yt - (2. * symm.ysym - ybeg)) * (yt - (2. * symm.ysym - ybeg))),
                    sqrt((xt0 - (2. * symm.xsym - xbeg)) * (xt0 - (2. * symm.xsym - xbeg)) +
                        (yt0 - (2. * symm.ysym - ybeg)) * (yt0 - (2. * symm.ysym - ybeg))));

                dend = min(sqrt((xt - (2. * symm.xsym - xend)) * (xt - (2. * symm.xsym - xend)) +
                    (yt - (2. * symm.ysym - yend)) * (yt - (2. * symm.ysym - yend))),
                    sqrt((xt0 - (2. * symm.xsym - xend)) * (xt0 - (2. * symm.xsym - xend)) +
                        (yt0 - (2. * symm.ysym - yend)) * (yt0 - (2. * symm.ysym - yend))));


                if (dbeg <= dend && dbeg <= tol * max(elm_list[i].a, z))
                {
                    xt = 2. * symm.xsym - xbeg;
                    yt = 2. * symm.ysym - ybeg;
                    return;
                }

                if (dend <= dbeg && dend <= tol * max(elm_list[i].a, z)) {
                    xt = 2. * symm.xsym - xend;
                    yt = 2. * symm.ysym - yend;
                    return;
                }

                id = cross(2. * symm.xsym - xbeg, 2. * symm.ysym - ybeg, 2. * symm.xsym - xend,
                    2. * symm.ysym - yend, xt0, yt0, xt, yt);

                if (id == 1)
                    return;
            }
        }
    }//end for label 300
    
    return;
}






float label_200(float xt, float yt, float xt0, float yt0, float& angle,bool flagB)
{

    float dtt = sqrt((xt - xt0) * (xt - xt0) + (yt - yt0) * (yt - yt0));
    if (dtt == 0) dtt = tips[ni].dl;
    if (!flagB)
    {
        if (tips[ni].ityp == 1 || tips[ni].ityp == 3)
        {
            angle = atan2f((yt - yt0), (xt - xt0))- 
                atan2f(tips[ni].sintem, tips[ni].costem);
        }
        else if (tips[ni].ityp == -1 || tips[ni].ityp == -3)
        {
            angle = atan2f((yt0 - yt), (xt0 - xt)) - 
                atan2f(tips[ni].sintem, tips[ni].costem);
        }
    }
    else
    {
        angle = atan2f((yt - yt0), (xt - xt0)) - 
            atan2f(tips[ni].sintem, tips[ni].costem);
    }  

   
    return dtt;
}



void recalculate_boundary_element_m(float xt, float yt, float xt0, float yt0,int mm, 
    float dtt)
{
    int m = numbe-1;  //instead of numbe
    float  sigxx = 0.0, sigyy = 0.0, sigxy = 0.0;   
    elm_list[m].mat_no = mm;
         
    elm_list[m].xm = 0.5 * (xt + xt0);
    elm_list[m].ym = 0.5 * (yt + yt0);
    elm_list[m].a = dtt / 2.0;
    elm_list[m].sinbet = (yt - yt0) / dtt;
    elm_list[m].cosbet = (xt - xt0) / dtt;
    float cosb = elm_list[m].cosbet;
    float sinb = elm_list[m].sinbet;

    if (elm_list[m].mat_no == mat_lining)
    {
        s4.b0[m * 2] = 0;
        s4.b0[m * 2+1] = 0;
    }
    else
    {
        NewFractureCentralPoint(elm_list[m].xm, elm_list[m].ym, sigxx, sigyy, sigxy);

        float ss = (sigyy - sigxx) * sinb * cosb + sigxy * (cosb * cosb - sinb * sinb);
        float sn = sigxx * sinb * sinb - 2.0 * sigxy * sinb * cosb + sigyy * cosb * cosb;
        s4.b0[m * 2] = -ss;
        s4.b0[m * 2 + 1] = -sn;
        if ( sn < 0 && abs(ss) > (-sn * tanf(30.0 / 180.0 * pi)) )    //!additional stress applied to crack growth element
        {
            s4.b0[m * 2] = -(abs(ss) - (-sn * tanf(0.524))) * std::copysign(1, ss);
            s4.b0[m * 2 + 1] = 0;
        }
    }

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - elm_list[m].ym);
    float pyy = symm.pyy1 + g.sky * (y0 - elm_list[m].ym);
    float pxy = symm.pxy1;
     
    b_elm[m].force1 = 2.0 * elm_list[m].a * (-((pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb)));     //!old b0()
    b_elm[m].force2 = 2.0 * elm_list[m].a * (-(pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb));  // !old b0()

    return;
}






float dxi(float& angle, bool flagB)
{
    /* -- add one element at the tip ---- 
    flagB is used to differentiate between dxi and dxiB 
    flagB =1 then it performs like dxiB in the original code*/  
   
   
    float xb = 0.0, yb = 0.0, xe = 0.0, ye = 0.0;
    float xt = 0.0, yt = 0.0, xt0 = 0.0, yt0 = 0.0, dtt = 0.0;
   
    xb = tips[ni].xbe;
    yb = tips[ni].ybe;
    xe = tips[ni].xen;
    ye = tips[ni].yen;

    int mm = tips[ni].mat_no;
    float ang = angle;         

    tips[ni].costem = (xe - xb) / tips[ni].dl;
    tips[ni].sintem = (ye - yb) / tips[ni].dl;
    //normal dxi
    if (!flagB)
    {
        if (tips[ni].ityp == 1 || tips[ni].ityp == 3)
        {
            xe = xb + tips[ni].dl * (cosf(ang) * tips[ni].costem - sinf(ang) * tips[ni].sintem);
            ye = yb + tips[ni].dl * (sinf(ang) * tips[ni].costem + cosf(ang) * tips[ni].sintem);
            xt = xe;
            yt = ye;
            xt0 = xb;
            yt0 = yb;
        }

        else if (tips[ni].ityp == -1 || tips[ni].ityp == -3)
        {
            xb = xe - tips[ni].dl * (cosf(ang) * tips[ni].costem - sinf(ang) * tips[ni].sintem);
            yb = ye - tips[ni].dl * (sinf(ang) * tips[ni].costem + cosf(ang) * tips[ni].sintem);
            xt = xb;
            yt = yb;
            xt0 = xe;
            yt0 = ye;
        }
    }
    //dxiB was called and flagB=true
    else
    {
        if (tips[ni].ityp == 4)
        {
            xe = xb + tips[ni].dl * (tips[ni].costem * cosf(ang) - sinf(ang) * tips[ni].sintem);
            ye = yb + tips[ni].dl * (sinf(ang) * tips[ni].costem + cosf(ang) * tips[ni].sintem);
            xt = xe;
            yt = ye;
            xt0 = xb;
            yt0 = yb;
        }
    }

    if (symm.ksym != 0)                     //optimization room Sara!  ksymm==4 extra
    {

        if (symm.ksym == 1 || symm.ksym == 4)
        {
            if ((xt - symm.xsym) * (xt0 - symm.xsym) <= 0)
            {
                //element cross the symmetry line
                yt = yt - (xt - symm.xsym) / (xt - xt0) * (yt - yt0);
                xt = symm.xsym;
                dtt = label_200(xt, yt, xt0, yt0, angle, flagB);
                recalculate_boundary_element_m(xt, yt, xt0, yt0, mm, dtt);   //label 500
                return dtt;
            }
        }
        if (symm.ksym == 2 || symm.ksym == 4)
        {
            //element cross the symmetry line
            if ((yt - symm.ysym) * (yt0 - symm.ysym) <= 0)
            {
                xt = xt - (yt - symm.ysym) / (yt - yt0) * (xt - xt0);
                yt = symm.ysym;
                dtt = label_200(xt, yt, xt0, yt0, angle, flagB);
                recalculate_boundary_element_m(xt, yt, xt0, yt0, mm, dtt);   //label 500
                return dtt;
            }
        }
        if (symm.ksym == 3 || symm.ksym == 4)
        {
            if ((xt - symm.xsym) * (yt0 - symm.ysym) == -(xt0 - symm.xsym) * (yt - symm.ysym))
            {
                xt = symm.xsym;
                yt = symm.ysym;
                dtt = label_200(xt, yt, xt0, yt0, angle, flagB);
                recalculate_boundary_element_m(xt, yt, xt0, yt0, mm, dtt);   //label 500
                return dtt;
            }
        }

    }
    check_all_tips(xt, yt, xt0, yt0, flagB);    //label55 in the oc      

    dtt = label_200(xt, yt, xt0, yt0,angle, flagB);
    recalculate_boundary_element_m(xt, yt, xt0, yt0,mm, dtt);   //label 500

    return dtt;
}
