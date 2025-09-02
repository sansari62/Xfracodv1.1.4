#include<stdafx.h>

#include<ExcavationCracks.h>
#include <cstdlib>
#include <ctime>
#include "CommonPara.h"
#include<Failure.h>
#include<Rectangle_check.h>

using namespace CommonPara_h::comvar;


   


void check_symm_and_set_n_valid(float xp, float yp, bool flag, int& n_valid)
{  
    float dist = 0.0, xb = 0.0, yb = 0.0, xe = 0.0, ye = 0.0;
    float thr = 0;   // 0.0006;
    if (symm.ksym == 1 || symm.ksym == 4)
    {
        for (int l = 0 + prenumbe; l < numbe; ++l)
        {
            BoundaryElement& be = elm_list[l];
            if (be.kod == 7) continue;

            dist = sqrt(pow(xp - (2.0 * symm.xsym - be.xm), 2) + 
                pow(yp - be.ym, 2));
            xb = be.xm - be.a * be.cosbet;
            yb = be.ym - be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - (2.0 * symm.xsym - xb), 2) +
                pow(yp - yb, 2)));
            xe = xb + 2.0 * be.a * be.cosbet;
            ye = yb + 2.0 * be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - (2.0 * symm.xsym - xe), 2) + 
                pow(yp - ye, 2)));

            if (flag == 0)
            {
                if (dist <= 1.1 * be.a) {
                    n_valid = 0;
                    return;
                }
            }
            else
            {
                if (dist <= s15.aaa) {
                    n_valid = 0;
                    return;
                }
            }
        }
    }

    if (symm.ksym == 2 || symm.ksym == 4)
    {
        for (int l = 0 + prenumbe; l < numbe ; ++l)
        {
            BoundaryElement& be = elm_list[l];
            if (be.kod == 7) continue;

            dist = sqrt(pow(xp - be.xm, 2) + pow(yp - (2.0 * symm.ysym - be.ym), 2));
            xb = be.xm - be.a * be.cosbet;
            yb = be.ym - be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - xb, 2) + pow(yp - (2.0 * symm.ysym - yb), 2)));
            xe = xb + 2.0 * be.a * be.cosbet;
            ye = yb + 2.0 * be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - xe, 2) + pow(yp - (2.0 * symm.ysym - ye), 2)));

            if (flag == 0)
            {
                if (dist <= 1.1 * be.a) {
                    n_valid = 0;
                    return;
                }
            }
            else
            {
                if (dist <= s15.aaa) {
                    n_valid = 0;
                    return;
                }
            }
        }
    }

    if (symm.ksym == 3 || symm.ksym == 4)
    {
        for (int l = 0 + prenumbe; l < numbe ; ++l)
        {
            BoundaryElement& be = elm_list[l];
            if (be.kod == 7) continue;

            dist = sqrt(pow(xp - (2.0 * symm.xsym - be.xm), 2) + pow(yp - (2.0 * symm.ysym - be.ym),2));
            xb = be.xm - be.a * be.cosbet;
            yb = be.ym - be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - (2.0 * symm.xsym - xb), 2) + pow(yp - (2.0 * symm.ysym - yb),2)));
            xe = xb + 2.0 * be.a * be.cosbet;
            ye = yb + 2.0 * be.a * be.sinbet;
            dist = min(dist, sqrt(pow(xp - (2.0 * symm.xsym - xe), 2) + pow(yp - (2.0 * symm.ysym - ye),2)));

            if (flag == 0)
            {
                if (dist <= 1.1 * be.a) {
                    n_valid = 0;
                    return;
                }
            }
            else
            {
                if (dist <= s15.aaa) {
                    n_valid = 0;
                    return;
                }
            }
        }
    }
}





void  check_point_in_rock(float xp, float yp, bool flag, int& n_valid)
{
    /*
     -- check if a point is in the rock mass ---
    flag is used to differentiate between this function and function check_point_rock1
    since the most of both codes are the same  so flag = 0 
    is this one and flag = 1 is check_point_in_rock1 
    */

    int k = 0;
    float x0, y0, dd0, dd1, dd2, xc = 0, yc = 0, dist = 0,
        xb, yb, xe, ye, sinb, cosb;
    const float epsilon = 1e-6;
    
    for (int l = 0; l < ntunnel; ++l)         
    {
        x0 = tunnl.x_cent[l];
        y0 = tunnl.y_cent[l];
        dd0 = tunnl.diameter[l];
        if (sqrt(pow(xp - x0, 2) + pow(yp - y0, 2)) < 0.5 * dd0)
        {
            n_valid = 0;
            return;
        }
    }

    for (int l = 0; l < nellipse; ++l)     
    {
        x0 = ellip_list[l].x_cent;
        y0 = ellip_list[l].y_cent;
        dd1 = 0.5 * ellip_list[l].diameter1;
        dd2 = 0.5 * ellip_list[l].diameter2;
        if (sqrt(pow(xp - x0, 2) / pow(dd1, 2) + 
            pow(yp - y0, 2) / pow(dd2, 2)) < 1.0)
        {
             n_valid = 0;
            return ;
        }
    }
    /*Rectangle1 rect = check_rectangle(false);
    Point p1;
    p1.x = xp;
    p1.y = yp;
    int p1_state = point_inside_rectangle(p1, rect);
    if (p1_state == 0)
    {
        n_valid = 0;
        return;
    }*/
    optional<Rectangle1> rect = check_rectangle(false);
    if(rect)
    {
        Point p1;
        p1.x = xp;
        p1.y = yp;
        int p1_state = point_inside_rectangle(p1, *rect);
        if (p1_state == 0)
        {
            n_valid = 0;
            return;
        }
    }
    float dist0 = 1e6;
    int n = -1;
    //**find the element closest to the point**  
    //  Sara optimization if symm = 0 just continue
    for (int m = 0 + prenumbe; m < numbe ; ++m)
    {
        if (elm_list[m].kod == 5 || elm_list[m].kod == 7) continue;

        dist = sqrt(pow(xp - elm_list[m].xm, 2) + pow(yp - elm_list[m].ym, 2));

        if (dist < dist0)
        {
            n = m;
            k = 0;
            dist0 = dist;
        }
        if (symm.ksym == 1 || symm.ksym == 4)
        {
            xc = 2.0 * symm.xsym - elm_list[m].xm;
            yc = elm_list[m].ym;
            dist = sqrt(pow(xp - xc, 2) + pow(yp - yc, 2));
            if (dist < dist0) 
            {
                n = m;
                k = 1;
                dist0 = dist;
            }
        }

        if (symm.ksym == 2 || symm.ksym == 4)
        {
            xc = elm_list[m].xm;
            yc = 2.0 * symm.ysym - elm_list[m].ym;
            dist = sqrt(pow(xp - xc, 2) + pow(yp - yc, 2));
            if (dist < dist0)
            {
                n = m;
                k = 2;
                dist0 = dist;
            }
        }

        if (symm.ksym == 3 || symm.ksym == 4)
        {
            xc = 2.0 * symm.xsym - elm_list[m].xm;
            yc = 2.0 * symm.ysym - elm_list[m].ym;
            dist = sqrt(pow(xp - xc, 2) + pow(yp - yc, 2));
            if (dist < dist0) 
            {
                n = m;
                k = 3;
                dist0 = dist;
            }
        }
    }

    if (n != -1)   //  n=-1 means no boundary, only fractures
    {

        float xpprime = (xp - elm_list[n].xm) * elm_list[n].cosbet + (yp - elm_list[n].ym) *
            elm_list[n].sinbet;
        float ypprime = -(xp - elm_list[n].xm) * elm_list[n].sinbet + (yp - elm_list[n].ym) *
            elm_list[n].cosbet;

        switch (k)
        {
        case 1:
            xc = 2.0 * symm.xsym - elm_list[n].xm;
            yc = elm_list[n].ym;
            sinb = -elm_list[n].sinbet;
            cosb = elm_list[n].cosbet;
            ypprime = -(xp - xc) * sinb + (yp - yc) * cosb;
            break;
        case 2:
            xc = elm_list[n].xm;
            yc = 2.0 * symm.ysym - elm_list[n].ym;
            sinb = elm_list[n].sinbet;
            cosb = -elm_list[n].cosbet;
            ypprime = -(xp - xc) * sinb + (yp - yc) * cosb;
            break;
        case 3:
            xc = 2.0 * symm.xsym - elm_list[n].xm;
            yc = 2.0 * symm.ysym - elm_list[n].ym;
            sinb = -elm_list[n].sinbet;
            cosb = -elm_list[n].cosbet;
            ypprime = -(xp - xc) * sinb + (yp - yc) * cosb;            
        }

        //!if the closest element is the interface, we need to check the pair element too
        if (elm_list[n].kod == 6)
        {
            int nn = 0;
            if (elm_list[n - 1].xm == elm_list[n].xm && elm_list[n - 1].ym == elm_list[n].ym)
                nn = n - 1;
            if (elm_list[n + 1].xm == elm_list[n].xm && elm_list[n + 1].ym == elm_list[n].ym)
                nn = n + 1;

            ypprime = min(ypprime, -(xp - elm_list[nn].xm) * elm_list[nn].sinbet +
                (yp - elm_list[nn].ym) * elm_list[nn].cosbet);
        }

        if (ypprime > 0)
        {
            n_valid = 0;
            return;
        }
    }
    float thr = 0;//0.0006;
  
    for (int l = 0 + prenumbe; l < numbe ; ++l)
    {
        BoundaryElement& be = elm_list[l];
        if (be.kod == 7) continue;

        dist = sqrt(pow(xp - be.xm, 2) + pow(yp - be.ym, 2));
        xb = be.xm - be.a * be.cosbet;
        yb = be.ym - be.a * be.sinbet;
        dist = min(dist, sqrt(pow(xp - xb, 2) + pow(yp - yb, 2)));
        xe = xb + 2.0 * be.a * be.cosbet;
        ye = yb + 2.0 * be.a * be.sinbet;
        dist = min(dist, sqrt(pow(xp - xe, 2) + pow(yp - ye, 2)));
        if (flag == 0)
        {
            if ((dist  <= 1.1 * be.a))
            {
                n_valid = 0;
                return;
            }
                
        }
        else
        {
            if (dist  <= s15.aaa) {
                n_valid = 0;
                return;
            }
        }
    }
     if(symm.ksym != 0)
        check_symm_and_set_n_valid(xp, yp, flag, n_valid);

    return ;
}






void ExcavationCracks() 
{
    float dist = 0;
    float  xb = 0, yb = 0, xe = 0, ye = 0, xp = 0, yp = 0;       
    int  n_valid;
    for (auto it = valid.begin(); it != valid.end(); ) {
             
        xp = (*it).first;
        yp = (*it).second;
        n_valid = 1;
        check_point_in_rock(xp, yp, 0, n_valid);
        if (n_valid == 0)
        {
            it = valid.erase(it);
            continue;
        }         

        dist = 1e6;
        for (int m = 0; m < numbe; ++m)
        {
            if (elm_list[m].kod == 5) continue;
			dist = min(dist, std::sqrt(std::pow(xp - elm_list[m].xm, 2) + std::pow(yp - elm_list[m].ym, 2)));
			xb = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
			yb = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;

			dist = min(dist, std::sqrt(std::pow(xp - xb, 2) + std::pow(yp - yb, 2)));
			xe = xb + 2.0 * elm_list[m].a * elm_list[m].cosbet;
			ye = yb + 2.0 * elm_list[m].a * elm_list[m].sinbet;
			dist = min(dist, std::sqrt(std::pow(xp - xe, 2) + std::pow(yp - ye, 2)));
           
        }

        if (dist > exca.d_wall) {
            ++it; 
            continue;
        }
        float randn;
        srand(time(nullptr));      
        randn = static_cast<float>(std::rand()) / RAND_MAX;
        if (randn > exca.rand_e) {
            ++it;
            continue;
        }
           
        srand(time(nullptr)); // Seed for random_comvar::number
        randn = static_cast<float>(std::rand()) / RAND_MAX;
        float alpha = randn * (2.0 * pi);

        float r = s15.a_ini;
        int legal;

        legal = CheckNewElement( r, xp, yp, cosf(alpha), sinf(alpha),0);
        float fos = 1;
        int im = 4;  // fracture initiation mode 0 - elastic, 1 - tensile; 2 - shear; 4 - blast
        if (legal == 1) 
            failure(xp, yp, r, alpha, im);        
        ++it;      
    }
}










