#include <Initiation.h>
#include <CommonPara.h>
#include <Failure.h>
#include <Work.h>
#include <random>
#include <Source.h>
#include <ExcavationCracks.h>
#include<Failure.h>

using namespace CommonPara_h::comvar;



bool for_i_loop(int k, float x, float y,int m, int numbe0)
{   
    bool jmp1500 = false;
    float xbeg1, ybeg1, xend1, yend1;
    float xbeg2, ybeg2, xend2, yend2, dist;

    // If fracture intersects this point, skip the point   
    // there is the possibility to write seperate fuction for 
    // if parts but need to think about the efficiency, Sara!
    for (int i = 0; i < numbe0; ++i)
    {       
        if (elm_list[i].kod == 5) 
        {
             xbeg1 = elm_list[i].xm - elm_list[i].a * elm_list[i].cosbet;
             ybeg1 = elm_list[i].ym - elm_list[i].a * elm_list[i].sinbet;
             xend1 = elm_list[i].xm + elm_list[i].a * elm_list[i].cosbet;
             yend1 = elm_list[i].ym + elm_list[i].a * elm_list[i].sinbet;

             dist = min(sqrt(std::pow(x - xbeg1,2) + std::pow(y - ybeg1,2)),
                sqrt(std::pow(x - xend1,2)  + std::pow(y - yend1,2)));

            if (dist < min(elm_list[m].a, elm_list[i].a) / 100.0) 
            {
                jmp1500 = true;
                return jmp1500;
            }
            if (symm.ksym == 1 || symm.ksym == 4)
                {
                    xend2 = 2.0 * symm.xsym - xbeg1;
                    yend2 = ybeg1;
                    xbeg2 = 2.0 * symm.xsym - xend1;
                    ybeg2 = yend1;

                    dist = min(sqrt(std::pow(x - xbeg2, 2) + std::pow(y - ybeg2, 2)),
                        sqrt(std::pow(x - xend2, 2) + std::pow(y - yend2, 2)));

                    if (dist < min(elm_list[m].a, elm_list[i].a) / 100.0)
                    {
                        jmp1500 = true;
                        return jmp1500;
                    }
                }

                if (symm.ksym == 2 || symm.ksym == 4)
                {
                    xend2 = xbeg1;
                    yend2 = 2.0 * symm.ysym - ybeg1;
                    xbeg2 = xend1;
                    ybeg2 = 2.0 * symm.ysym - yend1;

                    dist = min(sqrt(std::pow(x - xbeg2, 2) + std::pow(y - ybeg2, 2)),
                        sqrt(std::pow(x - xend2, 2) + std::pow(y - yend2, 2)));

                    if (dist < min(elm_list[m].a, elm_list[i].a) / 100.0)
                    {
                        jmp1500 = true;
                        return jmp1500;
                    }
                   
                }
                //Sara ! use else for optim
                if (symm.ksym == 3)
                {
                    xbeg2 = 2 * symm.xsym - xbeg1;
                    ybeg2 = 2 * symm.ysym - ybeg1;
                    xend2 = 2 * symm.xsym - xend1;
                    yend2 = 2 * symm.ysym - yend1;

                    dist = min(sqrt(std::pow(x - xbeg2, 2) + std::pow(y - ybeg2, 2)),
                        sqrt(std::pow(x - xend2, 2) + std::pow(y - yend2, 2)));

                    if (dist < min(elm_list[m].a, elm_list[i].a) / 100.0)
                    {
                        jmp1500 = true;
                        return jmp1500;
                    }
                }
            }        
    }
    return jmp1500;
}

// Continue with the code after label 301



void reassigning_boundary_values(int ID,int m,int j,int k, float beta, float x, float y,
    float& st, float& bet)
{
    BoundaryElement& be = elm_list[m];  //be is an alias 
    if (j == numbe)
        j--;

    float ss_m =be.sigma_s ; float sn_m =be.sigma_n ; float ustem_m =be.us;

    float untem_m =be.un; float us_m =be.us_neg; float un_m =be.un_neg;

    float ss_n = elm_list[j].sigma_s ;
    float sn_n = elm_list[j].sigma_n ; 
    float ustem_n = elm_list[j].us ;

    float untem_n = elm_list[j].un ; 
    float us_n = elm_list[j].us_neg ; 
    float un_n = elm_list[j].un_neg ;
    
    float us_a = 0, un_a = 0, ang = 0;

   
    if (ID == 1 || ID == 2)
    {
        ss_n = -ss_n;
        us_n = -us_n;
    }

    if (k == 1) 
    {
        us_a = (be.a * (us_n * cosf(beta/2) - un_n * sinf(beta/2)) + elm_list[j].a *
            (us_m * cosf(beta/2) + un_m * sinf(beta/2))) / (be.a + elm_list[j].a);

        un_a = (be.a * (us_n * sinf(beta/2) + un_n * cosf(beta/2)) + elm_list[j].a *
            (-us_m * sinf(beta/2) + un_m * cosf(beta/2))) / (be.a + elm_list[j].a);

        ang = static_cast<float>(atan2(be.sinbet,be.cosbet)) + beta / 2.0;
    }

    else //if (k == 2)
    {
        us_a = (elm_list[j].a * (us_m * cosf(beta/2) - un_m * sinf(beta/2)) +
            elm_list[j].a * (us_n * cosf(beta/2) + un_n * sinf(beta/2))) / (be.a + elm_list[j].a);

        un_a = (elm_list[j].a * (us_m * sinf(beta/2) + un_m * cosf(beta/2)) + 
            elm_list[j].a *
            (-us_n * sinf(beta/2) + un_n * cosf(beta/2))) / (be.a + elm_list[j].a);

        ang = static_cast<float>(atan2(be.sinbet,be.cosbet)) - beta / 2.0;
    }

    
    float cosb = cosf(ang);
    float sinb = sinf(ang);
    float ss, sn;

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 -be.ym);
    float pyy = symm.pyy1 + g.sky * (y0 -be.ym);
    float pxy = symm.pxy1;

    float st0 = pxx * cosb * cosb - 2.0 * pxy * sinb * cosb + pyy * sinb * sinb;
    float ss0 = (pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb);
    float sn0 = pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb;

    int mm = be.mat_no;

    if (k == 1) 
    {
        ss = (ss_m + ss_n) / 2.0 + (sn_m - sn_n) / 2.0 * tanf(beta/2);
        sn = (sn_m + sn_n) / 2.0 + (-ss_m + ss_n) / 2.0 * tanf(beta/2);

        st = rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr) * 
            ((us_n * cosf(beta/2) - un_n * sinf(beta/2)) - 
            (us_m * cosf(beta/2) + un_m * sinf(beta/2))) 
            / ((be.a + elm_list[j].a) * cosf(beta/2));
    }
    else 
    {
        ss = (ss_n + ss_m) / 2.0 + (sn_n - sn_m) / 2.0 * tanf(beta/2);
        sn = (sn_n + sn_m) / 2.0 + (-ss_n + ss_m) / 2.0 * tanf(beta/2);

        st = rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr) * ((us_m * cos(beta / 2) - un_m *
            sinf(beta/2)) - (us_n * cosf(beta/2) + un_n *
                sinf(beta/2))) / ((be.a + elm_list[j].a) * cosf(beta/2));
    }

    st += st0 + (sn - sn0) * rock1[mm].pr / (1 - rock1[mm].pr);

    bet = (sn == st)? pi / 2.0: 0.5 * atanf(2.0 * ss / (st - sn));   

    float sig1 = st * cosf(bet) * cosf(bet) + 2 * ss * sinf(bet) * cosf(bet) + sn * sinf(bet) * sinf(bet);
    float sig2 = st * sinf(bet) * sinf(bet) - 2 * ss * sinf(bet) * cosf(bet) + sn * cosf(bet) * cosf(bet);
    float sig12 = (sig1 - sig2) / 2.0;

    
    bet = (k == 1) ? atan2f(be.sinbet, be.cosbet) + beta / 2.0 + bet : 
        atan2f(be.sinbet, be.cosbet) - beta / 2.0 + bet;   
   

    if (sig1 > sig2)
    {
        bet = bet + pi / 2.0;
        float temp = sig1;
        sig1 = sig2;
        sig2 = temp;
    }

    ///////final part 

    float r;
    //Use deault fracture initiation element size if it was not given
    if (s15.a_ini == 0) 
    {
        r = min(be.a, elm_list[j].a) * 2.0;
    }
    // Use given fracture initiation element size
    else
    {
        r = s15.a_ini;
    }

    // Failure criterion for tensile failure
    float sss = rock1[mm].rst;
    if (sss == 0) 
    {
        sss = 1e-6; // avoid division by zero
    }
    float failt = sig2 / rock1[mm].rst;
    float alphat = bet;
    float dcheck = sqrt(powf((r / 2 * cosf(alphat) - x), 2) +
        powf((r / 2 * sinf(alphat) - y), 2));

    if (dcheck > r)
    {
        alphat = bet - pi;
    }

    // Failure criterion for shear failure
    float shear_angle = (pi / 4 + rock1[mm].rphi / 2);
    if (shear_angle > pi / 3.0)
    {
        shear_angle = pi / 3.0;
    }

    float alphas1 = bet - pi / 2 + shear_angle;
    float alphas2 = bet - pi / 2 - shear_angle;
    float s1 = (2.0 * rock1[mm].rcoh * cosf(rock1[mm].rphi) + (-sig2) *
        (1.0 + sinf(rock1[mm].rphi))) / (1.0 - sinf(rock1[mm].rphi));
    sn = sig1 * sinf(alphas1 - bet) * sinf(alphas1 - bet) + sig2 * cosf(alphas1 - bet) 
        * cosf(alphas1 - bet);

    if (s1 == 0) 
    {
        s1 = 1e-6;
    }
    float fails = abs(sig1) / s1;

    // Check if out of window range
    if (x > s5u.xmax || x < s5u.xmin || y > s5u.ymax || y < s5u.ymin)
    {
        failt = 0;
        fails = 0;        
    }

    //label 205

    int seed = int((1 + x + y) * 1000);          //random failure
    std::mt19937 gen(seed);        
    
    //random_number(rand);                 //random failure
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    ///float rand = std::generate_canonical<float, std::numeric_limits<float>::
        ///digits>(gen);   
    float randn = distribution(gen);
    randn = randn * s15.i_rand * s15.l_rand;          //i_rand = 1 or 0;
       
    if (failt >= s15.f_ini0 && randn <= std::pow((failt - s15.f_ini0) / (1.0001 - s15.f_ini0),2))
    {
        Sum_Failure(x, y, r, alphat, 1, failt, 1); //last "1" is location, ie.on boundary
        return;            
    }

    if ((symm.ksym == 2 || symm.ksym == 4) && std::abs(y - symm.ysym) <= 1e-5)
        return;
    if ((symm.ksym == 3) && (std::abs(x - symm.xsym) <= 1e-5 && std::abs(y - symm.ysym)
        <= 1e-5))
        return;

    float set, zet;
    float disp;
    float temp1 = x + 0.75 * r * cosf(alphas1);
    float temp2 = y + 0.75 * r * sinf(alphas1);
    point(temp1, temp2, sig1, sig2, bet, sig12, set, disp, zet);
    s1 = (2.0 * rock1[mm].rcoh * cosf(rock1[mm].rphi) + (-sig2) *
        (1. + sinf(rock1[mm].rphi))) / (1.0 - sinf(rock1[mm].rphi));

    if (s1 == 0) s1 = 1e-6;
    float sig_sum1 = std::abs(sig1) / s1;

    temp1 = x + 0.75 * r * cosf(alphas2);
    temp2 = y + 0.75 * r * sinf(alphas2);
    point(temp1, temp2, sig1, sig2, bet, sig12, set, disp, zet);

    s1 = (2.0 * rock1[mm].rcoh * cosf(rock1[mm].rphi) + (-sig2) * 
        (1. + sinf(rock1[mm].rphi))) / (1.0 - sinf(rock1[mm].rphi));
    if (s1 == 0) s1 = 1e-6;
    float sig_sum2 = std::abs(sig1) / s1;

    if (fails >= s15.f_ini0 && randn <= std::pow((fails - s15.f_ini0) /
        (1.0001 - s15.f_ini0),2)) 
    {
        if (sig_sum1 > sig_sum2)
            // last "1" is location, ie. on boundary
            Sum_Failure(x, y, r, alphas1, 2, fails, 1);        
        else //if (sig_sum1 <= sig_sum2)
            Sum_Failure(x, y, r, alphas2, 2, fails, 1); 
        return;
    }

    return;
}




std::pair<bool, bool> for_symm_similarity(int j, int m, float xt, float yt, float xbeg1, float ybeg1,
    float xend1, float yend1,int k, float& beta, float& x, float& y,int & ID,int symcase)
{
    //symcase check the symm case =3 forsymm =3,=2 for symm 2 or 4,=1: for symm 1 or 4
    float xt1, yt1;
    bool jmp_1500 = false;
    bool jmp_250 = false;   

    if (k == 2)
    {
        xt1 = xend1;
        yt1 = yend1;
    }
    else
    {
        xt1 = xbeg1;
        yt1 = ybeg1;
    }

    float dist = sqrt((xt - xt1) * (xt - xt1) + (yt - yt1) * (yt - yt1));
    float min_val = 0.0;

    if (symcase == 3)
         min_val = min(elm_list[m].a, elm_list[j].a) / 100;
    else
         min_val = min(elm_list[m].a, elm_list[j].a) / 1000.0;

    if (dist <= min_val)     //min used instead of Amin1
    {
        x = xt;
        y = yt;
        float temp = sqrt((xend1 + xbeg1) * (xend1 + xbeg1) + (yend1 + ybeg1) * 
            (yend1 + ybeg1));

        float sinbet_j = (yend1 - ybeg1) / temp;
        float cosbet_j = (xend1 - xbeg1) / temp;

        float sinb = sinbet_j * elm_list[m].cosbet - cosbet_j * elm_list[m].sinbet;
        float cosb = cosbet_j * elm_list[m].cosbet + sinbet_j * elm_list[m].sinbet;
        beta = (k == 1) ?(atan2(sinb, cosb)) :
           (atan2(-sinb, cosb));
        if (symcase != 3)
        {
            ID = symcase;            
            if (beta > 30.0 * pi / 180.0 || beta < -30.0 * pi / 180.0)
            {
                jmp_1500 = true;
                return make_pair(jmp_1500, jmp_250);
            }
           
        }        
        else
        {
            ID = 3;            
            if (beta < -30.0 / pi / 180. || beta > 30. * pi / 180.) //acute angle or sharp angle, stress not reliable
            {
                jmp_1500 = true;
                return make_pair(jmp_1500, jmp_250);
            }               
        }      
        
        jmp_250 = true;
        return make_pair(jmp_1500, jmp_250);
    }
    return make_pair(jmp_1500, jmp_250);
}


std::pair<bool, bool> symm_similar_internal(int j, int m, float xt, float yt, float xbeg1, float ybeg1,
    float xend1, float yend1, int k, float& beta, float& x, float& y, int& ID, int symcase)
{
    //symcase check the symm case =3 forsymm =3,=2 for symm 2 or 4,=1: for symm 1 or 4
    float xt1, yt1;
    bool jmp_1500 = false;
    bool jmp_250 = false;

    if (k == 2)
    {
        xt1 = xend1;
        yt1 = yend1;
    }
    else
    {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
        xt1 = xbeg1;
        yt1 = ybeg1;
    }

    float dist = sqrt((xt - xt1) * (xt - xt1) + (yt - yt1) * (yt - yt1));
    float min_val = 0.0;

    if (symcase == 3)
        min_val = min(elm_list[m].a, elm_list[j].a) / 100;
    else
        min_val = min(elm_list[m].a, elm_list[j].a) / 1000.0;

    if (dist <= min_val)     //min used instead of Amin1
    {
        x = xt;
        y = yt;
        float temp = sqrt((xend1 + xbeg1) * (xend1 + xbeg1) + (yend1 + ybeg1) *
            (yend1 + ybeg1));

        float sinbet_j = (yend1 - ybeg1) / temp;
        float cosbet_j = (xend1 - xbeg1) / temp;

        float sinb = sinbet_j * elm_list[m].cosbet - cosbet_j * elm_list[m].sinbet;
        float cosb = cosbet_j * elm_list[m].cosbet + sinbet_j * elm_list[m].sinbet;
        beta = (k == 1) ? atan2f(sinb, cosb) :
            atan2f(-sinb, cosb);
        
        ID = symcase;
        if (beta > 30.0 * pi / 180.0 || beta < -30.0 * pi / 180)
        {
            jmp_1500 = true;
            return make_pair(jmp_1500, jmp_250);
        }       

        jmp_250 = true;
        return make_pair(jmp_1500, jmp_250);
    }
    return make_pair(jmp_1500, jmp_250);
}




std::pair<bool,bool> for_j_loop(int k, int m, float xt, float yt, float & x, float& y,
                                int& ID, float& beta, int& jj, int numbe0,bool intern_cal)
    
    /*this function called by both initiation and internal functions
    flag intern_call differentiatie these calls*/
{
    //beta = 0;
    //jj = 0;
    //ID = 0;
    float xt1, yt1;
    int mm = elm_list[m].mat_no;
    std::pair<bool, bool> jmp;
    jmp.first = false;
    jmp.second = false;
   
    bool jmp_250 = false;
    bool jmp_1500 = false;

    float xbeg2, ybeg2, xend2, yend2, dist;
    float min_val = 100.0;
    float x180 = 180.0;
    if (intern_cal == true)
    {
        min_val = 1000.0;
        x180 = 180;
    }
    for (int j = 0; j < numbe0; ++j)
    {
        
        jj = j;
        if (elm_list[j].kod == 5 || elm_list[j].mat_no != mm)
            continue;        

        float xbeg1 = elm_list[j].xm - elm_list[j].a * elm_list[j].cosbet;
        float ybeg1 = elm_list[j].ym - elm_list[j].a * elm_list[j].sinbet;
        float xend1 = elm_list[j].xm + elm_list[j].a * elm_list[j].cosbet;
        float yend1 = elm_list[j].ym + elm_list[j].a * elm_list[j].sinbet;

        if (k == 2) 
        {
            xt1 = xend1;
            yt1 = yend1;
        }
        else
        {
            xt1 = xbeg1;
            yt1 = ybeg1;
        }        
        
        dist = sqrt((xt - xt1) * (xt - xt1) + (yt - yt1) * (yt - yt1));       
           
        if (dist <= (min(elm_list[m].a, elm_list[j].a) / min_val))   
        {
            x = xt;
            y = yt;
            ID = 0;
            float sinb = elm_list[j].sinbet * elm_list[m].cosbet - elm_list[j].cosbet * 
                elm_list[m].sinbet;
            float cosb = elm_list[j].cosbet * elm_list[m].cosbet + elm_list[j].sinbet *
                elm_list[m].sinbet;

            beta = (k == 1) ? atan2f(sinb, cosb):
                atan2f(-sinb, cosb);
           
            //acute angle or sharp angle, stress not reliable
            if (beta < ( - 30.0 * pi / x180) || beta > (30.0 * pi / 180.0) )
            {
                jmp_1500 = true;
                return make_pair(jmp_1500, jmp_250);
            }            
            jmp_250 = true;
            return make_pair(jmp_1500, jmp_250);
        }
    
        //label_201:
        
            if (symm.ksym == 3)
            {
                xbeg2 = 2 * symm.xsym - xbeg1;
                ybeg2 = 2 * symm.ysym - ybeg1;
                xend2 = 2 * symm.xsym - xend1;
                yend2 = 2 * symm.ysym - yend1;

                if(intern_cal)
                    jmp = symm_similar_internal(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                        ID,3);   //check again 
                else 
                    jmp = for_symm_similarity(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                        ID, 3);   
               
                if (jmp.first == true || jmp.second == true)
                    return jmp;
            }

            if (symm.ksym == 1 || symm.ksym == 4)
            {
                xend2 = 2.0 * symm.xsym - xbeg1;
                yend2 = ybeg1;
                xbeg2 = 2.0 * symm.xsym - xend1;
                ybeg2 = yend1;

                if (intern_cal)
                    jmp = symm_similar_internal(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                    ID,1);   //check again
                else
                    jmp = for_symm_similarity(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                    ID, 1);
                
                if (jmp.first == true || jmp.second == true)
                    return jmp;
            }

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                xend2 = xbeg1;
                yend2 = 2.0 * symm.ysym - ybeg1;
                xbeg2 = xend1;
                ybeg2 = 2.0 * symm.ysym - yend1;

                if (intern_cal)
                    jmp = symm_similar_internal(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                    ID,2); 
                else
                    jmp = for_symm_similarity(j, m, xt, yt, xbeg2, ybeg2, xend2, yend2, k, beta, x, y,
                    ID, 2);


                if (jmp.first == true || jmp.second == true)
                    return jmp;
            }
        
    }
    //jj = numbe;
    return make_pair(jmp_1500,jmp_250);
}




void InitiationB()
{
    float xtem[m0] = { 0 }, ytem[m0] = { 0 };
    float x = 0, st = 0;
    float y = 0;
    float beta = 0.0;
    int numbe0 = numbe;
    int ID;
    int jj = 0;
    std::pair<bool, bool> jmp = make_pair(false,false);
   
    bool jmp1500 = false;
    float randn = 0;
    float bet = 0;
    for (int m = 0; m < numbe0; ++m) 
    {
        if (elm_list[m].kod == 5)
            continue;       
        jmp1500 = false;

        float xbeg = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
        float ybeg = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;
        float xend = elm_list[m].xm + elm_list[m].a * elm_list[m].cosbet;
        float yend = elm_list[m].ym + elm_list[m].a * elm_list[m].sinbet;

        int mm = elm_list[m].mat_no;
        for (int k = 1; k <= 2; ++k)
        {
            jmp1500 == false;
            float xt, yt;
            if (k == 1)
            {
                xt = xend;
                yt = yend;
            }
            else
            {
                xt = xbeg;
                yt = ybeg;
            }
            //jj is the last j when loop finishes
            jmp = for_j_loop(k, m, xt, yt, x , y ,ID,beta,jj,numbe0, false);
            //if jmp to 1500 true continue with the next iterate
            if (jmp.first == true)
                continue;

            if (k == 1)
            {
                xtem[m] = x;
                ytem[m] = y;
            }
            // Avoid repeating one point  Sara ! optimization
            else if (k == 2)
            {
                jmp1500 = false;
                for (int i = 0; i < numbe0; ++i)
                {
                    if (std::abs(x - xtem[i]) < s15.aa / 1000 && 
                        std::abs(y - ytem[i]) < s15.aa / 1000.0)
                    {
                        jmp1500 = true;
                        break;
                    }
                }
                if (jmp1500 == true)
                    continue;
            }
            
            bool jmp_1500 = for_i_loop(k, x, y, m, numbe0);
            if (jmp_1500 == true)
                continue; 
            
            reassigning_boundary_values(ID, m, jj, k, beta, x, y,st, bet);            
        }    
    }
}




float call_point_and_cal_sig_sum(int index, float xp, float yp, float alphas1, int mm,
     float r)
{
    float set, zet, bet;
    float disp, sig1, sig2, sig12;
    float sig_sum1 = 0;
    float s1 = 0;
    if (index == 1)
    {
        point(xp + 0.5 * r * cosf(alphas1), yp + 0.5 * r * sinf(alphas1), sig1, sig2,
            bet, sig12, set, disp, zet);
        s1 = (2.0 * rock1[mm].rcoh * cosf(rock1[mm].rphi) + (-sig2) *
            (1. + sinf(rock1[mm].rphi))) / (1.0 - sinf(rock1[mm].rphi));
        if (s1 == 0) s1 = 1e-6;
        sig_sum1 = abs(sig1) / s1;
    }
    else
    {
        point(xp - 0.5 * r * cosf(alphas1), yp - 0.5 * r * sinf(alphas1), sig1, sig2, bet, sig12, set, disp, zet);
        s1 = (2.0 * rock1[mm].rcoh * cosf(rock1[mm].rphi) + (-sig2) *
            (1. + sinf(rock1[mm].rphi))) / (1.0 - sinf(rock1[mm].rphi));
        if (s1 == 0) s1 = 1e-6;
        sig_sum1 = abs(sig1) / s1;
    }
    return sig_sum1;

}



//********* Fracture initiation in rock ******************************

void InitiationR()
{      
    float xp = 0, yp = 0, sss,  alphat,
        alphas1, alphas2, ss, sn,  r;
    int mm, n_valid;
    float sig_sum1 = 0, sig_sum2 = 0;
    float set, zet, bet;
    float disp, sig1, sig2, sig12;   
    double failt, fails;

    for (auto it = valid1.begin(); it != valid1.end(); ) {

        xp = (*it).first;
        yp = (*it).second;        
        n_valid = 1;
        check_point_in_rock(xp, yp, 1, n_valid);
        if (n_valid == 0) {
            it = valid1.erase(it);
            continue;
        }

        int material = check_material_id(xp, yp);
        mm = material;

        point(xp, yp, sig1, sig2, bet, sig12, set, disp, zet);
        //random failure
        int seed = int((1 + xp + yp) * 1000);
        std::mt19937 gen(seed);
        std::uniform_real_distribution<float> distribution(0.0, 1.0);
        float randn = distribution(gen);
        randn = randn * s15.i_rand * s15.l_rand;

        //tensile failure
        sss = rock1[mm].rst;
        if (sss == 0) sss = 1e-6;     //avoid division by zero;
        failt = sig2 / rock1[mm].rst;
        alphat = bet;

        alphas1 = bet - pi / 2.0 + (pi / 4.0 + rock1[mm].rphi / 2.0);
        alphas2 = bet - pi / 2.0 - (pi / 4.0 + rock1[mm].rphi / 2.0);
        ss = -0.5 * (sig1 - sig2) * sinf(2 * (alphas2 - bet));
        sn = sig1 * sinf(alphas2 - bet) * sinf(alphas2 - bet) + sig2 * cosf(alphas2 - bet)
            * cosf(alphas2 - bet);

        sss = max(-sn, 0.0) * tanf(rock1[mm].rphi) + rock1[mm].rcoh;
        if (sss == 0) sss = 1e-6;
        fails = abs(ss) / sss;

        // use default fracture initiation element size if it was not given
        r = (s15.a_ini == 0) ? s15.aa : s15.a_ini;       // use given fracture initiation element size

        if (failt >= s15.f_ini0 && randn <= powf((failt - s15.f_ini0) /
            (1.0001 - s15.f_ini0), 2))
        {
            Sum_Failure(xp, yp, r, alphat, 1, failt, 2); // last "2" is location, ie. intact rock
            ktipgrow = true;
            ++it;
            continue;
        }

        sig_sum1 = call_point_and_cal_sig_sum(1, xp, yp, alphas1, mm, r);
        sig_sum1 += call_point_and_cal_sig_sum(2, xp, yp, alphas1, mm, r);
        sig_sum2 = call_point_and_cal_sig_sum(1, xp, yp, alphas2, mm, r);
        sig_sum2 += call_point_and_cal_sig_sum(2, xp, yp, alphas2, mm, r);

        if (fails >= s15.f_ini0 && randn <= std::powf((fails - s15.f_ini0) / (1.0001 - s15.f_ini0), 2))
        {
            if (sig_sum1 > sig_sum2)
                Sum_Failure(xp, yp, r, alphas1, 2, fails, 2); //last "2" is location, ie.intack rock
            else  //if (sig_sum1 <= sig_sum2) 
                Sum_Failure(xp, yp, r, alphas2, 2, fails, 2);
            ktipgrow = true;
        }
        ++it;   
    }   
  }




void initiation()
 {
      mf = 0;      // Number of failure points is set to zero initially
      if (s15.i_bound == 1)
          InitiationB();

      if (s15.i_intern == 1)
          InitiationR();
     
      Choose_Failure();           
  }
