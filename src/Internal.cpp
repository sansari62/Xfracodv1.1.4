#include<stdafx.h>

#include <Internal.h>
#include <CommonPara.h>
#include<ExcavationCracks.h>
#include<WinInterface.h>
#include<Initiation.h>
#include<Mainb.h>
#include<Failure.h>


using namespace comvar;
using namespace winvar;




void compute_disp_boundary_and_frac_surfaces(int & npoint, float  sig1, float sig2, float bet, float sig12,
    float set, int mm, stringstream& buffer)
{
    /*--- displacements on boundary and fracture surfaces */
    float ss = 0, sn = 0, ustem = 0, untem = 0, usneg = 0, unneg = 0;
    float uspos = 0, unpos = 0, uxpos = 0, uypos = 0, uxneg = 0,
        uyneg = 0, ux = 0, uy = 0;
    float disp = 0, zet = 0;
    int index = npoint;
    for (int m = 0; m < numbe; ++m) 
    {
        BoundaryElement& be = elm_list[m];
        if (be.kod == 7)
            continue;

        if (dispwin.ID_win == 1) 
        {
            if (be.xm > dispwin.xur || be.xm < dispwin.xll || be.ym > dispwin.yur || be.ym < dispwin.yll)
                continue;
        }
        else    // Id_win=2 
        {
            if (sqrt(pow(be.xm - dispwin.xc0, 2) + pow(be.ym - dispwin.yc0, 2)) > dispwin.radium)
                continue;
        }

        
        float disp1 = 0.0, zet1 = 0.0;
        be.bound(m, ss, sn, ustem, untem, usneg, unneg);

         uspos = usneg - s4.d0[2 * m ];
         unpos = unneg - s4.d0[2 * m + 1];
         uxpos = uspos * be.cosbet - unpos * be.sinbet;
         uypos = uspos * be.sinbet + unpos * be.cosbet;
         uxneg = usneg * be.cosbet - unneg * be.sinbet;
         uyneg = usneg * be.sinbet + unneg * be.cosbet;

         ux = uxneg;
         uy = uyneg;
         disp = sqrt(ux * ux + uy * uy);
         zet = (ux > 0.0) ? atanf(uy / ux) : ((ux < 0.0) ? pi + atanf(uy / ux) :
            1.57 * copysign(1, uy ));     //isign Sara!  if(ux==0.) zet=1.57*isign(1,INT(uy))
       

        // Write data to file and store in stress array
       // write_data_to_file( be.xm, be.ym, sig1, sig2, bet, sig12, set, disp, zet, mm, buffer);
               
        Stress st(be.xm, be.ym, 0, 0, 0, 0, 0, disp, zet, mm);   //new stress input defined
        //Stress st(be.xm, be.ym, sig1, sig2, bet, sig12, set, disp, zet, mm);//just for comparison
        stress[npoint] = st;
        npoint++;
        
        //maybe we can put all of following in one func and improve it Sara!  #2
        if (be.kod == 5)
        {
            disp1 = sqrt(pow(uxpos, 2) + pow(uypos, 2));
            zet1 = (uxpos > 0.) ? atanf(uypos / uxpos) :
                ((uxpos < 0.) ? pi + atanf(uypos / uxpos) :
                1.57 * copysign(1, uypos));

           // write_data_to_file( be.xm, be.ym, sig1, sig2, bet, sig12, set, disp1, zet1, mm, buffer);
                        
            Stress st(be.xm, be.ym, 0, 0, 0, 0, 0, disp1, zet1, mm);   //new stress input defined
            //Stress st(be.xm, be.ym, sig1, sig2, bet, sig12, set, disp1, zet1, mm);
            stress[npoint] = st;
            npoint++;
            
        }

        if (symm.ksym != 0)
        {

            if (symm.ksym == 1 || symm.ksym == 4)
            {
                float xval = 2.0 * symm.xsym - be.xm;
               // write_data_to_file(xval, be.ym, sig1, sig2, bet, sig12, set, disp, pi - zet, mm, buffer);
                Stress st(xval, be.ym, 0, 0, 0, 0, 0, disp, pi - zet, mm);   //new stress input defined
                //Stress st(xval, be.ym, sig1, sig2, bet, sig12, set, disp, pi - zet, mm);
                stress[npoint] = st;
                npoint++;

                if (be.kod == 5)
                {
                   // write_data_to_file(xval, be.ym, sig1, sig2, bet, sig12, set, disp1, pi - zet1, mm, buffer);

                    Stress st(xval, be.ym, 0, 0, 0, 0, 0, disp1, pi - zet1, mm);   //new stress input defined
                    //Stress st(xval, be.ym, sig1, sig2, bet, sig12, set, disp1, pi - zet1, mm);
                    stress[npoint] = st;
                    npoint++;
                }
            }

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                float yval = 2.0 * symm.ysym - be.ym;
               // write_data_to_file(be.xm, yval, sig1, sig2, bet, sig12, set, disp, -zet, mm, buffer);
                Stress st(be.xm, yval, 0, 0, 0, 0, 0, disp, -zet, mm);   //new stress input defined
                //Stress st(be.xm, yval, sig1, sig2, bet, sig12, set, disp,  - zet, mm);
                stress[npoint] = st;
                npoint++;

                if (be.kod == 5)
                {
                    //write_data_to_file(be.xm, yval, sig1, sig2, bet, sig12, set, disp1, -zet1, mm, buffer);
                    Stress st(be.xm, yval, 0, 0, 0, 0, 0, disp1, -zet1, mm);   //new stress input defined
                    //Stress st(be.xm, yval, sig1, sig2, bet, sig12, set, disp1,- zet1, mm);
                    stress[npoint] = st;
                    npoint++;
                }
            }
            if (symm.ksym == 3 || symm.ksym == 4)
            {
                float xval = 2.0 * symm.xsym - be.xm;
                float yval = 2.0 * symm.ysym - be.ym;
               // write_data_to_file(xval, yval, sig1, sig2, bet, sig12, set, disp, pi + zet, mm, buffer);
                Stress st(xval, yval, 0, 0, 0, 0, 0, disp, pi + zet, mm);   //new stress input defined
                //Stress st(xval, yval,sig1, sig2, bet, sig12, set, disp, pi + zet, mm);
                stress[npoint] = st;
                npoint++;

                if (be.kod == 5)
                {
                  // write_data_to_file(xval, yval, sig1, sig2, bet, sig12, set, disp1, pi + zet1, mm, buffer);
                    Stress st(xval, yval, 0, 0, 0, 0, 0, disp1, pi + zet1, mm);   //new stress input defined
                    //Stress st(xval, yval, sig1, sig2, bet, sig12, set, disp1, pi + zet1, mm);
                    stress[npoint] = st;
                    npoint++;
                }
            }
        }
    }
    return;
}





void  for_j_loop(int mm, float& sigxx, float& sigyy, float& sigxy, float& ux, float& uy, float xp, float yp)
{
    
    for (int j = 0; j < numbe; ++j)
    {
        BoundaryElement& be = elm_list[j];

        if (mm != be.mat_no)
            continue;
        int js = 2 * j;
        int jn = js + 1;
                  
        s2us.reset();            //reset function is initl();
        float xj = be.xm;
        float yj = be.ym;
        float aj = be.a;

        float cosbj = be.cosbet;
        float sinbj = be.sinbet;
        coeff(xp, yp, xj, yj, aj, cosbj, sinbj, +1, mm);

        switch (symm.ksym)
        {
        case 0:
            break;
        case 1:
            xj = 2.0 * symm.xsym - be.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            break;
        case 2:
            yj = 2.0 * symm.ysym - be.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            break;
        case 3:
            xj = 2.0 * symm.xsym - be.xm;
            yj = 2.0 * symm.ysym - be.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, -sinbj, +1, mm);
            break;
        case 4:
            xj = 2.0 * symm.xsym - be.xm;
            coeff(xp, yp, xj, yj, aj, cosbj, -sinbj, -1, mm);
            xj = be.xm;
            yj = 2.0 * symm.ysym - be.ym;
            coeff(xp, yp, xj, yj, aj, -cosbj, sinbj, -1, mm);
            xj = 2.0 * symm.xsym - be.xm;
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





void compute_stress_displ_at_specified_points(int& npoint, stringstream& buffer)
{
    /* compute displacements and stresses at specified points in body.
    */

    int material = 0, mm = 0;
    float sig1 = 0.0, sig2 = 0.0, sig12 = 0.0;
    float  bet = 0.0, set = 0.0;
   
    float xp, yp, y0, pxx = 0.0, pxy = 0.0, pyy = 0.0;
    float  ux = 0.0, uy = 0.0;
    float sigxx = 0.0, sigxy = 0.0, sigyy = 0.0;
    float disp = 0.0, zet = 0.0;
    int n_valid;      

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
        y0 = g.y_surf;    
        pxx = symm.pxx1 + g.skx * (y0 - yp);
        pyy = symm.pyy1 + g.sky * (y0 - yp);
        pxy = symm.pxy1;       

        mm = j_material;
        if(multi_region)
           mm = check_material_id(xp, yp);

        ux = 0.0, uy = 0.0, sigxx = pxx, sigyy = pyy, sigxy = pxy;

        if (mm == mat_lining)
        {
            sigxx = 0;
            sigyy = 0;
            sigxy = 0;
        }
        for_j_loop(mm, sigxx, sigyy, sigxy, ux, uy, xp, yp);
        bet = (sigxx == sigyy) ? pi / 2.0 : 0.5 * atanf(2.0 * sigxy / (sigxx - sigyy));

        float cosf1 = cosf(bet);
        float cosf2 = cosf1 * cosf1;
        float sinf1 = sinf(bet);
        float sinf2 = sinf1 * sinf1;
        sig1 = sigxx * cosf2 + 2 * sigxy * sinf1 * cosf1 + sigyy * sinf2;
        sig2 = sigxx * sinf2 - 2 * sigxy * sinf1 * cosf1 + sigyy * cosf2;        

        sig12 = (sig1 - sig2) / 2.0;
        set = (sig12 > 0) ? bet + pi / 4.0 : bet - pi / 4.0;
        disp = sqrt(ux * ux + uy * uy);
        zet = ux > 0. ? atanf(uy / ux) :
            (ux < 0. ? pi + atanf(uy / ux) :
                1.57 * copysign(1, (int)uy));    //isignSara!
        /*if ((sig1 == 0 || sig2 == 0) && (sig12 == 0) && (disp == 0))
        {
            ++it;
            continue;
        }*/
        Stress st(xp, yp, sig1, sig2, bet, sig12, set, disp, zet, mm);
        stress[npoint++] = st;
        ++it;         
    }

   // win_exchange.w_npoints = npoint;
   // write_data_to_file2(0, npoint, buffer);
    compute_disp_boundary_and_frac_surfaces(npoint, sig1, sig2, bet, sig12, set, mm, buffer);
    return;
}






void final_stress_assignment(int & npoint, float x, float y, float sig1, float sig2, float bet, float sig12, float set,
    float disp, float zet, int mm, stringstream& buffer)
{

    //write_data_to_file(x, y, sig1, sig2, bet, sig12, set, disp, zet, mm, buffer);
   
    Stress st(x, y, sig1, sig2, bet, sig12, set, disp, zet,mm);   //new stress input defined
    stress[npoint] = st;
    npoint++;

      if (symm.ksym == 1 || symm.ksym == 4)
        {
            float xval = 2.0 * symm.xsym - x;
           // write_data_to_file(xval, y, sig1, sig2, bet, sig12, set, disp, -zet, mm, buffer);

            Stress st(xval, y, sig1, sig2, pi - bet, sig12, pi - set, disp, pi - zet ,mm);   //new stress input defined
            stress[npoint] = st;
            npoint++;

        }

        if (symm.ksym == 2 || symm. ksym == 4)
        {
            float yval = 2.0 * symm.ysym - y;
           // write_data_to_file(x, yval, sig1, sig2, bet, sig12, set, disp, -zet, mm, buffer);

            Stress st(x, yval, sig1, sig2, -bet, sig12, -set, disp, -zet, mm);   //new stress input defined
            stress[npoint] = st;
            npoint++;

        }
        if (symm.ksym == 3 || symm.ksym == 4)
        {
            float xval = 2.0 * symm.xsym - x;
            float yval = 2.0 * symm.ysym - y;
          //  write_data_to_file(xval, yval, sig1, sig2, bet, sig12, set, disp, pi + zet, mm, buffer);

            Stress st(xval, yval, sig1, sig2, pi + bet, sig12, pi + set, disp, pi + zet, mm);   //new stress input defined
            stress[npoint] = st;
            npoint++;
        }
   
    return;
}





void reassigning_boundary_values2(int& npoint, int ID, int m, int j, int k, float beta, float x,
    float y, stringstream& buffer)
{    
    if (j == numbe)
       j--;
    BE& be = b_elm[m];  //be is an alias for boundary elemnt m
    BE& bj = b_elm[j];  //bj is an alias for boundary elemnt j


    BoundaryElement& be1 = elm_list[m];
    BoundaryElement& bj1 = elm_list[j];


    float ss_m = be.sigma_s;
    float sn_m = be.sigma_n;
    float ustem_m = be.us;

    float untem_m = be.un;
    float us_m = be.us_neg; 
    float un_m = be.un_neg;


    float ss_n = bj.sigma_s; 
    float sn_n = bj.sigma_n; 
    float ustem_n = bj.us;

    float untem_n = bj.un;
    float us_n = bj.us_neg;
    float un_n = bj.un_neg;

    float us_a =0, un_a = 0, ang = 0;

    if (be1.kod > 10)
    {// correction for gravity and introduced Ks and Kn
        ss_m = be.sigma_s + s4.d0[m * 2] * (rock1[1].e / 1e4);    
        sn_m = be.sigma_n + s4.d0[m * 2 + 1] * (rock1[1].e / 1e4);
    }

    if (bj1.kod > 10)
    {
        ss_n = bj.sigma_s + s4.d0[j * 2] * (rock1[1].e / 1e4);
        sn_n = bj.sigma_n + s4.d0[j * 2 +1] * (rock1[1].e / 1e4);
    }

    if (ID == 1 || ID == 2)
    {
        ss_n = -ss_n;
        us_n = -us_n;
    }

    if (k == 1)
    {  //for boundary element m
        us_a = (be1.a * (us_n * cosf(beta/2) - un_n * sinf(beta/2)) + bj1.a *
            (us_m * cosf(beta/2) + un_m * sinf(beta / 2))) / (be1.a + bj1.a);

        un_a = (be1.a * (us_n * sinf(beta/2) + un_n * cosf(beta/2)) + bj1.a *
            (-us_m * sinf(beta/2) + un_m * cosf(beta/2))) / (be1.a + bj1.a);

        ang = atan2f(be1.sinbet, be1.cosbet) + beta / 2.0;
    }

    else //if (k == 2)
    {  // //for boundary element j
        us_a = (bj1.a * (us_m * cosf(beta/2) - un_m * sinf(beta/2)) + bj1.a *
            (us_n * cos(beta/2) + un_n * sinf(beta / 2))) / (be1.a + bj1.a);

        un_a = (bj1.a * (us_m * sinf(beta/2) + un_m * cosf(beta/2)) + bj1.a *
            (-us_n * sinf(beta/2) + un_n * cosf(beta/2))) / (be1.a + bj1.a);

        ang = atan2f(be1.sinbet, be1.cosbet) - beta / 2.0;
    }


    float cosb = cosf(ang);
    float sinb = sinf(ang);
    float ss, sn;

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - y);
    float pyy = symm.pyy1 + g.sky * (y0 - y);
    float pxy = symm.pxy1;

    float st0 = pxx * cosb * cosb - 2.0 * pxy * sinb * cosb + pyy * sinb * sinb;
    float sn0 = pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb;

    int mm = bj1.mat_no;
    float st = 0;
    if (k == 1)
    {
        ss = (ss_m + ss_n) / 2.0 + (sn_m - sn_n) / 2.0 * tanf(beta/2);
        sn = (sn_m + sn_n) / 2.0 + (-ss_m + ss_n) / 2.0 * tanf(beta/2);

        st = rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr) *
            ((us_n * cosf(beta/2) - un_n * sinf(beta/2)) - (us_m * cosf(beta/2) +
                un_m * sinf(beta/2))) / ((be1.a + bj1.a) * cosf(beta/2));
    }
    else
    {
        ss = (ss_n + ss_m) / 2.0 + (sn_n - sn_m) / 2.0 * tanf(beta/2);
        sn = (sn_n + sn_m) / 2.0 + (-ss_n + ss_m) / 2.0 * tanf(beta/2);

        st = rock1[mm].e / (1 - rock1[mm].pr * rock1[mm].pr) * ((us_m * cosf(beta/2) - un_m *
            sinf(beta/2)) - (us_n * cosf(beta/2) + un_n *
                sinf(beta / 2))) / ((be1.a + bj1.a) * cosf(beta/2));
    }

    st += st0  + (sn - sn0) * rock1[mm].pr / (1 - rock1[mm].pr);

    float bet;
    bet = (sn == st) ? pi / 2.0 : 0.5 * atanf(2.0 * ss / (st - sn));


    //float sig1 = st * cosf(bet) * cosf(bet) + 2 * ss * sinf(bet) * cosf(bet) + sn * sinf(bet) * sinf(bet);
    //float sig2 = st * sinf(bet) * sinf(bet) - 2 * ss * sinf(bet) * cosf(bet) + sn * cosf(bet) * cosf(bet);
    float sin_bet = std::sin(bet);
    float cos_bet = std::cos(bet);
    float sin2 = sin_bet * sin_bet;  // sin(bet)^2
    float cos2 = cos_bet * cos_bet;  // cos(bet)^2
    float sin_cos = sin_bet * cos_bet;  // sin(bet) * cos(bet)

    float sig1 = std::fma(st, cos2, std::fma(2.0f * ss, sin_cos, sn * sin2));
    float sig2 = std::fma(st, sin2, std::fma(-2.0f * ss, sin_cos, sn * cos2));
    float sig12 = (sig1 - sig2) / 2.0;


    bet = (k == 1) ?atan2f(be1.sinbet, be1.cosbet)+ beta / 2.0 + bet :
        atan2f(be1.sinbet, be1.cosbet) - beta / 2.0 + bet; 

    float set = (sig12 > 0) ? bet + pi / 4.0 : bet - pi / 4.0;
   
    float alpha = (k == 1) ? atan2f(be1.sinbet, be1.cosbet) + beta / 2.0 :
       atan2f(be1.sinbet, be1.cosbet) - beta / 2.0;

    float ux = us_a * cosf(alpha) - un_a * sinf(alpha);
    float uy = us_a * sinf(alpha) + un_a * cosf(alpha);
   /* if (abs(ux) < 1e-7f)
        ux = 0;
    if (abs(uy) < 1e-7f)
        uy = 0;*/
    float disp = sqrt(ux * ux + uy * uy);
    float zet = atan2f(uy, ux);

    final_stress_assignment(npoint, x, y, sig1, sig2, bet, sig12, set, disp, zet, mm, buffer);
    return;
}






void compute_stress_on_boundary_surfaces(int& npoint, stringstream& buffer)
{
    std::vector<float> xtem(numbe,0), ytem(numbe,0);
    std::pair<bool, bool> jmp = make_pair(false,false);
   
    float x = 0;
    float y = 0;
    float beta = 0;
    int ID = 0;
    bool jmp_1500 = false;
    float st = 0,bet = 0;
    int jj = 0;
    float xt, yt;
    int index = npoint;
    for (int m = 0; m < numbe; ++m) 
    {
        BoundaryElement& be = elm_list[m];
        if (be.kod == 7 || be.kod == 5)
            continue;             

        if (dispwin.ID_win == 1)
        {
            if (be.xm > dispwin.xur || be.xm < dispwin.xll || be.ym > dispwin.yur ||
                be.ym < dispwin.yll)
                continue;
        }
        else // (ID_win == 2) 
        {
            if (sqrt(pow(be.xm - dispwin.xc0, 2) + pow(be.ym - dispwin.yc0, 2)) > dispwin.radium)
                continue;
        }       

        float xbeg = be.xm - be.a * be.cosbet;
        float ybeg = be.ym - be.a * be.sinbet;
        float xend = be.xm + be.a * be.cosbet;
        float yend = be.ym + be.a * be.sinbet;

        for (int k = 1; k <= 2; ++k) 
        {
            jmp_1500 = false;
            
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
            
            jmp = for_j_loop(k, m, xt, yt, x, y, ID, beta,jj,numbe,true);
            //if the label jmp to 1500 gets true need to continue with the next iterate
            if (jmp.first == true)
                continue;

            if (k == 1)
            {
                xtem[m] = x;                                                                    
                ytem[m] = y;
            }
            // Avoid repeating one point                                                 
            else  //if (k == 2)
            {
                jmp_1500 = false;
                for (int i = 0; i < numbe; ++i)
                {
                    if (std::abs(x - xtem[i]) < s15.aa / 1000 && (std::abs(y - ytem[i]) < s15.aa / 1000.0))
                    {
                        jmp_1500 = true;
                        break;
                    }
                }
                if (jmp_1500 == true)
                    continue;   
            }                                                             
            //if fracture intersects this point, skip the point             
             bool jmp1500 = for_i_loop(k, x, y, m, numbe);
            if (jmp1500 == true)
                continue;                        
            reassigning_boundary_values2(npoint, ID, m, jj, k, beta, x, y, buffer);
             
        }                                                   
    }
    return;
}






void save_buffer_to_file(ofstream& file4, stringstream& buffer)
{
    file4 << buffer.str();
}





void internal(int id , int& npoint)
{
    /* internal grid point stresses and displacements */    
    wstring filename = stress_dir + L"/Stress" + std::to_wstring(mcyc) + L".dat";
    std::ofstream file4(filename);
    auto old_flags = file4.flags();
    file4.setf(ios::fixed, ios::floatfield);  // Fixed-point format for floats
    stringstream buffer;

    buffer << "xp        yp        sxx        syy        sxy      disp      bet     set     zet      matregion" << std::endl;
    // compute displacements and stresses at specified points in body.
    npoint = 0;     
    compute_stress_displ_at_specified_points(npoint, buffer);    
    compute_stress_on_boundary_surfaces(npoint, buffer);
    int perc = 2;
    buffer.precision(perc);  
    
        for (size_t i = 0; i < npoint; ++i) 
        {
            Stress& st = stress[i];
            buffer << fixed << setprecision(perc+1) << setw(7) << st.w_xp << "  "
                << setw(7) << st.w_yp << "  "
                << setw(10) << scientific << setprecision(2) << st.w_sig1 << "  "
                << setw(10) << scientific << setprecision(2) << st.w_sig2 << "  "
                << setw(10) << scientific << setprecision(2) << st.w_sig12 << "  "
                << scientific << setprecision(2) << st.w_disp << "   "
                << setprecision(3) << fixed << st.w_bet << "   "               
                <<fixed << setprecision(perc+1) << st.w_set << "  "               
                << scientific << setprecision(perc) << st.w_zet << "   "
                << st.w_mat << endl;          
        }  
    save_buffer_to_file(file4, buffer);    
    buffer.flags(old_flags);
    file4.close();    
    return;
}