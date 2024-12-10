#include <common.h>
#include <CommonPara.h>
#include <WinInterface.h>


using namespace WinInterface_h::winvar;
using namespace CommonPara_h::comvar;


 bool check_symmetry_so( float xt1, float yt1, float xt2,
    float yt2, float psize,float xp, float yp)
{
    float jmp_100 = false;

    float xmin0 = min(xt1, xt2);
    float xmax0 = max(xt1, xt2);
    float ymin0 = min(yt1, yt2);
    float ymax0 = max(yt1, yt2);

    if (ymax0 < yp - psize / 2.0 || ymin0 > yp + psize / 2.0 ||
        xmax0 < xp - psize / 2.0 || xmin0 > xp + psize / 2.0 )
    {       
        return  jmp_100;
    }

    jmp_100 = true;

    return jmp_100;
}





 void label_100(int m,float psize,float& xlength, float& ylength, float & rlength, 
     float& tlength,float& px1, float& px2_reverse, float& py1, float& py2_reverse,
     float& pr1, float& pr2_reverse,float& pt1, float& pt2_reverse)
 {
     float dist0 = 1e-10;
     float cosr, sinr, sindiff, cosdiff, dist, prr, pt;
     float ge = 9.81;

     xlength += 2.0 * elm_list[m].a * std::abs(elm_list[m].cosbet);
     ylength += 2.0 * elm_list[m].a * std::abs(elm_list[m].sinbet);   

     float aperture = joint[m].aperture0 - s4.d0[m * 2];
     if (aperture <= joint[m].aperture_r)
         aperture = joint[m].aperture_r;

     //seems another function could be defined for following part  Sara!

     float px = pow(aperture, 3) / (12.0 * perm.viscosity * psize) *
         abs(elm_list[m].cosbet) * ge * perm.density;
     if (px == 0) 
         px = 1e-20;
     px1 += 2.0 * elm_list[m].a * std::abs(elm_list[m].cosbet) * px;
     px2_reverse += 2.0 * elm_list[m].a * std::abs(elm_list[m].cosbet) / px;


     float py = std::pow(aperture, 3) / (12.0 * perm.viscosity * psize)
         * std::abs(elm_list[m].sinbet) * ge * perm.density;
     if (py == 0)
         py = 1e-20;

     py1 += 2.0 * elm_list[m].a * std::abs(elm_list[m].cosbet) * py;
     py2_reverse += 2.0 * elm_list[m].a * std::abs(elm_list[m].sinbet) / py;

     //int ntunnel = ntunnel;
     for (int l = 0; l < ntunnel; ++l) 
     {
         dist = std::sqrt(pow((elm_list[m].xm - tunnl.x_cent[l]), 2) +
             pow((elm_list[m].ym - tunnl.y_cent[l]), 2));

         if (dist > dist0)
         {
             dist0 = dist;
             cosr = (elm_list[m].xm - tunnl.x_cent[l]) / dist;
             sinr = (elm_list[m].ym - tunnl.y_cent[l]) / dist;
         }
     }

     if (ntunnel == 0) 
     {
         dist = sqrt(elm_list[m].xm * elm_list[m].xm + elm_list[m].ym * elm_list[m].ym);
         cosr = elm_list[m].xm / dist;
         sinr = elm_list[m].ym / dist;
     }

     sindiff = elm_list[m].sinbet * cosr - elm_list[m].cosbet * sinr;
     cosdiff = elm_list[m].cosbet * cosr + elm_list[m].sinbet * sinr;
     rlength += 2.0 * elm_list[m].a * std::abs(cosdiff);
     tlength += 2.0 * elm_list[m].a * std::abs(sindiff);

     prr = std::pow(aperture, 3) / (12.0 * perm.viscosity * psize) *
         std::abs(cosdiff) * ge * perm.density;
     if (prr == 0)
         prr = 1e-20;
     pr1 += 2.0 * elm_list[m].a * std::abs(cosdiff) * prr;
     pr2_reverse += 2.0 * elm_list[m].a * std::abs(cosdiff) / prr;

     pt = std::pow(aperture, 3) / (12.0 * perm.viscosity * psize) *
         std::abs(sindiff) * ge * perm.density;

     if (pt == 0)
         pt = 1e-20;
     pt1 += 2.0 * elm_list[m].a * std::abs(sindiff) * pt;
     pt2_reverse += 2.0 * elm_list[m].a * std::abs(sindiff) / pt;
     return;
 }



 float  comput_rx_ry_perm(float rx, float& px1, float xlength,
     float& px2_reverse)
 {
     /* since this piece of code is 4 times repeated in calling function,
     it's written as a seprate function,
     rx could be ry, rr, rt and for px1 once it's px1, py1,pr1 and pt1
     for px2 and px2_reverse as well, each time with different
     args and values is called*/


     float permx, px2 = 0;

     px1 = px1 / xlength;
     if (px2_reverse == 0)
         px2 = 0;
     else
         px2 = 1 / px2_reverse * xlength;

     float permx_tem = std::sqrt(px1 * px2) + perm.perm0;
     if (rx <= 1.0)
         permx = pow(perm.perm0, (1 - rx)) * pow((permx_tem + perm.perm0) , rx);
     else if (rx > 1.0 && rx < 2.0)
         permx = permx_tem + pow(perm.perm0, (1 - (rx - 1))) *
         pow((permx_tem + perm.perm0), (rx - 1));
     else
         permx = 2.0 * permx_tem + perm.perm0;

     return permx;
 }

 void frm_lbl80_890(float xlength, float ylength, float rlength, float tlength,
     float& px1, float& px2_reverse, float& py1, float& py2_reverse,
     float& pr1, float& pr2_reverse, float& pt1, float& pt2_reverse,
     float xp, float yp, float psize, int np)
 {
     float rx = xlength / psize;
     float ry = ylength / psize;
     float rr = rlength / psize;
     float rt = tlength / psize;

     //float px2 = 0, py2 = 0, pr2 = 0, pt2 = 0;
     float permx = comput_rx_ry_perm(rx, px1, xlength, px2_reverse);
     float permy = comput_rx_ry_perm(ry, py1, ylength, py2_reverse);
     float permr = comput_rx_ry_perm(rr, pr1, rlength, pr2_reverse);
     float permt = comput_rx_ry_perm(rt, pt1, tlength, pt2_reverse);



     // Assuming "np" is a valid index
     permeability[np].w_xp = xp;
     permeability[np].w_yp = yp;
     permeability[np].w_permx = permx;
     permeability[np].w_permy = permy;
     permeability[np].w_permr = permr;
     permeability[np].w_permseta = permt;
 }


 void resetting_parameters(float& px1, float& py1, float& pr1, float& pt1, float& px2_reverse, 
     float& py2_reverse, float& pr2_reverse,float& pt2_reverse, float& xlength, float& ylength,
     float& rlength, float& tlength)
 {
     px1 = 0;
     py1 = 0;
     px2_reverse = 0;
     py2_reverse = 0;
     pr1 = 0;
     pt1 = 0;
     pr2_reverse = 0;
     pt2_reverse = 0;

     xlength = 1e-15;
     ylength = 1e-15;
     rlength = 1e-15;
     tlength = 1e-15;

     return;
 }




 void chk_no_frac_element_in_each_cell(float xp, float yp, float psize, float& px1, float& py1,
     float& pr1, float& pt1, float& px2_reverse, float& py2_reverse, float& pr2_reverse,
     float& pt2_reverse, float& xlength, float& ylength,float& rlength, float& tlength, int np)
 {
     /*check how many fracture element in each cell*/

     for (int m = 0; m < numbe; ++m)
     {
         if (elm_list[m].kod != 5) continue;        // It includes the boundary elements here  

         float xt1 = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
         float yt1 = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;

         float xt0 = elm_list[m].xm;
         float yt0 = elm_list[m].ym;

         float xt2 = elm_list[m].xm + elm_list[m].a * elm_list[m].cosbet;
         float yt2 = elm_list[m].ym + elm_list[m].a * elm_list[m].sinbet;

         float xmin0 = min(xt0, xt2);
         float xmax0 = max(xt0, xt2);
         float ymin0 = min(yt0, yt2);
         float ymax0 = max(yt0, yt2);

         float xt12 = 0, yt12 = 0, xt22 = 0, yt22 = 0;    //used these temp vars to prevent recomputing
         bool jmp_100 = false;

         if (ymax0 < yp - psize / 2.0 || ymin0 > yp + psize / 2.0 ||
             xmax0 < xp - psize / 2.0 || xmin0 > xp + psize / 2.0)
         {
             if (symm.ksym == 1 || symm.ksym == 4)
             {


                 xt12 = 2.0 * symm.xsym - xt1;
                 yt12 = yt1;
                 //xt0 = 2.0 * symm.xsym - elm_list[m].xm;
                 //yt0 = elm_list[m].ym;
                 xt22 = 2.0 * symm.xsym - xt2;
                 yt22 = yt2;
                 jmp_100 = check_symmetry_so(xt12, yt12, xt22, yt22, psize, xp, yp);


                 if (jmp_100 == true) break;
             }

             if (symm.ksym == 2 || symm.ksym == 4)
             {
                 xt12 = xt1;
                 yt12 = 2.0 * symm.ysym - yt1;
                 //xt0 = elm_list[m].xm;
                // yt0 = 2.0 * symm.ysym - elm_list[m].ym;
                 xt22 = xt2;
                 yt22 = 2.0 * symm.ysym - yt2;
                 jmp_100 = check_symmetry_so(xt12, yt12, xt22, yt22, psize, xp, yp);
                 if (jmp_100 == true) break;

             }

             if (symm.ksym == 3 || symm.ksym == 4)
             {
                 xt1 = 2.0 * symm.xsym - xt1;
                 yt12 = 2.0 * symm.ysym - yt1;
                 //xt0 = 2.0 * symm.xsym - elm_list[m].xm;
                 //yt0 = 2.0 * symm.ysym - elm_list[m].ym;
                 xt22 = 2.0 * symm.xsym - xt2;
                 yt22 = 2.0 * symm.ysym - yt2;
                 if (jmp_100 == true) break;

             }
             continue;
         }  //end if

         label_100(m, psize, xlength, ylength, rlength, tlength,px1, px2_reverse,
             py1, py2_reverse, pr1, pr2_reverse, pt1, pt2_reverse);

     }  //end for
     frm_lbl80_890(xlength, ylength, rlength, tlength,
         px1, px2_reverse, py1, py2_reverse,
         pr1, pr2_reverse, pt1, pt2_reverse,
         xp, yp, psize, np);
    
     return;
 }

 




 void frm_lbl_890_900(float psize, int& npointp)
 {
     int kk = 0;
     float xp = 0, yp = 0;
     int np;

     for (int i = 0; i < numbe; ++i) 
     {
         if (symm.ksym == 0) kk = 1;
         else
         if (symm.ksym == 1 || symm.ksym == 2 || symm.ksym == 3) kk = 2;
         else
         if (symm.ksym == 4) kk = 4;

         for (int k = 1; k <= kk; ++k)
         {
             switch (k)
             {
             case(1):
                 xp = elm_list[i].xm;
                 yp = elm_list[i].ym;
                 break;
             case(2):
                 switch(symm.ksym)
                 {
                 case(1):
                     xp = 2.0 * symm.xsym - elm_list[i].xm;
                     yp = elm_list[i].ym;
                     break;

                 case(2):
                     xp = elm_list[i].xm;
                     yp = 2.0 * symm.ysym - elm_list[i].ym;
                     break;

                 case(3):
                     xp = 2.0 * symm.xsym - elm_list[i].xm;
                     yp = 2.0 * symm.ysym - elm_list[i].ym;
                     break;

                 case(4):
                     xp = 2.0 * symm.xsym - elm_list[i].xm;
                     yp = elm_list[i].ym;
                     break;
                 
                 }
             case(3):
                 if (symm.ksym == 4)
                 {
                     xp = 2.0 * symm.xsym - elm_list[i].xm;
                     yp = 2.0 * symm.ysym - elm_list[i].ym;
                 }
                 break;

             case(4):
                 if (symm.ksym == 4)
                 {
                 xp = elm_list[i].xm;
                 yp = 2.0 * symm.ysym - elm_list[i].ym;
                 }
                 break;
             }

             npointp = npointp + 1;
             np = npointp;

             float px1 = 0, py1 = 0, px2_reverse = 0, py2_reverse = 0, pr1 = 0, pt1 = 0,
                 pr2_reverse = 0, pt2_reverse = 0, xlength = 1e-15, ylength = 1e-15,
                 rlength = 1e-15, tlength = 1e-15;

             //not sure which verison is most effeicent ina sprate function or here Sara!
             /*resetting_parameters(  px1,   py1,   pr1,   pt1,   px2_reverse,
                   py2_reverse,   pr2_reverse,   pt2_reverse,   xlength,   ylength,
                   rlength,   tlength);*/

             float xp1 = xp - psize / 2.0;
             float xp2 = xp + psize / 2.0;


             chk_no_frac_element_in_each_cell( xp,  yp,  psize, px1,  py1,
                  pr1,  pt1,  px2_reverse,  py2_reverse,  pr2_reverse,
                  pt2_reverse,  xlength,  ylength,  rlength,  tlength, np);

         }
     }
 }




void permeability_plot(int id, int& npointp) 
{
    npointp = 0;
    float dx = 0, dy = 0, ddr = 0, dda = 0;
    int num1 = 0, num2 = 0;

    if (dispwin.ID_win == 1)
    {
        dx = (dispwin.xur - dispwin.xll) / dispwin.numx;
        dy = (dispwin.yur - dispwin.yll) / dispwin.numy;
        num1 = dispwin.numx;
        num2 = dispwin.numy;
    }

    else //if (ID_win == 2) 
    {
        ddr = dispwin.radium / dispwin.numr;
        dda = 2 * pi / dispwin.numa;
        num1 = dispwin.numr;
        num2 = dispwin.numa;
    }

    // Permeability cell size
    float psize = 2.1 * max(0.5 * (dx + dy), ddr);
    bool jmp_flag = false;   //used for controling the loop instead of jmp to 890 in the org code

    for (int i1 = 0; i1 <= num1; ++i1)
    {
        for (int i2 = 0; i2 <= num2; ++i2)
        {
            float xp=0, yp=0 ;

            if (dispwin.ID_win == 1)
            {
                xp = dispwin.xll + i1 * dx;
                yp = dispwin.yll + i2 * dy;
            }

            else // (ID_win == 2)
            {
                xp = dispwin.xc0 + i1 * ddr * cosf(i2 * dda);
                yp = dispwin.yc0 + i1 * ddr * sinf(i2 * dda);
            }           

            for (int l = 0; l < ntunnel; ++l)
            {
                float x0 = tunnl.x_cent[l];
                float y0 = tunnl.y_cent[l];
                float dd0 = tunnl.diameter[l];

                if (sqrt((xp - x0) * (xp - x0) + (yp - y0) * (yp - y0)) < 0.5 * dd0)

                {
                    jmp_flag = true;
                    break;
                }
            }

            if (jmp_flag) continue;                //instead of jmp to 890 which is end of loop

            npointp++;
            npli = npointp;
            //resetting_parameters(px1, py1, pr1, pt1, px2_reverse,
            //    py2_reverse, pr2_reverse, pt2_reverse, xlength, ylength,
            //    //rlength, tlength);

            float px1 = 0, py1 = 0, px2_reverse = 0, py2_reverse = 0, pr1 = 0, pt1 = 0,
                pr2_reverse = 0, pt2_reverse = 0, xlength = 1e-15, ylength = 1e-15,
                rlength = 1e-15, tlength = 1e-15;

            // Check how many fracture elements are in each cell
            float xp1 = xp - psize / 2.0;
            float xp2 = xp + psize / 2.0;

            chk_no_frac_element_in_each_cell(xp, yp, psize, px1, py1,
                pr1, pt1, px2_reverse, py2_reverse, pr2_reverse,
                pt2_reverse, xlength, ylength, rlength, tlength, npli);
            
        }
    }
    frm_lbl_890_900( psize, npointp);
            
            //Sara!   data should be written in file
             /* write(7) npointp
                do 99 i = 1, npointp
                    write(7) permeability(i) % w_xp, permeability(i) % w_yp, permeability(i) % w_permx, permeability(i) % w_permy, &
                    permeability(i) % w_permr, permeability(i) % w_permseta
                    99 continue

                    WRITE(57, 250)
                    250 FORMAT('  xp,       yp,      x-perm (m/s),   y-perm (m/s),   r-perm (m/s),    tan-perm (m/s) ')
                    do 199 i = 1, npointp
                        write(57, 299) permeability(i) % w_xp, permeability(i) % w_yp, permeability(i) % w_permx, permeability(i) % w_permy, &
                        permeability(i) % w_permr, permeability(i) % w_permseta
                        299 format(2f10.4, 1x, 4e16.4)
                        199 continue
                        REWIND(57)*/

            win_exchange.w_npointp = npointp;

        }
    

