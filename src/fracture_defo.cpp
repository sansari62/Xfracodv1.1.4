#include<stdafx.h>
#include <CommonPara.h>
#include <fracture_defo.h>
#include <Failure.h>

using namespace CommonPara_h::comvar;



void write_to_file(stringstream& buffer, int jpoint)
{
    int index = 0;
    while (index < jpoint)
    {
       // ////old one 
       // buffer << std::setw(6) << std::fixed << std::setprecision(3) << (wjoint[index].w_xp + wjoint[index + 1].w_xp) / 2 << std::setw(9) <<
       //     std::fixed << std::setprecision(3) << (wjoint[index].w_yp + wjoint[index + 1].w_yp) / 2 << std::setw(13) << std::scientific <<
       //     std::setprecision(3) << wjoint[index].w_ds <<
       //     // std::setw(8) << std::fixed << std::setprecision(3) <<
       //     // wjoint[index].w_bet << 
       //     std::setw(13) << std::scientific << std::setprecision(3) << wjoint[index].w_dn <<
       //     //std::setw(8) << std::fixed << std::setprecision(3) << wjoint[index].w_set << 
       //     std::setw(13) << std::scientific << std::setprecision(3) << wjoint[index].w_aperture << endl;
       // //<< std::setw(8) << std::fixed << std::setprecision(3) <<
       //// wjoint[index].w_zet << endl;
       // ///////
        int m = wjoint[index].m_indx;
        //float tx = (wjoint[index].w_xp + wjoint[index + 1].w_xp) / 2;
        //float ty = (wjoint[index].w_yp + wjoint[index + 1].w_yp) / 2;
        float ta= elm_list[m].a ;
        float cosm = elm_list[m].cosbet;
        float sinm = elm_list[m].sinbet;
        float tx = elm_list[m].xm;
        float ty = elm_list[m].ym;
        float xb = tx - ta * cosm;
        float yb = ty - ta * sinm;
        float xe = tx + ta * cosm;
        float ye = ty + ta * sinm;
        //char deli = ';';
        string deli = ",";
        float aperture = joint[m].aperture0 - 2*s4.d0[2 * m + 1];
        if (aperture < joint[m].aperture_r)
            aperture = joint[m].aperture_r;

       /* buffer << std::setw(6) << std::fixed << std::setprecision(3) <<xb << deli << std::setw(9) <<
            std::fixed << std::setprecision(3) << yb << deli << std::setw(9) <<
            std::fixed << std::setprecision(3)<<xe << deli << std::setw(9) <<
            std::fixed << std::setprecision(3)<<ye << deli <<std::setw(13) << std::scientific <<
            std::setprecision(3) << 2* wjoint[index].w_ds << deli <<
           std::setw(13) << std::scientific << std::setprecision(3) << 2*wjoint[index].w_dn << deli <<
             std::setw(13) <<std::scientific << std::setprecision(3) << aperture << endl;*/


        buffer << std::fixed << std::setprecision(6) << xb << deli << std::fixed << std::setprecision(6) << yb << deli << 
            std::fixed << std::setprecision(6) << xe << deli << std::fixed <<
            std::setprecision(6) << ye << deli << std::scientific << std::setprecision(6) << 2 * wjoint[index].w_ds << deli <<
            std::scientific <<std::setprecision(4) << 2 * wjoint[index].w_dn << deli<< std::setprecision(4) <<aperture << endl;
           
        index+=2;
    }

}





Wjoint compute_join_attr_and_add_them(int m, int& jpoint, int round)
{
    float bet, xp = 0, yp = 0, ds = 0, dn = 0, set1 = 0, zet1 = 0, aperture = 0, bet1 = 0;

    if(round == 1)
        bet = static_cast<float> (atan2(elm_list[m].sinbet, elm_list[m].cosbet));
    else
        bet = static_cast<float>(atan2(elm_list[m].sinbet, elm_list[m].cosbet)) + pi;


     xp = elm_list[m].xm + 0.3 * elm_list[m].a * cosf(bet + pi / 2) + 0.2 * elm_list[m].a * cosf(bet);
     yp = elm_list[m].ym + 0.3 * elm_list[m].a * sinf(bet + pi / 2) + 0.2 * elm_list[m].a * sinf(bet);
     ds = abs(s4.d0[2 * m]);

     if (s4.d0[2 * m ] > 0)
         bet1 = bet + pi;
     else if (s4.d0[2 * m] < 0)
         bet1 = bet;


     dn = abs(s4.d0[2 * m + 1]);

    if (s4.d0[2 * m + 1] > 0)
        set1 = bet - pi / 2;
    else if (s4.d0[2 * m + 1] <  0)
        set1 = bet + pi / 2;


     aperture = joint[m].aperture0 - s4.d0[2 * m + 1];
    if (aperture < joint[m].aperture_r)
        aperture = joint[m].aperture_r;

    zet1 = bet + pi/2;           
    Wjoint jointex;
    jointex.assign_val(xp, yp, ds/2, bet1, dn/2, set1, aperture/2, zet1,m);
    wjoint[jpoint] = jointex;
    jpoint++;
           
    return jointex;
}





void frac_defo_Circ_win(int& jpoint) {

    Wjoint joint1, joint2, jointex;
    for (int m = 0; m < numbe; ++m)
    {
        if (elm_list[m].kod != 5) continue;

        if (sqrt(pow(elm_list[m].xm - dispwin.xc0, 2) +
            pow(elm_list[m].ym - dispwin.yc0, 2)) > dispwin.radium)
            continue;

        int mm = elm_list[m].mat_no;
        if(multi_region)
            mm = check_material_id(elm_list[m].xm, elm_list[m].ym);
        
        joint1 = compute_join_attr_and_add_them(m, jpoint, 1);    //return value for the first joint element

        joint2 = compute_join_attr_and_add_them(m, jpoint, 2);    //return value for the second joint element

        if (symm.ksym != 0)
        {
            if (symm.ksym == 1 || symm.ksym == 4)
            {

                jointex.assign_val(2. * symm.xsym - joint2.w_xp, joint2.w_yp, joint1.w_ds, pi - joint1.w_bet, joint1.w_dn,
                    pi - joint1.w_set, joint1.w_aperture, pi - joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;             


                jointex.assign_val(2. * symm.xsym - joint2.w_xp, joint2.w_yp, joint2.w_ds, pi - joint2.w_bet,
                    joint2.w_dn, pi - joint2.w_set, joint2.w_aperture, pi - joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                jointex.assign_val(joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint1.w_ds, -joint1.w_bet, joint1.w_dn,
                    -joint1.w_set, joint1.w_aperture, -joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;

                jointex.assign_val(joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint2.w_ds, -joint2.w_bet, joint2.w_dn,
                    -joint2.w_set, joint2.w_aperture, -joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }

            if (symm.ksym == 3 || symm.ksym == 4)
            {
                jointex.assign_val(2. * symm.xsym - joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint1.w_ds, pi + joint1.w_bet, joint1.w_dn,
                    pi + joint1.w_set, joint1.w_aperture, pi + joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;


                jointex.assign_val(2. * symm.xsym - joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint2.w_ds, pi + joint2.w_bet, joint2.w_dn,
                    pi + joint2.w_set, joint2.w_aperture, pi + joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }
        }

    }
}






void frac_defo_Rec_win(int& jpoint) {

    Wjoint joint1, joint2, jointex;
    for (int m = 0; m < numbe; ++m)
    {
        if (elm_list[m].kod != 5) continue;
        
        if (elm_list[m].xm > dispwin.xur || elm_list[m].xm < dispwin.xll ||
            elm_list[m].ym > dispwin.yur || elm_list[m].ym < dispwin.yll) continue;     
        
        int mm = elm_list[m].mat_no;
        if (multi_region)
            mm = check_material_id(elm_list[m].xm, elm_list[m].ym);
        joint1 = compute_join_attr_and_add_them(m, jpoint, 1);    // value for the first joint element

        joint2 = compute_join_attr_and_add_them(m, jpoint, 2);    // value for the second joint element

        if (symm.ksym != 0)
        {
            if (symm.ksym == 1 || symm.ksym == 4)
            {

                jointex.assign_val(2. * symm.xsym - joint2.w_xp, joint2.w_yp, joint1.w_ds, pi - joint1.w_bet, joint1.w_dn,
                    pi - joint1.w_set, joint1.w_aperture, pi - joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;    

                jointex.assign_val(2. * symm.xsym - joint2.w_xp, joint2.w_yp, joint2.w_ds, pi - joint2.w_bet,
                    joint2.w_dn, pi - joint2.w_set, joint2.w_aperture, pi - joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }

            if (symm.ksym == 2 || symm.ksym == 4)
            {
                jointex.assign_val(joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint1.w_ds, -joint1.w_bet, joint1.w_dn,
                    -joint1.w_set, joint1.w_aperture, -joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;

                jointex.assign_val(joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint2.w_ds, -joint2.w_bet, joint2.w_dn,
                    -joint2.w_set, joint2.w_aperture, -joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }

            if (symm.ksym == 3 || symm.ksym == 4)
            {
                jointex.assign_val(2. * symm.xsym - joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint1.w_ds, pi + joint1.w_bet, joint1.w_dn,
                    pi + joint1.w_set, joint1.w_aperture, pi + joint1.w_zet,m);

                wjoint[jpoint] = jointex;
                jpoint++;


                jointex.assign_val(2. * symm.xsym - joint2.w_xp, 2. * symm.ysym - joint2.w_yp, joint2.w_ds, pi + joint2.w_bet, joint2.w_dn,
                    pi + joint2.w_set, joint2.w_aperture, pi + joint2.w_zet,m);
                wjoint[jpoint] = jointex;
                jpoint++;

            }
        }

    }
}





void fracture_defo(int id, int& jpoint) {

    wstring filename = fd_dir + L"/Frac_deform" + std::to_wstring(state) + L".csv";
    std::ofstream outfile(filename);

    //outfile << "xp         yp         ds          dn          aperture \n";
   // outfile << "x1         y1       x2       y2       ds           dn          aperture \n";
    //outfile << "x1;y1;x2;y2;ds;dn;aperture\n";
    outfile << "x1"<<","<<"y1"<<","<<"x2"<< "," << "y2" << ","<<"ds" << "," << "dn" << "," << "aperture"<<"\n";
    
   
    jpoint = 0;
    stringstream buffer;
    //------- displ. on fracture surfaces----
    if (dispwin.ID_win == 1)
    {
        frac_defo_Rec_win(jpoint);
    }
    else
    {
        frac_defo_Circ_win(jpoint);
    }

    write_to_file(buffer,jpoint);
    outfile << buffer.str();    
    outfile.close();
    return;
}
