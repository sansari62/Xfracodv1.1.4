#include<stdafx.h>
#include <Mainb.h>
#include "CommonPara.h"

using namespace CommonPara_h::comvar;

//void check_row_for_nan(int row) {
//    for (int j = 0; j < 613; ++j) {
//        if (_isnan(comvar::s4.c[row][j])) {
//            printf("NaN at [%d][%d] = %g\n", row, j, comvar::s4.c[row][j]);
//        }
//    }
//}

bool check_diagonal_for_zero()
{
    int n = 2 * numbe;
    for (int j = 0; j < n; ++j) {
        if (fabs(comvar::s4.c[j][j]) < 1e-12) {
            printf("Near-zero at diagonal [%d][%d] = %.3e\n", j, j, comvar::s4.c[j][j]);
            comvar::s4.c[j][j] = 1e-5;
            return true;
        }
    }
    return false;
}


std::pair<int, int>  partial_pivoting(int j,int n) {
    int best = j;
    double max_val = fabs(s4.c[j][j]);
    for (int k = j + 2; k < n; k += 2)
    {
        if (fabs(s4.c[k][j]) > max_val)
        {
            max_val = fabs(s4.c[k][j]);
            best = k;
        }
    }
    if (max_val < 1e-12)
    {
        std::cerr << "Pivot failure: column " << j << " has near-zero values.\n";
        return make_pair(j, -1); // Skip or handle singularity
    }

    // Swap row pair (j,j+1) with (best,best+1)
    if (best != j)
    {
        //for (int col = 0; col < n; ++col)
        //{
            std::swap(s4.c[j], s4.c[best]);
            std::swap(s4.c[j + 1], s4.c[best + 1]);
        //}
        std::swap(s4.b[j], s4.b[best]);
        std::swap(s4.b[j + 1], s4.b[best + 1]);
        return make_pair(j, best);
    }
    return make_pair(j, -1);
}



void solve(int n, int mode)
{
    int nb = n - 1;    
    std::pair<int, int> swp;
    swp.first = -1;
    swp.second = -1;
    if (ni == 147)
    {
        cout << "c1erst in befgin of solve: c[613][613]= " << s4.c[613][613] << endl;
    }
    for (int j = 0; j < nb; ++j)
    {
        /*if (ni == 147)
        {
            if (j == 613 || j == 612)
                cout << "c= " << s4.c[613][613] << " j=" << j << endl;
        }*/
        for (int jj = j + 1; jj < n; ++jj)
        {
            if (s4.c[jj][j] == 0)
                continue;
            //if (fabs(s4.c[j][j]) <= 1e-12)
            //{
            //    s4.c[j][j] = 1e-5;
            //    /*printf("Near-zero at diagonal [%d][%d] = %.3e\n", j, j, comvar::s4.c[j][j]);
            //    swp = partial_pivoting(j, n);

            //    if (swp.second== -1) {
            //        cout << "partial pivoting is not possible" << n << j << endl;
            //        return;
            //    }*/
            //}         
            float xx = s4.c[jj][j] / s4.c[j][j];            
            for (int i = j; i < n; ++i)
            {
                /*if ( ni == 147)
                {
                    if (jj == 613)
                    {
                        if (i == 613)
                            cout << "BEFORE: s4.c[613][613] = " << s4.c[613][613]
                            << ", subtracting s4.c[" << j << "][613] * xx = "
                            << s4.c[j][613] << " * " << xx << " = " << s4.c[j][613] * xx << endl;
                    }
                }*/

                s4.c[jj][i] -= s4.c[j][i] * xx;

                /*if (ni == 147)
                {
                    if (jj == 613)
                    {
                        if (i == 613)
                            cout << "AFTER: s4.c[613][613] = " << s4.c[613][613] << endl;
                    }
                }*/
            }
            s4.b[jj] -= s4.b[j] * xx;
        }
    }
    /*if (swp.first != -1)
    {
        int k = swp.first;
        int best = swp.first;
        std::swap(s4.c[k], s4.c[best]);
        std::swap(s4.c[k + 1], s4.c[best + 1]);
        
        std::swap(s4.b[k], s4.b[best]);
        std::swap(s4.b[k + 1], s4.b[best + 1]);
    }*/
    if (fabs(s4.c[n - 1][n - 1]) <= 1e-12)
       {
            std::cerr << "Warning: small or zero c[" <<n-1 << "][" << n-1 << "] corrected to 1e-6\n";
            s4.c[n - 1][n - 1] = 1e-6;
    }
       
    s4.d[n - 1] = s4.b[n - 1] / s4.c[n - 1][n - 1];
    if (elm_list[int(n / 2) - 1].kod == 5 && mode == 0)
    {
        if (std::abs(s4.d[n - 1]) > d_max)
            s4.d[n - 1] = int(s4.d[n - 1] / abs(s4.d[n - 1])) * d_max;
    }

    for (int j = 1; j <= nb; ++j)
    {
        int jj = n - 1 - j;
        double sum = 0.0;
        int l = jj + 1;        

        for (int i = l; i < n; ++i)
        {
           //
            sum += s4.c[jj][i] * s4.d[i];
           // double a = s4.c[jj][i];
           // double b = s4.d[i];
            /*if (ni == 200) {
                double c = 0.0;
                if (!std::isfinite(a) || !std::isfinite(b))
                    std::cout << "Invalid value: c[" << jj << "][" << i << "] = " << a << ", d[" << i << "] = " << b << std::endl;
                 c += a * b;
                 cout << "c=" << c << endl;
            }*/
            //sum += a * b;
        }
        /*if (ni == 200) {
            if (s4.c[jj][jj] == 0)
                std::cout << "zero value: c[" << jj << "]" << endl;
        }*/

        s4.d[jj] = (s4.b[jj] - sum) / s4.c[jj][jj];


        if (elm_list[int(jj / 2)].kod == 5 && mode == 0)
        {
            if (std::abs(s4.d[jj]) > d_max)
                s4.d[jj] = int(s4.d[jj] / abs(s4.d[jj])) * d_max;
        }
    }

    return;
}





void coeff(float xi, float yi, float xj, float yj, float aj, float cosbj, float sinbj,
    int msym, int material)
{
    int mm = material;
    float con = 1.0 / (4. * pi * (1. - rock1[mm].pr));
    float cons = rock1[mm].e / (1. + rock1[mm].pr);
    float pr1 = 1. - 2. * rock1[mm].pr;
    float pr2 = 2. * (1. - rock1[mm].pr);

    float cos2b = cosbj * cosbj - sinbj * sinbj;
    float sin2b = 2. * sinbj * cosbj;
    float cosb2 = cosbj * cosbj;
    double sinb2 = sinbj * sinbj;

    float xb = (xi - xj) * cosbj + (yi - yj) * sinbj;
    float yb = -(xi - xj) * sinbj + (yi - yj) * cosbj;

    float r1s = (xb - aj) * (xb - aj) + yb * yb;
    float r2s = (xb + aj) * (xb + aj) + yb * yb;
    float fl1 = 0.5 * log(r1s);
    float fl2 = 0.5 * log(r2s);

    float fb2 = con * (fl1 - fl2);
    float fb3;

    if (yb == 0.0)
    {
        fb3 = 0.;
        if (abs(xb) < aj)
        {
            fb3 = con * pi;
        }
    }
    else
    {
        fb3 = -con * (atanf((xb + aj) / yb) - atanf((xb - aj) / yb));
    }

    float fb4 = con * (yb / r1s - yb / r2s);
    float fb5 = con * ((xb - aj) / r1s - (xb + aj) / r2s);

    float fb6 = con * (((xb - aj) * (xb - aj) - yb * yb) / (r1s * r1s) - ((xb + aj) * (xb + aj) - yb * yb) / (r2s * r2s));
    float fb7 = 2.0 * con * yb * ((xb - aj) / (r1s * r1s) - (xb + aj) / (r2s * r2s));
    
    float uxds = -pr1 * sinbj * fb2 + pr2 * cosbj * fb3 + yb * (sinbj * fb4 - cosbj * fb5);
    float uxdn = -pr1 * cosbj * fb2 - pr2 * sinbj * fb3 - yb * (cosbj * fb4 + sinbj * fb5);
    float uyds = pr1 * cosbj * fb2 + pr2 * sinbj * fb3 - yb * (cosbj * fb4 + sinbj * fb5);
    float uydn = -pr1 * sinbj * fb2 + pr2 * cosbj * fb3 - yb * (sinbj * fb4 - cosbj * fb5);

    float sxxds = cons * (2. * cosb2 * fb4 + sin2b * fb5 + yb * (cos2b * fb6 - sin2b * fb7));
    float sxxdn = cons * (-fb5 + yb * (sin2b * fb6 + cos2b * fb7));
    float syyds = cons * (2. * sinb2 * fb4 - sin2b * fb5 - yb * (cos2b * fb6 - sin2b * fb7));
    float syydn = cons * (-fb5 - yb * (sin2b * fb6 + cos2b * fb7));
    float sxyds = cons * (sin2b * fb4 - cos2b * fb5 + yb * (sin2b * fb6 + cos2b * fb7));
    float sxydn = cons * (-yb * (cos2b * fb6 - sin2b * fb7));

    s2us.uxs += msym * uxds;
    s2us.uxn += uxdn;
    s2us.uys += msym * uyds;
    s2us.uyn += uydn;

    s2us.sxxs += msym * sxxds;
    s2us.sxxn += sxxdn;
    s2us.syys += msym * syyds;
    s2us.syyn += syydn;
    s2us.sxys += msym * sxyds;
    s2us.sxyn += sxydn;
    return;
}






void interface_matrix(int i, int j)
{
    /* set up the interface matrix */
    matx.A_is_js = 0;
    matx.A_is_jn = 0;
    matx.A_in_js = 0;
    matx.A_in_jn = 0;

    matx.B_is_js = 0;
    matx.B_is_jn = 0;
    matx.B_in_js = 0;
    matx.B_in_jn = 0;

    int is = 2 * i;
    int in = is + 1;
    BoundaryElement& be = elm_list[i];
    BoundaryElement& bj = elm_list[j];

    float xi = be.xm;
    float yi = be.ym;
    float cosbi = be.cosbet;
    float sinbi = be.sinbet;
    int mm = be.mat_no; // to confirm this later

    int js = 2 * j;
    int jn = js + 1;

    s2us.reset();
    float xj = bj.xm;
    float yj = bj.ym;
    float cosbj = bj.cosbet;
    float sinbj = bj.sinbet;
    int mmj = bj.mat_no;
    float aj = bj.a;

    if (mm != mmj) return; // avoid mixing 2 different interfaces
    // -----------------------------------------------------------
    coeff(xi, yi, xj, yj, aj, cosbj, sinbj, +1, mm);
    switch (symm.ksym + 1)
    {
    case 1:
        break;
    case 2:
        xj = 2.0 * symm.xsym - bj.xm;
        coeff(xi, yi, xj, yj, aj, cosbj, -sinbj, -1, mm);
        break;
    case 3:
        yj = 2.0 * symm.ysym - bj.ym;
        coeff(xi, yi, xj, yj, aj, -cosbj, sinbj, -1, mm);
        break;
    case 4:
        xj = 2.0 * symm.xsym - bj.xm;
        yj = 2.0 * symm.ysym - bj.ym;
        coeff(xi, yi, xj, yj, aj, -cosbj, -sinbj, +1, mm);
        break;
    case 5:
        xj = 2.0 * symm.xsym - bj.xm;
        coeff(xi, yi, xj, yj, aj, cosbj, -sinbj, -1, mm);
        xj = bj.xm;
        yj = 2.0 * symm.ysym - bj.ym;
        coeff(xi, yi, xj, yj, aj, -cosbj, sinbj, -1, mm);
        xj = 2.0 * symm.xsym - bj.xm;
        coeff(xi, yi, xj, yj, aj, -cosbj, -sinbj, +1, mm);
        break;
    }
    matx.A_is_js = (s2us.syys - s2us.sxxs) * sinbi * cosbi + s2us.sxys * (cosbi * cosbi - sinbi * sinbi);
    matx.A_is_jn = (s2us.syyn - s2us.sxxn) * sinbi * cosbi + s2us.sxyn * (cosbi * cosbi - sinbi * sinbi);
    matx.A_in_js = s2us.sxxs * sinbi * sinbi - 2.0 * s2us.sxys * sinbi * cosbi + s2us.syys * cosbi * cosbi;
    matx.A_in_jn = s2us.sxxn * sinbi * sinbi - 2.0 * s2us.sxyn * sinbi * cosbi + s2us.syyn * cosbi * cosbi;

    matx.B_is_js = s2us.uxs * cosbi + s2us.uys * sinbi;
    matx.B_is_jn = s2us.uxn * cosbi + s2us.uyn * sinbi;
    matx.B_in_js = -s2us.uxs * sinbi + s2us.uys * cosbi;
    matx.B_in_jn = -s2us.uxn * sinbi + s2us.uyn * cosbi;
    return;
}






void label_500(int i, int j, int mm, int mmj, int is, int in, int js, int jn, int mode)
{
    /* this function implements all code after label 500 in mainb
    as before  label500 there is a jump to end of j loop means it's not part of
    normal run of mainb */


    // !------------zero matrix, (1, m, n, k) -> (n, k)
    if (elm_list[i].kod == 6)
    {
        //-------- - read existing influence coefficients-----------------------       
        // Sara! optimization can be done, for case 2 is always similar for 1 depends only on mm       
        if (i <= numbe_old && j <= numbe_old && mode > 0)
        {
            //othogonal line components--or -Non-othogonal line components, same region -
            if (i == j || mm == mmj)
            {
                switch (b_elm[i].ipair)
                {
                case 1:
                    s4.c[is][js] = s4.c_s[is][js];
                    s4.c[is][jn] = s4.c_s[is][jn];
                    s4.c[in][js] = s4.c_s[in][js];
                    s4.c[in][jn] = s4.c_s[in][jn];
                    break;
                case 2:
                    s4.c[is][js] = s4.c_d[is][js];
                    s4.c[is][jn] = s4.c_d[is][jn];
                    s4.c[in][js] = s4.c_d[in][js];
                    s4.c[in][jn] = s4.c_d[in][jn];
                    break;
                }
            }
            //----Non-othogonal line components, different region -- to me is the same !sara
            else
            {
                switch (b_elm[i].ipair)
                {

                case 1:
                    s4.c[is][js] = -s4.c_s[is][js];
                    s4.c[is][jn] = -s4.c_s[is][jn];
                    s4.c[in][js] = -s4.c_s[in][js];
                    s4.c[in][jn] = -s4.c_s[in][jn];
                    break;
                case 2:
                    s4.c[is][js] = s4.c_d[is][js];
                    s4.c[is][jn] = s4.c_d[is][jn];
                    s4.c[in][js] = s4.c_d[in][js];
                    s4.c[in][jn] = s4.c_d[in][jn];
                    break;

                }
            }
            return;  //goto 300                                    
        }
        // !----othogonal line components------------------------------------------
        if (i == j)
        {
            interface_matrix(i, i);
            //Shear / normal stress condition to balance the equation  

            switch (b_elm[i].ipair)
            {
            case 1:
                s4.c[is][js] = matx.A_is_js;
                s4.c[is][jn] = matx.A_is_jn;
                s4.c[in][js] = matx.A_in_js;
                s4.c[in][jn] = matx.A_in_jn;
                break;
            case 2:
                s4.c[is][js] = matx.B_is_js;
                s4.c[is][jn] = matx.B_is_jn;
                s4.c[in][js] = matx.B_in_js;
                s4.c[in][jn] = matx.B_in_jn;
                break;
            }
        }
        //----Non-othogonal line components, same region -----------------------------
        else
        {
            //to me is redundant and better to combine the conditiion with the previous
            if (mm == mmj)
            {
                interface_matrix(i, j);
                switch (b_elm[i].ipair)
                {
                case 1:
                    s4.c[is][js] = matx.A_is_js;
                    s4.c[is][jn] = matx.A_is_jn;
                    s4.c[in][js] = matx.A_in_js;
                    s4.c[in][jn] = matx.A_in_jn;
                    break;
                case 2:
                    s4.c[is][js] = matx.B_is_js;
                    s4.c[is][jn] = matx.B_is_jn;
                    s4.c[in][js] = matx.B_in_js;
                    s4.c[in][jn] = matx.B_in_jn;
                    break;
                }
            }
            //!----Non-othogonal line components, different region ----
            else   //mm!=mmj
            {
                switch (b_elm[i].ipair)
                {
                case 1:
                    interface_matrix(i + 1, j);
                    s4.c[is][js] = -matx.A_is_js;
                    s4.c[is][jn] = -matx.A_is_jn;
                    s4.c[in][js] = -matx.A_in_js;
                    s4.c[in][jn] = -matx.A_in_jn;
                    break;
                case 2:
                    interface_matrix(i - 1, j);
                    s4.c[is][js] = matx.B_is_js;
                    s4.c[is][jn] = matx.B_is_jn;
                    s4.c[in][js] = matx.B_in_js;
                    s4.c[in][jn] = matx.B_in_jn;
                    break;
                }
            }
        }
        //save the coeff.all time, fictious element c may not be used for mode = 0, by will be used in Bound() for work1

        s4.c_s[is][js] = matx.A_is_js;
        s4.c_s[is][jn] = matx.A_is_jn;
        s4.c_s[in][js] = matx.A_in_js;
        s4.c_s[in][jn] = matx.A_in_jn;

        s4.c_d[is][js] = matx.B_is_js;
        s4.c_d[is][jn] = matx.B_is_jn;
        s4.c_d[in][js] = matx.B_in_js;
        s4.c_d[in][jn] = matx.B_in_jn;
    } //if codi 

    return;
}






void mainb_ini(int mode, int i_st, int i_en, int j_st, int j_en)
{
    int in = 0, is = 0, mm = 0, mmj = 0, jn = 0, js = 0, n = 0;
    float xi = 0, yi = 0, cosbi = 0, sinbi = 0, xj, yj, cosbj = 0, sinbj = 0, aj = 0;
    double ks = 0.0, kn = 0.0;
    float ph = 0.0, pd = 0.0;
    double S_is_js = 0, S_is_jn = 0, S_in_js = 0, S_in_jn = 0, d_is_js = 0, d_is_jn = 0, d_in_js = 0, d_in_jn = 0;
    double sinbi2 = 0;
    double cosbi2 = 0;
    double sincosbi = 0;
    for (int i = i_st; i < i_en; ++i)
    {
        BoundaryElement elm_i = elm_list[i];
        is = 2 * i;
        in = is + 1;
        xi = elm_i.xm;
        yi = elm_i.ym;
        cosbi = elm_i.cosbet;
        sinbi = elm_i.sinbet;
        mm = elm_i.mat_no;
        sinbi2 = sinbi * sinbi;
        cosbi2 = cosbi * cosbi;
        sincosbi = sinbi * cosbi;

        if (mm == 0)
        {
            MessageBox(nullptr, L"material not defined_i mainb.", L"Error", MB_OK);
            file2 << "program stopped due to undefined material number for element = " << i << endl;
            logfile << "program stopped due to undefined material number for element = " << i << endl;
            exit(EXIT_FAILURE);
        }
        for (int j = j_st; j < j_en; ++j)
        {
            BoundaryElement elm_j = elm_list[j];
            mmj = elm_j.mat_no;
            js = 2 * j;
            jn = js + 1;
            if (mm != mmj || elm_list[i].kod == 6)
            {
                label_500(i, j, mm, mmj, is, in, js, jn, mode);
                continue;
            }
            //label 2111-------------------------
            s2us.reset();
            xj = elm_j.xm;
            yj = elm_j.ym;
            cosbj = elm_j.cosbet;
            sinbj = elm_j.sinbet;
            aj = elm_j.a;
            //-------------------------------
            coeff(xi, yi, xj, yj, aj, cosbj, sinbj, +1, mm);

            switch (symm.ksym + 1)
            {
            case 1:
                break;
            case 2:
                xj = 2.0 * symm.xsym - elm_j.xm;
                coeff(xi, yi, xj, yj, aj, cosbj, -sinbj, -1, mm);
                break;
            case 3:
                yj = 2.0 * symm.ysym - elm_j.ym;
                coeff(xi, yi, xj, yj, aj, -cosbj, sinbj, -1, mm);
                break;
            case 4:
                xj = 2.0 * symm.xsym - elm_j.xm;
                yj = 2.0 * symm.ysym - elm_j.ym;
                coeff(xi, yi, xj, yj, aj, -cosbj, -sinbj, +1, mm);
                break;
            case 5:
                xj = 2.0 * symm.xsym - elm_j.xm;
                coeff(xi, yi, xj, yj, aj, cosbj, -sinbj, -1, mm);
                xj = elm_j.xm;
                yj = 2.0 * symm.ysym - elm_j.ym;
                coeff(xi, yi, xj, yj, aj, -cosbj, sinbj, -1, mm);
                xj = 2.0 * symm.xsym - elm_j.xm;
                coeff(xi, yi, xj, yj, aj, -cosbj, -sinbj, +1, mm);
                //label 240
            }
            S_is_js = (s2us.syys - s2us.sxxs) * sincosbi + s2us.sxys * (cosbi2 - sinbi2);
            S_is_jn = (s2us.syyn - s2us.sxxn) * sincosbi + s2us.sxyn * (cosbi2 - sinbi2);
            S_in_js = s2us.sxxs * sinbi2 - 2.0 * s2us.sxys * sincosbi + s2us.syys * cosbi2;
            S_in_jn = s2us.sxxn * sinbi2 - 2.0 * s2us.sxyn * sincosbi + s2us.syyn * cosbi2;
            d_is_js = s2us.uxs * cosbi + s2us.uys * sinbi;
            d_is_jn = s2us.uxn * cosbi + s2us.uyn * sinbi;
            d_in_js = -s2us.uxs * sinbi + s2us.uys * cosbi;
            d_in_jn = -s2us.uxn * sinbi + s2us.uyn * cosbi;          

            if (i == j && elm_i.kod != 5)
            {
                S_is_js += comvar::k_num;
                S_in_jn += comvar::k_num;
            }
            switch (elm_i.kod)
            {
            case 1:
            case 5:
            case 11:
            case 15:
                s4.c[is][js] = S_is_js;
                s4.c[is][jn] = S_is_jn;
                s4.c[in][js] = S_in_js;
                s4.c[in][jn] = S_in_jn;
                break;
            case 2:
            case 7:
            case 12:
            case 17:
                s4.c[is][js] = d_is_js;
                s4.c[is][jn] = d_is_jn;
                s4.c[in][js] = d_in_js;
                s4.c[in][jn] = d_in_jn;
                break;

            case 3:
            case 13:
                s4.c[is][js] = d_is_js;
                s4.c[is][jn] = d_is_jn;
                s4.c[in][js] = S_in_js;
                s4.c[in][jn] = S_in_jn;
                break;
            case 4:
            case 14:
                s4.c[is][js] = S_is_js;
                s4.c[is][jn] = S_is_jn;
                s4.c[in][js] = d_in_js;
                s4.c[in][jn] = d_in_jn;
                break;
            case 6:
            case 16:
                label_500(i, j, mm, mmj, is, in, js, jn, mode);
                continue;
            }
            //label2010
            //---------------save influence coefficients, save the coeff. all time ------
            s4.c_s[is][js] = S_is_js;
            s4.c_s[is][jn] = S_is_jn;
            s4.c_s[in][js] = S_in_js;
            s4.c_s[in][jn] = S_in_jn;

            s4.c_d[is][js] = d_is_js;
            s4.c_d[is][jn] = d_is_jn;
            s4.c_d[in][js] = d_in_js;
            s4.c_d[in][jn] = d_in_jn;
            //--------end of coeff.saving------------------------------------------
                //---- - gravity and infinity problem-----------------------------
            //label2020
            if (elm_i.kod > 10 && i == j)
            {
                ks = rock1[1].e / 1e4;
                kn = rock1[1].e / 1e4;
                switch (elm_i.kod)
                {
                case 11:
                    s4.c[is][js] = S_is_js + ks;
                    s4.c[is][jn] = S_is_jn;
                    s4.c[in][js] = S_in_js;
                    s4.c[in][jn] = S_in_jn + kn;
                    break;
                case 12:
                    s4.c[is][js] = d_is_js;
                    s4.c[is][jn] = d_is_jn;
                    s4.c[in][js] = d_in_js;
                    s4.c[in][jn] = d_in_jn;
                    break;
                case 13:
                    s4.c[is][js] = d_is_js;
                    s4.c[is][jn] = d_is_jn;
                    s4.c[in][js] = S_in_js;
                    s4.c[in][jn] = S_in_jn + kn;
                    break;
                case 14:
                    s4.c[is][js] = S_is_js + ks;
                    s4.c[is][jn] = S_is_jn;
                    s4.c[in][js] = d_in_js;
                    s4.c[in][jn] = d_in_jn;
                }
            }
            //label 2025
            if (elm_i.kod != 5 || i != j) continue;  //2080 is goto 300 continue loop j
            if (i == numbe - 1 && mode > 0)
            {
                ks = 10e22;
                kn = 10e22;
                ph = 30.0 * 3.1416 / 180.0;
                pd = 0;
            }
            else
            {
                ks = b_elm[i].aks;
                kn = b_elm[i].akn;
                ph = b_elm[i].phi;
                pd = b_elm[i].phid;
            }
            if (i == numbe - 1 && mode == 1)
            {
                    //second round
                    switch (b_elm[i].jstate + 1)
                    {
                    case 1:
                    case 3:
                        s4.c[is][js] += ks;
                        s4.c[in][jn] += kn;
                        break;
                    case 2:
                        s4.c[is][js] += ks;
                    }
                }               
           
            else if (i == numbe - 1 && mode == 2)
            {
                //third run
                switch (b_elm[i].jstate + 1)
                {
                case 1:
                    s4.c[is][js] += ks;
                    s4.c[in][jn] += kn;
                    break;
                case 2:
                    s4.c[in][jn] += kn;
                    s4.c[is][jn] = s4.c[is][jn];
                    break;
                case 3:
                    s4.c[in][jn] += kn;
                    s4.c[is][jn] = s4.c[is][jn] + kn * tanf(ph) * b_elm[i].jslipd;
                }
            }
            else
            {
                switch (b_elm[i].jstate + 1)
                {
                case 1:
                    s4.c[is][js] += ks;
                    s4.c[in][jn] += kn;
                    break;
                case 3:
                    s4.c[in][jn] += kn;
                    s4.c[is][jn] = s4.c[is][jn] + kn * tanf(ph) * b_elm[i].jslipd;
                    s4.c[in][js] = s4.c[in][js] + kn * tanf(pd) * (-b_elm[i].jslipd);
                }
            }

            //label2080 goto 300
               //next iteration of j loop
            //label 500 written as a function   
        }

    }
}





void mainb_work0(int mode)
{ 
    mainb_ini(mode, 0, numbe, 0, numbe);
    //  solve system of algebric equations.
    int n = 2 * numbe;
    //check_diagonal_for_zero();
    solve(n, mode);
}






void mainb_work1_ini()
{
    int in = 0, is = 0, mm = 0, mmj = 0, jn = 0, js = 0, n = 0;
    float xi = 0, yi = 0, cosbi = 0, sinbi = 0, xj, yj, cosbj = 0, sinbj = 0, aj = 0;
    double ks = 0.0, kn = 0.0;
    float ph = 0.0, pd = 0.0;
    double S_is_js = 0, S_is_jn = 0, S_in_js = 0, S_in_jn = 0, d_is_js = 0,
        d_is_jn = 0, d_in_js = 0, d_in_jn = 0;
    double sinbi2 = 0;
    double cosbi2 = 0;
    double sincosbi = 0;

    for (int i = 0; i < numbe_old; ++i)
    {
        BoundaryElement &elm_i = elm_list[i];
        is = 2 * i;
        in = is + 1;
        xi = elm_i.xm;
        yi = elm_i.ym;
        cosbi = elm_i.cosbet;
        sinbi = elm_i.sinbet;
        mm = elm_i.mat_no;
        sinbi2 = sinbi * sinbi;
        cosbi2 = cosbi * cosbi;
        sincosbi = sinbi * cosbi;

        if (mm == 0)
        {
            MessageBox(nullptr, L"material not defined_i mainb.", L"Error", MB_OK);
            file2 << "program stopped due to undefined material number for element = " << i << endl;
            logfile << "program stopped due to undefined material number for element = " << i << endl;
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < numbe_old; ++j)
        {
            BoundaryElement &elm_j = elm_list[j];
            mmj = elm_j.mat_no;
            js = 2 * j;
            jn = js + 1;
            if (mm != mmj || elm_list[i].kod == 6)
            {
                label_500(i, j, mm, mmj, is, in, js, jn, 1);
                continue;
            }
            //--------------- - read existing influence coefficients----------
            switch (elm_i.kod)
            {
            case 1:
            case 5:
            case 11:
            case 15:
                s4.c[is][js] = s4.c_s[is][js];
                s4.c[is][jn] = s4.c_s[is][jn];
                s4.c[in][js] = s4.c_s[in][js];
                s4.c[in][jn] = s4.c_s[in][jn];
                break;       //goto label2025;

            case 2:
            case 7:
            case 12:
            case 17:
                s4.c[is][js] = s4.c_d[is][js];
                s4.c[is][jn] = s4.c_d[is][jn];
                s4.c[in][js] = s4.c_d[in][js];
                s4.c[in][jn] = s4.c_d[in][jn];
                break;

            case 3:
            case 13:
                s4.c[is][js] = s4.c_d[is][js];
                s4.c[is][jn] = s4.c_d[is][jn];
                s4.c[in][js] = s4.c_s[in][js];
                s4.c[in][jn] = s4.c_s[in][jn];
                break;

            case 4:
            case 14:
                s4.c[is][js] = s4.c_s[is][js];
                s4.c[is][jn] = s4.c_s[is][jn];
                s4.c[in][js] = s4.c_d[in][js];
                s4.c[in][jn] = s4.c_d[in][jn];
                break;
            case 6:
            case 16:
                label_500(i, j, mm, mmj, is, in, js, jn, 1);
                continue;
            }

            if (elm_i.kod != 5 || i != j) continue;  //2080 is goto 300 continue loop j

            ks = b_elm[i].aks;
            kn = b_elm[i].akn;
            ph = b_elm[i].phi;
            pd = b_elm[i].phid;
            switch (b_elm[i].jstate + 1)
            {
            case 1:
                s4.c[is][js] += ks;
                s4.c[in][jn] += kn;
                break;
            case 3:
                s4.c[in][jn] += kn;
                s4.c[is][jn] = s4.c[is][jn] + kn * tanf(ph) * b_elm[i].jslipd;
                s4.c[in][js] = s4.c[in][js] + kn * tanf(pd) * (-b_elm[i].jslipd);
            }
        }
    }
    return;
}






void mainb_work1(int mode)
{
    mainb_work1_ini();
    mainb_ini(mode, 0, numbe, numbe - 1, numbe);
    mainb_ini(mode, numbe - 1, numbe, 0, numbe);
    int n = 2 * numbe;
    //check_diagonal_for_zero();
    /*if (fabs(comvar::s4.c[n-1][n-1]) < 1e-12 || comvar::s4.c[n - 1][n - 1]==0) {
        printf("Near-zero at diagonal [%d][%d] = %.3e\n", n-1, n-1, comvar::s4.c[n - 1][n - 1]);
        comvar::s4.c[n - 1][n - 1] = 1e-5;
      }*/
    solve(n, mode);
}
