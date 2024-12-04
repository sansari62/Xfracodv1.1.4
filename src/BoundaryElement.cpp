#include "BoundaryElement.h"
#include <CommonPara.h>

using namespace comvar;

float val = 0.0;

BoundaryElement::BoundaryElement():mat_no(0),xm(val), ym(val), kod(0), sigma_s(val), sigma_n(val), a(val), sinbet(val), cosbet(val),
ipair(0), jstate(0), jslipd(0), us(val), un(val), forcen(val), forces(val), force1(val), force2(val), jmode(val), phi(val), phid(val),
 aks(val), akn(val), us_neg(0),un_neg(0),ss_old(val),sn_old(val),coh(val) {}



BoundaryElement::BoundaryElement(float x, float y, float am, float cos, float sin, int kd, int mat)
{
	xm = x;
	ym = y;
	a = am;
	sinbet = sin;
	cosbet = cos;
	mat_no = mat;
	kod = kd;

}




void BoundaryElement::read_from_file(ifstream& f,int m)
{
	f >> xm  >> ym  >> cosbet  >> sinbet  >> a  >> kod  >> ipair  >> jstate  >> jslipd  >> jmode  >>
		force1  >> force2  >> tips[m].costem >> tips[m].sintem >> us >> un >> forces >> forcen >>
		aks  >> akn  >> phi  >> phid  >> coh  >> mat_no  >> joint[m].aperture0  >> joint[m].aperture_r  >>
		sigma_s  >> sigma_n  >> us_neg  >> un_neg  >> ss_old  >> sn_old ;

}




void BoundaryElement::save_to_file(ofstream& f, int m)
{
	f << xm << " " << ym << " " << cosbet << " " << sinbet << " " 
		<< a << " " << kod << " " << ipair << " " << jstate << " "
		<< jslipd << " " << jmode << " " <<
		force1 << " " << force2 << " " << tips[m].costem << " " << tips[m].sintem
		<< " " << us << " " << un << " " << forces << " " << forcen
		<< " " <<		aks << " " << akn << " " << phi << " " << phid
		<< " " << coh << " "  << mat_no << " " << joint[m].aperture0
		<< " " << joint[m].aperture_r << " " <<
		sigma_s << " " << sigma_n << " " << us_neg << " " << un_neg <<
		" " << ss_old << " " << sn_old << std::endl;
}



void BoundaryElement::bound(int i, float& ss, float& sn, float& ustem, float& untem, float& usneg, float& unneg)
{
    /*compute boundary displacements and stresses.
    boundary stress / displacement using d0(m) at any given
    iteration steps*/

    int mm = mat_no;
    int is = i * 2;
    int in = is + 1;
    

    float y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - ym);
    float pyy = symm.pyy1 + g.sky * (y0 - ym);
    float pxy = symm.pxy1;
    
    usneg = 0.0;
    unneg = 0.0;

    ss = (pyy - pxx) * sinbet * cosbet + pxy * (cosbet * cosbet - sinbet * sinbet);
    sn = pxx * sinbet * sinbet - 2.0 * pxy * sinbet * cosbet + pyy * cosbet * cosbet;

    if (mm == mat_lining)
    {
        ss = 0;
        sn = 0;
    }

    for (int j = 0; j < numbe; j++)
    {
        int js = 2 * j;
        int jn = js + 1;      
        if (mm != elm_list[j].mat_no) continue;

        ss += s4.c_s[is][js] * s4.d0[js] + s4.c_s[is][jn] * s4.d0[jn];
        sn += s4.c_s[in][js] * s4.d0[js] + s4.c_s[in][jn] * s4.d0[jn];
        usneg += s4.c_d[is][js] * s4.d0[js] + s4.c_d[is][jn] * s4.d0[jn];
        unneg += s4.c_d[in][js] * s4.d0[js] + s4.c_d[in][jn] * s4.d0[jn];
    }
    
    if (kod == 5 || kod == 15)
    {
        ustem = s4.d0[is];
        untem = s4.d0[in];
    }
    else   //for 1 to 7 and 11 to 17 material comvar::number
    {
        ustem = usneg;
        untem = unneg;
    }
    return;
}



void  MatrixB(int m)
{
    
    int ms = 2 * m;
    int mn = ms + 1;
    

    float y0;

    y0 = g.y_surf;
    float pxx = symm.pxx1 + g.skx * (y0 - elm_list[m].ym);
    float pyy = symm.pyy1 + g.sky * (y0 - elm_list[m].ym);
    float pxy = symm.pxy1;

    float cosb = elm_list[m].cosbet;
    float sinb = elm_list[m].sinbet;
    float ss = (pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb);
    float sn = pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb;

    switch (elm_list[m].kod)
    {
    //Sara optimization
    case 1:
    case 5:
    case 11:
    case 15:
        s4.b0[ms] = s4.b1[ms] - ss;
        s4.b0[mn] = s4.b1[mn] - sn - watercm.jwater[m] * watercm.pwater[m];
        if (elm_list[m].mat_no == mat_lining)
        {
            s4.b0[ms] = s4.b1[ms];
            s4.b0[mn] = s4.b1[mn] - watercm.jwater[m] * watercm.pwater[m];
        }
        break;
    case 2:
    case 7:
    case 12:
    case 17:
        s4.b0[ms] = s4.b1[ms];
        s4.b0[mn] = s4.b1[mn];
        break;

    case 3:
    case 13:
        s4.b0[mn] = s4.b1[mn] - sn - watercm.jwater[m] * watercm.pwater[m];
        if (elm_list[m].mat_no == mat_lining) 
            s4.b0[mn] = s4.b1[mn] - watercm.jwater[m] * watercm.pwater[m];
        s4.b0[ms] = s4.b1[ms];
        break;

    case 4:
    case 14:
        s4.b0[ms] = s4.b1[ms] - ss;
        if (elm_list[m].mat_no == mat_lining)  
            s4.b0[ms] = s4.b1[ms];
        s4.b0[mn] = s4.b1[mn];
        break;

    case 6:
    case 16:
        if (elm_list[m].mat_no == mat_lining || elm_list[m + 1].mat_no == mat_lining)
        {
            if (elm_list[m].ipair == 1)
            {
                s4.b0[ms] = s4.b1[ms] - ss;
                s4.b0[mn] = s4.b1[ms] - sn;
            }
            else
            {
                if (elm_list[m].ipair == 2)
                {
                    s4.b0[ms] = 0;
                    s4.b0[mn] = 0;
                }
            }
        }
        else
        {
            s4.b0[ms] = 0;
            s4.b0[mn] = 0;
        }
        break;

    }
    if (elm_list[m].mat_no == mat_lining)    //Force1 etc is not function of s4.b0, but pxx etc.
    {

        elm_list[m].force1 = 0;
        elm_list[m].force2 = 0;
    }
    else
    {
        elm_list[m].force1 = 2.0 * elm_list[m].a * (-ss);
        elm_list[m].force2 = 2.0 * elm_list[m].a * (-sn);
    }
    return;
}

