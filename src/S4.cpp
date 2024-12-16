#include "S4.h"
#include<CommonPara.h>

using namespace  CommonPara_h::comvar;




S4::S4(int length) :
    c(length / 2, vector<float>(length / 2, 0.0)), c_s(length / 2, vector<float>(length / 2, 0.0)),
    c_d(length / 2, vector<float>(length / 2, 0.0)),
    d0(length / 2, 0.0), d(length / 2, 0.0), b(length / 2, 0.0), b1(length / 2 / 2, 0.0),
    df0(length / 2, 0.0), df(length / 2, 0.0), b0(length / 2, 0.0), b0_old(length / 2, 0.0) {

    c.reserve(length/2);
    c_s.reserve(length/2);
    c_d.reserve(length/2);
    d0.reserve(length/2);  
    d.reserve(length/2);
    b.reserve(length/2);
    b1.reserve(length/2);
    b0.reserve(length/2);
    b0_old.reserve(length/2);
    df0.reserve(length/2);
    df.reserve(length/2);
}




void S4::limit_d()
{
    float ds_ratio = 0, dn_ratio = 0, d_ratio = 0;
    if (ktipgrow == false) return;

    float d0s_max = 0;
    float d0n_max = 0;
    float ds_max = 0;
    float dn_max = 0;

    for (int k = 0; k < numbe; k++)
    {
        int tk = k * 2;
        if (d0s_max < abs(d0[tk])) 
            d0s_max = abs(d0[tk]);

        if (d0n_max < abs(d0[tk+1]))    
            d0n_max = abs(d0[tk+1]);

        if (ds_max < abs(d[tk ])) 
            ds_max = abs(d[tk]);

        if (dn_max < abs(d[tk+1])) 
            dn_max = abs(d[tk+1]);
    }

    if (d0s_max != 0)
        ds_ratio = ds_max / d0s_max;
    else
        ds_ratio = 0;

    if (d0n_max != 0)
        dn_ratio = dn_max / d0n_max;
    else
        dn_ratio = 0;

    d_ratio = max(ds_ratio, dn_ratio);    //max instead of amax1

    if (d_ratio > 1.0)
    {
        for (int k = 0; k < numbe; k++)
        {
            d[k * 2] = d[k * 2] / d_ratio / 2.0;
            d[k * 2 +1] = d[k * 2 +1] / d_ratio / 2.0;
        }
    }
}




void S4::save_to_file(ofstream& f)
{
    for (int m = 0; m < numbe * 2; ++m)
    {
        f << b[m] << " " << d[m] << " " << b0[m] << " " << b1[m] << " " << d0[m] 
            << " " << df0[m] << " " << df[m] << std::endl;
        for (int n = 0; n < numbe * 2; ++n)
        {
            f << c[m][n] << " " << c_s[m][n] << " " << c_d[m][n] << std::endl;
        }
    }
    return;
}


