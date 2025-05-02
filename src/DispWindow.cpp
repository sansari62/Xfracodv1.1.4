#include<stdafx.h>

#include "DispWindow.h"
#include <CommonPara.h>
using namespace CommonPara_h::comvar;



DispWindow::DispWindow(): xll(0), xur(0), yll(0), yur(0), numx(0),numy(0), xc0(0), yc0(0), 
radium(0), numr(0), numa(0), ID_win(0) {}




void DispWindow::save_to_file(ofstream & f)
{
    f << xll << " " << xur << " " << yll << " " << yur << " " << numx << " " << numy <<
        " " << xc0 << " " << yc0 << radium << " " << numr << " " << numa << " " <<
        ID_win<< std::endl;
}




void CheckRange() {

    s15.aa = 0;
    s15.aaa = 1000000.0;
    float xl = 100000;
    float xh = -100000;
    float yl = 100000;
    float yh = -100000;

    for (int i = 0; i < numbe; ++i) 
    {
        s15.aa = (s15.aa < elm_list[i].a) ? elm_list[i].a : s15.aa;
        s15.aaa = (s15.aaa > elm_list[i].a) ? elm_list[i].a : s15.aaa;
           

        float xbeg = elm_list[i].xm - elm_list[i].a * elm_list[i].cosbet;
        float xend = elm_list[i].xm + elm_list[i].a * elm_list[i].cosbet;
        float ybeg = elm_list[i].ym - elm_list[i].a * elm_list[i].sinbet;
        float yend = elm_list[i].ym + elm_list[i].a * elm_list[i].sinbet;

        xl = min( xl, min(xbeg, xend));
        xh = max( xh, max(xbeg, xend));
        yl = min( yl, min(ybeg, yend));
        yh = max( yh, max(ybeg, yend));

        if (symm.ksym == 0)
            continue;
        else
        {
            if (symm.ksym == 1 || symm.ksym == 4)
            {
                xl = min(xl, min(2.0 * symm.xsym - xbeg, 2.0 * symm.xsym - xend));
                xh = max(xh, max(2.0 * symm.xsym - xbeg, 2.0 * symm.xsym - xend));
            }
            if (symm.ksym == 2 || symm.ksym == 4)
            {
                yl = min(yl, min(2.0 * symm.ysym - ybeg, 2.0 * symm.ysym - yend));
                yh = max(yh, max(2.0 * symm.ysym - ybeg, 2.0 * symm.ysym - yend));
            }
            if (symm.ksym == 3 || symm.ksym == 4)
            {
                xl = min(xl, min(2.0 * symm.xsym - xbeg, 2.0 * symm.xsym - xend));
                xh = max(xh, max(2.0 * symm.xsym - xbeg, 2.0 * symm.xsym - xend));
                yl = min(yl, min(2.0 * symm.ysym - ybeg, 2.0 * symm.ysym - yend));
                yh = max(yh, max(2.0 * symm.ysym - ybeg, 2.0 * symm.ysym - yend));
            }
        }
    }

    if (s5u.xmin == -100000 && s5u.xmax == 100000 &&
        s5u.ymin == -100000 && s5u.ymax == 100000) 
    {
        s5u.xmin = xl - (xh - xl) * 0.5; // xll
        s5u.ymin = yl - (yh - yl) * 0.5; // yll
        s5u.xmax = xh + (xh - xl) * 0.5; // xur
        s5u.ymax = yh + (yh - yl) * 0.5; // yur
    }
    return;
}

