#include<stdafx.h>
#include <GeologicalForm.h>
#include "CommonPara.h"

using namespace CommonPara_h::comvar;
using namespace std;




GeologicalForm::GeologicalForm():mat_no(1), elem_no(1), bound_type(1) {}
GeologicalForm::GeologicalForm(int mat, int eleno, int kode) : mat_no(mat), elem_no(eleno), bound_type(kode) {}
int GeologicalForm::getMatno(){ return mat_no;}
int GeologicalForm::getElem_no() { return elem_no; }
int GeologicalForm::getBoundtype() { return bound_type; }





void GeologicalForm::chk_potential_crack_growth(float xbeg, float ybeg,float xend,float yend,int numbe0,int itype)
{   
    float x1 = 0, y1 = 0,x2 = 0, y2 = 0;
   
    if (itype == 0) return;
    else
    {
        if (itype == -1 || itype == 2)
        {
            for (int m = numbe0 ; m < numbe; m++)
            {
                x1 = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
                y1 = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;
                if ((abs(xbeg - x1) <= 1e-5) && (abs(ybeg - y1) <= 1e-5))
                {                   
                    tips[no].xbe = elm_list[m].xm - 3.0 * elm_list[m].a * elm_list[m].cosbet;
                    tips[no].ybe = elm_list[m].ym - 3.0 * elm_list[m].a * elm_list[m].sinbet;
                    tips[no].xen = tips[no].xbe + 2. * elm_list[m].a * elm_list[m].cosbet;
                    tips[no].yen = tips[no].ybe + 2. * elm_list[m].a * elm_list[m].sinbet;
                    tips[no].dl = 2. * elm_list[m].a;
                    tips[no].ityp = -1;
                    tips[no].costem = elm_list[m].cosbet;     
                    tips[no].sintem = elm_list[m].sinbet;
                    tips[no].mpointer = m;
                    tips[no].mat_no = elm_list[m].mat_no;
                    no++;
                    break;    //Sara ! not sure yet
                }
            }
        }
        if (itype == 1 || itype == 2)
        {
            for (int m = numbe0 ; m < numbe; m++)  //Sara not sure about index
            {
                x2 = elm_list[m].xm + elm_list[m].a * elm_list[m].cosbet;
                y2 = elm_list[m].ym + elm_list[m].a * elm_list[m].sinbet;

                if ((abs(xend - x2) <= 1e-5) && (abs(yend - y2) <= 1e-5)) 
                {                    
                    tips[no].xbe = xend;
                    tips[no].ybe = yend;
                    tips[no].xen = tips[no].xbe + 2. * elm_list[m].a * elm_list[m].cosbet;
                    tips[no].yen = tips[no].ybe + 2. * elm_list[m].a * elm_list[m].sinbet;
                    tips[no].dl = 2. * elm_list[m].a;

                    tips[no].ityp = 1;
                    tips[no].costem = elm_list[m].cosbet;
                    tips[no].sintem = elm_list[m].sinbet;
                    tips[no].mpointer = m;
                    tips[no].mat_no = elm_list[m].mat_no;
                    no++;
                    break; // notsure Sara
                }
            }
        }
    }        
     return;
}




bool approximatelyEqual(float a1, float a2) {
    return std::round(a1 * 1e5) == std::round(a2 * 1e5);
}




//check dupliction of newly added element
bool  isNewElementUnique(const BoundaryElement& newelem)
{
    for (int ii = 0; ii < numbe ; ii++)
    {
        BoundaryElement & be = elm_list[ii];
        if (be.xm == newelem.xm &&
            be.ym == newelem.ym &&
            //be.a == newelem.a &&
            approximatelyEqual(be.a,newelem.a) &&
            be.sinbet == newelem.sinbet &&
            be.cosbet == newelem.cosbet) {
            return false;
        }
    }
    return true;
}




// check if element m cross any of the saved elements or newly added elements 
void GeologicalForm::cross_current_or_save_element(float xb1, float yb1, float xe1, float ye1, float xcross,
    float ycross, int m, int ii, int &k)
    {   
    elm_list[m].xm = 0.5 * (xb1 + xcross);
    elm_list[m].ym = 0.5 * (yb1 + ycross);
    float a0 = elm_list[m].a;
    
    elm_list[m].a = 0.5 * sqrt(pow(xb1 - xcross, 2) + pow(yb1 - ycross, 2));
    b_elm[m].force1 = b_elm[m].force1 * elm_list[m].a / a0;
    b_elm[m].force2 = b_elm[m].force2 * elm_list[m].a / a0;
    BoundaryElement newelement;
    newelement.xm = 0.5 * (xe1 + xcross);
    newelement.ym = 0.5 * (ye1 + ycross);
    newelement.a = 0.5 * sqrt(pow(xe1 - xcross, 2) + pow(ye1 - ycross, 2));

    newelement.sinbet = elm_list[m].sinbet;
    newelement.cosbet = elm_list[m].cosbet;   
    bool unique = isNewElementUnique(newelement);
    float length = 2 * newelement.a;
    if (!unique || newelement.a / ( elm_list[m].a) < 0.125)
    {
        k--;
        return;
    }
    elm_list[ii] = newelement;   
    b_elm[ii].force1 = b_elm[m].force1 * elm_list[ii].a / a0;
    b_elm[ii].force2 = b_elm[m].force2 * elm_list[ii].a / a0;
    elm_list[ii].kod = elm_list[m].kod;
   
    elm_list[ii].mat_no = elm_list[m].mat_no;
    b_elm[ii].aks = b_elm[m].aks;
    b_elm[ii].akn = b_elm[m].akn;
    b_elm[ii].phi = b_elm[m].phi;
    b_elm[ii].phid = b_elm[m].phid;
    
    b_elm[ii].coh = b_elm[m].coh;
    joint[ii].aperture0 = joint[m].aperture0;
    joint[ii].aperture_r = joint[m].aperture_r;
    s4.b0[2 * ii ] = s4.b0[2 * m ];
    s4.b0[2 * ii +1] = s4.b0[2 * m +1];
    return ;
}





// all newlly added element checked for crossing with other elements
void GeologicalForm::ctl_cross_elements(int& k,int m, int numbe0,int num)
{
    
    float x1 = 0.0, y1 = 0.0, xb1, yb1, xe1, ye1, xb2, xe2, x2,
        y2, yb2, ye2, xcross = 0, ycross = 0;
    float tan1 = 0, tan2 = 0;
    
    //from here newly added
    if (elm_list[m].cosbet == 0)
    {
        tan1 = 10e20 * elm_list[m].sinbet;
    }
    else
    {
        tan1 = elm_list[m].sinbet / elm_list[m].cosbet;
    }
    x1 = elm_list[m].xm;
    y1 = elm_list[m].ym;
    xb1 = elm_list[m].xm - elm_list[m].a * elm_list[m].cosbet;
    yb1 = elm_list[m].ym - elm_list[m].a * elm_list[m].sinbet;
    xe1 = elm_list[m].xm + elm_list[m].a * elm_list[m].cosbet;
    ye1 = elm_list[m].ym + elm_list[m].a * elm_list[m].sinbet;

    //element m compare with all elements before it     
    for (int j = 0; j < numbe - 1; j++)
    {
        if (elm_list[j].cosbet == 0)
        {
            tan2 = 10e20 * elm_list[j].sinbet;
        }
        else 
        {
            tan2 = elm_list[j].sinbet / elm_list[j].cosbet;
        }

        if (tan1 == tan2) 
        {
            continue;   
        }
        x2 = elm_list[j].xm;
        y2 = elm_list[j].ym;
        xb2 = elm_list[j].xm - elm_list[j].a * elm_list[j].cosbet;;
        yb2 = elm_list[j].ym - elm_list[j].a * elm_list[j].sinbet;
        xe2 = elm_list[j].xm + elm_list[j].a * elm_list[j].cosbet;
        ye2 = elm_list[j].ym + elm_list[j].a * elm_list[j].sinbet;

        xcross = (y2 - y1 + tan1 * x1 - tan2 * x2) / (tan1 - tan2);
        ycross = (tan1 * tan2 * (x1 - x2) + tan1 * y2 - tan2 * y1) / (tan1 - tan2);

        //------for saved elements
        if ((xcross > max(xb1, xe1) || xcross < min(xb1, xe1)) ||
            (ycross > max(yb1, ye1) || ycross < min(yb1, ye1)) ||
            (xcross > max(xb2, xe2) || xcross < min(xb2, xe2)) ||
            (ycross > max(yb2, ye2) || ycross < min(yb2, ye2))) 
        {
           continue;
        }
        float threshold = 1e-5;
        if (!((abs(xcross - xb2) <= threshold && abs(ycross - yb2) <= threshold) ||
            (abs(xcross - xe2) <= threshold && abs(ycross - ye2) <= threshold)) )
        {          
            k++;               //one new element will be added            
            int jj = numbe0 + num + k - 1 ;   // element_no instead of num, in the orig code num is one of input parametes to elem
            cross_current_or_save_element(xb2, yb2, xe2, ye2, xcross, ycross, j, jj,k);
        }

        if ((abs(xcross - xb1) <= threshold && abs(ycross - yb1) <= threshold) ||
            (abs(xcross - xe1) <= threshold && abs(ycross - ye1) <= threshold))
        {
            continue;       //goto 200
        }
        k++;         //one new element will be added 
        int ii = numbe0 + num + k - 1;
        cross_current_or_save_element(xb1, yb1, xe1, ye1, xcross, ycross, m, ii,k);       
    }       
    return;
  }





int GeologicalForm::def_boundary_elements_for_Geoform(int num, float xbeg, float ybeg, float xend, float yend,
    float bvs, float bvn, float gradsy,float gradny, int  itype, int jmat)
{
    /*
    numbe is the global variable shows the total no of elements, so
    for adding new elements for a new form we need to insert it after the last element in the list 
    kode = 5 means joints
    itype = -1 tip elem at beg point; =+1, tip elem at end point;
!    itype = 0  no tip elem; =+2, tip elem at both points.
    ieven is for uneven distribution of elements, not used
    */

    float x1=0, y1=0, xb1=0, xe1=0, yb1=0, ye1=0, xb2=0, xe2=0, x2=0, y2=0,yb2=0,ye2=0, xcross = 0,ycross = 0;
    int mm = mat_no;
    float st = sqrt(pow((xend - xbeg), 2) + pow((yend - ybeg), 2));// an estimation of length of fracture 
    
    float angd = comvar::pi / num;  
    float ang0 = 0.5 * angd;
    float xs = xbeg;
    float ys = ybeg;
 
    float xd = (xend - xbeg) / num;
    float yd = (yend - ybeg) / num;
    int k = 0;  
    int numbe0 = numbe;   // numbe0 is the current number of elements or the old no.
    int m = 0;

    float sw = st / num;    
    float sinb = yd / sw;
    float cosb = xd / sw;
    float tmp_a = 0.5 * sw;    // half of element length
    bool unique = true;
    BoundaryElement newelement;  
  
    /* we add element_no or num elements to the current list of elements
    for the new geological form
    for each element first the coordination is computed 
    */

   for (int ne = 0; ne < num; ne++) 
    {   
        xs += xd;
        ys += yd;
       
        newelement.xm = xs - 0.5 * xd;  // the x and y is  coordinates of the middle of element
        newelement.ym = ys - 0.5 * yd;  //optimize it later
        newelement.a = tmp_a;
        newelement.sinbet = sinb;
        newelement.cosbet = cosb; 
        if (numbe > 0)
        {
            unique = isNewElementUnique(newelement);
            if (!unique)   continue;
        }        
        numbe++;
        if (numbe >= m0 - 1)
        {
            MessageBox(nullptr, L"Maximum BE limit exceeded!", L"Message!", MB_OK);
            exit(0);
            return 1;
        }
        m = numbe - 1;   //the position of new element to add 
        newelement.kod = bound_type;     // for all elems are equal
        newelement.mat_no = mat_no; 
        elm_list[m] = newelement;

        //for fractures
        if (bound_type == 5) 
        {
            BE newelm(s8[jmat].aks0, s8[jmat].akn0, s8[jmat].phi0, s8[jmat].phid0, s8[jmat].coh0);
            b_elm[m] = newelm;            
            joint[m].aperture0 = s8[jmat].apert0;
            joint[m].aperture_r = s8[jmat].apert_r;
            watercm.pwater[m] = 0;
            watercm.jwater[m] = 0;
        }
        int ms = 2 * m;
        int mn = ms  + 1;       

        s4.b1[ms] = bvs + elm_list[m].ym * gradsy;
        s4.b1[mn] = bvn + elm_list[m].ym * gradny;

        //adjust stress boundary values to account for initial stresses.
        MatrixB(m);
        ctl_cross_elements(k, m, numbe0,num); // there are more than one element    
    }
    if (k > 0 && numbe - numbe0 < num)
    {
        int t = numbe;
        for (int l = 0; l < k; l++)
        {
            int t1 = numbe0 + num + l;
            elm_list[t] = elm_list[t1];
            t++;
        }
    }
    numbe = numbe + k;  // k new elements added 
    if (numbe >= m0 - 1)
    {
        MessageBox(nullptr, L"Maximum BE limit exceeded!", L"Message!", MB_OK);
        exit(0);
        return 1;
    }
    chk_potential_crack_growth(xbeg, ybeg, xend, yend, numbe0,itype);
    numbe_old = numbe;      
    return numbe;
}


 