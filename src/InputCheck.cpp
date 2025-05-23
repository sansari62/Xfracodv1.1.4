#include <stdafx.h>
#include <CommonPara.h>
#include<Source.h>
#include<rstor_chek.h>
#include<Rectangle_check.h>


using namespace CommonPara_h::comvar;



const int len = 400;
int ncr[len] = {0}, itip[len] = {0};
float xcr[len][len] = {0.0}, ycr[len][len] = { 0.0 }, acr[len][len] = {0.0};


static void Interface(int num, float xbeg, float ybeg, float xend, float yend, int mat1, int mat2)
{
    float st, angd, ang0, a0, xs, ys, xd, yd, sw;
    int  m, n, mn, ms, nn, ns;
    float pxx, pyy, pxy;
    float cosb, sinb, ss, sn;

    // Define locations, size, orientations of multiregion interfaces
    st = sqrt((xend - xbeg) * (xend - xbeg) + (yend - ybeg) * (yend - ybeg));
    angd = pi / num;
    ang0 = 0.5 * angd;
    a0 = st * 0.5;

    xs = xbeg;
    ys = ybeg;
    xd = (xend - xbeg) / num;
    yd = (yend - ybeg) / num;
    sw = sqrt(xd * xd + yd * yd);
    int numbe0 = numbe;
    for (int ne = 0; ne < num; ++ne)
    {
        numbe += 2;
        m = numbe - 2;
        n = numbe - 1;
        BoundaryElement& be = elm_list[m];

        sw = st / num;
        xd = (xend - xbeg) * sw / st;
        yd = (yend - ybeg) * sw / st;
        xs += xd;
        ys += yd;

        be.xm = xs - 0.5 * xd;
        be.ym = ys - 0.5 * yd;
        be.a = 0.5 * sw;
        be.sinbet = yd / sw;
        be.cosbet = xd / sw;
        be.kod = 6;
        b_elm[m].ipair = 1;
        be.mat_no = mat1;

        elm_list[n].xm = be.xm;
        elm_list[n].ym = be.ym;
        elm_list[n].a = be.a;
        elm_list[n].sinbet = -be.sinbet;
        elm_list[n].cosbet = -be.cosbet;
        elm_list[n].kod = 6;
        b_elm[n].ipair = 2;
        elm_list[n].mat_no = mat2;
        //************************* Avoid double definition of elements********************************
        bool jmp_label = false;
        for (int ii = 0; ii < m - 1; ++ii)
        {
            if (be.xm == elm_list[ii].xm && be.ym == elm_list[ii].ym &&
                be.a == elm_list[ii].a && be.sinbet == elm_list[ii].sinbet && be.cosbet == elm_list[ii].cosbet)
            {
                numbe -= 2;    // I think those two newly added elements should be removed and the size of vector elm_list controlled by vec.size               
                jmp_label = true;
                break;
            }
        }
        if (jmp_label) continue;   //jmp to label110
        ms = 2 * m;
        mn = ms + 1;
        s4.b1[ms] = 0;
        s4.b1[mn] = 0;

        ns = 2 * n;
        nn = ns + 1;
        s4.b1[ns] = 0;
        s4.b1[nn] = 0;

        //------------------ Adjust stress boundary values to account for initial stresses-------------------------
        float y0 = g.y_surf;
        pxx = symm.pxx1 + g.skx * (y0 - be.ym);
        pyy = symm.pyy1 + g.sky * (y0 - be.ym);
        pxy = symm.pxy1;

        cosb = be.cosbet;
        sinb = be.sinbet;
        ss = (pyy - pxx) * sinb * cosb + pxy * (cosb * cosb - sinb * sinb);
        sn = pxx * sinb * sinb - 2.0 * pxy * sinb * cosb + pyy * cosb * cosb;

        // Zero at interface - it is not real boundary value; for concrete lining insitu stresses should not be included
        if (be.mat_no == mat_lining || elm_list[n].mat_no == mat_lining)
        {
            s4.b0[ms] = -ss;                //interface element
            s4.b0[mn] = -sn;
        }
        else
        {
            s4.b0[ms] = 0;                  //interface element
            s4.b0[mn] = 0;
        }
        s4.b0[ns] = 0;                      //interface element
        s4.b0[nn] = 0;

        // Calculate forces         !Force1 etc is not function of b0, but pxx etc.
        //still work done by interface if it is between lining and rock
        if (be.mat_no == mat_lining)
        {
            b_elm[m].force1 = 0;
            b_elm[m].force2 = 0;
        }
        //still work done by interface if it is between lining and rock
        else
        {
            b_elm[m].force1 = 2.0 * be.a * (-ss);
            b_elm[m].force2 = 2.0 * be.a * (-sn);
        }

        if (elm_list[n].mat_no == mat_lining)
        {
            b_elm[n].force1 = 0;
            b_elm[n].force2 = 0;
        }
        else
        {
            b_elm[n].force1 = 2.0 * elm_list[n].a * (-ss);
            b_elm[n].force2 = 2.0 * elm_list[n].a * (-sn);
        }
    }
    //label110
    numbe_old = numbe;
    return;
}





void reorder_fractures(fstream& file25)
{          

    float xbeg = 0, ybeg = 0, xend = 0, yend = 0, xtem = 0, ytem = 0;
    float d1, d2;
    int itl, itr;
    for (int i = 0; i < nf; i++)
    {
        Fracture& f = frac_list[i];        // f is alias for frac_list[i]
        xbeg = f.get_xbeg();
        ybeg = f.get_ybeg();
        xend = f.get_xend();
        yend = f.get_yend();
       
        if (ncr[i] < 2) continue;
        for (int nn = 0; nn < ncr[i]; nn++)
        {
            for (int n = 0; n < ncr[i] - 1; n++)
            {
                d1 = sqrt(pow(xcr[i][n] - xbeg, 2) + pow(ycr[i][n] - ybeg, 2));
                d2 = sqrt(pow(xcr[i][n + 1] - xbeg, 2) + pow(ycr[i][n + 1] - ybeg, 2));
                if (d2 < d1)
                {
                    xtem = xcr[i][n + 1];
                    ytem = ycr[i][n + 1];
                    xcr[i][n + 1] = xcr[i][n];
                    ycr[i][n + 1] = ycr[i][n];
                    xcr[i][n] = xtem;
                    ycr[i][n] = ytem;
                }
            }
        }
    }   //end i loop 


    int k = 0;

    for (int i = 0; i < nf; ++i) 
    {
        Fracture& f = frac_list[i];        
        for (int n = 0; n < ncr[i] + 1; ++n) 
        {
            itl = 0;
            itr = 0;

            if (n == 0)
            {
                xbeg = f.get_xbeg();
                ybeg = f.get_ybeg();
                itl = 1;
            }
            else
            {
                xbeg = xcr[i][n-1];    // change index y of xcr and ycr from n-1 to n because of warnings
                ybeg = ycr[i][n-1];
            }
                                 

            float xend = xcr[i][n];
            float yend = ycr[i][n];
            //I removed +1 in the condition !Sara float check later 
            if (n == ncr[i] )
            {
                xend = f.get_xend();
                yend = f.get_yend();
                itr = 1;
            }

            float tem = (frac_list[i].elem_no * sqrt(pow(xend - xbeg, 2) + pow(yend - ybeg, 2)) /
                sqrt(pow(f.get_xend() - f.get_xbeg(), 2) + pow(f.get_yend() - f.get_ybeg(), 2)));
            int num = (int)tem;

            if (tem - num > 0)  num++;
            if ((xbeg == xend && ybeg == yend) || num == 0 )
                continue;
           
            int kode = 5;
            float bvs = 0;
            float bvn = 0;
            int jmat = f.jmat;  //matf[i]; joint_id is joint_mat    Sara!
            int material = f.mat_no;
            int itype = 0;

            if (itl == 0) 
            {

                if (itr == 0)  itype = 0;
                else
                    if (itr == 1) itype = 1;
            }
            else
               {
                    if (itr == 0)  itype = -1;
                    else
                        if (itr == 1)  itype = 2;
                }
            k++;           

            try
            {
                file25 << num << " " << xbeg << " " << ybeg << " " << xend << " " << yend << " " << kode << " "
                    << itype << " " << jmat << " " << material << std::endl;
            }
            catch (std::ifstream::failure e)
            {
                std::cerr << "Exception opening/reading file:in fracture_reordering1\n";

            }

        }

    }   //end for 

    file25.seekg(0, std::ios::beg);
    nf = k;

    int  mat_no = 0, joint_mat = 0, elem_no = 0, bound_type = 5;
    if (!file25)
    {
        std::cerr << "Failed to open file for reading." << std::endl;
        return;
    }
    //reassign the values from file 25 to fracture objects    

    try
    {
        for (int i = 0; i < nf; i++)
        {
            Fracture& f = frac_list[i];       
            file25 >> elem_no>> xbeg >> ybeg >> xend >>
                yend >> bound_type >> itip[i] >> joint_mat >> mat_no;
            f.frac_reassign_values(elem_no,xbeg,ybeg,xend,yend,bound_type,mat_no,joint_mat);                         

        }
    }
    catch (std::ifstream::failure e)
    {
        std::cerr << "Exception opening/reading file:in fracture_reordering2\n";

    }
    return;
}




void check_fracture_cross()
{   

    float  xb1, xe1, xb2, xe2, yb1, yb2, ye1, ye2, xcross, ycross;   
    float d1, db1, de1, d2, db2, de2;
    float dtol = s5u.dtol;
    float tan1 = 0, tan2 = 0;    
   
    for (int i = 0; i < nf; i++)
    {
        Fracture& f = frac_list[i];       
        xb1 = f.get_xbeg();
        yb1 = f.get_ybeg();
        xe1 = f.get_xend();
        ye1 = f.get_yend();
        
        if (xe1 == xb1)
            tan1 = 10e20;
        else
            tan1 = (ye1 - yb1) / (xe1 - xb1);

        for (int j = 0; j < nf; j++)
        {
            if (j == i)
                continue;
            Fracture& f2 = frac_list[j];
            xb2 = f2.get_xbeg();
            yb2 = f2.get_ybeg();
            xe2 = f2.get_xend();
            ye2 = f2.get_yend();        
               

            if (xe2 == xb2)
                tan2 = 10e20;
            else
                tan2 = (ye2 - yb2) / (xe2 - xb2);
            if (tan1 == tan2)
                continue;


            if (xe2 == xb2)
                xcross = xe2;
            else
                xcross = (yb2 - yb1 + tan1 * xb1 - tan2 * xb2) / (tan1 - tan2);
            ycross = (tan1 * tan2 * (xb1 - xb2) + tan1 * yb2 - tan2 * yb1) / (tan1 - tan2);

            d1 = std::sqrt(std::pow(xe1 - xb1, 2) + std::pow(ye1 - yb1, 2));
            db1 = std::sqrt(std::pow(xcross - xb1, 2) + std::pow(ycross - yb1, 2));
            de1 = std::sqrt(std::pow(xcross - xe1, 2) + std::pow(ycross - ye1, 2));
            d2 = std::sqrt(std::pow(xe2 - xb2, 2) + std::pow(ye2 - yb2, 2));
            db2 = std::sqrt(std::pow(xcross - xb2, 2) + std::pow(ycross - yb2, 2));
            de2 = std::sqrt(std::pow(xcross - xe2, 2) + std::pow(ycross - ye2, 2));

            //distance less than half element length
            if (db1 < dtol)
            {
                f.take_xbeg(xcross);
                f.take_ybeg(ycross);
                xb1 = xcross;
                yb1 = ycross;
            }

            //distance less than half element length
            if (de1 < dtol)
            {
                f.take_xend(xcross);
                f.take_yend(ycross);
                xe1 = xcross;
                ye1 = ycross;
            }

            if (db2 < dtol)
            {
                f2.take_xbeg(xcross);
                f2.take_ybeg(ycross);
                xb2 = xcross;
                yb2 = ycross;
            }
            if (de2 < dtol)
            {
                f2.take_xend(xcross);
                f2.take_yend(ycross);
                xe2 = xcross;
                ye2 = ycross;

            }

            if (xcross < (min(xb1, xe1) - dtol) ||
                ycross < (min(yb1, ye1) - dtol) ||
                xcross >(max(xb1, xe1) + dtol)  ||
                ycross >(max(yb1, ye1) + dtol)  ||
                xcross < (min(xb2, xe2) - dtol) ||
                ycross < (min(yb2, ye2) - dtol) ||
                xcross >(max(xb2, xe2) + dtol)  ||
                ycross >(max(yb2, ye2) + dtol)) {
                continue;
            }
           
            xcr[i][ncr[i]] = xcross;
            ycr[i][ncr[i]] = ycross;
            ncr[i]++;
        } // end j loop
    }   //end i loop
   
    return;
}  //end check_fracture_cross()





void reset_boundary_bvs_bvn(float& bvs, float& bvn, int i, int kode)
{
    switch (kode)
    {
    case 1:
    case 11:
        bvn = bund_list[i].get_bvn();
        bvs = bund_list[i].get_bvs();
        break;
    case 2:
    case 12:
        bvn = bund_list[i].get_dn();
        bvs = bund_list[i].get_ds();
        break;
    case 3:
    case 13:
        bvn = bund_list[i].get_bvn();
        bvs = bund_list[i].get_ds();
        break;
    case 4:
    case 14:
        bvn = bund_list[i].get_dn();
        bvs = bund_list[i].get_bvs();
        break;
    case 7:
        bvn = bund_list[i].get_dn();
        bvs = bund_list[i].get_bvs();
    }
    return;
}




void populate_boudary_i(fstream& file25)
{
    float bvs, bvn, xb, yb, xe, ye, gsy, gny;
    int elm_no, kode, material;
    file25.seekg(0, std::ios::beg);

    for (int i = 0; i < nb; ++i)
    {                
        //read every boundary data from file to reorder it 
        file25 >> elm_no >>xb >>  yb >> xe >> ye >>
            kode>> bvs >> bvn >> gsy >> gny >> material ;
        
        Edge edge(material,elm_no, kode , xb, yb, xe, ye, 0, 0, 0, 0, gsy, gny);
        bund_list[i] = edge;
        Edge& b = bund_list[i];

        switch (kode)
        {
        case 1:
        case 11:
            b.take_bvn(bvn);
            b.take_bvs(bvs);
            break;
        case 2:
        case 12:
            b.take_ds(bvs);
            b.take_dn(bvn);
            break;
        case 3:
        case 13:
            b.take_bvn(bvn);
            b.take_ds(bvs);
            break;
        case 4:
        case 14:
            b.take_dn(bvn);
            b.take_bvs(bvs);
            break;

        case 7:
            b.take_dn(bvn);             //gost element
            b.take_bvs(bvs);
        }            
                              
    }
    file25.seekg(0, std::ios::beg);
    return;
}




void reorder_boundaries(fstream & file25)
{  
    file25.clear(); // Clear any error flags
    file25.seekp(0, std::ios::beg);
    for (int i = 0; i < nb; ++i)
    {
        Edge& obj = bund_list[i];   //obj is alias for bund_list[i]

        float xbeg = obj.get_xbeg();
        float ybeg = obj.get_ybeg();             

        if (ncr[i] < 2)
            continue;

        for (int nn = 0; nn < ncr[i]; ++nn) 
        {
            for (int n = 0; n < ncr[i] - 1; ++n)
            {
                float d1 = std::sqrt(std::pow(xcr[i][n] - xbeg, 2) + std::pow(ycr[i][n] - ybeg, 2));
                float d2 = std::sqrt(std::pow(xcr[i][n + 1] - xbeg, 2) + std::pow(ycr[i][n + 1] - ybeg, 2));

                if (d2 < d1)
                {
                    float xtem = xcr[i][n + 1];
                    float ytem = ycr[i][n + 1];
                    xcr[i][n + 1] = xcr[i][n];
                    ycr[i][n + 1] = ycr[i][n];
                    xcr[i][n] = xtem;
                    ycr[i][n] = ytem;
                }
            }
        }    
    }

    int k = 0;
    
    for (int i = 0; i < nb; ++i)
    {
        for (int n = 0; n < ncr[i] +1; ++n)
        {
            float xbeg, ybeg,xend, yend;
            Edge & b = bund_list[i];
            xbeg = xcr[i][n];  //Sara changed index
            ybeg = ycr[i][n];

            if (n == 0)
            {
                xbeg = b.get_xbeg();
                ybeg = b.get_ybeg();
            }
          
            xend = (n == ncr[i]) ? b.get_xend() : xcr[i][n];
            yend = (n == ncr[i]) ? b.get_yend() : ycr[i][n];
            float tem = (b.elem_no * std::sqrt(std::pow(xend - xbeg, 2) + std::pow(yend - ybeg, 2))) /
                std::sqrt(std::pow(b.get_xend() - b.get_xbeg(), 2) + std::pow(b.get_yend() - b.get_ybeg(), 2));

            int num = int(tem);
            //if (tem - num > 0)
            if (tem - int(tem) > 0.0009)
                num++;

            if ( (xbeg == xend && ybeg == yend) || (num == 0) )
               continue;
             

            int kode = b.bound_type;
            float bvs, bvn, gradsy, gradny;
            int material;

            reset_boundary_bvs_bvn(bvs,bvn,i,kode);
           
            gradsy = b.get_gsy();
            gradny = b.get_gny();
            material = b.mat_no;    //mbe

            k++;
            file25 << num << " " << xbeg << " " << ybeg << " " << xend << " " << yend << " " << kode <<
                " " << bvs << " " << bvn << " " << gradsy <<
                " " << gradny << " " << material << std::endl;

        }
    }
   // file25.close();
    nb = k;
    populate_boudary_i(file25);

    return;
}




void check_cross_boundaries()
{
    float xb1, yb1, xe1, ye1, xb2, yb2, xe2, ye2;
    float tan1, tan2;
    float& dtol = s5u.dtol;

    for (int i = 0; i < nb; i++)
    {
        Edge& obj = bund_list[i];   //obj is alias for bund_list[i]

        xb1 = obj.get_xbeg();
        yb1 = obj.get_ybeg();
        xe1 = obj.get_xend();
        ye1 = obj.get_yend();

        if (xe1 == xb1)
            tan1 = 10e20 * (ye1 - yb1);
        else
            tan1 = (ye1 - yb1) / (xe1 - xb1);

        ncr[i] = 0;

        for (int j = 0; j < nf; j++)
        {
            Fracture& f = frac_list[j];       
            xb2 = f.get_xbeg();
            yb2 = f.get_ybeg();
            xe2 = f.get_xend();
            ye2 = f.get_yend();

            if ((std::abs(xb2 - xb1) < dtol && std::abs(yb2 - yb1) < dtol) ||
                (std::abs(xb2 - xe1) < dtol && std::abs(yb2 - ye1) < dtol)) {
                if (itip[j] == 2)
                    itip[j] = 1;
                else if (itip[j] == -1)
                    itip[j] = 0;
                continue;
            }

            if ((std::abs(xe2 - xb1) < dtol && std::abs(ye2 - yb1) < dtol) ||
                (std::abs(xe2 - xe1) < dtol && std::abs(ye2 - ye1) < dtol)) {
                if (itip[j] == 2)
                    itip[j] = -1;

                else if (itip[j] == 1)
                    itip[j] = 0;
                continue;
            }

            if (xe2 == xb2)
                tan2 = 10e20 * (ye2 - yb2);
            else
                tan2 = (ye2 - yb2) / (xe2 - xb2);
                

            if (tan1 == tan2)  continue;               

            float xcross = (yb2 - yb1 + tan1 * xb1 - tan2 * xb2) / (tan1 - tan2);
            float ycross = (tan1 * tan2 * (xb1 - xb2) + tan1 * yb2 - tan2 * yb1) / (tan1 - tan2);

            float d2  =  sqrt((xe2 - xb2) * (xe2 - xb2) + (ye2 - yb2) * (ye2 - yb2));
            float db2 = sqrt((xcross - xb2) * (xcross - xb2) + (ycross - yb2) * (ycross - yb2));
            float de2 = sqrt((xcross - xe2) * (xcross - xe2) + (ycross - ye2) * (ycross - ye2));

            int & nfe = frac_list[j].elem_no;      //nfe is an alias 

            if ( db2 < (d2 / nfe / 2.) )
            {
                frac_list[j].take_xy_beg(xcross, ycross);      //assign xcross and ycross to fbeginx and fbegy
                
                if (itip[j] == 2)
                    itip[j] = 1;
                else if (itip[j] == -1)
                    itip[j] = 0;
            }

            if (de2 < (d2 / nfe / 2.))
            {
                frac_list[j].take_xy_end(xcross, ycross);     //assign xcross and ycross to fendx and fendy
               
                if (itip[j] == 2)
                    itip[j] = -1;
                else if (itip[j] == 1)
                    itip[j] = 0;
            }

            const float epsilon = 1e-6;


            if (xcross < min(xb1, xe1) - dtol - epsilon || ycross < min(yb1, ye1) - dtol - epsilon ||
                xcross > max(xb1, xe1) + dtol + epsilon || ycross > max(yb1, ye1) + dtol + epsilon ||
                xcross < min(xb2, xe2) - dtol - epsilon || ycross < min(yb2, ye2) - dtol - epsilon ||
                xcross > max(xb2, xe2) + dtol + epsilon || ycross > max(yb2, ye2) + dtol + epsilon)
                continue;

            
            xcr[i][ncr[i]] = xcross;
            ycr[i][ncr[i]] = ycross;
            ncr[i]++;

             //since fj changed we reassign all vars and used in computing ang.   
            // boundary i is not changed so xb1,yb1,xe1,ye1 no need to reassign again     //test Sara!
            xb2 = f.get_xbeg();
            yb2 = f.get_ybeg();
            xe2 = f.get_xend();
            ye2 = f.get_yend();

            float ang0 = atan2f((ye1 - yb1), (xe1 - xb1));
            float ange = atan2f((ye2 - yb1), (xe2 - xb1));
            float angb = atan2f((yb2 - yb1), (xb2 - xb1));        //diff bet frac j and bound i  coordinates

            if (angb > ang0 && ange < ang0) 
            {
                float tem = (nfe * sqrt((xcross - xe2) * (xcross - xe2) +
                    (ycross - ye2) * (ycross - ye2)) /
                    sqrt((xe2 - xb2) * (xe2 - xb2) +
                        (ye2 - yb2) * (ye2 - yb2)));

                nfe = (int) tem;
                if (tem - nfe > 0)
                    nfe++;
                f.take_xbeg(xcross);
                f.take_ybeg(ycross);
                if (itip[j] == 2)
                    itip[j] = 1;
                else if (itip[j] == -1)
                    itip[j] = 0;
            }

            else if (angb < ang0 && ange > ang0)
            {
                float tem = (nfe * sqrt((xcross - xb2) * (xcross - xb2) +
                    (ycross - yb2) * (ycross - yb2)) /
                    sqrt((xe2 - xb2) * (xe2 - xb2) +
                        (ye2 - yb2) * (ye2 - yb2)));

                nfe = (int) tem;
                if (tem - nfe > 0)
                    nfe++;
                f.take_xend(xcross);
                f.take_yend(ycross);
                if (itip[j] == 2)
                    itip[j] = -1;
                else if (itip[j] == 1)
                    itip[j] = 0;
            }
        }
    }

    return;
}



void populate_arc_i(int k,fstream & file25)
{
    na = k;
    float bvs, bvn;
    float  xc, yc, arcbeg, arcend, gsy = 0.0, gny = 0.0, r;
    int elm_no, kode, material;
    string lineData;
    file25.seekg(0, std::ios::beg);
    
    for (int i = 0; i < na; ++i)
    {
        //read every boundary data from file to reorder it 
        getline(file25, lineData);
        std::stringstream ss(lineData);
        ss >> elm_no >> xc >> yc >> r >> arcbeg >> arcend>>
            kode >> bvs >> bvn >> gsy >> gny >> material;

        Arch arc( xc, yc, r, arcbeg, arcend, 0, 0, 0,0, gsy, gny,  elm_no, material, kode);
        arc_list[i] = arc;
        Arch& a = arc_list[i];

        switch (kode)
        {
        case 1:
        case 11:
           a.take_sn(bvn);
            a.take_ss(bvs);
            break;

        case 2:
        case 12:
           a.take_ds(bvs);
            a.take_dn(bvn);
            break;

        case 3:
        case 13:
            a.take_sn(bvn);
            a.take_ds(bvs);
            break;

        case 4:
        case 14:
            a.take_dn(bvn);
            a.take_ss(bvs);
            break;

        case 7:
            a.take_dn(bvn);             //gost element
            a.take_ss(bvs);
            break;
        }

    }

    return;
}



void reset_arc_bvs_bvn( float& bvs, float& bvn, int i, int kode)
{
    Arch& arc = arc_list[i];
    switch (kode)
    {
    case 1:
    case 11:
        bvn = arc.get_arcsn();
        bvs = arc.get_arcss();
        break;
    case 2:
    case 12:
        bvn = arc.get_dn();
        bvs = arc.get_ds();
        break;
    case 3:
    case 13:
        bvn = arc.get_arcsn();
        bvs = arc.get_ds();
        break;
    case 4:
    case 14:
        bvn = arc.get_dn();
        bvs = arc.get_arcss();
    }
    return;
}




void arc_reordering(fstream& file25)
{    
    file25.clear(); // Clear any error flags
    file25.seekp(0, std::ios::beg);
    //reordering all defined arcs 
    for (int i = 0; i < na; ++i) 
    {
        if (ncr[i] < 2)  //Sara! 14.10 change 2 to 1
            continue;

        for (int nn = 0; nn < ncr[i]; ++nn) 
        {
            for (int n = 0; n < ncr[i] - 1; ++n) 
            {
                float d1 = acr[i][n];
                float d2 = acr[i][n + 1];
                if (d2 <= d1) 
                {
                    float tem = acr[i][n + 1];
                    acr[i][n + 1] = acr[i][n];
                    acr[i][n] = tem;
                }
            }
        }
    }

    int k = 0;
    for (int i = 0; i < na; ++i)
    {
        Arch& arc = arc_list[i];
        for (int n = 0; n < ncr[i] + 1; ++n)
        {
            float angb = (n == 0)? arc.get_arcbeg(): acr[i][n - 1];
            
            float ange = (n == ncr[i])? arc.get_arcend(): acr[i][n];
           
            float tem = arc.elem_no * (ange - angb) / ( arc.get_arcend() - arc.get_arcbeg());
            int num = (int)tem;   
            if (tem - int(tem) > 0.0009)
                num++;

            if (angb == ange || num == 0)
                continue;
            int kode = arc.bound_type;
            float bvn, bvs;
            reset_arc_bvs_bvn(bvs, bvn, i, kode);
                      
            float gradsy = arc.get_gsy();
            float gradny = arc.get_gny();
            int material = arc.mat_no;
            int itype = 0;
            int jmat = 0;
            k++;
            file25 << num << "  " << arc.get_arcx() << "  " << arc.get_arcy() << "  " << arc.get_arcr() << "  " << angb << "  " << ange << "  " << kode << "  " << bvs << 
                "  " << bvn <<
                "  " << gradsy << "  " << gradny << "  " << material << std::endl;
        }   
    }
    //file25.close();
    populate_arc_i(k,file25);

    return;
}




void final_wrap_up_for_Archs()
{
    for (int i = 0; i < na; ++i)
    {
        Arch& arc = arc_list[i];        //alias for each arc
        int nume = arc.elem_no;
        float xcen = arc.get_arcx();
        float ycen = arc.get_arcy();
        float diam = 2.0 * arc.get_arcr();
        float ang1 = arc.get_arcbeg();
        float ang2 = arc.get_arcend();
        int kode = arc.bound_type;
        float bvn, bvs;
        reset_arc_bvs_bvn(bvs, bvn, i, kode);
        float gradsy = arc.get_gsy();
        float gradny = arc.get_gny();
        //new above added//not sure about this just for checking test11
        for (int m = 0; m < nume; ++m)
        {
            float& dtol = s5u.dtol;       //alias for dtol
            float seta1 = ang1 + ((ang2 - ang1) / static_cast<float>(nume)) * static_cast<float>(m );
            float seta2 = ang1 + ((ang2 - ang1) / static_cast<float>(nume)) * static_cast<float>(m+1);
            float xbeg = xcen + 0.5 * diam * cosf(seta1);
            float ybeg = ycen + 0.5 * diam * sinf(seta1);
            float xend = xcen + 0.5 * diam * cosf(seta2);
            float yend = ycen + 0.5 * diam * sinf(seta2);

            if (symm.ksym == 0)
            {
                if (cosf(seta1) > 0)
                {
                    xbeg -= dtol / 1000.0;
                    ybeg -= dtol / 1000.0;
                }
                else if (cosf(seta1) <= 0)    //two equal in coditions Sara! < > else if or if
                {
                    xbeg += dtol / 1000.0;
                    ybeg += dtol / 1000.0;
                }
                if (sinf(seta2) > 0)
                {
                    xend -= dtol / 1000.0;
                    yend -= dtol / 1000.0;
                }
                else if (sinf(seta2) <= 0) {
                    xend += dtol / 1000.0;
                    yend += dtol / 1000.0;
                }
            }

            int num = 1;
            int jmat = 1; // Note: jmat here is meaningless
            int itype = 0; // Note: itype here is meaningless
            arc.def_boundary_elements_for_Geoform(num, xbeg, ybeg, xend, yend, bvs, bvn, gradsy, gradny, itype, jmat);  //need to change it probably Sara!
        }
    }
}



void final_wrap_up_for_Fracs()
{
    for (int i = 0; i < nf; ++i) 
    {
        int itype = itip[i];
        Fracture& f = frac_list[i];
        float xb = f.get_xbeg();
        float yb = f.get_ybeg();
        float xe = f.get_xend();
        float ye = f.get_yend();

        if (symm.ksym == 1 || symm.ksym == 4) 
        {
            if (xb == symm.xsym)
            {
                if (itype == 2 || itype == 1)
                    itip[i] = 1;
                else if (itype == -1)
                    itip[i] = 0;
            }
            else
            {
                if (xe == symm.xsym)
                {
                    if (itype == 2 || itype == -1)
                        itip[i] = -1;
                    else if (itype == 1)
                        itip[i] = 0;
                    continue;
                }
            }
        }
       
        if (symm.ksym == 2 || symm.ksym == 4) 
        {
            if (yb == symm.ysym)
            {
                if (itype == 2 || itype == 1)
                    itip[i] = 1;
                else if (itype == -1)
                    itip[i] = 0;
                
            }
            else if (ye == symm.ysym)
            {
                if (itype == 2 || itype == -1)
                    itip[i] = -1;
                else if (itype == 1)
                    itip[i] = 0;
            }
            
        }

        if (symm.ksym == 3 || symm.ksym == 4) 
        {
            if (xb == symm.xsym && yb == symm.ysym)
            {
                if (itype == 2 || itype == 1)
                    itip[i] = 1;
                if (itype == -1)
                    itip[i] = 0;
            }
            else if (xe == symm.xsym && ye == symm.ysym)
            {
                if (itype == 2 || itype == -1)
                    itip[i] = -1;
                else if (itype == 1)
                    itip[i] = 0;
            }       
            
        }       
        f.bound_type = 5;
        f.def_boundary_elements_for_Geoform(f.elem_no, xb, yb, xe, ye,  0.0, 0.0, 0.0, 0.0, itip[i],f.jmat);   //need to change it Sara!
    }
}




void final_wrap_up()
{        
        //iterate on all boundaries and model each as a set of elements
        for (int i = 0; i < nb; ++i) 
        {
            Edge & obj = bund_list[i];
            int kode = obj.bound_type;
            float bvn, bvs;
            reset_boundary_bvs_bvn(bvs, bvn, i, kode);
            float gradsy = obj.get_gsy();
            float gradny = obj.get_gny();

            obj.def_boundary_elements_for_Geoform( obj.elem_no, obj.get_xbeg(), obj.get_ybeg(), obj.get_xend(), obj.get_yend(), bvs,
                bvn, gradsy, gradny, 0, 0);      //do we need to send the attr of the object to the method of that class while we access to all of it Sara!
        }

        //final check for archs
        final_wrap_up_for_Archs();
        //final check for Fracs
        final_wrap_up_for_Fracs();
        
        //------------------------------Linear Interface----------------------------------
        //iterate on all linear interfaces defined in input
            for (int i = 0; i < npli; ++i)
            {
                Edge_interface& ei = lin_intrfce_list[i];  //alias for linear interface i
                float xbeg = ei.get_xbeg();
                float ybeg = ei.get_ybeg();
                float xend = ei.get_xend();
                float yend = ei.get_yend();
                int mat2 = ei.get_pos_matno();
                int mat1 = ei.get_neg_matno();
                int num = ei.get_elemno();
                Interface(num, xbeg, ybeg, xend, yend, mat1, mat2);
            }

        // -------------------------------Arch interface -----------------------------
        for (int i = 0; i < nq; ++i)
        {
            Arch_interface& ai  = arc_intrface_list[i];  //alias for linear interface i
            int nume = ai.get_elenum();
            float xcen = ai.get_xcent();
            float ycen = ai.get_ycent();
            float diam = ai.get_diameter();
            float ang1 = ai.get_beg_ang();
            float ang2 = ai.get_end_ang();
            int mat2 =ai.get_pos_mat();
            int mat1 = ai.get_neg_mat();
           
            for (int m = 0; m < nume; ++m)
            {
                float seta1 = ang1 + ((ang2 - ang1) / static_cast<float>(nume)) * static_cast<float>(m);
                float seta2 = ang1 + ((ang2 - ang1) / static_cast<float>(nume)) * static_cast<float>(m+1);
                float xbeg = xcen +  diam * cosf(seta1);
                float ybeg = ycen + diam * sinf(seta1);
                float xend = xcen +  diam * cosf(seta2);
                float yend = ycen +  diam * sinf(seta2);
                Interface(1, xbeg, ybeg, xend, yend, mat1, mat2);  // interface function defined in source.cpp
            }
        }

        return;        
}





void cross_arc_with_boundaries(int i, float  xc, float yc, float angb, float ange, float r)
{
    /*iterate over all boundaries in input function 
     defined before and
    // check the possibility of crossing the arcs */


    for (int j = 0; j < nb; ++j)
    {
        Edge& obj = bund_list[j];   

        float xb = obj.get_xbeg();
        float yb = obj.get_ybeg();
        float xe = obj.get_xend();
        float ye = obj.get_yend();


        float db = sqrt(pow(xc - xb, 2) + pow(yc - yb, 2));
        float ang = static_cast<float>(atan2(yb - yc, xb - xc));

        if (!(db > 1.05 * r || db < 0.95 * r) && (!(ang < (angb - 1e-5) || ang >(ange + 1e-5))))
        {
            
            acr[i][ncr[i]] = ang;
            ncr[i]++;
            obj.take_xbeg(xc + cosf(ang) * r);
            obj.take_ybeg(yc + sinf(ang) * r);
        }

        //label750:
        float de = sqrt(powf(xc - xe, 2) + powf(yc - ye, 2));
        if (de > 1.05 * r || de < 0.95 * r)
            continue;

        ang = static_cast<float>(atan2(ye - yc, xe - xc));
        if (ang < (angb - 1e-5) || ang >(ange + 1e-5))
            continue;

       
        acr[i][ncr[i]] = ang;
        ncr[i]++;
        obj.take_xend(xc + cosf(ang) * r);
        obj.take_yend(yc + sinf(ang) * r);
    }

    return;               
}



void cross_arc_with_fractures(int i, float  xc, float yc, float angb, float ange, float r)
{
    //check the potential crossing between each arc and all of fractures

    for (int j = 0; j < nf; ++j)
    {
        Fracture& f = frac_list[j];        // f is alias for frac_list[i]
        float xb = f.get_xbeg();
        float yb = f.get_ybeg();
        float xe = f.get_xend();
        float ye = f.get_yend();

        float db = sqrt((xc - xb)* (xc - xb) + (yc - yb)* (yc - yb));
        float ang = static_cast<float>(atan2(yb - yc, xb - xc));
        if (!(db > 1.05 * r || db < 0.95 * r || ang < (angb - 1e-5) || ang >(ange + 1e-5) ) )
        {
            acr[i][ncr[i]] = ang;
            ncr[i]++;
            f.take_xbeg(xc + cosf(ang) * r);
            f.take_ybeg(yc + sinf(ang) * r);

            if (itip[j] == 2)
                itip[j] = 1;
            else if (itip[j] == -1)
                itip[j] = 0;
           
        }
        float de = sqrt((xc - xe)* (xc - xe) + (yc - ye)* (yc - ye));
        if (de > 1.05 * r || de < 0.95 * r)
            continue;

        ang = static_cast<float>(atan2(ye - yc, xe - xc));
        if (ang < (angb - 1e-5) || ang >(ange + 1e-5))
            continue;

        acr[i][ncr[i]] = ang;
        ncr[i]++;
        frac_list[j].take_xend(xc + cosf(ang) * r);
        frac_list[j].take_yend(yc + sinf(ang) * r);

        if (itip[j] == 2)
            itip[j] = -1;
        else if (itip[j] == 1)
            itip[j] = 0;
    }
    return;
}



     
void check_cross_arcs(fstream& file25)
{

    //because of this location diff in test12 need to keep it
    pi = 4.0 * atan(1.0); //3.14159      
    for (int i = 0; i < na; ++i)
    {
        Arch& arc = arc_list[i];        //alias for arc
        float xc = arc.get_arcx();
        float yc = arc.get_arcy();
        float r = arc.get_arcr();

        if ((arc.get_arcend() - arc.get_arcbeg()) == 2.0 * pi)
        {
            arc.take_arcbeg(-pi);
            arc.take_arcend(pi);
        }
        float angb = arc.get_arcbeg();
        float ange = arc.get_arcend();

        ncr[i] = 0;
        //check the potential crossing between each arc and all of fractures
        cross_arc_with_fractures(i, xc, yc, angb, ange, r);
        //check the potential crossing between each arc and all of boundaries
        cross_arc_with_boundaries(i, xc, yc, angb, ange,r);
    }

    arc_reordering(file25);
    return;
}





void inputcheck()
{
    
   /*---------------------------------------------------------------------- -
         check crossing of fracturs, boundaries and arcs
    ---------------------------------------------------------------------- -*/   

    std::fstream file25("temp001.dat", std::ios::in | std::ios::out | std::ios::trunc);
    if (!file25.is_open()) {
        std::cerr << "Error opening file temp001!" << std::endl;
        return;
    }
    check_fracture_cross();
    reorder_fractures(file25);
    check_cross_boundaries();
    reorder_boundaries(file25);
    check_cross_arcs(file25);
    if (restor_flg)
    {
        if (na>0)
            check_boreholes();
        if(nb>3)
            check_rectangle();
    }       
    final_wrap_up();
    return;
}
