#include<stdafx.h>
#include<Remov_preexisting_BEs.h>


int  if_tip_element(int m)
{
    for (int k = 0; k < no; ++k)
    {
        Tip& t = tips[k];
        if (t.mpointer == m)
            return k;
    }
    return -1;
}



void fix_tip_pointer1(int new_numbe, int k,int direction)
{    
     tips[k].mpointer = new_numbe;
            if (direction == 0)
            {               
                tips[k].xbe = elm_list[new_numbe].xm - 3.0 * elm_list[new_numbe].a * elm_list[new_numbe].cosbet;
                tips[k].ybe = elm_list[new_numbe].ym - 3.0 * elm_list[new_numbe].a * elm_list[new_numbe].sinbet;
                tips[k].xen = tips[k].xbe + 2. * elm_list[new_numbe].a * elm_list[new_numbe].cosbet;
                tips[k].yen = tips[k].ybe + 2. * elm_list[new_numbe].a * elm_list[new_numbe].sinbet;
                tips[k].dl = 2. * elm_list[new_numbe].a;
                tips[k].costem = elm_list[new_numbe].cosbet;
                tips[k].sintem = elm_list[new_numbe].sinbet;
            }
              
            else //direction=1 right tip element
            {               
                tips[k].xbe = elm_list[new_numbe].xm + elm_list[new_numbe].a * elm_list[new_numbe].cosbet;
                tips[k].ybe = elm_list[new_numbe].ym + elm_list[new_numbe].a * elm_list[new_numbe].sinbet;
                tips[k].xen = tips[k].xbe + 2. * elm_list[new_numbe].a * elm_list[new_numbe].cosbet;
                tips[k].yen = tips[k].ybe + 2. * elm_list[new_numbe].a * elm_list[new_numbe].sinbet;
                tips[k].dl = 2. * elm_list[new_numbe].a;
                tips[k].ityp = 1;
                tips[k].costem = elm_list[new_numbe].cosbet;
                tips[k].sintem = elm_list[new_numbe].sinbet;
            }
    return;
}


//void setting_newelement(BoundaryElement& newelement, int m, int new_numbe, int p1_stat, int p2_stat)
//{
//    newelement.kod = elm_list[m].kod;
//    newelement.mat_no = elm_list[m].mat_no;
//    elm_list[new_numbe] = newelement;
//    b_elm[new_numbe].force1 = b_elm[m].force1 * elm_list[new_numbe].a / elm_list[m].a;
//    b_elm[new_numbe].force2 = b_elm[m].force2 * elm_list[new_numbe].a / elm_list[m].a;
//    b_elm[new_numbe].aks = b_elm[m].aks;
//    b_elm[new_numbe].akn = b_elm[m].akn;
//    b_elm[new_numbe].phi = b_elm[m].phi;
//    b_elm[new_numbe].phid = b_elm[m].phid;
//
//    b_elm[new_numbe].coh = b_elm[m].coh;
//    joint[new_numbe].aperture0 = joint[m].aperture0;
//    joint[new_numbe].aperture_r = joint[m].aperture_r;
//    s4.b0[2 * new_numbe] = s4.b0[2 * m];
//    s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
//   // fix_tip_pointer1(m, new_numbe, p1_stat, p2_stat);
//}


bool check_segmnt_valid(const Point& a, const Point& b,
    int m, BoundaryElement& newelem)
{
    const BoundaryElement& orig_elem = elm_list[m];
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float len = std::sqrt(dx * dx + dy * dy);
    bool valid = true;
    if (len / orig_elem.a < 0.125f)
    {
        //int merge = 0;
        //specialLabel_200(merge, m);
        valid = false; // skip degenerate
    }
    else
    {
        // BoundaryElement newelem;
        newelem.xm = 0.5 * (a.x + b.x);
        newelem.ym = 0.5 * (a.y + b.y);
        newelem.a = len / 2.0;
        newelem.cosbet = dx / len;
        newelem.sinbet = dy / len;
        newelem.kod = orig_elem.kod;
        newelem.mat_no = orig_elem.mat_no;
       // if (!isNewElementUnique(newelem)) valid = false;
    }
    return valid;
}



bool isALeftOfB(const Point& A, const Point& B) {
    if (A.x < B.x)
        return true;
    if (A.x > B.x)
        return false;

    // Same x-coordinate
    if (A.x < 0) {
        // Left side: higher Y comes first
        return A.y > B.y;
    }
    else {
        // Right side (or center): lower Y comes first
        return A.y < B.y;
    }
}


//bool isALeftOfB(const Point& A, const Point& B) {
//    if (A.x < B.x)
//        return true;
//    if (A.x > B.x)
//        return false;
//    // If X equal, compare Y
//    return A.y < B.y;
//}

void addClippedElement2(int m, int new_numbe, BoundaryElement& new_elem, int first, std::vector<new_element_para>& vec) {
    if (first == 1)
    {
        elm_list[new_numbe] = new_elem;
        float ratio = new_elem.a / elm_list[m].a;
        b_elm[new_numbe] = b_elm[m];
        b_elm[new_numbe].force1 *= ratio;
        b_elm[new_numbe].force2 *= ratio;

        joint[new_numbe] = joint[m];
        s4.b0[2 * new_numbe] = s4.b0[2 * m];
        s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
       
    }
    else
    {
        new_element_para t;
        t.b01 = s4.b0[2 * m];
        t.b02 = s4.b0[2 * m + 1];
        t.j = joint[m];
        t.new_el = new_elem;
        t.ratio = new_elem.a / elm_list[m].a;
        t.be1 = b_elm[m];       
        vec.push_back(t);
    }
    return;
}



void one_point_intersection(Point outside_pt, Point clipped, int m, int& new_numbe, int tip_index,
    std::vector<new_element_para>& vec) {

    BoundaryElement new_elem;
    bool segm = false;
    bool direct;
    if (isALeftOfB(clipped, outside_pt))
    {
        segm = check_segmnt_valid(clipped, outside_pt, m, new_elem);
        direct = 1;//right
    }
    else
    {
        segm = check_segmnt_valid(outside_pt, clipped, m, new_elem);
        direct = 0; //left
    }
    //one segment remained
    if (segm) {
        addClippedElement2(m, new_numbe, new_elem, 1,vec);
        if (tip_index != -1)
        {
            if (direct == 0 && tips[tip_index].ityp == -1)
                fix_tip_pointer1(new_numbe, tip_index, direct);

            else if (direct == 1 && tips[tip_index].ityp == 1)
                fix_tip_pointer1(new_numbe, tip_index, direct);
            else
                tips[tip_index].ityp = 0;
        }
        new_numbe++;
    }
    else
    {
        //the remaining segment after cross in not eligible
        if (tip_index != -1)
        {
            tips[tip_index].ityp = 0;
        }
    }
}


void two_intersection_case(Point p1, Point p2, Point ip1, Point ip2, int m, int& new_numbe, int tip_index, 
    std::vector<new_element_para>& vec)
{
    if ((ip1.x - p1.x) * (p2.x - p1.x) + (ip1.y - p1.y) * (p2.y - p1.y) >
        (ip2.x - p1.x) * (p2.x - p1.x) + (ip2.y - p1.y) * (p2.y - p1.y)) {
        std::swap(ip1, ip2);
    }

    BoundaryElement new_elem1, new_elem2;
    bool segm1 = check_segmnt_valid(p1, ip1, m, new_elem1);
    bool segm2 = check_segmnt_valid(ip2, p2, m, new_elem2);
    if (segm1 && segm2)
    {
        float a = elm_list[m].a;
        addClippedElement2(m, new_numbe, new_elem1, 1,vec);
        // Create [ip2, p2]
        second_clipped = true;
        addClippedElement2(m, new_numbe, new_elem2, 2,vec);
        new_element_para& last = vec.back();
        last.ratio = new_elem2.a / a;
       
        if (tip_index != -1)
        {
            if (tips[tip_index].ityp == -1)
                fix_tip_pointer1(new_numbe, tip_index, 0);
            else if (tips[tip_index].ityp == 1)
                last.tip_indx = tip_index;

        }
        new_numbe++;
    }
    else if (segm1) {
        addClippedElement2(m, new_numbe, new_elem1, 1,vec);
        if (tip_index != -1)
        {
            if (tips[tip_index].ityp == -1)
                tips[tip_index].mpointer = new_numbe;
        }
        new_numbe++;
    }
    else  if (segm2) {
        addClippedElement2(m, new_numbe, new_elem2, 1,vec);
        if (tip_index != -1)
        {
            if (tips[tip_index].ityp == 1)
                tips[tip_index].mpointer = new_numbe;
        }
        new_numbe++;
    }
    else
    {
        if (tip_index != -1)
        {
            tips[tip_index].ityp = 0;
        }
    }
    return;
}
