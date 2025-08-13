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


bool check_segmnt_valid(const Point& a, const Point& b,
    int m, BoundaryElement& newelem)
{
    const BoundaryElement& orig_elem = elm_list[m];
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float len = std::sqrt(dx * dx + dy * dy);
    float new_a = len / 2;
    bool valid = true;
    if (new_a / orig_elem.a < 0.125f)
    {        
        valid = false; // skip degenerate
    }
    else
    {
        newelem.xm = 0.5 * (a.x + b.x);
        newelem.ym = 0.5 * (a.y + b.y);
        newelem.a = len / 2.0;
        newelem.cosbet = dx / len;
        newelem.sinbet = dy / len;
        newelem.kod = orig_elem.kod;
        newelem.mat_no = orig_elem.mat_no;
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



void  update_s4(int m, int new_numbe)
{
    int is = 2 * new_numbe;
    int in =  is+1;
    int it = 2 * m;
    int itt = it + 1;

    s4.b0[is] = s4.b0[it];
    s4.b0[in] = s4.b0[itt];
    s4.d0[is] = s4.d0[it];
    s4.d0[in] = s4.d0[itt];
    s4.b0_old[is] = s4.b0_old[it];
    s4.b0_old[in] = s4.b0_old[itt];

    s4.b[is] = s4.b[it];
    s4.b[in] = s4.b[itt];
    s4.b1[is] = s4.b1[it];
    s4.b1[in] = s4.b1[itt];
    s4.d[is] = s4.d[it];
    s4.d[in] = s4.d[itt];       

    s4.c[is] = s4.c[it];
    s4.c[in] = s4.c[itt];

    s4.c_d[is] = s4.c_d[it];
    s4.c_d[in] = s4.c_d[itt];

    s4.c_s[is] = s4.c_s[it];
    s4.c_s[in] = s4.c_s[itt];   

    s4.df0[is] = s4.df0[it];
    s4.df0[in] = s4.df0[itt];

}



void addClippedElement2(int m, int new_numbe, BoundaryElement& new_elem, int first, std::vector<new_element_para>& vec) {
    if (first == 1)
    {
        elm_list[new_numbe] = new_elem;
        float ratio = new_elem.a / elm_list[m].a;
        b_elm[new_numbe] = b_elm[m];
        b_elm[new_numbe].force1 *= ratio;
        b_elm[new_numbe].force2 *= ratio;

        joint[new_numbe] = joint[m];       
        update_s4(m,new_numbe);
       
    }
    else
    {
        new_element_para t;      
        t.j = joint[m];
        t.new_el = new_elem;
        t.ratio = new_elem.a / elm_list[m].a;
        t.be1 = b_elm[m];    
        t.s4_indx = new_numbe;
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
                fix_tip_pointer1(new_numbe, tip_index, 0);
                //tips[tip_index].mpointer = new_numbe;
            else
                tips[tip_index].ityp = 0;
        }
        new_numbe++;
    }
    else  if (segm2) {
        addClippedElement2(m, new_numbe, new_elem2, 1,vec);
        if (tip_index != -1)
        {
            if (tips[tip_index].ityp == 1)
                //tips[tip_index].mpointer = new_numbe;
                fix_tip_pointer1(new_numbe, tip_index, 1);
            else
                tips[tip_index].ityp = 0;

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
