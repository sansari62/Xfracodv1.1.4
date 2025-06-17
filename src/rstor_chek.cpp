#include<stdafx.h>
#include <CommonPara.h>
#include<GeologicalForm.h>
#include <optional>
#include<Tip.h>

using namespace CommonPara_h::comvar;


struct Point {
    float x, y;
};


void fix_tip_pointer(int m,int new_numbe)
{
    for (int k = 0; k < no; ++k)
    {
        Tip& t = tips[k];
        if (t.mpointer == m)
            t.mpointer = new_numbe;
    }
    return;
}



int pointInSector(const Point& p, float xc, float yc, float R, float ang1_rad, 
    float ang2_rad) {
    float dx = p.x - xc;
    float dy = p.y - yc;
    float R2 = R * R;
    const float EPS = 1e-6f;

    float dist2 = dx * dx + dy * dy;
    if (dist2 > R * R + EPS)   
        return 1;
   // if (dist2 < R * R - EPS) //→ inside
//return 0;
    //if (fabs(dist2 - R * R) < EPS) //→ on arc
      //  return 2;

    // Compute angle from center to point in radians
    float angle_rad = std::atan2(dy, dx);
    if (angle_rad < 0.0f)
        angle_rad += 2.0f * pi;

    // Normalize angles to [0, 2pi)
    float ang1 = fmodf(ang1_rad + 2.0f * pi, 2.0f * pi);
    float ang2 = fmodf(ang2_rad + 2.0f * pi, 2.0f * pi);

    bool angle_inside;   

    // Standard angular inclusion test
    if (std::abs(ang1 - ang2) < EPS) {
        angle_inside = true;  // full circle
    }
    else if (ang1 < ang2) {
        angle_inside = angle_rad >= ang1 && angle_rad <= ang2;
    }
    else {
        angle_inside = angle_rad >= ang1 || angle_rad <= ang2;
    }
    if (!angle_inside)
        return 1; // angle not inside → outside

    if (std::abs(dist2 - R2) < EPS)
        return 2; // on arc    

    return 0; // inside
}





bool computeIntersectionWithCircleSector(
    const Point& p1, const Point& p2, float xc, float yc, float R,
    float ang1_deg, float ang2_deg, Point& ipOut)
{
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float fx = p1.x - xc;
    float fy = p1.y - yc;

    float a = dx * dx + dy * dy;
    float b = 2 * (fx * dx + fy * dy);
    float c = fx * fx + fy * fy - R * R;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < -1e-6f) return false;  // Definitely no intersection
    if (discriminant < 0) discriminant = 0;   // Treat near-zero as tangent
    //if (discriminant < 0) return false;

    discriminant = std::sqrt(discriminant);
    float t1 = (-b - discriminant) / (2 * a);
    float t2 = (-b + discriminant) / (2 * a);

    for (float t : {t1, t2}) {
        if (t < 0 || t > 1) continue;
        Point ip{ p1.x + t * dx, p1.y + t * dy };
        if (pointInSector(ip, xc, yc, R, ang1_deg, ang2_deg)== 2) {
            ipOut = ip;
            return true;
        }
    }

    return false;
}





void addClippedElement(const Point& a, const Point& b,
    int m, int& new_numbe)
{
    const BoundaryElement& orig_elem = elm_list[m];
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f) return; // skip degenerate

    BoundaryElement newelem;
    newelem.xm = 0.5f * (a.x + b.x);
    newelem.ym = 0.5f * (a.y + b.y);
    newelem.a = len / 2.0f;
    newelem.cosbet = dx / len;
    newelem.sinbet = dy / len;
    newelem.kod = orig_elem.kod;
    newelem.mat_no = orig_elem.mat_no;

    if (!isNewElementUnique(newelem)) return;

    elm_list[new_numbe] = newelem;
    float ratio = newelem.a / orig_elem.a;

    b_elm[new_numbe] = b_elm[m];
    b_elm[new_numbe].force1 *= ratio;
    b_elm[new_numbe].force2 *= ratio;

    joint[new_numbe] = joint[m];
    s4.b0[2 * new_numbe] = s4.b0[2 * m];
    s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];

    fix_tip_pointer(m, new_numbe);
    new_numbe++;

}




void findallintersect(const Point& p1, const Point& p2, float xc, float yc, float R,
    float ang1, float ang2,int& new_numbe, int m)
{
    // Step 1: Find all intersections
    Point ip1, ip2;
    bool found1 = computeIntersectionWithCircleSector(p1, p2, xc, yc, R, ang1, ang2, ip1);
    bool found2 = computeIntersectionWithCircleSector(p2, p1, xc, yc, R, ang1, ang2, ip2);

    if (found1 && found2) {
        // Check order along the segment
        if ((ip1.x - p1.x) * (p2.x - p1.x) + (ip1.y - p1.y) * (p2.y - p1.y) >
            (ip2.x - p1.x) * (p2.x - p1.x) + (ip2.y - p1.y) * (p2.y - p1.y)) {
            std::swap(ip1, ip2);
        }

        // Create [p1, ip1]
        addClippedElement(p1, ip1, m, new_numbe);

        // Create [ip2, p2]
        addClippedElement(ip2, p2, m, new_numbe);
    }
    return; 

}



void clipBoundaryElements(    
    float xc, float yc, float R, float ang1_deg, float ang2_deg)
{
    int new_numbe = 0;
    for (int m = 0; m < numbe; ++m) {
        auto& elem = elm_list[m];
        // Reconstruct endpoints
        float dx = elem.a * elem.cosbet;
        float dy = elem.a * elem.sinbet;

        Point p1{ elem.xm - dx, elem.ym - dy };
        Point p2{ elem.xm + dx, elem.ym + dy };

        int p1_status = pointInSector(p1, xc, yc, R, ang1_deg, ang2_deg);
        int p2_status = pointInSector(p2, xc, yc, R, ang1_deg, ang2_deg);
        bool remove =
            (p1_status == 0 && p2_status == 0) ||                     // both inside
            (p1_status == 2 && p2_status == 2) ||                     // both on arc
            (p1_status == 0 && p2_status == 2) ||                     // one inside, one on arc
            (p1_status == 2 && p2_status == 0);                       // one inside, one on arc

        if (remove) {
            // Entire element is inside → skip
            int merge = 0;
            specialLabel_200(merge,m);
            continue;
        }
        else if (p1_status ==1 && p2_status ==1) {

            Point intersection;
            if (computeIntersectionWithCircleSector(p1, p2, xc, yc, R, ang1_deg,
                ang2_deg, intersection))
            {
                std::cout << "Element crosses the arc. Should clip\n";
                findallintersect(p1, p2, xc, yc, R, ang1_deg, ang2_deg,new_numbe, m);
            }
            // Entirely outside → keep as-is
            else
            {
                fix_tip_pointer(m, new_numbe);
                elm_list[new_numbe] = elem;
                b_elm[new_numbe] = b_elm[m];
                new_numbe++;
            }
        }
        else {
            if (p1_status == 1 && p2_status == 2 || p1_status == 2 && p2_status == 1)
            {
                fix_tip_pointer(m, new_numbe);
                elm_list[new_numbe] = elem;
                b_elm[new_numbe] = b_elm[m];
                new_numbe++;
            }
            else
            {

                // Partially inside → clip and recompute center           
                Point intersection;
                if (!computeIntersectionWithCircleSector(p1, p2, xc, yc, R, ang1_deg, ang2_deg, intersection))
                {
                    fix_tip_pointer(m, new_numbe);
                    elm_list[new_numbe] = elem;
                    b_elm[new_numbe] = b_elm[m];
                    new_numbe++;
                    continue;
                }

                if (p1_status == 0) {
                    p1 = intersection;
                }
                else {
                    p2 = intersection;
                }
                // Recompute center and direction
                float cx = 0.5f * (p1.x + p2.x);
                float cy = 0.5f * (p1.y + p2.y);
                float dx = p2.x - p1.x;
                float dy = p2.y - p1.y;
                float len = std::sqrt(dx * dx + dy * dy);
                if (len < 1e-6f) {
                    // Entire element is inside → skip
                    int merge = 0;
                    specialLabel_200(merge, m);
                    continue;
                }
                //add new element
                BoundaryElement newelement;
                newelement.xm = cx;
                newelement.ym = cy;
                newelement.a = len / 2;
                newelement.sinbet = dy / len;
                newelement.cosbet = dx / len;
                bool unique = isNewElementUnique(newelement);
                if (!unique)
                {
                    new_numbe--;
                    continue;
                }
                newelement.kod = elm_list[m].kod;
                newelement.mat_no = elm_list[m].mat_no;
                elm_list[new_numbe] = newelement;
                b_elm[new_numbe].force1 = b_elm[m].force1 * elm_list[new_numbe].a / elem.a;
                b_elm[new_numbe].force2 = b_elm[m].force2 * elm_list[new_numbe].a / elem.a;
                b_elm[new_numbe].aks = b_elm[m].aks;
                b_elm[new_numbe].akn = b_elm[m].akn;
                b_elm[new_numbe].phi = b_elm[m].phi;
                b_elm[new_numbe].phid = b_elm[m].phid;

                b_elm[new_numbe].coh = b_elm[m].coh;
                joint[new_numbe].aperture0 = joint[m].aperture0;
                joint[new_numbe].aperture_r = joint[m].aperture_r;
                s4.b0[2 * new_numbe] = s4.b0[2 * m];
                s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
                fix_tip_pointer(m, new_numbe);
                new_numbe++;
            }

        }
    }
    numbe = new_numbe;
    arrangetip();
    return;
}


void check_boreholes()
{
    for (int i = 0; i < na; ++i)
    {
        Arch& arc = arc_list[i];
        float xc = arc.get_arcx();
        float yc = arc.get_arcy();
        float r = arc.get_arcr();
        clipBoundaryElements(xc, yc, r,arc.get_arcbeg(), arc.get_arcend());
    }
}