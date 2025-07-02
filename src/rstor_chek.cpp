#include<stdafx.h>
#include<rstor_chek.h>
#include <CommonPara.h>
#include<GeologicalForm.h>
#include <optional>
#include<Tip.h>
#include <algorithm>

using namespace CommonPara_h::comvar;

bool second_clipped = false;

struct new_element_para
{
    BoundaryElement new_el;
    BE be1;
    float ratio = 0;
    Joint j;
    float b01 = 0, b02 = 0;
};
std::vector<new_element_para>  new_elements;



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
    const float EPS = 1e-5f;
    const float TWO_PI = 2.0f * pi;


    float dist2 = dx * dx + dy * dy;
    if (dist2 - R * R > EPS)
        return 1;
   
    // Compute angle from center to point in radians
    float angle_rad = std::atan2(dy, dx);
    if (angle_rad < 0.0f)
        angle_rad += 2.0f * pi;

    // Normalize angles to [0, 2pi)
    float ang1 = fmodf(ang1_rad + 2.0f * pi, 2.0f * pi);
    float ang2 = fmodf(ang2_rad + 2.0f * pi, 2.0f * pi);

    bool angle_inside;   

    // Standard angular inclusion test
    //if (std::abs(ang1 - ang2) < EPS) {
    //    angle_inside = true;  // full circle
    //}
    float delta = fmodf((ang2 - ang1 + TWO_PI), TWO_PI);
    bool full_circle = std::abs(delta) < EPS || std::abs(delta - TWO_PI) < EPS ||  std::abs(ang1 - ang2) < EPS;


    if (full_circle) {
        angle_inside = true;
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


bool computeSegmentIntersection(
    const Point& p1, const Point& p2,
    float qx1, float qy1, float qx2, float qy2,
    Point& intersection)
{
    float s1_x = p2.x - p1.x;
    float s1_y = p2.y - p1.y;
    float s2_x = qx2 - qx1;
    float s2_y = qy2 - qy1;

    float denom = (-s2_x * s1_y + s1_x * s2_y);
    if (std::abs(denom) < 1e-8f)
        return false; // Parallel

    float s = (-s1_y * (p1.x - qx1) + s1_x * (p1.y - qy1)) / denom;
    float t = (s2_x * (p1.y - qy1) - s2_y * (p1.x - qx1)) / denom;

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        // Intersection detected
        intersection.x = p1.x + (t * s1_x);
        intersection.y = p1.y + (t * s1_y);
        return true;
    }
    else return false;

}


void addClippedElement(const Point& a, const Point& b,
    int m, int& new_numbe, int first)
{
    const BoundaryElement& orig_elem = elm_list[m];
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f)
    {
        int merge = 0;
        specialLabel_200(merge, m);
        return; // skip degenerate
    }

    BoundaryElement newelem;
    newelem.xm = 0.5f * (a.x + b.x);
    newelem.ym = 0.5f * (a.y + b.y);
    newelem.a = len / 2.0f;
    newelem.cosbet = dx / len;
    newelem.sinbet = dy / len;
    newelem.kod = orig_elem.kod;
    newelem.mat_no = orig_elem.mat_no;   
    if (!isNewElementUnique(newelem)) return;
    if (first == 1)

    {
        if (new_numbe > m)
        {
            new_element_para t;
            t.b01 = s4.b0[2 * m];
            t.b02 = s4.b0[2 * m + 1];
            t.j = joint[m];
            t.new_el = newelem;
            t.ratio = newelem.a / orig_elem.a;
            t.be1 = b_elm[m];
            new_elements.push_back(t);
        }
        else
        {
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
    }
    else
    {
        new_element_para t;
        t.b01 = s4.b0[2 * m];
        t.b02 = s4.b0[2 * m + 1];
        t.j = joint[m];
        t.new_el = newelem;
        t.ratio = newelem.a / orig_elem.a;
        t.be1 = b_elm[m];
        new_elements.push_back(t);
    }

}



std::vector<Point> computeIntersectionWithCircleSector(
    const Point& p1, const Point& p2, float xc, float yc, float R,
    float ang1_deg, float ang2_deg,int p1_status,int p2_status)
{
    std::vector<Point> intersections;
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float fx = p1.x - xc;
    float fy = p1.y - yc;

    float a = dx * dx + dy * dy;
    float b = 2 * (fx * dx + fy * dy);
    float c = fx * fx + fy * fy - R * R;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < -1e-6f) return intersections;  // Definitely no intersection
    if (discriminant < 0) discriminant = 0;   // Treat near-zero as tangent
    //if (discriminant < 0) return false;

    discriminant = std::sqrt(discriminant);
    float t1 = (-b - discriminant) / (2 * a);
    float t2 = (-b + discriminant) / (2 * a);

    for (float t : {t1, t2}) {
        if (t < 0 || t > 1) continue;
        Point ip{ p1.x + t * dx, p1.y + t * dy };
        int pos = pointInSector(ip, xc, yc, R, ang1_deg, ang2_deg);
        if (pos == 0 || pos == 2) {
        //if (pointInSector(ip, xc, yc, R, ang1_deg, ang2_deg)== 2) {
            intersections.push_back(ip);            
        }
    }
    if (intersections.empty() && (p1_status == 0 || p2_status == 0))
    {
        // Convert angles to Cartesian points on circle perimeter
        Point arc_pt1 = { xc + R * cos(ang1_deg), yc + R * sin(ang1_deg) };
        Point arc_pt2 = { xc + R * cos(ang2_deg), yc + R * sin(ang2_deg) };

        // Check intersection with radial lines
        Point radial_inter1, radial_inter2;
        bool intersects_radial1 = computeSegmentIntersection(p1, p2, xc, yc, arc_pt1.x, arc_pt1.y, radial_inter1);
        bool intersects_radial2 = computeSegmentIntersection(p1, p2, xc, yc, arc_pt2.x, arc_pt2.y, radial_inter2);

        if (intersects_radial1) intersections.push_back(radial_inter1);
        if (intersects_radial2) intersections.push_back(radial_inter2);
    }

    return intersections;
}



void findallintersect(const Point& p1, const Point& p2, float xc, float yc, float R,
    float ang1, float ang2,int& new_numbe, int m,int p1_stat, int p2_stat)
{
    // Step 1: Find all intersections
    Point ip1, ip2;   
    std::vector<Point> intersections = computeIntersectionWithCircleSector(p1, p2, xc, yc, R, ang1,ang2,p1_stat,p2_stat);

    if (intersections.empty()) {        
            auto& elem = elm_list[m];
            fix_tip_pointer(m, new_numbe);
            elm_list[new_numbe] = elem;
            b_elm[new_numbe] = b_elm[m];
            joint[new_numbe].aperture0 = joint[m].aperture0;
            joint[new_numbe].aperture_r = joint[m].aperture_r;
            s4.b0[2 * new_numbe] = s4.b0[2 * m];
            s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
            new_numbe++; 
            return;
    }
    else if (intersections.size() == 1)
    {
        Point outside_pt = (p1_stat == 1) ? p1 : p2;
        addClippedElement(outside_pt, intersections[0], m, new_numbe, 1);
    }
    else
    {
        struct IntersectionData {
            float t;
            Point pt;
        };
        std::vector<IntersectionData> sortedPoints;

        float dx = p2.x - p1.x;
        float dy = p2.y - p1.y;

        for (const Point& ip : intersections) {
            float t;
            if (std::abs(dx) >= std::abs(dy))
                t = (ip.x - p1.x) / dx;
            else
                t = (ip.y - p1.y) / dy;

            sortedPoints.push_back({ t, ip });
        }

        // Sort by t
        std::sort(sortedPoints.begin(), sortedPoints.end(),
            [](const IntersectionData& a, const IntersectionData& b) {
                return a.t < b.t;
            });
        
            addClippedElement(p1, sortedPoints[0].pt, m, new_numbe, 1);
            second_clipped = true;
            addClippedElement(sortedPoints[1].pt, p2, m, new_numbe, 2);
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
        else
        {
            findallintersect(p1, p2, xc, yc, R, ang1_deg, ang2_deg, new_numbe, m,p1_status,p2_status);

        }         
            
    }
    if (second_clipped = true)
    {
        for (const auto& e : new_elements)
        {
            elm_list[new_numbe] = e.new_el;
            b_elm[new_numbe] = e.be1;
            joint[new_numbe] = e.j;
            s4.b0[2 * new_numbe] = e.b01;
            s4.b0[2 * new_numbe + 1] = e.b02;
            b_elm[new_numbe].force1 *=  e.ratio ;
            b_elm[new_numbe].force2 *=  e.ratio;
            new_numbe++;
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