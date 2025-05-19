#include<stdafx.h>

#include <CommonPara.h>
#include<GeologicalForm.h>
#include <optional>
#include<Tip.h>

using namespace CommonPara_h::comvar;


void remove_elements_inside_borehole(
     float xc, float yc, float R)
{
    int new_numbe = 0;
    for (int i = 0; i < numbe; ++i) {
        BoundaryElement elm = comvar::elm_list[i];
        float dx = elm.xm - xc;
        float dy = elm.ym - yc;
        float dist2 = dx * dx + dy * dy;

        if (dist2 >= R * R) {
            // Keep the element
            elm_list[new_numbe++] = elm_list[i];
        }
        // else: skip the element (it is "removed")
    }
    numbe = new_numbe;  // Update number of active elements
    numbe_old = numbe;
}



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


struct Point {
    float x, y;
};

std::optional<Point> computeIntersectionWithCircle(
    const Point& p1, const Point& p2, float xc, float yc, float R)
{
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float fx = p1.x - xc;
    float fy = p1.y - yc;

    float a = dx * dx + dy * dy;
    float b = 2 * (fx * dx + fy * dy);
    float c = fx * fx + fy * fy - R * R;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return std::nullopt;

    discriminant = std::sqrt(discriminant);
    float t1 = (-b - discriminant) / (2 * a);
    float t2 = (-b + discriminant) / (2 * a);

    float t = -1;
    if (t1 >= 0 && t1 <= 1) t = t1;
    else if (t2 >= 0 && t2 <= 1) t = t2;

    if (t < 0) return std::nullopt;

    return Point{ p1.x + t * dx, p1.y + t * dy };
}

bool pointInCircle(const Point& p, float xc, float yc, float R) {
    float dx = p.x - xc, dy = p.y - yc;
    return dx * dx + dy * dy < R * R;
}


void clipBoundaryElements(    
    float xc, float yc, float R)
{
    int new_numbe = 0;
    for (int m = 0; m < numbe; ++m) {
        auto& elem = elm_list[m];
        // Reconstruct endpoints
        float dx = elem.a * elem.cosbet;
        float dy = elem.a * elem.sinbet;

        Point p1{ elem.xm - dx, elem.ym - dy };
        Point p2{ elem.xm + dx, elem.ym + dy };

        bool p1_in = pointInCircle(p1, xc, yc, R);
        bool p2_in = pointInCircle(p2, xc, yc, R);

        if (p1_in && p2_in) {
            // Entire element is inside → skip
            int merge = 0;
            specialLabel_200(merge,m);
            continue;
        }
        else if (!p1_in && !p2_in) {
            // Entirely outside → keep as-is
            fix_tip_pointer(m, new_numbe);
            elm_list[new_numbe++] = elem;
        }
        else {
            // Partially inside → clip and recompute center
            auto intersection = computeIntersectionWithCircle(p1, p2, xc, yc, R);
            if (!intersection) continue; // edge case: no intersection

            Point new_p = intersection.value();

            if (p1_in) p1 = new_p;
            else       p2 = new_p;

            // Recompute center and direction
            float cx = 0.5f * (p1.x + p2.x);
            float cy = 0.5f * (p1.y + p2.y);
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            float len = std::sqrt(dx * dx + dy * dy);
            if (len == 0) continue; // degenerate element
            //add new element
            BoundaryElement newelement;
            newelement.xm =cx;
            newelement.ym = cy;
            newelement.a = len/2;
            newelement.sinbet = dy / len;
            newelement.cosbet = dx / len;
            bool unique = isNewElementUnique(newelement);
            if (!unique)
            {
                new_numbe--;
                continue;
            }    
            newelement.kod= elm_list[m].kod;
            newelement.mat_no= elm_list[m].mat_no;
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
        //remove_elements_inside_borehole(xc, yc, r);
        clipBoundaryElements(xc, yc, r);
    }
}