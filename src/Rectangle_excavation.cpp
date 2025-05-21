#define NOMINMAX        // must go before windows.h

#include<stdafx.h>

#include "CommonPara.h"
#include <set>
#include <map>
#include <optional>
#include <algorithm>
#include<Tip.h>
#include<GeologicalForm.h>


using namespace CommonPara_h::comvar;


struct Point {
    float x, y;

    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }
};

//struct Edge {
//    Point p1, p2;
//};


bool are_floats_equal(float a, float b, float tol = 1e-6f) {
    return std::fabs(a - b) < tol;
}

struct Rectangle1 {
    Point corners[4];
};


void fix_tip_pointer1(int m, int new_numbe)
{
    for (int k = 0; k < no; ++k)
    {
        Tip& t = tips[k];
        if (t.mpointer == m)
            t.mpointer = new_numbe;
    }
    return;
}


bool point_inside_rectangle(const Point& p, const Rectangle1& rect) {
    float min_x = std::min(
        std::min(rect.corners[0].x, rect.corners[1].x),
        std::min(rect.corners[2].x, rect.corners[3].x)
    );

    float max_x = std::max(
        std::max(rect.corners[0].x, rect.corners[1].x),
        std::max(rect.corners[2].x, rect.corners[3].x)
    );

    float min_y = std::min(
        std::min(rect.corners[0].y, rect.corners[1].y),
        std::min(rect.corners[2].y, rect.corners[3].y)
    );

    float max_y = std::max(
        std::max(rect.corners[0].y, rect.corners[1].y),
        std::max(rect.corners[2].y, rect.corners[3].y)
    );
    return (p.x >= min_x && p.x <= max_x && p.y >= min_y && p.y <= max_y);
}



std::optional<Point> get_intersection(Point A, Point B, Point C, Point D) {
    
    float a1 = B.y - A.y;
    float b1 = A.x - B.x;
    float c1 = a1 * A.x + b1 * A.y;

    float a2 = D.y - C.y;
    float b2 = C.x - D.x;
    float c2 = a2 * C.x + b2 * C.y;

    float det = a1 * b2 - a2 * b1;
    if (std::fabs(det) < 1e-6f) return std::nullopt;

    float x = (b2 * c1 - b1 * c2) / det;
    float y = (a1 * c2 - a2 * c1) / det;

    auto between = [](float a, float b, float c) {
        return (std::min(a, b) - 1e-6f <= c && c <= std::max(a, b) + 1e-6f);
        };

    if (between(A.x, B.x, x) && between(A.y, B.y, y) &&
        between(C.x, D.x, x) && between(C.y, D.y, y)) {
        return Point{ x, y };
    }
    return std::nullopt;
}

std::vector<Point> get_all_intersections(const Point& p1, const Point& p2, const Rectangle1& rect) {
    std::vector<Point> intersections;
    for (int i = 0; i < 4; ++i) {
        Point r1 = rect.corners[i];
        Point r2 = rect.corners[(i + 1) % 4];
        auto ipt = get_intersection(p1, p2, r1, r2);
        if (ipt.has_value()) {
            intersections.push_back(ipt.value());
        }
    }
    return intersections;
}

// ---------- Rectangle Validation ----------

bool validate_and_extract_rectangle( Rectangle1& out_rect) {
    if (nb != 4) return false;

    std::set<Point> unique_points;
    std::map<Point, int> point_count;

    for (size_t i = 0; i < nb; ++i) {
        Edge& edge = bund_list[i]; 
        Point p1,p2;
        p1.x = edge.get_xbeg();
        p1.y = edge.get_ybeg();
        unique_points.insert(p1);
        p2.x = edge.get_xend();
        p2.y = edge.get_yend();
        unique_points.insert(p2);
        point_count[p1]++;
        point_count[p2]++;
    }

    if (unique_points.size() != 4) return false;

    for (const auto& [pt, count] : point_count) {
        if (count != 2) return false;
    }

    // Sort corners clockwise
    float cx = 0, cy = 0;
    std::vector<Point> corners;
    for (const auto& pt : unique_points) {
        corners.push_back(pt);
        cx += pt.x;
        cy += pt.y;
    }
    cx /= 4.0f; cy /= 4.0f;

    std::sort(corners.begin(), corners.end(), [cx, cy](const Point& a, const Point& b) {
        return std::atan2(a.y - cy, a.x - cx) < std::atan2(b.y - cy, b.x - cx);
        });

    for (int i = 0; i < 4; ++i)
        out_rect.corners[i] = corners[i];

    return true;
}

// ---------- Main Clipping Logic ----------

void process_elements(const Rectangle1& rect) {

    std::vector<BoundaryElement> updated;

    int new_numbe = 0;
    for (int m = 0; m < numbe; ++m) {
        auto& elem = elm_list[m]; 
        Point p1, p2;
        p1.x = elem.xm - 0.5f * elem.a * elem.cosbet;
        p1.y = elem.ym - 0.5f * elem.a * elem.sinbet;
        p2.x = elem.xm + 0.5f * elem.a * elem.cosbet;
        p2.y = elem.ym + 0.5f * elem.a * elem.sinbet;

        bool p1_inside = point_inside_rectangle(p1, rect);
        bool p2_inside = point_inside_rectangle(p2, rect);

        if (p1_inside && p2_inside) {
            // Entire element is inside → skip
            int merge = 0;
            specialLabel_200(merge, m);
            continue;
        } 

        auto intersections = get_all_intersections(p1, p2, rect);

        if (intersections.empty()) {
            //updated.push_back(elem); // Fully outside
            fix_tip_pointer1(m, new_numbe);
            elm_list[new_numbe] = elem;
            b_elm[new_numbe] = b_elm[m];
            new_numbe++;
        }
        else if (intersections.size() == 1) {
            Point clipped = intersections[0];
            Point outside_pt = p1_inside ? p2 : p1;

            float dx = outside_pt.x - clipped.x;
            float dy = outside_pt.y - clipped.y;
            float len = std::sqrt(dx * dx + dy * dy);
            if (len == 0) continue; // degenerate element

            float new_xm = 0.5f * (outside_pt.x + clipped.x);
            float new_ym = 0.5f * (outside_pt.y + clipped.y);
            
           // updated.push_back({ new_xm, new_ym, new_a, new_sinbet, new_cosbet });
            BoundaryElement newelement;
            newelement.xm = new_xm;
            newelement.ym = new_ym;
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
            new_numbe++;
        }
        //else {
        //    // In rare case of two intersections: skip or clip both ends
        //    // Optional: could be split into two new elements if needed
        //}
    }
    numbe = new_numbe;
    arrangetip();
    return ;
}



void check_rectangle()
{
    Rectangle1 rect1;
   if( validate_and_extract_rectangle(rect1))
        process_elements(rect1);
   else
       MessageBox(nullptr, L"Check edge definition in input!", L"Error!", MB_OK);

}