#define NOMINMAX        // must go before windows.h

#include<stdafx.h>
#include "CommonPara.h"
#include <set>
#include <map>
#include <optional>
#include <algorithm>
#include<Tip.h>
#include<GeologicalForm.h>
#include<rstor_chek.h>


using namespace CommonPara_h::comvar;


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



int point_inside_rectangle(const Point& point, const Rectangle1& rect, float eps = 1e-6f) {
    float min_x = std::min(rect.corners[0].x, rect.corners[2].x);
    float max_x = std::max(rect.corners[0].x, rect.corners[2].x);
    float min_y = std::min(rect.corners[0].y, rect.corners[2].y);
    float max_y = std::max(rect.corners[0].y, rect.corners[2].y);

    bool inside_x = (point.x > min_x + eps && point.x < max_x - eps);
    bool inside_y = (point.y > min_y + eps && point.y < max_y - eps);
    bool on_x = (point.x >= min_x - eps && point.x <= max_x + eps);
    bool on_y = (point.y >= min_y - eps && point.y <= max_y + eps);

    if (inside_x && inside_y)
        return 0;  // strictly inside
    if (on_x && on_y)
        return 2;  // on the border
    return 1;      // outside
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
        const float EPS = 1e-5f;
        return (std::min(a, b) - EPS <= c && c <= std::max(a, b) + EPS);
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
    if (nb < 3) return false;

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




void setting_newelement(BoundaryElement& newelement, int m, int new_numbe)
{
    newelement.kod = elm_list[m].kod;
    newelement.mat_no = elm_list[m].mat_no;
    elm_list[new_numbe] = newelement;
    b_elm[new_numbe].force1 = b_elm[m].force1 * elm_list[new_numbe].a / elm_list[m].a;
    b_elm[new_numbe].force2 = b_elm[m].force2 * elm_list[new_numbe].a / elm_list[m].a;
    b_elm[new_numbe].aks = b_elm[m].aks;
    b_elm[new_numbe].akn = b_elm[m].akn;
    b_elm[new_numbe].phi = b_elm[m].phi;
    b_elm[new_numbe].phid = b_elm[m].phid;

    b_elm[new_numbe].coh = b_elm[m].coh;
    joint[new_numbe].aperture0 = joint[m].aperture0;
    joint[new_numbe].aperture_r = joint[m].aperture_r;
    s4.b0[2 * new_numbe] = s4.b0[2 * m];
    s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
    fix_tip_pointer1(m, new_numbe);
}
// ---------- Main Clipping Logic ----------

void process_elements(const Rectangle1& rect) {

    int new_numbe = 0;
    for (int m = 0; m < numbe; ++m) {
        auto& elem = elm_list[m]; 
        Point p1, p2;
        p1.x = elem.xm - 0.5 * elem.a * elem.cosbet;
        p1.y = elem.ym - 0.5 * elem.a * elem.sinbet;
        p2.x = elem.xm + 0.5 * elem.a * elem.cosbet;
        p2.y = elem.ym + 0.5 * elem.a * elem.sinbet;

        int p1_state = point_inside_rectangle(p1, rect);
        int p2_state = point_inside_rectangle(p2, rect);

        if ((p1_state == 0 && p2_state == 0) ||
            (p1_state == 2 && p2_state == 2) ||
            (p1_state == 0 && p2_state == 2) ||
            (p1_state == 2 && p2_state == 0)) {
            int merge = 0;
            specialLabel_200(merge, m);            
            continue;
        }

        auto intersections = get_all_intersections(p1, p2, rect);

        if (intersections.empty()) {
             // Fully outside            
            fix_tip_pointer1(m, new_numbe);
            elm_list[new_numbe] = elem;
            b_elm[new_numbe] = b_elm[m];
            //probablu true
            joint[new_numbe].aperture0 = joint[m].aperture0;
            joint[new_numbe].aperture_r = joint[m].aperture_r;
            s4.b0[2 * new_numbe] = s4.b0[2 * m];
            s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
            new_numbe++;
        }
        else if (intersections.size() == 1) {
            if (p1_state == 1 && p2_state == 1)
            {
                cout << "new condit";
            }
            Point clipped = intersections[0];
            Point outside_pt = (p1_state == 1) ? p1 : p2;

            float dx = outside_pt.x - clipped.x;
            float dy = outside_pt.y - clipped.y;
            float len = std::sqrt(dx * dx + dy * dy);
            if (len < 1e-6f) {
                int merge = 0;
                specialLabel_200(merge, m);
                continue;
            }  // skip degenerate or tiny element
            float new_xm = 0.5 * (outside_pt.x + clipped.x);
            float new_ym = 0.5 * (outside_pt.y + clipped.y);
            
            BoundaryElement newelement;
            newelement.xm = new_xm;
            newelement.ym = new_ym;
            newelement.a = len / 2;
            newelement.sinbet = dy / len;
            newelement.cosbet = dx / len;
            bool unique = isNewElementUnique(newelement);
            if (!unique)
            {
                //new_numbe--;
                continue;
            }
            setting_newelement(newelement, m, new_numbe);
            new_numbe++;
        }
        else if (intersections.size() == 2 && p1_state == 1 && p2_state == 1) {
             Point ip1, ip2;
             ip1 = intersections[0];
             ip2 = intersections[1];

            if ((ip1.x - p1.x) * (p2.x - p1.x) + (ip1.y - p1.y) * (p2.y - p1.y) >
                (ip2.x - p1.x) * (p2.x - p1.x) + (ip2.y - p1.y) * (p2.y - p1.y)) {
                std::swap(ip1, ip2);
            }
            addClippedElement(p1, ip1, m, new_numbe);
            // Create [ip2, p2]
            addClippedElement(ip2, p2, m, new_numbe);






            //float d1 = std::hypot(p1.x - clip1.x, p1.y - clip1.y);
            //float d2 = std::hypot(p2.x - clip2.x, p2.y - clip2.y);

            //// Determine order based on proximity
            //Point out1_start = (d1 < d2) ? p1 : p2;
            //Point out1_end = (d1 < d2) ? clip1 : clip2;

            //Point out2_start = (d1 < d2) ? clip2 : clip1;
            //Point out2_end = (d1 < d2) ? p2 : p1;

            //auto try_add_segment = [&](const Point& start, const Point& end) {
            //    float dx = end.x - start.x;
            //    float dy = end.y - start.y;
            //    float len = std::sqrt(dx * dx + dy * dy);
            //    if (len < 1e-6f) {
            //        // Entire element is inside → skip
            //        int merge = 0;
            //        specialLabel_200(merge, m);
            //        return;
            //    }

            //    float new_xm = 0.5 * (start.x + end.x);
            //    float new_ym = 0.5 * (start.y + end.y);

            //    BoundaryElement newelement;
            //    newelement.xm = new_xm;
            //    newelement.ym = new_ym;
            //    newelement.a = len / 2;
            //    newelement.sinbet = dy / len;
            //    newelement.cosbet = dx / len;

            //    bool unique = isNewElementUnique(newelement);
            //    if (!unique) return;
            //    setting_newelement(newelement, m, new_numbe);

            //    /*newelement.kod = elem.kod;
            //    newelement.mat_no = elem.mat_no;

            //    elm_list[new_numbe] = newelement;

            //    float scale = newelement.a / elem.a;
            //    b_elm[new_numbe].force1 = b_elm[m].force1 * scale;
            //    b_elm[new_numbe].force2 = b_elm[m].force2 * scale;
            //    b_elm[new_numbe].aks = b_elm[m].aks;
            //    b_elm[new_numbe].akn = b_elm[m].akn;
            //    b_elm[new_numbe].phi = b_elm[m].phi;
            //    b_elm[new_numbe].phid = b_elm[m].phid;
            //    b_elm[new_numbe].coh = b_elm[m].coh;

            //    joint[new_numbe].aperture0 = joint[m].aperture0;
            //    joint[new_numbe].aperture_r = joint[m].aperture_r;

            //    s4.b0[2 * new_numbe] = s4.b0[2 * m];
            //    s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];*/

            //    //fix_tip_pointer1(m, new_numbe);
            //    new_numbe++;
            //    };

            //// Try adding two outside segments
            //try_add_segment(out1_start, out1_end);
            //try_add_segment(out2_start, out2_end);
        }
        else
        {
            cout << "new condition needed";
        }

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
   /*else
   {
       MessageBox(nullptr, L"Check edge definition in input!", L"Error!", MB_OK);
       exit(0);
   }*/

}