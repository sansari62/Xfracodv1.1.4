//#define NOMINMAX        // must go before windows.h
#include<stdafx.h>
#include<Rectangle_check.h>
#include "CommonPara.h"
#include <set>
#include <map>
#include <optional>
#include <algorithm>
#include<GeologicalForm.h>
#include<Remov_preexisting_BEs.h>


using namespace CommonPara_h::comvar;

std::vector<new_element_para>  new_elements1;





int point_inside_rectangle(const Point& point, const Rectangle1& rect) {
    float eps = 1e-6f;
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







// ---------- Main Clipping Logic ----------

void process_elements(const Rectangle1& rect) {

    int new_numbe = 0;
    for (int m = 0; m < numbe; ++m) {
        auto& elem = elm_list[m];
        // Reconstruct endpoints
        float dx = elem.a * elem.cosbet;
        float dy = elem.a * elem.sinbet;
        int tip_index = if_tip_element(m);

        Point p1{ elem.xm - dx, elem.ym - dy };
        Point p2{ elem.xm + dx, elem.ym + dy };

        int p1_state = point_inside_rectangle(p1, rect);
        int p2_state = point_inside_rectangle(p2, rect);

        if ((p1_state == 0 && p2_state == 0) ||
            (p1_state == 2 && p2_state == 2) ||
            (p1_state == 0 && p2_state == 2) ||
            (p1_state == 2 && p2_state == 0)) {
            if (tip_index != -1)
                tips[tip_index].ityp = 0;
            continue;
        }

        auto intersections = get_all_intersections(p1, p2, rect);
        
        if (intersections.empty()) {
            // Fully outside            
            if (tip_index != -1)
                tips[tip_index].mpointer = new_numbe;
            elm_list[new_numbe] = elem;
            b_elm[new_numbe] = b_elm[m];
            joint[new_numbe] = joint[m];
            update_s4(m, new_numbe);
            //s4.b0[2 * new_numbe] = s4.b0[2 * m];
            //s4.b0[2 * new_numbe + 1] = s4.b0[2 * m + 1];
            new_numbe++;
        }
        else if (intersections.size() == 1) {
            if (p1_state == 1 && p2_state == 1)
            {
                two_intersection_case(p1, p2, intersections[0], intersections[0], m, new_numbe, tip_index, new_elements1);
                
            }
            else
            {
                Point outside_pt = (p1_state == 1) ? p1 : p2;
                one_point_intersection(outside_pt, intersections[0], m, new_numbe, tip_index, new_elements1);
            }
        }

        //&&p1_state == 1 && p2_state == 1
        else if (intersections.size() == 2) {
            Point ip1, ip2;
            ip1 = intersections[0];
            ip2 = intersections[1];

            if ((ip1.x - p1.x) * (p2.x - p1.x) + (ip1.y - p1.y) * (p2.y - p1.y) >
                (ip2.x - p1.x) * (p2.x - p1.x) + (ip2.y - p1.y) * (p2.y - p1.y)) {
                std::swap(ip1, ip2);
            }
            two_intersection_case(p1, p2, ip1, ip2, m, new_numbe, tip_index, new_elements1);

        }
    }
    if (second_clipped = true)
    {
        for (const auto& e : new_elements1)
        {
            elm_list[new_numbe] = e.new_el;
            b_elm[new_numbe] = e.be1;
            joint[new_numbe] = e.j;
           // s4.b0[2 * new_numbe] = e.b01;
           // s4.b0[2 * new_numbe + 1] = e.b02;
            update_s4(e.s4_indx,new_numbe);
            b_elm[new_numbe].force1 *= e.ratio;
            b_elm[new_numbe].force2 *= e.ratio;
            if (e.tip_indx != -1)
                fix_tip_pointer1(new_numbe, e.tip_indx, 1);                
            new_numbe++;
        }
    }
   
    numbe = new_numbe;
    arrangetip();
    return ;
}



std::optional<Rectangle1> check_rectangle(bool flag)
{
    Rectangle1 rect1;
    if (flag)
    {
        if (validate_and_extract_rectangle(rect1))
        {
            rect_exca = true;
            process_elements(rect1);
        }
        return rect1;
    }
    else        
    {
        if (validate_and_extract_rectangle(rect1))
        {
            return rect1;     // return valid rectangle
        }
    }    
    return std::nullopt;
}