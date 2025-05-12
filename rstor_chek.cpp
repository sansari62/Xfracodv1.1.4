//#include <stdafx.h>
#include <CommonPara.h>
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
}
