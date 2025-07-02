#pragma once
#ifndef rstor_chek_h
#define rstor_chek_h

extern bool second_clipped;
struct Point {
    float x, y;

    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }
};
void check_boreholes();
void addClippedElement(const Point& a, const Point& b,
    int m, int& new_numbe,int);

#endif // !rstor_chek_h