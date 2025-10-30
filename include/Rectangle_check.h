#pragma once
#ifndef Rectangle_check_h
#define Rectangle_check_h
#include<rstor_chek.h>
#include <optional>

struct Rectangle1 {
    Point corners[4];
};

int point_inside_rectangle(const Point& point, const Rectangle1& rect);
std::optional<Rectangle1> check_rectangle(bool);

#endif // !Rectangle_check_h

