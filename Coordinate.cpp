#include "Coordinate.hpp"

Coordinate::Coordinate()
{
    x = 0;
    y = 0;
}

Coordinate::Coordinate(int x1, int y1)
{
    x = x1;
    y = y1;
}

Coordinate::Coordinate(const Coordinate &c)
{
    x = c.x;
    y = c.y;
}

//https://stackoverflow.com/questions/19871647/how-do-i-insert-objects-into-stl-set
//https://stackoverflow.com/questions/34047772/stdset-custom-comparator-for-2d-points
bool Coordinate::operator<(const Coordinate &right) const
{
    //    return x > right.x;
    if (x < right.x)
        return true;
    if (x > right.x)
        return false;
    if (y < right.y)
        return true;
    return false;
}