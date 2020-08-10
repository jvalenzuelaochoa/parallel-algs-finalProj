class Coordinate
{
public:
    double x;
    double y;

    Coordinate();
    Coordinate(int x1, int y1);
    Coordinate(const Coordinate &c);

    bool operator<(const Coordinate &right) const;
};