#include <vector>
#include <iostream>
#include <fstream>
#ifndef TOP_LEVEL
    #include "Coordinate.hpp"
#endif
using namespace std;

// Parse a file contatining a string of multiple 2D points
vector<Coordinate> parseFileForCoordinates(ifstream *inp)
{

    vector<Coordinate> v;

    string sNum;
    bool x_next = true;
    double x_temp = 0;
    Coordinate *c_temp;
    while (!inp->eof())
    {
        *inp >> sNum;
        for (int i = 0, len = sNum.size(); i < len; i++)
        {
            // check whether parsing character is punctuation or not
            if (sNum[i] == ' ' || sNum[i] == ',' || sNum[i] == '(' || sNum[i] == ')')
            {
                sNum.erase(i--, 1);
                len = sNum.size();
            }
        }
        if (x_next)
        {
            x_temp = stod(sNum);
        }
        else
        {
            c_temp = new Coordinate(x_temp, stod(sNum));
            v.push_back(*c_temp);
        }
        x_next = !x_next;
    }

    printf("Number of elements in file: %d\n", static_cast<int>(v.size()));
    return v;
}

// Generate output file
void storePolygon(vector<Coordinate> v, string filename)
{

    ofstream myfile;
    myfile.open(filename, ofstream::trunc);

    int i = 0;
    while (true)
    {
        myfile << '(' << v[i].x << ", " << v[i].y << ')';

        if (++i >= (int)(v.size()))
            break;

        myfile << " ";
    }
    myfile.close();
}

void displayCoordinate(Coordinate c)
{
    cout << '(' << c.x << ", " << c.y << ')';
}

// Generate output file
void displayCoordinateVec(vector<Coordinate> v)
{

    int i = 0;
    while (true)
    {
        displayCoordinate(v[i]);

        if (++i >= (int)(v.size()))
            break;

        cout << " -> ";
    }
    cout << endl;
}

int ccw(Coordinate a, Coordinate b, Coordinate c){
    float area = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
    if (area < 0 ) return -1;
    if (area > 0 ) return  1;
    return 0;
}
