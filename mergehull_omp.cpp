#define TOP_LEVEL
#include <iostream>
#include <omp.h>
#include <cmath>
#include <vector>
#include <set>
#include <fstream>
#include "Coordinate.hpp"
#include "common.cpp"
using namespace std;

bool debug = false;

vector<Coordinate> joinHulls(vector<Coordinate> Hl, vector<Coordinate> Hr)
{

  // n1 -> number of points in polygon a
  // n2 -> number of points in polygon b
  int n1 = Hl.size(), n2 = Hr.size();

  int ia = 0, ib = 0;
  for (int i = 1; i < n1; i++)
    if (Hl[i].x > Hl[ia].x)
      ia = i;

  // ib -> leftmost point of b
  for (int i = 1; i < n2; i++)
    if (Hr[i].x < Hr[ib].x)
      ib = i;

  // finding the upper tangent
  int inda = ia, indb = ib;
  bool done = 0;
  while (!done)
  {
    done = 1;
    if (n1 > 1)
    {
      while (ccw(Hr[indb], Hl[inda], Hl[(inda + 1) % n1]) >= 0)
        inda = (inda + 1) % n1;
    }

    if (n2 > 1)
    {
      while (ccw(Hl[inda], Hr[indb], Hr[(n2 + indb - 1) % n2]) <= 0)
      {
        indb = (n2 + indb - 1) % n2;
        done = 0;
      }
    }
    cout << 1 << endl;
  }

  int uppera = inda, upperb = indb;
  inda = ia, indb = ib;
  done = 0;
  int g = 0;
  while (!done) //finding the lower tangent
  {
    done = 1;
    if (n2 > 1)
    {
      while (ccw(Hl[inda], Hr[indb], Hr[(indb + 1) % n2]) >= 0)
        indb = (indb + 1) % n2;
    }

    if (n1 > 1)
    {
      while (ccw(Hr[indb], Hl[inda], Hl[(n1 + inda - 1) % n1]) <= 0)
      {
        inda = (n1 + inda - 1) % n1;
        done = 0;
      }
    }
    cout << 2 << endl;
  }

  int lowera = inda, lowerb = indb;
  vector<Coordinate> ret;

  if (debug)
  {
    cout << "Hl : ";
    displayCoordinateVec(Hl);
    cout << "Hr : ";
    displayCoordinateVec(Hr);
  }

  //ret contains the convex hull after merging the two convex hulls
  //with the points sorted in anti-clockwise order
  int ind = uppera;
  ret.push_back(Hl[uppera]);
  while (ind != lowera)
  {
    ind = (ind + 1) % n1;
    ret.push_back(Hl[ind]);
  }

  ind = lowerb;
  ret.push_back(Hr[lowerb]);
  while (ind != upperb)
  {
    ind = (ind + 1) % n2;
    ret.push_back(Hr[ind]);
  }
  return ret;
}

vector<Coordinate> mergeHull(vector<Coordinate> P)
{
  if (debug)
    displayCoordinateVec(P);
  if (P.size() < 3)
    return P;

  vector<Coordinate> Hl;
  vector<Coordinate> Hr;
  // https://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-parallel.html
  // #pragma omp parallel num_threads(2)
  //   {
  //     int num = omp_get_thread_num();
  //     if (num == 0)
  // https://stackoverflow.com/questions/50549611/slicing-a-vector-in-c`
  Hl = mergeHull(vector<Coordinate>(P.begin(), P.begin() + P.size() / 2));
  // if (num == 1)
  Hr = mergeHull(vector<Coordinate>(P.begin() + P.size() / 2, P.end()));
  // }

  return joinHulls(Hl, Hr);
}

int main(int argc, char **argv)
{

  ifstream inp;

  if (argc == 3)
    debug = true;

  if (argc == 2)
  {
    inp.open(argv[1]);
  }
  else
  {
    inp.open("pre_sorted_x.txt");
  }
  if (!inp)
    return (EXIT_FAILURE);

  vector<Coordinate> v = parseFileForCoordinates(&inp);
  inp.close();

  const int ARRAY_SIZE = static_cast<int>(v.size());

  vector<Coordinate> hull = mergeHull(v);
  cout << "posthull" << endl;

  printf("Number of elements in hull: %ld\n", hull.size());

  storePolygon(hull, "polygon.txt");

  return 0;
}
