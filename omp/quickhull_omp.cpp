#define TOP_LEVEL
#include <iostream>
#include <omp.h>
#include <cmath>
#include <vector>
#include <set>
#include <fstream>
#include "../common/Coordinate.hpp"
#include "../common/common.cpp"
#include <time.h>
using namespace std;

bool debug = false;

//reference https://stackoverflow.com/questions/28258590/using-openmp-to-get-the-index-of-minimum-element-parallelly
typedef std::pair<unsigned int, double> IndexValuePair;

IndexValuePair myMin(IndexValuePair a, IndexValuePair b)
{
  return a.second < b.second ? a : b;
}

IndexValuePair myMax(IndexValuePair a, IndexValuePair b)
{
  return a.second > b.second ? a : b;
}

// function from https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
// Returns the side of point p with respect to line
// joining points p1 and p2.
int findSide(Coordinate p1, Coordinate p2, Coordinate p)
{
  double val = (p.y - p1.y) * (p2.x - p1.x) -
               (p2.y - p1.y) * (p.x - p1.x);

  if (val > 0)
    return 1; //left
  if (val < 0)
    return -1; //right
  return 0;
}

// function from https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
// returns a value proportional to the distance
// between the point p and the line joining the
// points p1 and p2
double lineDist(Coordinate p1, Coordinate p2, Coordinate p)
{
  return abs((p.y - p1.y) * (p2.x - p1.x) -
             (p2.y - p1.y) * (p.x - p1.x));
}

// function from https://www.geeksforgeeks.org/perpendicular-distance-between-a-point-and-a-line-in-2-d/
// Function to find distance
double shortest_distance(Coordinate p1, Coordinate p2, Coordinate p)
{
  double x1, y1, a, b, c;

  x1 = p.x;
  y1 = p.y;

  a = p1.y - p2.y;
  b = p2.x - p1.x;
  c = (p1.x - p2.x) * p1.y + (p2.y - p1.y) * p1.x;

  double d = fabs((a * x1 + b * y1 + c)) /
             (sqrt(a * a + b * b));
  //printf("Perpendicular distance is %f\n", d);
  return d;
}

// function from https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
// End points of line L are p1 and p2.  side can have value
// 1 or -1 specifying each of the parts made by the line L
vector<Coordinate> subHull(vector<Coordinate> P, Coordinate p1, Coordinate p2)
{
  int ind = -1;
  double max_dist = 0;
  int numthreads = 2;
  int totalnumthreads, tid;

  // array compaction to get Pprime
  int i, n, chunk, d;
  n = P.size();
  chunk = 1;
  int c[n], p[n];

#pragma omp parallel private(i) num_threads(n)
  {

#pragma omp for schedule(static, chunk)
    for (i = 0; i < n; i++)
    {
      if (findSide(p1, p2, P[i]) > 0)
      {
        c[i] = 1;
      }
      else
      {
        c[i] = 0;
      }
    }
  }

  if (debug)
  {
    for (i = 0; i < n; i++)
    {
      printf("C array= %d\n", c[i]);
    }
    printf("n is= %d\n", n);
  }

#pragma omp parallel private(i) num_threads(n)
  {
#pragma omp for schedule(static, 1)
    for (i = 0; i < n; i++)
    {
      p[i] = c[i];
    }
  }

  int val;

  for (d = 1; d < n; d = 2 * d)
  {
#pragma omp parallel private(i, val) num_threads(n)
    {
#pragma omp for schedule(static, chunk)
      for (i = 1; i < n; i++)
      {
        if (i >= d)
          val = p[i - d];
      }

#pragma omp barrier
#pragma omp for schedule(static, chunk)

      for (i = 1; i < n; i++)
      {
        if (i >= d)
          p[i] = p[i] + val;
      }
    }
  }

  if (debug)
  {
    for (i = 0; i < n; i++)
    {
      printf("Prefix sum of c= %d\n", p[i]);
    }
  }

  int sum_c = 0; //sum of p

#pragma omp parallel for reduction(+ \
                                   : sum_c)

  for (i = 0; i < n; i++)
    sum_c = sum_c + c[i];

  if (debug)
    printf("sum_c = %d\n", sum_c);
  vector<Coordinate> Pprime(sum_c);

#pragma omp parallel private(i) num_threads(n)
  {
#pragma omp for schedule(static, chunk)
    for (i = 0; i < n; i++)
    {
      if (c[i] == 1)
      {
        Pprime[p[i] - 1] = P[i];
      }
    }
  }

  if (Pprime.size() < 2)
  {
    Pprime.push_back(p1);
    return Pprime;
  }
  else
  {
    for (int i = 0; i < Pprime.size(); i++)
    {
      double temp = shortest_distance(p1, p2, Pprime[i]);
      if (temp > max_dist)
      {
        ind = i;
        max_dist = temp;
      }
    }
    vector<Coordinate> Hl;
    vector<Coordinate> Hr;

#pragma omp parallel private(tid)
    {
      // Recur for the two parts divided by a[ind]
      /* Obtain and print thread id */
      tid = omp_get_thread_num();

      Hl = subHull(Pprime, p1, Pprime[ind]);
      Hr = subHull(Pprime, Pprime[ind], p2);

      /* Only master thread does this */
      if (tid == 0)
      {
        totalnumthreads = omp_get_num_threads();
        //printf("Number of threads = %d\n", totalnumthreads);
      }
    }

    //ref: https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
    for (int i = 0; i < Hr.size(); i++)
    {
      Hl.push_back(Hr[i]);
    }

    return Hl;
  }
}

vector<Coordinate> quickHull(vector<Coordinate> P)
{

  // Finding the point with minimum and
  // maximum x-coordinate
  int min_x = 0, max_x = 0;
  for (int i = 1; i < P.size(); i++)
  {
    if (P[i].x < P[min_x].x)
      min_x = i;
    if (P[i].x > P[max_x].x)
      max_x = i;
  }

  if (debug)
  {
    cout << "The seq min_x is: " << min_x << ", the max_x is " << max_x << endl;
  }

  //use reduction to find min and max

  int i, n, chunk, s;
  n = P.size();
  chunk = 1;
  double arrmin[n], arrmax[n];
  int xmin, xmax;

//omp block begins
#pragma omp parallel private(i) num_threads(n)
  {
#pragma omp for schedule(static, chunk)

    for (i = 0; i < n; i++)
    {
      arrmin[i] = P[i].x;
      arrmax[i] = P[i].x;
    }
  }
  //omp block ends

  IndexValuePair minValueIndex(0, INT_MAX);

#pragma omp declare reduction(minPair:IndexValuePair              \
                              : omp_out = myMin(omp_out, omp_in)) \
    initializer(omp_priv = IndexValuePair(0, INT_MAX))

#pragma omp parallel for reduction(minPair \
                                   : minValueIndex)
  for (i = 0; i < n; i++)
  {
    if (arrmin[i] < minValueIndex.second)
    {
      minValueIndex.first = i;
      minValueIndex.second = arrmin[i];
    }
  }

  if (debug)
  {
    cout << "minimum value = " << minValueIndex.second << endl;
    cout << "min index = " << minValueIndex.first << endl;
  }

  IndexValuePair maxValueIndex(0, INT_MIN);

#pragma omp declare reduction(maxPair:IndexValuePair              \
                              : omp_out = myMax(omp_out, omp_in)) \
    initializer(omp_priv = IndexValuePair(0, INT_MIN))

#pragma omp parallel for reduction(maxPair \
                                   : maxValueIndex)
  for (i = 0; i < n; i++)
  {
    if (arrmax[i] > maxValueIndex.second)
    {
      maxValueIndex.first = i;
      maxValueIndex.second = arrmax[i];
    }
  }

  if (debug)
  {
    cout << "max value = " << maxValueIndex.second << endl;
    cout << "max index = " << maxValueIndex.first << endl;
  }

  min_x = minValueIndex.first;
  max_x = maxValueIndex.first;

  // Recur for the two parts divided by a[ind]
  vector<Coordinate> Hl = subHull(P, P[min_x], P[max_x]);
  vector<Coordinate> Hr = subHull(P, P[max_x], P[min_x]);
  //ref: https://stackoverflow.com/questions/201718/concatenating-two-stdvectors

  for (int i = 0; i < Hr.size(); i++)
  {
    Hl.push_back(Hr[i]);
  }
  return Hl;
}

int main(int argc, char **argv)
{

  ifstream inp;

  if (argc == 3)
    debug = true;

  if (argc >= 2)
  {
    inp.open(argv[1]);
  }
  else
  {
    inp.open("pre_sorted.txt");
  }
  if (!inp)
    return (EXIT_FAILURE);

  vector<Coordinate> v = parseFileForCoordinates(&inp);
  inp.close();

  const int ARRAY_SIZE = static_cast<int>(v.size());

  clock_t start, end;
  start = clock();

  vector<Coordinate> hull = quickHull(v);

  end = clock();
  // Calculating total time taken by the program.
  double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
  printf("Time taken by program is :%f\n", time_taken);

  printf("Number of elements in hull: %ld\n", hull.size());

  storePolygon(hull,"polygon.txt");

  return 0;
}