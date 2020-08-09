#include <iostream>
#include <omp.h>
#include<cmath>
#include <vector>
#include<set>
#include <fstream>

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>


#include <thrust/count.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/partition.h>
#include <thrust/reverse.h>
#include <thrust/pair.h>

using namespace std;
#define BLOCK_SIZE 1024
#define CudaChainPoint  float2


// Make predicate easy to write
typedef thrust::tuple<float, float, int> FloatTuple3;

// Predicate
struct is_interior_tuple {
  __host__ __device__

  bool operator()(const FloatTuple3 &p) {
    return thrust::get<2>(p) > 0;
  }
};

/////////////////////////////////////////////////////////
// reference: https://github.com/sina-masnadi/CudaChain/blob/master/src/CudaChain.cu
// isLeft():tests if a point is Left|On|Right of an infinite line.
// Input: three points P0, P1, and P2
// Return:
//	>0 for P2 left of the line through P0 and P1
// =0 for P2 on the line
// <0 for P2 right of the line
// See: Algorithm 1 on Area of Triangles
__host__ __device__

int isLeft(CudaChainPoint *P0, CudaChainPoint *P1, float &P2_x, float &P2_y) {
  ///      float val = (p.y - p1.y) * (p2.x - p1.x) - 
  //      (p2.y - p1.y) * (p.x - p1.x); 

  float val =  (P1->x - P0->x) * (P2_y - P0->y) - (P2_x - P0->x) * (P1->y - P0->y);
  if (val >0) return 1;
  if ( fabs(P2_x - P0->x) < 0.0001 && fabs(P2_y - P0->y) <0.0001) return 1;
  return 0;

}

// function from https://www.geeksforgeeks.org/perpendicular-distance-between-a-point-and-a-line-in-2-d/
// Function to find distance
__host__ __device__

float shortest_distance(CudaChainPoint *P0, CudaChainPoint *P1, float &P2_x, float &P2_y)
{
  // P0 : p1
  // P1 : p2
  // P2_x : p
  float x1, y1, a, b, c;


  x1 = P2_x;
  y1 = P2_y;
  
  a = P0->y-P1->y;
  b = P1->x-P0->x;
  c = (P0->x-P1->x)*P0->y + (P1->y-P0->y)*P0->x;

    
  float d = fabs((a * x1 + b * y1 + c)) /  
    (sqrt(a * a + b * b)); 
  //printf("Perpendicular distance is %f\n", d); 
  return d; 
}

// Preprocess by discarding interior points, i.e., 1st round of discarding
__global__
void kernelIsLeft(float *h_extreme_x, float *h_extreme_y,
		  float *v_x, float *v_y, int *flag, int n) {
  __shared__ float2 s_extreme[2]; // Stored in shared memory

  if (threadIdx.x == 0) {
    for (int t = 0; t < 2; t++) {
      s_extreme[t].x = h_extreme_x[t];
      s_extreme[t].y = h_extreme_y[t];
    }
  }
  __syncthreads();

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n) // Check
    {
      flag[i] = isLeft(&s_extreme[0], &s_extreme[1], v_x[i], v_y[i]);
    }
}

// Preprocess by discarding interior points, i.e., 1st round of discarding
__global__
void kernelCalcDist(float *h_extreme_x, float *h_extreme_y,
		    float *v_x, float *v_y, float *dist, int *flag, int n) {
  __shared__ float2 s_extreme[2]; // Stored in shared memory

  if (threadIdx.x == 0) {
    for (int t = 0; t < 2; t++) {
      s_extreme[t].x = h_extreme_x[t];
      s_extreme[t].y = h_extreme_y[t];
    }
  }
  __syncthreads();

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n && flag[i] ==1) // Check
    {
      dist[i] = shortest_distance(&s_extreme[0], &s_extreme[1], v_x[i], v_y[i]);
    }
}

class Coordinate {
public:
  float x;
  float y;

  Coordinate()
  {
    x = 0;
    y = 0;
  }

  Coordinate(float x1, float y1)
  {
    x = x1;
    y = y1;
  }

  Coordinate(const Coordinate &c)
  {
    x = c.x;
    y = c.y;
  }
  
  //https://stackoverflow.com/questions/19871647/how-do-i-insert-objects-into-stl-set
  //https://stackoverflow.com/questions/34047772/stdset-custom-comparator-for-2d-points
  bool operator< (const Coordinate &right) const
  {
    //    return x > right.x;
    if (x < right.x) return true;
    if (x > right.x) return false;
    if (y < right.y) return true;
    return false;
  }

};

int ccw(Coordinate a, Coordinate b, Coordinate c){
  float area = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
  if (area < 0 ) return -1;
  if (area > 0 ) return  1;
  return 0;

}
// function from https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
// Returns the side of point p with respect to line 
// joining points p1 and p2. 
int findSide(Coordinate p1, Coordinate p2, Coordinate p) 
{ 
  float val = (p.y - p1.y) * (p2.x - p1.x) - 
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
float lineDist(Coordinate p1, Coordinate p2, Coordinate p) 
{ 
  return abs ((p.y - p1.y) * (p2.x - p1.x) - 
	      (p2.y - p1.y) * (p.x - p1.x)); 

}



// function from https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
// End points of line L are p1 and p2.  side can have value 
// 1 or -1 specifying each of the parts made by the line L 
/*vector<Coordinate> subHull(vector<Coordinate> P, Coordinate p1, Coordinate p2) 
  { 
  int ind = -1; 
  float max_dist = 0; 

  vector<Coordinate> Pprime;

  for (int i=0; i< P.size() ; i++) 
  { 
  if (findSide(p1, p2, P[i]) >0)
  { 
  Pprime.push_back(P[i]); 
  } 
  }
  if(Pprime.size() <2)
  {
  Pprime.push_back(p1);
  return Pprime;
  }
  else
  {
  for (int i=0; i< Pprime.size() ; i++) 
  {
  float temp = shortest_distance(p1, p2, Pprime[i]); 
  if ( temp > max_dist)
  { 
  ind = i;
  max_dist = temp;
  } 
  }
  // Recur for the two parts divided by a[ind]
  vector<Coordinate> Hl = subHull(Pprime, p1, Pprime[ind]);
  vector<Coordinate> Hr = subHull(Pprime, Pprime[ind], p2);
  //ref: https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
  for (int i=0; i< Hr.size() ; i++) 
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
  for (int i=1; i< P.size(); i++) 
  { 
  if (P[i].x < P[min_x].x) 
  min_x = i; 
  if (P[i].x > P[max_x].x) 
  max_x = i; 
  } 
  cout << "The min_x is: " <<min_x<<", the max_x is "<<max_x << endl;


  int data[6] = {1, 0, 2, 2, 1, 3};
  thrust::pair<int *, int *> result = thrust::minmax_element(thrust::host, data, data + 6);
  cout << "The thrust result.first is: " <<result.first<<endl;
  cout << "The thrust result.second is: " <<result.second<<endl;
  cout << "The thrust result.second - data is : " <<result.second - data<<endl;
  // Recur for the two parts divided by a[ind]
  vector<Coordinate> Hl = subHull(P, P[min_x], P[max_x]);
  vector<Coordinate> Hr = subHull(P, P[max_x], P[min_x]);
  //ref: https://stackoverflow.com/questions/201718/concatenating-two-stdvectors

  for (int i=0; i< Hr.size() ; i++) 
  {
  Hl.push_back(Hr[i]);
  }
  return Hl;

  
  } */

int main(int argc, char **argv)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    fprintf(stderr, "error: no devices supporting CUDA.\n");
    exit(EXIT_FAILURE);
  }
  int dev = 0;
  cudaSetDevice(dev);

  cudaDeviceProp devProps;
  if (cudaGetDeviceProperties(&devProps, dev) == 0)
    {
      printf("Number of multiprocessors:     %d\n", devProps.multiProcessorCount);
      printf("Number of threads per Block:   %d\n", devProps.maxThreadsPerBlock);
    }

  ifstream inp;

  if (argc == 2)
    {
      inp.open(argv[1]);
    }
  else
    {
      inp.open("pre_sorted.txt");
    }
  if (!inp)
    return(EXIT_FAILURE);

  vector<Coordinate> v;
  string sNum;
  bool x_next = true;
  float x_temp = 0;
  Coordinate *c_temp;
  while(!inp.eof())
    {
      inp >> sNum;
      for (int i = 0, len = sNum.size(); i < len; i++)
        {
	  // check whether parsing character is punctuation or not
	  if (sNum[i]== ' ' || sNum[i] == ',' || sNum[i] == '(' || sNum[i] == ')')
            {
	      sNum.erase(i--, 1);
	      len = sNum.size();
            }
        }
      if(x_next)
        {
	  x_temp = stof(sNum);
        }
      else
        {
	  c_temp = new Coordinate(x_temp, stof(sNum));
	  v.push_back(*c_temp);
        }
      x_next = !x_next;
    }

  printf("Number of elements in file: %d\n", static_cast<int>(v.size()));

  const int n = static_cast<int>(v.size());

  thrust::host_vector<float> hx(n);
  thrust::host_vector<float> hy(n);
    
  for (int i=0; i< n; i++) 
    {
      hx[i] = v[i].x;
      hy[i] = v[i].y;
     
    }

  thrust::device_vector<float> dx = hx;
  thrust::device_vector<float> dy = hy;

    
  thrust::device_vector<float> dist(n);
  thrust::device_vector<int> head(n);
  thrust::device_vector<int> keys(n);
  thrust::device_vector<int> first_pts(n);
  thrust::device_vector<int> flag(n);

     
  int min_x = 0, max_x = 0; 
  for (int i=1; i< n; i++) 
    { 
      if (hx[i] < hx[min_x]) 
	min_x = i; 
      if (hx[i] > hx[max_x]) 
	max_x = i; 
    } 
  cout << "The min_x is: " <<min_x<<", the max_x is "<<max_x << endl;
  typedef thrust::device_vector<float>::iterator FloatIter;

  // Find the four extreme points, i.e., min x, max x, min y, and max y
  thrust::pair<FloatIter, FloatIter> extremex = thrust::minmax_element(dx.begin(), dx.end());
    

  // One method to find min / max (Correct)
  thrust::device_vector<float>::iterator minx = extremex.first;
  thrust::device_vector<float>::iterator maxx = extremex.second;
    

  //typedef thrust::device_vector<float>::iterator floatIter;
  //floatIter minx = thrust::min_element(dx.begin(),dx.end());
  //floatIter maxx = thrust::max_element(dx.begin(),dx.end());
    
  //thrust::pair<float *, float *> result = thrust::minmax_element(thrust::host, dx, dx+n);
  cout << "minx = : " <<*minx<<endl;
  cout << "maxx = : " <<*maxx<<endl;
  cout << "minx index = : " <<minx - dx.begin()<<endl;
  cout << "maxx index = : " <<maxx - dx.begin()<<endl;

  // Store the four extreme points temporarily
  thrust::device_vector<float> d_extreme_x(2);
  thrust::device_vector<float> d_extreme_y(2);

  thrust::device_vector<float> dpstack(0);
  thrust::device_vector<float> dhull(0);
  thrust::device_vector<int> ps(0); //psize
  


  // max-min: p0x p0y p1x p1y
  dpstack.push_back(*maxx);
  dpstack.push_back(dy[maxx - dx.begin()]);
  dpstack.push_back(*minx);
  dpstack.push_back(dy[minx - dx.begin()]);
  ps.push_back(n);
  // min - max: p0x p0y p1x p1y
  dpstack.push_back(*minx);
  dpstack.push_back(dy[minx - dx.begin()]);
  dpstack.push_back(*maxx);
  dpstack.push_back(dy[maxx - dx.begin()]);
  ps.push_back(n);

  


  // Get the pointers to the arrays, to be used as launch arguments
  float *d_extreme_x_ptr = thrust::raw_pointer_cast(&d_extreme_x[0]);
  float *d_extreme_y_ptr = thrust::raw_pointer_cast(&d_extreme_y[0]);
  float *d_x_ptr = thrust::raw_pointer_cast(&dx[0]);
  float *d_y_ptr = thrust::raw_pointer_cast(&dy[0]);
  int *flag_ptr = thrust::raw_pointer_cast(&flag[0]);
  float * dist_ptr = thrust::raw_pointer_cast(&dist[0]);
  // Defining a zip_iterator type can be a little cumbersome ...
  typedef thrust::device_vector<float>::iterator FloatIterator;
  typedef thrust::device_vector<int>::iterator IntIterator;
  typedef thrust::tuple<FloatIterator, FloatIterator, IntIterator> FloatIteratorTuple;
  typedef thrust::zip_iterator<FloatIteratorTuple> Float3Iterator;

  bool debug = true;
  int di = 0;

  clock_t start, end;
  start = clock(); 

  while(dpstack.size() >0)
    //while(di<2)
    {
      if(debug){
	di++;
	std::cout << "/////DEBUG CYCLE = " << di << std::endl;

	std::cout << "ps.back() = " << ps.back()<< std::endl;}
      
      //p1 p2   // min - max: p0x p0y p1x p1y
      // p1y p1x p0y p0x


      d_extreme_y[1] = dpstack.back();
      dpstack.pop_back();
      d_extreme_x[1] = dpstack.back();
      dpstack.pop_back();
      d_extreme_y[0] = dpstack.back();
      dpstack.pop_back();
      d_extreme_x[0] = dpstack.back();
      dpstack.pop_back();
      if(debug){
	for (int t = 0; t < 2; t++) {
	  std::cout<<"d_extreme x "<< t <<" = " <<        d_extreme_x[t] <<std::endl;
	  std::cout<<"d_extreme y "<< t <<" = " <<        d_extreme_y[t] <<std::endl;
	}
      }
      
      kernelIsLeft <<< (n + 1023) / 1024, 1024 >>>(d_extreme_x_ptr, d_extreme_y_ptr,
							   d_x_ptr, d_y_ptr, flag_ptr, n);
      ps.pop_back();
      if(debug){
	for(int i = 0; i < n; i++)
	  std::cout << "flag[" << i << "] = " << flag[i] << std::endl;}

      // create some zip_iterators
      Float3Iterator P_first = thrust::make_zip_iterator(thrust::make_tuple(dx.begin(), dy.begin(), flag.begin()));
      Float3Iterator P_last = thrust::make_zip_iterator(thrust::make_tuple(dx.end(), dy.end(), flag.end()));

      // pass the zip_iterators into Partion()
      // Partion
      Float3Iterator first_of_R = thrust::partition(P_first, P_last,
						    is_interior_tuple());                   // Find Interior
  
      Float3Iterator first_of_L = P_first;

      
      if(debug){
	std::cout << "first of R - first of L  = " << first_of_R - first_of_L << std::endl;}

      int s = first_of_R - first_of_L;

      if((first_of_R - first_of_L) <2)
	{
	  dhull.push_back(d_extreme_y[0]);
	  dhull.push_back(d_extreme_x[0]);
	  if(debug){
	    std::cout << "pushed point x = " <<d_extreme_x[0]
		      << "pushed point y = " <<d_extreme_y[0] << std::endl;
	  }
	  if((first_of_R - first_of_L) > 0)
	    {
	      dhull.push_back(dy[0]);
	      dhull.push_back(dx[0]);
	      if(debug){
		std::cout << "pushed point x = " <<dx[0]
			  << "pushed point y = " <<dy[0] << std::endl;
	      }


	    }
			  
	}
      else
	{
	  
	  kernelCalcDist <<< (n + 1023) / 1024, 1024 >>>(d_extreme_x_ptr, d_extreme_y_ptr,
							 d_x_ptr, d_y_ptr, dist_ptr, flag_ptr, s);
	  if(debug){
	    for(int i = 0; i < n; i++)
	      std::cout << "dist[" << i << "] = " << dist[i] << std::endl;}


	  FloatIterator maxdist = thrust::max_element(dist.begin(),dist.begin()+s);
	  if(debug){
	    cout << "maxdist = : " <<*maxdist<<endl;
	    cout << "maxdist index = : " <<maxdist - dist.begin()<<endl;}

  	  if(debug){
	    for(int i = 0; i < n; i++)
	      std::cout << "partitioned flag[" << i << "] = " << flag[i] << std::endl;}

	  //p1 p2   // min - max: p0x p0y p1x p1y
	  //push back Hr
	  dpstack.push_back(dx[maxdist - dist.begin()]);
	  dpstack.push_back(dy[maxdist - dist.begin()]);
	  dpstack.push_back( d_extreme_x[1]);
	  dpstack.push_back( d_extreme_y[1]);
	  ps.push_back(s);

	  //Hl
	  dpstack.push_back( d_extreme_x[0]);
	  dpstack.push_back( d_extreme_y[0]);	  
	  dpstack.push_back(dx[maxdist - dist.begin()]);
	  dpstack.push_back(dy[maxdist - dist.begin()]);
	  ps.push_back(s);	  
  
	}

    }

  if(debug){ 
    for(int i = 0; i < dhull.size(); i++)
      std::cout << "dull[" << i << "] = " << dhull[i] << std::endl;}
  /*    
	vector<Coordinate> hull; 

	clock_t start, end;
	start = clock(); 

	hull = quickHull(v);
  */
  end = clock();
  // Calculating total time taken by the program. 
  float time_taken = float(end - start) / float(CLOCKS_PER_SEC); 
  cout << "Time taken by program is : " 
       << time_taken << " sec " << endl; 
    
  ofstream myfile;
  printf("Number of elements in hull: %ld\n", dhull.size());
  myfile.open("qhull.txt", ofstream::trunc);


  while(dhull.size()>0)
    {
      myfile << '(' << dhull.back();
      dhull.pop_back();
      myfile << ','  << dhull.back() << ')';
      dhull.pop_back();
      myfile << " ";
    }
  myfile.close();

  return 0;
}


