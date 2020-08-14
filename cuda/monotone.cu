#define TOP_LEVEL
#include <iostream>
#include <omp.h>
#include <cmath>
#include <vector>
#include <set>
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

#include "../common/Coordinate.hpp"
#include "../common/common.cpp"

using namespace std;
#define BLOCK_SIZE 1024

// Make predicate easy to write
typedef thrust::tuple<float, float, int> FloatTuple3;
typedef thrust::tuple<float, float, int, int,int,int> FloatTuple6;
// Predicate
struct is_interior_tuple {
  __host__ __device__

  bool operator()(const FloatTuple3 &p) {
    return thrust::get<2>(p) > 0;
  }
};

// Predicate
struct is_interior_6_tuple {
  __host__ __device__

  bool operator()(const FloatTuple6 &p) {
    return thrust::get<2>(p) > 0;
  }
};

//reference: https://stackoverflow.com/questions/38923671/thrust-cuda-find-maximum-per-each-groupsegment
struct my_max_func
{

  template <typename T1, typename T2>
  __host__ __device__
  T1 operator()(const T1 t1, const T2 t2){
    T1 res;
    if (thrust::get<0>(t1) > thrust::get<0>(t2)){
      thrust::get<0>(res) = thrust::get<0>(t1);
      thrust::get<1>(res) = thrust::get<1>(t1);}
    else {
      thrust::get<0>(res) = thrust::get<0>(t2);
      thrust::get<1>(res) = thrust::get<1>(t2);}
    return res;
  }
};



__host__ __device__
bool cmpf(float A, float B, float epsilon = 0.0001f)
{
  return (fabs(A - B) < epsilon);
}

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
int isLeft(float& p0x, float& p0y, float& p1x, float& p1y, float& p2x, float& p2y) {
  ///      float val = (p.y - p1.y) * (p2.x - p1.x) -
  //      (p2.y - p1.y) * (p.x - p1.x);
  bool e1, e2;

  e1 = cmpf(p2x, p1x) && cmpf(p2y,p1y);
  e2 = cmpf(p2x, p0x) && cmpf(p2y,p0y);


  if(e1 || e2 )
    return 1;
  else
    {
      float val =  (p1x - p0x) * (p2y - p0y) - (p2x - p0x) * (p1y - p0y);
      if (val >=0) return 1;
    }
  return 0;

}
//reference https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
__host__ __device__

float sign (float& p1x, float& p1y, float& p2x,float& p2y, float& p3x,float& p3y)
{
  return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
}

//pt is the point under test
//v1 is the first point
//v2 is the fartheset point
//v3 is the last point
__host__ __device__

int PointInTriangle (float& ptx, float& pty, float& v1x,float& v1y, float& v2x,float& v2y, float& v3x, float& v3y)
{
  float d1, d2, d3;

  bool has_neg, has_pos;
  bool e1, e2, e3;

  e1 = cmpf(ptx, v1x) && cmpf(pty,v1y);
  e2 = cmpf(ptx, v2x) && cmpf(pty,v2y);
  e3 = cmpf(ptx, v3x) && cmpf(pty,v3y);

  if(e1 || e2 || e3)
    return 0;

  d1 = sign(ptx, pty, v1x, v1y, v2x, v2y);
  d2 = sign(ptx,pty, v2x,v2y, v3x,v3y);
  d3 = sign(ptx,pty, v3x,v3y, v1x, v1y);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  if(!(has_neg && has_pos))
    return 1;
  else
    {
      if(v1x < v3x) //upside triangle
	{
	  //isLeft(float& p0x, float& p0y, float& p1x, float& p1y, float& p2x, float& p2y) p2 is the target point
	  if(isLeft( v1x, v1y,  v2x,v2y, ptx, pty)==1 ||isLeft( v2x, v2y,  v3x,v3y, ptx, pty)==1 )
	    return 0;
	  else
	    return 1;
	}
      else // upside down
	{
	  //isLeft(float& p0x, float& p0y, float& p1x, float& p1y, float& p2x, float& p2y) p2 is the target point
	  if(isLeft( v1x, v1y,  v2x,v2y, ptx, pty)==0 ||isLeft( v2x, v2y,  v3x,v3y, ptx, pty)==0 )
	    return 0;
	  else
	    return 1;
	  
	}
      return 0;
    }
  //    return !(has_neg && has_pos);
}


// function from https://www.geeksforgeeks.org/perpendicular-distance-between-a-point-and-a-line-in-2-d/
// Function to find distance
__host__ __device__

float shortest_distance(float &P0_x,float &P0_y,float &P1_x,float &P1_y, float &P2_x, float &P2_y)
{
  // P0 : p1
  // P1 : p2
  // P2_x : p
  float x1, y1, a, b, c;


  x1 = P2_x;
  y1 = P2_y;

  a = P0_y-P1_y;
  b = P1_x-P0_x;
  c = (P0_x-P1_x)*P0_y + (P1_y-P0_y)*P0_x;


  float d = fabs((a * x1 + b * y1 + c)) /
    (sqrt(a * a + b * b));
  //printf("Perpendicular distance is %f\n", d);
  return d;
}
// Preprocess by discarding interior points, i.e., 1st round of discarding
__global__
void kernelCalcDist(int *first_pts, int *last_pts, int *keys, float *v_x, float *v_y, float *dist, int n) {


  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n ) // Check
    {
      //get first and last point

      dist[i] = shortest_distance(v_x[first_pts[i]],v_y[first_pts[i]],v_x[last_pts[keys[i]]],v_y[last_pts[keys[i]]], v_x[i], v_y[i]);
    }
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
      flag[i] = isLeft(s_extreme[0].x, s_extreme[0].y,s_extreme[1].x, s_extreme[1].y, v_x[i], v_y[i]);
    }
}


__global__
void kernelIsFirst(int *head, int * first_pts,  int n) {

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n) // Check
    {
      if( head[i] == 1)
	first_pts[i] = i;
      else
	first_pts[i] = 0;
    }
}

__global__
void kernelUpdateHead(int *head, int *d_idxs_out, int n) {

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n) // Check
    {
      head[d_idxs_out[i]] = 1;

    }
}

__global__
void kernelLabelInterior(int *first_pts, int *last_pts, int *d_idxs_out, int *keys, int *flag, float *v_x, float *v_y, int n) {


  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n ) // Check
    {
      //if interior flag = 0 otherwise flag =1
      //float& ptx, float& pty, float& v1x,float& v1y, float& v2x,float& v2y, float& v3x, float& v3y)
      if(PointInTriangle (v_x[i], v_y[i],v_x[first_pts[i]],v_y[first_pts[i]],
			  v_x[d_idxs_out[keys[i]]],v_y[d_idxs_out[keys[i]]],
			  v_x[last_pts[keys[i]]],v_y[last_pts[keys[i]]]) == 1)
	flag[i] = 0;
      else
	flag[i] = 1;

    }


}

int main(int argc, char **argv)
{
  bool debug = false;
  //bool debug = true;
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


  const int constn = static_cast<int>(v.size());
  int n = v.size();
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

  if(debug)
    {
      cout << "The min_x is: " <<min_x<<", the max_x is "<<max_x << endl;
    }
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
  if(debug)
    {

      cout << "minx = : " <<*minx<<endl;
      cout << "maxx = : " <<*maxx<<endl;
      cout << "minx index = : " <<minx - dx.begin()<<endl;
      cout << "maxx index = : " <<maxx - dx.begin()<<endl;
    }
  //cout << "maxy = : " <<dy[maxx - dx.begin()]<<endl;
  float maxy = dy[maxx - dx.begin()];
  if(debug)
    {

      cout << "maxy = : " <<maxy<<endl;
    }
  // Store the four extreme points temporarily
  thrust::device_vector<float> d_extreme_x(2);
  thrust::device_vector<float> d_extreme_y(2);

  d_extreme_x[0] = *minx;
  d_extreme_x[1] = *maxx;

  d_extreme_y[0] = dy[minx - dx.begin()];
  d_extreme_y[1] = dy[maxx - dx.begin()];

  if(debug)
    {

      cout << "d_extreme_x[0] = : " <<d_extreme_x[0]<<endl;
      cout << "d_extreme_y[0] = : " <<d_extreme_y[0]<<endl;
      cout << "d_extreme_x[1] = : " <<d_extreme_x[1]<<endl;
      cout << "d_extreme_y[1] = : " <<d_extreme_y[1]<<endl;
    }



  // Get the pointers to the arrays, to be used as launch arguments
  float *d_extreme_x_ptr = thrust::raw_pointer_cast(&d_extreme_x[0]);
  float *d_extreme_y_ptr = thrust::raw_pointer_cast(&d_extreme_y[0]);
  float *d_x_ptr = thrust::raw_pointer_cast(&dx[0]);
  float *d_y_ptr = thrust::raw_pointer_cast(&dy[0]);
  int *flag_ptr = thrust::raw_pointer_cast(&flag[0]);
  int *first_pts_ptr = thrust::raw_pointer_cast(&first_pts[0]);
  int *head_ptr = thrust::raw_pointer_cast(&head[0]);
  float *dist_ptr = thrust::raw_pointer_cast(&dist[0]);
  int *keys_ptr = thrust::raw_pointer_cast(&keys[0]);

  // 1st discarding :  Block size can be 1024 in this kernel
  kernelIsLeft <<< (n + 1023) / 1024, 1024 >>>(d_extreme_x_ptr, d_extreme_y_ptr,
					       d_x_ptr, d_y_ptr, flag_ptr, n);
  if(debug)
    {

      for(int i = 0; i < n; i++)
	std::cout << "flag[" << i << "] = " << flag[i] << std::endl;
    }
  // Defining a zip_iterator type can be a little cumbersome ...
  // same size vector are dx, dy, flag, head, keys, first_pts
  typedef thrust::device_vector<float>::iterator FloatIterator;
  typedef thrust::device_vector<int>::iterator IntIterator;
  typedef thrust::tuple<FloatIterator, FloatIterator, IntIterator> FloatIteratorTuple;
  typedef thrust::zip_iterator<FloatIteratorTuple> Float3Iterator;

  // Defining a zip_iterator type can be a little cumbersome ...
  // same size vector are dx, dy, flag, head, keys, first_pts
  //typedef thrust::device_vector<float>::iterator FloatIterator;
  //typedef thrust::device_vector<int>::iterator IntIterator;
  typedef thrust::tuple<FloatIterator, FloatIterator, IntIterator, IntIterator,IntIterator,IntIterator> FloatIteratorTupleAll;
  typedef thrust::zip_iterator<FloatIteratorTupleAll> Float6Iterator;

  // create some zip_iterators
  Float3Iterator P_first = thrust::make_zip_iterator(thrust::make_tuple(dx.begin(), dy.begin(), flag.begin()));
  Float3Iterator P_last = thrust::make_zip_iterator(thrust::make_tuple(dx.end(), dy.end(), flag.end()));

  // pass the zip_iterators into Partion()
  // Partion
  Float3Iterator first_of_R = thrust::stable_partition(P_first, P_last,
						       is_interior_tuple());                   // Find Interior
  Float3Iterator first_of_L = P_first;
  if(debug)
    {

      std::cout<<" first_of_R - first_of_L = " << first_of_R -first_of_L << std::endl;
    }
  //head init
  head[0] = 1;
  //if all points in flag are 1 or   //if all points in flag are 0
  if((first_of_R - first_of_L) == head.size() || (first_of_R - first_of_L) ==0)
    head[head.size()-1] = 1;
  else
    //if some are zeros some are 1
    {
      	head[first_of_R - first_of_L-1] = 1;
	/*
	  if(dx[first_of_R - first_of_L] > dx[first_of_R - first_of_L - 1] )
	  head[first_of_R - first_of_L] = 1;
	  else
	  head[first_of_R - first_of_L - 1] = 1;
	*/
    }
  thrust::inclusive_scan(head.begin(), head.end(), keys.begin()); // in-place scan

  if(debug)
    {

      for(int i = 0; i < n; i++)
	{
	  std::cout << "partitioned flag[" << i << "] = " << flag[i] << std::endl;

	}
    }
  thrust::fill(flag.begin(), flag.end(), 0);
  if(debug)
    {

      for(int i = 0; i < n; i++)
	{
	  std::cout << "Zeroed flag[" << i << "] = " << flag[i] << std::endl;

	}

      for(int i = 0; i < n; i++)
	{
	  std::cout << "head[" << i << "] = " << head[i] << std::endl;
	}
      for(int i = 0; i < n; i++)
	{
	  std::cout << "keys[" << i << "] = " << keys[i] << std::endl;
	}
    }

  //subtract 1 from keys

  thrust::transform(keys.begin(),
		    keys.end(),
		    thrust::make_constant_iterator(1),
		    keys.begin(),
		    thrust::minus<int>());
  if(debug)
    {

      for(int i = 0; i < n; i++)
	{
	  std::cout << "subtracted keys[" << i << "] = " << keys[i] << std::endl;
	}
    }
  /*
    thrust::device_vector<int> gindex(n);
    thrust::sequence(gindex.begin(), gindex.end());
    for(int i = 0; i < n; i++)
    {
    std::cout << "global index[" << i << "] = " << gindex[i] << std::endl;
    }
  */
  //first_pts
  kernelIsFirst <<< (n + 1023) / 1024, 1024 >>>(head_ptr,first_pts_ptr, n);
  if(debug)
    {

      for(int i = 0; i < n; i++)
	{
	  std::cout << "first_pts[" << i << "] = " << first_pts[i] << std::endl;
	}
    }
  thrust::inclusive_scan_by_key(keys.begin(), keys.end(), first_pts.begin(), first_pts.begin()); // in-place
  if(debug)
    {

      for(int i = 0; i < n; i++)
	{
	  std::cout << "inclusive scanned first_pts[" << i << "] = " << first_pts[i] << std::endl;
	}
    }

  // Get the position of the first points in each sub-region
  FloatIteratorTuple pos_L = first_of_L.get_iterator_tuple();
  FloatIteratorTuple pos_R = first_of_R.get_iterator_tuple();
  FloatIteratorTuple pos_last = P_last.get_iterator_tuple();

  // Partly Sort for each sub-regions
  //  : ascending X and descending X
  FloatIterator first_of_L_x = thrust::get<0>(pos_L);
  FloatIterator first_of_R_x = thrust::get<0>(pos_R);
  FloatIterator last_of_R_x = thrust::get<0>(pos_last);
  thrust::stable_sort_by_key(first_of_L_x, first_of_R_x, first_of_L);
  thrust::stable_sort_by_key(first_of_R_x, last_of_R_x, first_of_R);


  if(debug)
    {
      std::cout << "dx[first_or_R -first_of_L] = " << dx[first_of_R -first_of_L -1] << std::endl;
    }
  float tempMaxy = dy[first_of_R -first_of_L - 1] ;
  if(debug)
    {

      std::cout << "after tempMaxy" << std::endl;
    }
  /*
  if(tempMaxy != maxy)
    {
        if(debug)
    {

      std::cout << "if tempMaxy != maxy" << std::endl;
    }

      dx[first_of_R -first_of_L-1] = *maxx+0.0001;

      std::cout << "dx[first_or_R -first_of_L] = " << dx[first_of_R -first_of_L -1] << std::endl;

      thrust::stable_sort_by_key(first_of_L_x, first_of_R_x, first_of_L);
      thrust::stable_sort_by_key(first_of_R_x, last_of_R_x, first_of_R);
      }*/
  // Sort in ascending order, and then reverse
  if(first_of_R - first_of_L !=0){
    thrust::reverse(thrust::get<0>(pos_R), thrust::get<0>(pos_last));
    thrust::reverse(thrust::get<1>(pos_R), thrust::get<1>(pos_last));
    thrust::reverse(thrust::get<2>(pos_R), thrust::get<2>(pos_last));
  }

  if(debug)
    {
      for(int i = 0; i < n; i++)
	std::cout << "sorted dx[" << i << "] = " << dx[i] << std::endl;
      //std::cout << "dy[24] = " << dy[24] << std::endl;
    }
  int headcount;

  clock_t start, end;
  start = clock();

  //recursive step
  while(true){

    headcount = thrust::count(head.begin(), head.end(), 1);
    if(debug)
      {
	std::cout << "n = " << n << std::endl;
	std::cout << "headcount = " << headcount << std::endl;
      }
    thrust::device_vector<int> allone(n);
    thrust::device_vector<int> last_pts(headcount);
    thrust::device_vector<int> reducedkeys(headcount);
    thrust::fill(allone.begin(), allone.end(), 1);
    typedef thrust::device_vector<int>::iterator IntIter;
    int *last_pts_ptr = thrust::raw_pointer_cast(&last_pts[0]);

    // Find the four extreme points, i.e., min x, max x, min y, and max y
    // thrust::pair<FloatIter, FloatIter> extremex = thrust::minmax_element(dx.begin(), dx.end());

    thrust::pair<IntIter, IntIter> new_end = thrust::reduce_by_key(keys.begin(), keys.end(), allone.begin(), reducedkeys.begin(), last_pts.begin());
    thrust::inclusive_scan(last_pts.begin(), last_pts.end(), last_pts.begin()); // in-place scan
    last_pts[headcount-1] = 0;
    if(debug)
      {
	for(int i = 0; i < last_pts.size(); i++)
	  std::cout << "last_pts[" << i << "] = " << last_pts[i] << std::endl;
	for(int i = 0; i < reducedkeys.size(); i++)
	  std::cout << "reducedkeys[" << i << "] = " << reducedkeys[i] << std::endl;
      }


    kernelCalcDist <<< (n + 1023) / 1024, 1024 >>>(first_pts_ptr, last_pts_ptr, keys_ptr, d_x_ptr, d_y_ptr, dist_ptr, n);
    if(debug)
      {
	for(int i = 0; i < n; i++)
	  std::cout << "dist[" << i << "] = " << dist[i] << std::endl;
      }

    //max by id ref: https://stackoverflow.com/questions/38923671/thrust-cuda-find-maximum-per-each-groupsegment
    // thrust method
    //thrust::device_vector<float> d_vals(h_vals, h_vals+vsize);
    //thrust::device_vector<int> d_keys(h_keys, h_keys+vsize);
    thrust::device_vector<int> d_keys_out(headcount);
    thrust::device_vector<float> d_vals_out(headcount);
    thrust::device_vector<int> d_idxs(n);
    thrust::device_vector<int> d_idxs_out(headcount);
    int *d_idxs_out_ptr = thrust::raw_pointer_cast(&d_idxs_out[0]);

    thrust::sequence(d_idxs.begin(), d_idxs.end());
    if(debug)
      {
	for(int i = 0; i < headcount; i++)
	  std::cout << "before max by id d_idxs_out[" << i << "] = " << d_idxs_out[i] << std::endl;
      }
    cudaDeviceSynchronize();
    //unsigned long long et = dtime_usec(0);

    //reference: https://stackoverflow.com/questions/38923671/thrust-cuda-find-maximum-per-each-groupsegment    
    //thrust::sort_by_key(keys.begin(), keys.begin()+n, thrust::make_zip_iterator(thrust::make_tuple(dist.begin(), d_idxs.begin())));
    thrust::reduce_by_key(keys.begin(), keys.begin()+n,
			  thrust::make_zip_iterator(thrust::make_tuple(dist.begin(),d_idxs.begin())), d_keys_out.begin(),
			  thrust::make_zip_iterator(thrust::make_tuple(d_vals_out.begin(), d_idxs_out.begin())),
			  thrust::equal_to<int>(), my_max_func());

    cudaDeviceSynchronize();
    //et = dtime_usec(et);
    //std::cout << "Thrust time: " << et/(float)USECPSEC << "s" << std::endl;
    if(debug)
      {
	for(int i = 0; i < headcount; i++)
	  std::cout << "max by id dist d_vals_out[" << i << "] = " << d_vals_out[i] << std::endl;
      }

    if(debug)
      {
	for(int i = 0; i < headcount; i++)
	  std::cout << "max by id d_keys_out[" << i << "] = " << d_keys_out[i] << std::endl;
      }
    if(debug)
      {
	for(int i = 0; i < headcount; i++)
	  std::cout << "max by id d_idxs_out[" << i << "] = " << d_idxs_out[i] << std::endl;
      }

    //line 8
    // detect interior points
    if(debug)
      {
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "before detected interior flag[" << i << "] = " << flag[i] << std::endl;

	  }
      }
    /*
      int tempi = PointInTriangle (dx[5], dy[5],dx[first_pts[5]],dy[first_pts[5]],
      dx[d_idxs_out[keys[5]]],dy[d_idxs_out[keys[5]]],
      dx[last_pts[keys[5]]],dy[last_pts[keys[5]]]);
      std::cout<< "temp i = " << tempi<<std::endl;

      std::cout<<dx[5]<< dy[5]<<dx[first_pts[5]]<<dy[first_pts[5]]<<
      dx[d_idxs_out[keys[5]]]<<dy[d_idxs_out[keys[5]]]<<
      dx[last_pts[keys[5]]]<<dy[last_pts[keys[5]]])<< std::endl;

    */

    if(debug){
      std::cout<<"PointIntriangle "<<dx[5]<<" " << dy[5]<<" " <<dx[first_pts[5]]<<" " <<dy[first_pts[5]]
	       <<" " <<dx[d_idxs_out[keys[5]]]<<" " <<dy[d_idxs_out[keys[5]]]
	       <<" " <<dx[last_pts[keys[5]]]<<" " <<dy[last_pts[keys[5]]]<< std::endl;
    }
    kernelLabelInterior <<< (n + 1023) / 1024, 1024 >>>(first_pts_ptr, last_pts_ptr,
							d_idxs_out_ptr, keys_ptr, flag_ptr,
							d_x_ptr, d_y_ptr, n);

    if(debug)
      {
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "detected interior flag[" << i << "] = " << flag[i] << std::endl;

	  }
      }

    //stable_partition and then resize()


    //update heads
    kernelUpdateHead <<< (headcount + 1023) / 1024, 1024 >>>(head_ptr, d_idxs_out_ptr, headcount);
    if(debug)
      {
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "head[" << i << "] = " << head[i] << std::endl;
	  }
      }
    //update keys and first_pts
    thrust::inclusive_scan(head.begin(), head.end(), keys.begin()); // in-place scan
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "keys[" << i << "] = " << keys[i] << std::endl;
	  }
      }
    //subtract 1 from keys

    thrust::transform(keys.begin(),
		      keys.end(),
		      thrust::make_constant_iterator(1),
		      keys.begin(),
		      thrust::minus<int>());
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "subtracted keys[" << i << "] = " << keys[i] << std::endl;
	  }
      }
    //first_pts
    kernelIsFirst <<< (n + 1023) / 1024, 1024 >>>(head_ptr,first_pts_ptr, n);
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }
      }
    thrust::inclusive_scan_by_key(keys.begin(), keys.end(), first_pts.begin(), first_pts.begin()); // in-place
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "inclusive scanned first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "before resize flag[" << i << "] = " << flag[i] << std::endl;

	  }
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "before resize first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "before resize head[" << i << "] = " << head[i] << std::endl;
	  }


	for(int i = 0; i < n; i++)
	  {
	    std::cout << "before resize keys[" << i << "] = " << keys[i] << std::endl;
	  }
      }
    // create some zip_iterators
    Float6Iterator stable_par_first = thrust::make_zip_iterator(thrust::make_tuple(dx.begin(), dy.begin(),
										   flag.begin(), head.begin(),
										   keys.begin(), first_pts.begin()));
    Float6Iterator stable_par_last =thrust::make_zip_iterator(thrust::make_tuple(dx.end(), dy.end(),
										 flag.end(), head.end(),
										 keys.end(), first_pts.end()));

    // pass the zip_iterators into Partion()
    // Partion
    Float6Iterator stable_par_R = thrust::stable_partition(stable_par_first, stable_par_last,
							   is_interior_6_tuple());                   // Find Interior
    Float6Iterator stable_par_L = stable_par_first;
    if(debug)
      {
	std::cout << "n = " << n << std::endl;
      }
    if(n == stable_par_R - stable_par_L)
      break;
    n = stable_par_R - stable_par_L;
    if(debug)
      {

	std::cout << " exteriorpointcount = stable_par_R - stable_par_L = " << n << std::endl;
      }
    //resize()
    /*
      dx.resize(exteriorpointcount);
      dy.resize(exteriorpointcount);
      flag.resize(exteriorpointcount);
      head.resize(exteriorpointcount);
      keys.resize(exteriorpointcount);
      first_pts.resize(exteriorpointcount);
      dist.resize(exteriorpointcount);
    */
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "after resize flag[" << i << "] = " << flag[i] << std::endl;

	  }
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "after resize first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "after resize head[" << i << "] = " << head[i] << std::endl;
	  }


	for(int i = 0; i < n; i++)
	  {
	    std::cout << "after resize keys[" << i << "] = " << keys[i] << std::endl;
	  }
	for(int i = 0; i < n; i++)
	  {
	    std::cout << "after resize dx[" << i << "] = " << dx[i] << std::endl;
	  }

      }
    //update first points
    kernelIsFirst <<< (n + 1023) / 1024, 1024 >>>(head_ptr,first_pts_ptr, n);

    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "2nd update first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }
      }
    thrust::inclusive_scan_by_key(keys.begin(), keys.end(), first_pts.begin(), first_pts.begin()); // in-place
    if(debug)
      {

	for(int i = 0; i < n; i++)
	  {
	    std::cout << "2nd update inclusive scanned first_pts[" << i << "] = " << first_pts[i] << std::endl;
	  }

      }
  }

  end = clock();
  // Calculating total time taken by the program.
  float time_taken = float(end - start) / float(CLOCKS_PER_SEC);
  cout << "Time taken by program is : "
       << time_taken << " sec " << endl;

  printf("Number of elements in hull: %d\n", n);
  ofstream myfile;
  myfile.open("polygon.txt", ofstream::trunc);

  int i = 0;
  while(true)
    {
      myfile << '(' << dx[i] << ','  << dy[i] << ')';

      if (++i >= n) break;

      myfile << " ";
    }
  myfile.close();

  return 0;
}


