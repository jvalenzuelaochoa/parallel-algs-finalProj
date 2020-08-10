#define TOP_LEVEL
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include <fstream>
#include "../common/Coordinate.hpp"
#include "../common/common.cpp"
using namespace std;

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

    vector<Coordinate> v = parseFileForCoordinates(&inp);
    inp.close();

    const int ARRAY_SIZE = static_cast<int>(v.size());

    clock_t start, end;
    start = clock();

    vector<Coordinate> graham_stack;
    graham_stack.push_back(v[0]);
    graham_stack.push_back(v[1]);

    for(int i = 2; i < ARRAY_SIZE; i++) {
        Coordinate next = v[i];
        Coordinate p = graham_stack.back();
        graham_stack.pop_back();

        while((int)(graham_stack.size()) != 0 && ccw(graham_stack.back(), p, next) <= 0)
        {
            p = graham_stack.back();
            graham_stack.pop_back();
        }
        graham_stack.push_back(p);
        graham_stack.push_back(v[i]);
    }

    Coordinate p = graham_stack.back();
    graham_stack.pop_back();

    if(ccw(graham_stack.back(), p, v[0]) > 0)
    {
        graham_stack.push_back(p);
    }

    end = clock();
    // Calculating total time taken by the program.
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    printf( "Time taken by program is :%f\n", time_taken);

    storePolygon(graham_stack,"polygon.txt");

    return 0;
}
