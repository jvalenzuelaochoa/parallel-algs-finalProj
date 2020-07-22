#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include <fstream>
using namespace std;


class Coordinate {
    public:
      int x;
      int y;

    Coordinate()
    {
        x = 0;
        y = 0;
    }

    Coordinate(int x1, int y1)
    {
       x = x1;
       y = y1;
    }

    Coordinate(const Coordinate &c)
    {
       x = c.x;
       y = c.y;
    }

  };

int ccw(Coordinate a, Coordinate b, Coordinate c){
    float area = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
    if (area < 0 ) return -1;
    if (area > 0 ) return  1;
    return 0;

}

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
    int x_temp = 0;
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
            x_temp = stoi(sNum);
        }
        else
        {
            c_temp = new Coordinate(x_temp, stoi(sNum));
            v.push_back(*c_temp);
        }
        x_next = !x_next;
    }

    printf("Number of elements in file: %d\n", static_cast<int>(v.size()));

    const int ARRAY_SIZE = static_cast<int>(v.size());


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

    ofstream myfile;
    myfile.open("polygon.txt", ofstream::trunc);

    int i = 0;
    while(true)
    {
        myfile << '(' << graham_stack[i].x << ','  << graham_stack[i].y << ')';

        if (++i >= (int)(graham_stack.size())) break;

        myfile << " ";
    }
    myfile.close();

    return 0;
}