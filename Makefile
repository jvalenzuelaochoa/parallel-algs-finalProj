FILE?=pre_sorted
POINTS?=points
HULL?=polygon
SORTEDHULL?=sortedpolygon
SIZE?=10
CPP_SOURCES = common/Coordinate.cpp

.PHONY: all test clean mergehull quickhull quickhull_cu seq

all: test

input:
	python tools/create_input.py --output $(POINTS) --size $(SIZE) --grid 500000

test: presort seq plot

plot:
	python tools/plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(HULL).txt

presort:
	python tools/sort_points.py --quiet --pointsfile $(POINTS).txt --output $(FILE)

aftersort:
	python tools/sort_points.py --quiet --pointsfile $(HULL).txt --output $(SORTEDHULL)
	cp sortedpolygon.txt polygon.txt

seq:
	./graham

quickhull:
	./quickhull "pre_sorted.txt"

quickhull_cu:
	./quickhull_cu

testquickhull: build_omp presort quickhull aftersort plot

quickhullnew: build_cu input presort quickhull_cu aftersort plot

qrun: input presort quickhull_cu aftersort plot

mergehull:
	./mergehull "pre_sorted_x.txt"

testmergehull: build_omp presort mergehull plot

mergehullnew: input build_omp presort mergehull

build: clean build_cu build_omp

build_cu:
	nvcc cuda/graham.cu $(CPP_SOURCES) -o graham
	nvcc cuda/quickhull.cu $(CPP_SOURCES) -o quickhull_cu

build_omp:
	g++ -m64 omp/quickhull_omp.cpp $(CPP_SOURCES) -fopenmp -o quickhull
	g++ -m64 omp/mergehull_omp.cpp $(CPP_SOURCES) -fopenmp -o mergehull

run: build input presort seq plot

clean:
	rm -rf *.exe *.exp *.lib polygon.txt

cleanall: clean
	rm -rf *.txt
