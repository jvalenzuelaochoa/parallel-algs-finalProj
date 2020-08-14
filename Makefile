FILE?=pre_sorted
POINTS?=points
HULL?=polygon
SORTEDHULL?=sortedpolygon
SIZE?=10
CPP_SOURCES = common/Coordinate.cpp

.PHONY: all test clean mergehull quickhull monotone graham

all: test

input:
	python tools/create_input.py --output $(POINTS) --size $(SIZE) --grid 500000
	python tools/sort_points.py --quiet --pointsfile $(POINTS).txt --output $(FILE)

test: presort graham plot

plot:
	python tools/plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(HULL).txt

aftersort:
	python tools/sort_points.py --quiet --pointsfile $(HULL).txt --output $(SORTEDHULL)
	cp sortedpolygon.txt polygon.txt

graham:
	./graham

quickhull:
	./quickhull "pre_sorted.txt"

monotone:
	./monotone

testquickhull: build_omp quickhull aftersort plot

monotonenew: build_cu input monotone aftersort plot

qrun: input presort monotone aftersort plot

mergehull:
	./mergehull "pre_sorted_x.txt"

testmergehull: build_omp mergehull plot

mergehullnew: input build_omp mergehull

build: clean build_cu build_omp

build_cu:
	nvcc cuda/graham.cu $(CPP_SOURCES) -o graham
	nvcc cuda/monotone.cu $(CPP_SOURCES) -o monotone

build_omp:
	g++ -m64 omp/quickhull_omp.cpp $(CPP_SOURCES) -fopenmp -o quickhull
	g++ -m64 omp/mergehull_omp.cpp $(CPP_SOURCES) -fopenmp -o mergehull

run: build input graham plot

clean:
	rm -rf *.exe *.exp *.lib polygon.txt

cleanall: clean
	rm -rf *.txt
