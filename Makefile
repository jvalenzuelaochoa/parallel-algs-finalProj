FILE?=pre_sorted
POINTS?=points
HULL?=polygon
SORTEDHULL?=sortedpolygon
SIZE?=10
CPP_SOURCES = Coordinate.cpp

all: test

input:
	python create_input.py --output $(POINTS) --size $(SIZE) --grid 500000

test: presort seq plot

plot:
	python plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(HULL).txt

presort:
	python sort_points.py --quiet --pointsfile $(POINTS).txt --output $(FILE)

aftersort:
	python sort_points.py --quiet --pointsfile $(HULL).txt --output $(SORTEDHULL)

seq:
	./graham

quickhull:
	./quickhull "pre_sorted.txt" debug

testquickhull: build presort quickhull aftersort
	cp sortedpolygon.txt polygon.txt
	make plot

quickhullnew: input presort
	./quickhull_cu
	make aftersort
	cp sortedpolygon.txt polygon.txt
	make plot

quickhullnew: input testquickhull

mergehull:
	./mergehull "pre_sorted_x.txt" debug

testmergehull: build presort mergehull plot

build: clean
	nvcc graham.cu $(CPP_SOURCES) -o graham
	nvcc quickhull.cu -o quickhull_cu
	g++ quickhull_omp.cpp $(CPP_SOURCES) -fopenmp -o quickhull
	g++ mergehull_omp.cpp $(CPP_SOURCES) -fopenmp -o mergehull

run: build input presort seq plot

clean:
	rm -rf *.exe *.exp *.lib

cleanall: clean
	rm -rf *.txt
