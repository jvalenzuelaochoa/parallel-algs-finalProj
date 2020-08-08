FILE?=sort_pointsed
POINTS?=points
HULL?=polygon
SORTEDHULL?=sortedpolygon
SIZE?=10
CPP_SOURCES = Coordinate.cpp

all: test

input:
	python create_input.py --output $(POINTS) --size $(SIZE)

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
	./quickhull "sort_pointsed.txt" debug

testquickhull: build input presort quickhull aftersort
	cp sortedpolygon.txt polygon.txt
	make plot

quickhullnew: input testquickhull

mergehull:
	./mergehull "sort_pointsed_x.txt" debug

testmergehull: build input presort mergehull aftersort
	cp sortedpolygon.txt polygon.txt
	make plot

build: clean
	nvcc graham.cu $(CPP_SOURCES) -o graham
	g++ quickhull_omp.cpp $(CPP_SOURCES) -fopenmp -o quickhull
	g++ mergehull_omp.cpp $(CPP_SOURCES) -fopenmp -o mergehull

run: build input presort seq plot

clean:
	rm -rf *.exe *.exp *.lib

cleanall: clean
	rm -rf *.txt
