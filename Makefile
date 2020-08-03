FILE?=pre_sorted
POINTS?=points
HULL?=polygon
SORTEDHULL?=sortedpolygon
SIZE?=100

input:
	python create_input.py --output $(POINTS) --size $(SIZE)

all: test

test: presort seq plot

plot:
	python plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(HULL).txt

plotsorted:
	python plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(SORTEDHULL).txt

presort:
	python pre_sort.py --quiet --pointsfile $(POINTS).txt --output $(FILE).txt

aftersort:
	python pre_sort.py --quiet --pointsfile $(HULL).txt --output $(SORTEDHULL).txt

seq:
	./graham

testquickhull: 
	g++ quickhull_omp.cpp -fopenmp
	./a.out
	make aftersort
	make plotsorted

quickhullnew: 
	make input
	make presort
	g++ quickhull_omp.cpp -fopenmp
	./a.out
	make aftersort
	make plotsorted


build: clean
	nvcc graham.cu -o graham

clean:
	rm -rf *.exe *.exp *.lib
