FILE?=pre_sorted
POINTS?=points
HULL?=polygon
SIZE?=1000

input:
	python create_input.py --output $(POINTS) --size $(SIZE)

all: test

test: presort plot

plot:
	python plot.py --quiet --pointsfile $(POINTS).txt --polygonfile $(HULL).txt

presort:
	python pre_sort.py --quiet --pointsfile $(POINTS).txt --output $(FILE).txt

seq: clean

clean:
	rm -rf *.exe *.exp *.lib