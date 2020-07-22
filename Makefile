OP_FILE?=output
FILE?=output
DEFAULT=inp
SIZE?=1000
MIN?=0

input:
	python create_input.py --output $(OP_FILE) --size $(SIZE)

all: hello

test:
	./q1 "output.txt"
	python verify_min.py --quiet --filename $(FILE).txt
	python verify_units.py --quiet --filename $(FILE).txt
	./q2 "output.txt"
	python verify_range.py --quiet --filename $(FILE).txt --output q2a
	python verify_range.py --quiet --filename $(FILE).txt --output q2b
	python verify_range.py --quiet --filename $(FILE).txt --output q2c --lump

hello:
	nvcc cuda_hello.cu -o hello
	./hello

hw3: clean q1
	python verify_min.py --quiet --filename $(DEFAULT).txt
	python verify_units.py --quiet --filename $(DEFAULT).txt

q1:
	nvcc q1.cu -o q1
	./q1

q2:
	nvcc q2.cu -o q2
	./q2

clean:
	rm -f *.exe *.lib *.exp