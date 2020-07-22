import argparse
import random

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--output', type=str, default='nums', help='file to be generated(.txt)')
parser.add_argument('--size', dest='size', default = 1000, type=int,
                    help='Size of list generated')
parser.add_argument('--min', dest='min', type=int, default=0,
                    help='Minimum value desired for list')

args = parser.parse_args()

size = args.size
min = args.min
num_list = []

o_file_str = args.output + '.txt'
o_file = open(o_file_str, 'w')

# Genrate random int list
for i in range(size):
    num_list.append(random.randint(min, min + size))

# Make sure there's at least one instance of min in the list
if not(min in num_list):
    num_list[random.randint(0, size-1)] = min

o_file.write(', '.join(map(str, num_list)))
o_file.close