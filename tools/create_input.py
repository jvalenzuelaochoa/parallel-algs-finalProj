import argparse
import random

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--output', type=str, default='nums', help='file to be generated(.txt)')
parser.add_argument('--size', dest='size', default = 1000, type=int,
                    help='Size of list generated')
parser.add_argument('--grid', dest='grid', type=int, default=2000,
                    help='desired grid size')
parser.add_argument('--tridimensional', action='store_true', default=False)

args = parser.parse_args()

size = args.size
grid = args.grid
if grid < size:
    grid = size

num_list = []

o_file_str = args.output + '.txt'
o_file = open(o_file_str, 'w')

xs = random.sample(range(-1*grid, grid), size)
ys = random.sample(range(-1*grid, grid), size)
zs = random.sample(range(-1*grid, grid), size)

# Genrate random int list
for x, y, z in zip(xs,ys,zs):
    if args.tridimensional:
        num_list.append((x, y, z))
    else:
        num_list.append((x, y))

o_file.write(' '.join(map(str, num_list)))
o_file.close
