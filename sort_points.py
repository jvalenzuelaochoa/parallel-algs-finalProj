import matplotlib.pyplot as plt
import argparse
import math

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--pointsfile', type=str, help='input file to be analyzed')
parser.add_argument('--output', type=str, help='output file to be analyzed')
parser.add_argument('--quiet', action='store_true', default=False)

args = parser.parse_args()

def parse_coordinates(file):
    file1 = open(file, 'r')
    Lines = file1.readlines()
    for line in Lines:
        str_arr = line.split(')')

    num_list = []
    for i in range(len(str_arr)):
        str_arr[i] = str_arr[i].replace(' ', '')
        str_arr[i] = str_arr[i].replace('(', '')
        coordinates = str_arr[i].split(',')
        if len(str_arr[i]) > 2:
            num_list.append((float(coordinates[0]), float(coordinates[1])))

    return num_list

# Helper Class to handle Coodinates
class Coordinate:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.angle = 0.0

    def __repr__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ' : ' + str(self.angle)

coord = parse_coordinates(args.pointsfile)

points = []

for point in coord:
    points.append(Coordinate(point[0], point[1]))

# get the lowest point in the grid and save it as a reference
points.sort(key=lambda p: p.y, reverse=False)
ref = points[0]

# Compute the angle between all points to the reference
for point in points:
    if point != ref:
        point.angle = math.atan2(point.y-ref.y, point.x-ref.x)

# resort elements by their angle to the reference
points.sort(key=lambda p: p.angle, reverse=False)

# store points in file
o_file_str = args.output + '.txt'
o_file = open(o_file_str, 'w')

num_list =[]
for point in points:
    num_list.append((point.x, point.y))

o_file.write(' '.join(map(str, num_list)))
o_file.close

# resort elements by their angle to the reference
points.sort(key=lambda p: p.x, reverse=False)

# store points in file
o_file_str = args.output + '_x.txt'
o_file = open(o_file_str, 'w')

num_list =[]
for point in points:
    num_list.append((point.x, point.y))

o_file.write(' '.join(map(str, num_list)))
o_file.close