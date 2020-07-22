import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--pointsfile', type=str, help='input file to be analyzed')
parser.add_argument('--polygonfile', type=str, help='input file to be analyzed')
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
            num_list.append([coordinates[0], coordinates[1]])

    return num_list

coord = parse_coordinates(args.polygonfile)

coord.append(coord[0]) #repeat the first point to create a 'closed loop'

xs, ys = zip(*coord) #create lists of x and y values

plt.figure()
# Add Polygon to plot
plt.plot(xs,ys)

# Add all points to plot
coord = parse_coordinates(args.pointsfile)
for (x,y) in coord:
    circle = circle = plt.Circle((x, y), radius=0.3, fc='y')
    plt.gca().add_patch(circle)

# Display plot
plt.show()