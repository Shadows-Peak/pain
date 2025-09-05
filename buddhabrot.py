import numpy as np
import math
from rtree import index
from PIL import Image

catch_distance = 10
pointEvals = 30
num_points = 200000000 # Recommended to be higher than the number of pixels
dimensions = 2

precision = 10000 # Subdivide each side into this many segments for pixels
resolution = [10000,10000]

# Define the range for each dimension
x_min, x_max = -2, 1
y_min, y_max = -1, 1

# Generate random x and y coordinates separately
x_coords = np.random.uniform(x_min, x_max, num_points)
y_coords = np.random.uniform(y_min, y_max, num_points)

# Combine them into an array of points
random_points = np.column_stack((x_coords, y_coords))



pointLength = len(random_points)

escapedPoints = []
allPointLocations = []

print("Calculating points...")

for point in random_points:
    pointLocations = []
    zFirst = complex(point[0], point[1])
    zNext = zFirst
    for i in range(pointEvals):
        try:
            zNext = zNext**2 + zFirst
        except OverflowError:
            escapedPoints.append(point)
            allPointLocations += pointLocations
            break
        if abs(zNext) > catch_distance or np.isnan(abs(zNext)):
            escapedPoints.append(point)
            allPointLocations += pointLocations
            break
        pointLocations.append(zNext)
    if len(allPointLocations) % (pointLength // 100) == 0 and len(allPointLocations) != 0:
        print(f"Point calculations: {int((len(allPointLocations) / (pointLength * pointEvals)) * 100)}% done")

print("Done calculating points.")
#print(allPointLocations)

xs = [z.real for z in allPointLocations]
ys = [z.imag for z in allPointLocations]
xs.sort()
ys.sort()
print("X range: ", (min(xs) if xs else None, max(xs) if xs else None))
print("Y range: ", (min(ys) if ys else None, max(ys) if ys else None))

min_x = min(xs) if xs else None
max_x = max(xs) if xs else None
min_y = min(ys) if ys else None
max_y = max(ys) if ys else None

gap_x = abs(max_x-min_x)
gap_y = abs(max_y-min_y)

pointXIter = resolution[0]/(x_max-x_min)
pointYIter = resolution[1]/(y_max-y_min)

catchBinIterX = resolution[0]/precision
catchBinIterY = resolution[1]/precision

img = Image.new("RGB", (resolution[0], resolution[1]), (0, 0, 0))  # RGB black

x_ranges = [
    (min_x + i * (gap_x / precision), min_x + (i + 1) * (gap_x / precision))
    for i in range(precision)
]

y_ranges = [
    (min_y + i * (gap_y / precision), min_y + (i + 1) * (gap_y / precision))
    for i in range(precision)
]

#print("X ranges: ", x_ranges)
#print("Y ranges: ", y_ranges)

print("Creating spatial index...")
idx = index.Index()
for i, point in enumerate(allPointLocations):
    idx.insert(i, (point.real, point.imag, point.real, point.imag))
    # Print percent done every 1%
    if len(allPointLocations) > 0 and i % (len(allPointLocations) // 100) == 0:
        print(f"Spatial index: {int((i / len(allPointLocations)) * 100)}% done")
print("Done creating spatial index.")

def filter_complex_points(points, x_range, y_range):
    return sum(1 for _ in idx.intersection((x_range[0], y_range[0], x_range[1], y_range[1])))

checks = len(allPointLocations)
full = resolution[0] * resolution[1]
current = 0
print("Total pixels to check: ", full)
print("Total points to check: ", checks)

percentCheckThresholdPixels = (full/100)
percentCheckThresholdChecks = (checks/100)

previousxBin = -1
previousyBin = -1
previousCount = -1

lastNum = -1

locationAndCount = {}
highestCount = -1
for x in range(resolution[0]):
    for y in range(resolution[1]):
        current += 1
        pixel = (x,y)
        
        xBin = math.floor(x / catchBinIterX) + 1
        yBin = math.floor(y / catchBinIterY) + 1

        if xBin == previousxBin and yBin == previousyBin:
            locationAndCount[(x,y)] = count
            continue
        previousxBin = xBin
        previousyBin = yBin

        # Check if the point is in the escaped points

        xToCheck = None
        yToCheck = None
        count = 0
        x_range = x_ranges[xBin-1]
        y_range = y_ranges[yBin-1]
        
        if math.floor((current % (percentCheckThresholdPixels))/10) == 0 and current != full and lastNum != f"{(current/full)*100:.0f}":
            lastNum = f"{(current/full)*100:.0f}"
            print(f"{(current/full)*100:.0f}% done")
        
        count = filter_complex_points(allPointLocations, x_range, y_range)

        brightnessIndex = round(255*(count/checks))
        previousCount = brightnessIndex
        #print(f"Pixel {pixel} has count {count} in bin {(xBin,yBin)} with x range {x_range} and y range {y_range}")
        locationAndCount[(x,y)] = count
        if count > highestCount:
            highestCount = count

print("Done Calculating Bins.")

for (x, y), count in locationAndCount.items():
    brightnessIndex = round(255*(count/highestCount))
    img.putpixel((x, y), (brightnessIndex, brightnessIndex, brightnessIndex))

print("Complete!")
rotated_img = img.transpose(Image.ROTATE_270)
rotated_img.save("buddhabrot.png")