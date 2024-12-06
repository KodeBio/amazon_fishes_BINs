# -*- coding: utf-8 -*-
# Code for assigning basin information to coordinates based on shapefile data

### ----------------------------
### 1. Import Required Packages
### ----------------------------
import shapefile
from shapely.geometry import Polygon, Point

print("Packages loaded successfully.")

### ----------------------------
### 2. Define Input and Output Files
### ----------------------------

coord_name = 'data/amz_coords.txt'  # Input file with coordinates
out_name = 'data/amz_basin.txt'     # Output file with basin assignments
shf_name = 'data/hybas_sa_lev04_v1c.shp'  # Shapefile containing basin polygons

### ----------------------------
### 3. Read Shapefile Data
### ----------------------------
print(">> Reading shapefile data...")

sf = shapefile.Reader(shf_name)
shapes = sf.shapes()

print("<< Shapefile data loaded successfully.")

### ----------------------------
### 4. Read Coordinate Data
### ----------------------------
print(">> Reading coordinates...")

coords = []
with open(coord_name, 'r') as f:
    next(f)  # Skip the header
    for coord in f:
        # Convert latitude and longitude to a shapely Point object
        point = Point(map(float, str.split(coord)[::-1]))  # Reverse order (lon, lat)
        coords.append(point)

print(f"<< Loaded {len(coords)} coordinates successfully.")

### ----------------------------
### 5. Assign Basins to Coordinates
### ----------------------------
print(">> Assigning basins to coordinates...")

with open(out_name, 'w') as fw:
    fw.write('lat\tlon\tname\n')
    
    for i, point in enumerate(coords):
        c = 0
        name = 'undefined'  # Default value if no match is found
        
        # Check if the point is within any basin polygon
        for j in range(len(shapes)):
            polygon = Polygon(shapes[j].points)
            if polygon.contains(point):
                name = sf.record(j)[0]  # Assign the basin name
                break
            c += 1
        
        # Write the result to the output file
        fw.write('{}\t{}\t{}\n'.format(point.y, point.x, name))

        if (i + 1) % 100 == 0 or i == len(coords) - 1:
            print(f"Processed {i + 1}/{len(coords)} coordinates.")

print(">> Basin assignment completed. Results saved to", out_name)
