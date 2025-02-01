import arcpy
import pandas as pd
import os
import numpy as np
import sqlite3
import math

def get_id_sjo(key):
    sr = arcpy.SpatialReference(4326)
    id = int(key)
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{3}.shp")
    cols = ['TARGET_FID', 'SHAPE@', 'FID_ID']
    numpy_array = arcpy.da.TableToNumPyArray(fp, cols)
    df = pd.DataFrame(numpy_array)
    df_filter = df[df['FID_ID'] == id]

    cell = df_filter.iloc[0, 1]
    x = cell.centroid.X
    y = cell.centroid.Y
    point = arcpy.Point(x, y)
    geom = arcpy.PointGeometry(point, sr)
    return geom
  
def extract_coordinates(shape):
    sr = arcpy.SpatialReference(4326)
    projected_shape = shape.projectAs(sr)
    point = arcpy.Point(projected_shape.centroid.X, projected_shape.centroid.Y)
    return point

def get_ids():
    output_file = os.path.join(project_folder, "output_file.db")

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)
    cursor = conn.cursor()

    # Retrieve FID and file paths from the database
    query = """
    SELECT FID, File_path, MainFileID
    FROM recordManager 
    WHERE Processed = (SELECT MAX(Processed) FROM recordManager)
    AND FTP != 'Invalid';
    """

    arr_FID = []
    filepaths = []
    MFID = []
    for row in cursor.execute(query):
        arr_FID.append(row[0])
        filepaths.append(row[1])
        MFID.append(row[2])

    conn.close()
    return arr_FID, MFID

# P! -----------



import time
sr = arcpy.SpatialReference(4326)
arcpy.env.addOutputsToMap = True
count = 1
for filename in os.listdir(directory):
    if filename.endswith('.shp'):
        key = filename.split("_")[1]
        if int(key) in num_arr:
            shapefile_path = os.path.join(directory, filename)
            desc = arcpy.Describe(shapefile_path)
            spatial_reference = desc.spatialReference

            # Load shapefile into a DataFrame
            cols = ['FID', 'Id', 'gridcode', 'SHAPE@']
            numpy_array = arcpy.da.TableToNumPyArray(shapefile_path, cols)
            df = pd.DataFrame(numpy_array)

            max_gridcode = df['gridcode'].max()
            min_gridcode = df['gridcode'].min()
            range_gridcode = round(((max_gridcode - min_gridcode) * 0.45) + min_gridcode)
            #print(f"{range_gridcode}, is range. {min_gridcode} is min")

            filtered_df = df[df['gridcode'] <= range_gridcode]
            filtered_df = filtered_df.reset_index()

            un_arr = filtered_df['gridcode'].unique()
            
            # Create an in-memory feature layer for selection
            arcpy.MakeFeatureLayer_management(shapefile_path, "temp_layer")

            # Construct the SQL query to select rows with gridcode in un_arr
            sql_query = f"gridcode IN ({','.join(map(str, un_arr))})"

            # Select features based on the SQL query
            arcpy.SelectLayerByAttribute_management("temp_layer", "NEW_SELECTION", sql_query)

            # Define the output path for the dissolved shapefile
            output_dissolve_path = os.path.join(gdb, "dissolved_output")

            # Perform the pairwise dissolve on the selected features
            arcpy.analysis.PairwiseDissolve("temp_layer", output_dissolve_path,  multi_part = "SINGLE_PART")
            if count == 1:
                break

# P2

def get_central_points(id_key):
    id_key = int(id_key)
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
    central_point_coords = None

    with arcpy.da.SearchCursor(fp, ["SHAPE@", "FID_ID"]) as cursor:
        for row in cursor:
            if row[1] == id_key:
                point_geom = row[0]  # Geometry from SHAPE@
                point_geom_sr_4326 = point_geom.projectAs(sr_4326)
                central_point_coords = (point_geom_sr_4326.centroid.X, point_geom_sr_4326.centroid.Y)
                break
    return central_point_coords


import math
import arcpy

sr_4326 = arcpy.SpatialReference(4326)

# Define the shapefile path and the coordinates of the central point (red dot)
shapefile_path = r"C:\Script_Arc\Script_Isolated\SJO\SpatialJoin_Output3.shp"

central_point_coords = None  # Initialize as None, will be updated when FID_ID 4 is found

# Iterate through the shapefile and find the row with FID_ID = 4
with arcpy.da.SearchCursor(shapefile_path, ["SHAPE@", "FID_ID"]) as cursor:
    for row in cursor:
        if row[1] == 4:
            # Get the centroid coordinates of the geometry (central point)
            point_geom = row[0]  # Geometry from SHAPE@
            point_geom_sr_4326 = point_geom.projectAs(sr_4326)
            central_point_coords = (point_geom_sr_4326.centroid.X, point_geom_sr_4326.centroid.Y)
            break  # Exit the loop once FID_ID 4 is found

print("Central point coordinates (in SR 4326):", central_point_coords)

# Load and process the dissolved_output shapefile
with arcpy.da.SearchCursor("dissolved_output", ["SHAPE@", "OID@"]) as cursor:
    left_points = []
    right_points = []
    above_points = []
    below_points = []

    for row in cursor:
        point_geom = row[0]
        point_oid = row[1]
        
        # Reproject the geometry to EPSG:4326
        point_geom_sr_4326 = point_geom.projectAs(sr_4326)
        
        # Get the coordinates in EPSG:4326
        point_coords = (point_geom_sr_4326.centroid.X, point_geom_sr_4326.centroid.Y)

        
        # Calculate the distance to the central point, considering both X and Y
        distance = math.sqrt((point_coords[0] - central_point_coords[0])**2 +
                             (point_coords[1] - central_point_coords[1])**2)
        
        # Determine the relative position: left, right, above, or below
        if point_coords[0] < central_point_coords[0]:
            left_points.append((point_oid, distance))
        else:
            right_points.append((point_oid, distance))
        
        if point_coords[1] < central_point_coords[1]:
            below_points.append((point_oid, distance))
        else:
            above_points.append((point_oid, distance))

# Sort the points by distance and select the two closest on each side
left_points = sorted(left_points, key=lambda x: x[1])
right_points = sorted(right_points, key=lambda x: x[1])
above_points = sorted(above_points, key=lambda x: x[1])
below_points = sorted(below_points, key=lambda x: x[1])

# Output the results
print("Two closest points on the left side:", left_points)
print("Two closest points on the right side:", right_points)
print("Two closest points above the central point:", above_points)
print("Two closest points below the central point:", below_points)


