import arcpy
import pandas as pd
import numpy as np
import os

# Function to extract X, Y coordinates
def extract_coordinates(shape):
    # Project the shape to WGS 84 (EPSG 4326) if needed
    projected_shape = shape.projectAs(sr)
    return projected_shape.centroid.X, projected_shape.centroid.Y


# Define the spatial reference
sr = arcpy.SpatialReference(4326)

for map in project.listMaps():
    for layer in map.listLayers():
        print(layer.name)
        if layer.name == "GRIT_Minor_Structures":
            fp2 = layer.dataSource

cols = ['Id', 'SHAPE@']
numpy_array = arcpy.da.TableToNumPyArray(fp2, cols)
df = pd.DataFrame(numpy_array)
sr = arcpy.SpatialReference(4326)

df['X'], df['Y'] = zip(*df['SHAPE@'].apply(extract_coordinates))
df.drop(columns=['SHAPE@'], inplace=True)

project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)

# Define the output Excel file path
output_excel_path = os.path.join(project_folder, "PointLayer.xlsx")

# Save the DataFrame to an Excel file
df.to_excel(output_excel_path, index=False)

print(f"DataFrame successfully exported to {output_excel_path}")

# Define the path to the Excel file
excel_path = os.path.join(project_folder, "PointLayer.xlsx")

# Read the Excel file into a DataFrame
df = pd.read_excel(excel_path)

# Set the geodatabase workspace
gdb = os.path.join(project_folder, 'Default.gdb')
arcpy.env.workspace = gdb

# Define the spatial reference (WGS 84)
sr = arcpy.SpatialReference(4326)

# Define the name for the new feature class
feature_class_name = "PolyPointLayer"

# Create a new PolyPoint feature class
arcpy.CreateFeatureclass_management(
    out_path=gdb,
    out_name=feature_class_name,
    geometry_type="MULTIPOINT",
    spatial_reference=sr
)

# Add the necessary fields to the feature class
arcpy.AddField_management(feature_class_name, "Id", "LONG")

# Insert the data into the new feature class
with arcpy.da.InsertCursor(feature_class_name, ["Id", "SHAPE@"]) as cursor:
    for index, row in df.iterrows():
        # Create a point geometry from X, Y
        point = arcpy.Point(row['X'], row['Y'])
        array = arcpy.Array([point])
        multipoint = arcpy.Multipoint(array, sr)
        
        # Insert the row into the feature class
        cursor.insertRow((row['Id'], multipoint))

print(f"PolyPoint shapefile created at {os.path.join(gdb, feature_class_name)}")


# ------------------------------------

# Paths to your input feature classes
project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)
gdb = os.path.join(project_folder, 'Default.gdb')
poly_point_fc = os.path.join(gdb, "PolyPointLayer_Buffer1")  # PolyPoint feature class
buffer_fc = "Joined_Roads_Buffer"  # The buffer feature class

# Output paths for the spatial join result and the two output shapefiles
spatial_join_output = os.path.join(gdb, "SpatialJoinOutput")
intersecting_output = os.path.join(gdb, "IntersectingPoints")
non_intersecting_output = os.path.join(gdb, "NonIntersectingPoints")

# Step 1: Perform a spatial join to identify intersecting and non-intersecting points
arcpy.analysis.SpatialJoin(
    target_features=poly_point_fc,
    join_features=buffer_fc,
    out_feature_class=spatial_join_output,
    join_type="KEEP_ALL",
    match_option="INTERSECT"
)

# Step 2: Create a layer for the spatial join result
spatial_join_layer = "spatial_join_layer"
arcpy.management.MakeFeatureLayer(spatial_join_output, spatial_join_layer)

# Step 3: Select points that intersect (where Join Count is greater than 0)
arcpy.management.SelectLayerByAttribute(
    spatial_join_layer,
    "NEW_SELECTION",
    '"Join_Count" > 0'
)

# Step 4: Export the intersecting points to a new shapefile
arcpy.management.CopyFeatures(spatial_join_layer, intersecting_output)

# Step 5: Select points that do not intersect (where Join Count is 0)
arcpy.management.SelectLayerByAttribute(
    spatial_join_layer,
    "NEW_SELECTION",
    '"Join_Count" = 0'
)

# Step 6: Export the non-intersecting points to another shapefile
arcpy.management.CopyFeatures(spatial_join_layer, non_intersecting_output)

print(f"Intersecting points saved to {intersecting_output}")
print(f"Non-intersecting points saved to {non_intersecting_output}")

# ------------------------------------------
sr = arcpy.SpatialReference(4326)
cols_df = ['OBJECTID', 'SHAPE@']
cols = ['Id', 'SHAPE@']
nip = os.path.join(gdb, "NonIntersectingPoints")
ip = os.path.join(gdb, "IntersectingPoints")
default = os.path.join(gdb, "PolyPointLayer")

numpy_array = arcpy.da.TableToNumPyArray(nip, cols)
df = pd.DataFrame(numpy_array)
df['X'], df['Y'] = zip(*df['SHAPE@'].apply(extract_coordinates))
# Define the output Excel file path
output_excel_path = os.path.join(project_folder, "NIP_PointLayer.xlsx")

# Save the DataFrame to an Excel file
df.to_excel(output_excel_path, index=False)
print(df.head())

numpy_array = arcpy.da.TableToNumPyArray(ip, cols)
df = pd.DataFrame(numpy_array)
df['X'], df['Y'] = zip(*df['SHAPE@'].apply(extract_coordinates))
# Define the output Excel file path
output_excel_path = os.path.join(project_folder, "IP_PointLayer.xlsx")

# Save the DataFrame to an Excel file
df.to_excel(output_excel_path, index=False)
print(df.head())

default = os.path.join(gdb, "PolyPointLayer")

numpy_array = arcpy.da.TableToNumPyArray(default, cols)
df = pd.DataFrame(numpy_array)
df['X'], df['Y'] = zip(*df['SHAPE@'].apply(extract_coordinates))
output_excel_path = os.path.join(project_folder, "default_PointLayer.xlsx")

# Save the DataFrame to an Excel file
df.to_excel(output_excel_path, index=False)
print(df.head())

default = os.path.join(project_folder, "IP_PointLayer.xlsx")
print('done')

# -------------------------------------------------------------
# Read the Excel file into a DataFrame
output_excel_path_ip = os.path.join(project_folder, "IP_PointLayer.xlsx")
df = pd.read_excel(output_excel_path_ip)

# Extract only the 'Id' column
df_id_only = df[['Id']]

print(df_id_only)

for map in project.listMaps():
    for layer in map.listLayers():
        if layer.name == "GRIT_Minor_Structures":
            fp2 = layer.dataSource
            break

cols = ['Id', 'SHAPE@']
numpy_array = arcpy.da.TableToNumPyArray(fp2, cols)
df = pd.DataFrame(numpy_array)
filtered_df = df[df['Id'].isin(df_id_only['Id'])]
print(filtered_df)

