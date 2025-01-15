import arcpy
import pandas as pd
import os
import numpy as np

def get_id_sjo(key):
    sr = arcpy.SpatialReference(4326)
    id = int(key)
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output119.shp")
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
project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)

key = 999
get_id_sjo(key)

folder = os.path.join(project_folder, 'ConvertedSHP')
file = os.path.join(folder, 'Group_999_FID_1996.shp')

cols = ['Id', 'gridcode', 'SHAPE@']
numpy_array = arcpy.da.TableToNumPyArray(file, cols)
df = pd.DataFrame(numpy_array)
print(df.head())

# ----

max_gridcode = df['gridcode'].max()
min_gridcode = df['gridcode'].min()
range_gridcode = round(((max_gridcode - min_gridcode) * 0.35) + min_gridcode)
filtered_df = df[df['gridcode'] <= range_gridcode]
filtered_df = filtered_df.reset_index()
print(filtered_df.head())

# ----

directory = os.path.join(project_folder, 'ConvertedSHP')
file_desc = os.path.join(project_folder, f'SJO/SpatialJoin_Output119.shp')

desc = arcpy.Describe(file_desc)
spatial_reference = desc.spatialReference

# Process each point to find the closest one
start_point = get_id_sjo(key)
min_dist = float('inf')
master_val = None

filtered_df['distance'] = 999

for index, row in filtered_df.iterrows():
    cell = row['SHAPE@']
    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
    projected_centroid = centroid_geom.projectAs(spatial_reference)

    point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
    end_point = arcpy.PointGeometry(point, spatial_reference)

    distance = start_point.distanceTo(end_point)
    filtered_df.at[index, 'distance'] = distance
    #print(f"Distance: {distance} and the id: {row['Id']}, {row['gridcode']}")

# ----

lowest_ele_id = 998

# Sort the DataFrame by 'Rank' (distance) and reset the index
sorted_df = filtered_df.sort_values(by='distance', ascending=True).reset_index(drop=True)

# Find the index of the row where 'Id' == lowest_ele_id
index_of_row_start = sorted_df[sorted_df['Id'] == lowest_ele_id].index[0]

# Get the 5 rows above and 5 rows below using iloc
start_index = max(index_of_row_start - 30, 0)  # Ensure we don't go out of bounds at the start
end_index = min(index_of_row_start + 31, len(sorted_df))  # Ensure we don't go out of bounds at the end

surrounding_rows = sorted_df.iloc[start_index:end_index].reset_index(drop=True)
surrounding_rows['Id_diff'] = abs(surrounding_rows['Id'] - lowest_ele_id)

# Find the row with the maximum difference
max_diff_row = surrounding_rows.loc[surrounding_rows['Id_diff'].idxmax()]

# Extract the 'Id' value from that row
furthest_id = max_diff_row['Id']
index_of_row_end = sorted_df[sorted_df['Id'] == furthest_id].index[0]

end_point = sorted_df.iloc[index_of_row_end]['SHAPE@']
start_point = sorted_df.iloc[index_of_row_start]['SHAPE@']
spoint = extract_coordinates(start_point)
epoint = extract_coordinates(end_point)
data_dict = {}
array = arcpy.Array([spoint, epoint])
line = arcpy.Polyline(array, sr)
data_dict[999] = line

output_fc = os.path.join(gdb,"testline")
#arcpy.CopyFeatures_management(line, output_fc)


id = 999

project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)
gdb = os.path.join(project_folder, 'Default.gdb')
line_fin = os.path.join(gdb, "PolylineLayer")
arcpy.env.workspace = gdb

spatial_reference = arcpy.SpatialReference(4326)

with arcpy.da.InsertCursor(line_fin, ["SHAPE@", "FID"]) as cursor:
    for key, val in data_dict.items():
        try:
            cursor.insertRow([val, key])
            print('Inserted.')
        except Exception as e:
            print(f"Error creating polyline for key {key} with value {val}: {e}")

print('done')
