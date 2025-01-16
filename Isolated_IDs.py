import arcpy
import pandas as pd
import os
import numpy as np
import sqlite3

def get_id_sjo(key):
    sr = arcpy.SpatialReference(4326)
    id = int(key)
    sfj = 14090
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{1}.shp")
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

# Start Here
project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)
gdb = os.path.join(project_folder, 'Default.gdb')
line_fin = os.path.join(gdb, "PolylineLayer")
arcpy.env.workspace = gdb

sr = arcpy.SpatialReference(4326)
project = arcpy.mp.ArcGISProject("CURRENT")
project_folder = os.path.dirname(project.filePath)
directory = os.path.join(project_folder, 'ConvertedSHP')
spatial_reference = 4326
num_arr, mainfileFID_arr = get_ids()
mainfileFID_arr = sorted(mainfileFID_arr)

for filename in os.listdir(directory):
    if filename.endswith('.shp'):
        key = filename.split("_")[1]
        if int(key) in num_arr:
            print(key)
            shapefile_path = os.path.join(directory, filename)
            desc = arcpy.Describe(shapefile_path)
            spatial_reference = desc.spatialReference
            
            # Load shapefile into a DataFrame
            cols = ['FID', 'Id', 'gridcode', 'SHAPE@']
            numpy_array = arcpy.da.TableToNumPyArray(shapefile_path, cols)
            df = pd.DataFrame(numpy_array)

            max_gridcode = df['gridcode'].max()
            min_gridcode = df['gridcode'].min()
            range_gridcode = round(((max_gridcode - min_gridcode) * 0.35) + min_gridcode)
            filtered_df = df[df['gridcode'] <= range_gridcode]
            filtered_df = filtered_df.reset_index()

            start_point = get_id_sjo(key)
            filtered_df['distance'] = 999

            for index, row in filtered_df.iterrows():
                cell = row['SHAPE@']
                centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
                projected_centroid = centroid_geom.projectAs(spatial_reference)

                point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
                end_point = arcpy.PointGeometry(point, spatial_reference)

                distance = start_point.distanceTo(end_point)
                filtered_df.at[index, 'distance'] = distance
                
                # ----
                
            # First, find the row with the minimum distance
            lowest_gridcode_df = filtered_df.nsmallest(5, 'gridcode')
            fid_of_min_value = lowest_gridcode_df.at[lowest_gridcode_df['distance'].idxmin(), 'Id']

            # ---- Conert to 25 lowest and closest to Point

            lowest_ele_id = fid_of_min_value

            # Sort the DataFrame by 'Rank' (distance) and reset the index
            sorted_df = filtered_df.sort_values(by='distance', ascending=True).reset_index(drop=True)

            # Find the index of the row where 'Id' == lowest_ele_id
            index_of_row_start = sorted_df[sorted_df['Id'] == lowest_ele_id].index[0]
  
            # Get the 5 rows above and 5 rows below using iloc
            start_index = max(index_of_row_start - 30, 0)  # Ensure we don't go out of bounds at the start
            end_index = min(index_of_row_start + 31, len(sorted_df))  # Ensure we don't go out of bounds at the end

            surrounding_rows = sorted_df.iloc[start_index:end_index].reset_index(drop=True)
            surrounding_rows['Id_diff'] = abs(surrounding_rows['Id'] - lowest_ele_id).reset_index(drop = True)

            # Find the row with the maximum difference
            max_diff_row = surrounding_rows.loc[surrounding_rows['Id_diff'].idxmax()]

            # Extract the 'Id' value from that row
            furthest_id = max_diff_row['Id']
            index_of_row_end = sorted_df[sorted_df['Id'] == furthest_id].index[0]

            # ----

            end_point = sorted_df.iloc[index_of_row_end]['SHAPE@']
            start_point = sorted_df.iloc[index_of_row_start]['SHAPE@']
            spoint = extract_coordinates(start_point)
            epoint = extract_coordinates(end_point)

            # ----
            data_dict = {}
            array = arcpy.Array([spoint, epoint])
            line = arcpy.Polyline(array, sr)
            data_dict[lowest_ele_id] = line

            output_fc = os.path.join(gdb,"testline")
            arcpy.CopyFeatures_management(line, output_fc)
            # ----
            
            data_dict = {}
            array = arcpy.Array([spoint, epoint])
            line = arcpy.Polyline(array, sr)
            data_dict[lowest_ele_id] = line

            with arcpy.da.InsertCursor(line_fin, ["SHAPE@", "FID"]) as cursor:
                for key, val in data_dict.items():
                    try:
                        cursor.insertRow([val, key])
                        print('Inserted.')
                    except Exception as e:
                        print(f"Error creating polyline for key {key} with value {val}: {e}")
            
            # ----

            
print(sorted_df)     
print('done')             



"""
start_point = sorted_df.iloc[index_of_row_start]['SHAPE@']
centroid_geom = arcpy.PointGeometry(start_point.centroid, start_point.spatialReference)
projected_centroid = centroid_geom.projectAs(spatial_reference)
point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
sp_point = arcpy.PointGeometry(point, spatial_reference)
spoint = extract_coordinates(start_point)

for index, row in surrounding_rows.iterrows():
    cell = row['SHAPE@']
    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
    projected_centroid = centroid_geom.projectAs(spatial_reference)

    point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
    end_point = arcpy.PointGeometry(point, spatial_reference)

    distance = sp_point.distanceTo(end_point)
    surrounding_rows.at[index, 'distance_end_point'] = distance

print(surrounding_rows[['Id', 'gridcode', 'distance', 'distance_end_point']])

# Sort by 'distance to the point centered around' (Column A) in ascending order
sorted_df = surrounding_rows.sort_values(by='distance', ascending=True)

# Further sort by 'distance to the point looking at' (Column B) in descending order
sorted_df = surrounding_rows.sort_values(by='distance_end_point', ascending=False)

print(sorted_df[['Id', 'gridcode', 'distance', 'distance_end_point']])

"""

"""
data_dict = {}
array = arcpy.Array([spoint, epoint])
line = arcpy.Polyline(array, sr)
data_dict[lowest_ele_id] = line

output_fc = os.path.join(gdb,"testline")
arcpy.CopyFeatures_management(line, output_fc)
"""

def isolated_points():
    sr = arcpy.SpatialReference(4326)
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
                print(f"{range_gridcode}, is range. {min_gridcode} is min")
                
                filtered_df = df[df['gridcode'] <= range_gridcode]
                filtered_df = filtered_df.reset_index()

                start_point = get_id_sjo(key)
                filtered_df['distance'] = 999

                for index, row in filtered_df.iterrows():
                    cell = row['SHAPE@']
                    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
                    projected_centroid = centroid_geom.projectAs(spatial_reference)

                    point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
                    end_point = arcpy.PointGeometry(point, spatial_reference)

                    distance = start_point.distanceTo(end_point)
                    filtered_df.at[index, 'distance'] = distance

                    # ----

                # First, find the row with the minimum distance
                lowest_gridcode_df = filtered_df.nsmallest(5, 'gridcode')
                fid_of_min_value = lowest_gridcode_df.at[lowest_gridcode_df['distance'].idxmin(), 'Id']

                # ---- Conert to 25 lowest and closest to Point

                lowest_ele_id = fid_of_min_value

                # Sort the DataFrame by 'Rank' (distance) and reset the index
                sorted_df = filtered_df.sort_values(by='distance', ascending=True).reset_index(drop=True)

                # Find the index of the row where 'Id' == lowest_ele_id
                index_of_row_start = sorted_df[sorted_df['Id'] == lowest_ele_id].index[0]

                # Get the 5 rows above and 5 rows below using iloc
                start_index = max(index_of_row_start - 40, 0)  # Ensure we don't go out of bounds at the start
                end_index = min(index_of_row_start + 41, len(sorted_df))  # Ensure we don't go out of bounds at the end

                surrounding_rows = sorted_df.iloc[start_index:end_index].reset_index(drop=True)
                surrounding_rows['Id_diff'] = abs(surrounding_rows['Id'] - lowest_ele_id).reset_index(drop = True)

                # Find the row with the maximum difference
                max_diff_row = surrounding_rows.loc[surrounding_rows['Id_diff'].idxmax()]

                # Extract the 'Id' value from that row
                furthest_id = max_diff_row['Id']
                index_of_row_end = sorted_df[sorted_df['Id'] == furthest_id].index[0]

                # ----

                end_point = sorted_df.iloc[index_of_row_end]['SHAPE@']
                start_point = sorted_df.iloc[index_of_row_start]['SHAPE@']
                spoint = extract_coordinates(start_point)
                epoint = extract_coordinates(end_point)

                # ----
                data_dict = {}
                array = arcpy.Array([spoint, epoint])
                line = arcpy.Polyline(array, sr)
                data_dict[lowest_ele_id] = line

                output_fc = os.path.join(gdb,"testline")
                arcpy.CopyFeatures_management(line, output_fc)
                # ----

                data_dict = {}
                array = arcpy.Array([spoint, epoint])
                line = arcpy.Polyline(array, sr)
                data_dict[lowest_ele_id] = line
                line_fin = os.path.join(gdb, "PolylineLayer")
                with arcpy.da.InsertCursor(line_fin, ["SHAPE@", "FID"]) as cursor:
                    for key, val in data_dict.items():
                        try:
                            cursor.insertRow([val, key])
                        except Exception as e:
                            print(f"Error creating polyline for key {key} with value {val}: {e}")

                # ----
