import os
import shutil
import arcpy
from arcpy import env
from arcpy.sa import *
import sqlite3
import pandas as pd
import urllib.request
from arcpy.sa import ExtractByMask
import time
import random
import sys
import io
from datetime import datetime
from pathlib import Path
import zipfile
import requests

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Start timing
        result = func(*args, **kwargs)  # Call the original function
        end_time = time.time()  # End timing
        duration = end_time - start_time  # Calculate how long it took
        print(f"Function '{func.__name__}' took {duration:.4f} seconds to run.")
        return result  # Return the result of the original function
    return wrapper


def check_gdb():
    """
    Creates File Geodatabase
    
    """
    check_file = os.path.join(project_folder, "Hold_DEMs")
    make_file = os.path.join(project_folder, "Hillshade_Converted")
    if not os.path.exists(check_file):
        os.makedirs(check_file)
    if not os.path.exists(make_file):
        os.makedirs(make_file)
    return check_file

def check_sqlite():
    """
    Creates SQLite File
    
    """
    # Define the path for the SQLite database file
    output_file = os.path.join(project_folder, "output_file.db")

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)

    cursor = conn.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS recordManager (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        FID INTEGER,
        File_name TEXT,
        File_path TEXT,
        FTP TEXT,
        Processed INTEGER,
        MainFileID INTEGER,
        X_LOCATION REAL,
        Y_LOCATION REAL,
        DateTimeAdded TEXT
    )
    """)
    conn.close()
    return    
    
def get_sp():
    """
    Get File ID's from SQLite Database
    
    """
    check_sqlite()
    db = os.path.join(project_folder, "Default.gdb")
    PL = os.path.join(db, "PointLayer")

    if os.path.exists(PL):
        os.remove(PL)

    # Define the path for the SQLite database file
    output_file = os.path.join(project_folder, "output_file.db")

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)
    cursor = conn.cursor()
    cursor.execute("SELECT max(FID), max(Processed) FROM recordManager")
    result = cursor.fetchone()
    
    # Check for results and update max_fid accordingly
    if result is None or result[0] is None:
        max_fid = 1
    else:
        max_fid = int(result[0]) + 1

    if result is None or result[1] is None:
        max_processed = 0  
    else:
        max_processed = int(result[1])  

    # Close the database connection
    conn.close()
    return  [max_fid, max_processed]


def holderFolder():
    """
    Create files for holding geospatial data
    
    """
    # Create the join folder
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    folder = "HolderFolder"
    dest = os.path.join(project_folder, folder)
    if not os.path.exists(dest):
        os.makedirs(dest)
    
    # Create folder to house the geodatabase
    lidar_gdb_path = os.path.join(project_folder, 'lidar.gdb')
    if not arcpy.Exists(lidar_gdb_path):
        arcpy.CreateFileGDB_management(project_folder, 'lidar.gdb')

    # URL to the raw ZIP file
    url = "https://github.com/MTS-DWSC/ND_Lidar_DEM/raw/main/ND_Index.zip"
    check_file = os.path.join(project_folder, 'ND_Index.shp')
    
    if not os.path.exists(check_file):
        response = requests.get(url, verify=False)
        if response.status_code == 200:
            # Extract the ZIP file
            with zipfile.ZipFile(io.BytesIO(response.content)) as z:
                z.extractall(project_folder)  
        else:
            print(f"Failed to download the file. Status code: {response.status_code}")
    else:
        pass

    # Grab data source from API
    grit_layer = os.path.join(lidar_gdb_path, "GRIT_Minor_Structures")
    if arcpy.Exists(grit_layer):
        return
    else:
        # URL of the REST service
        service_url = "https://dotsc.ugpti.ndsu.nodak.edu:6443/arcgis/rest/services/GRIT_all/grit20_bridges_all_feature/MapServer/0"
        
        # Create a temporary layer from the REST service
        temp_layer = "in_memory/GRIT_Minor_Structures"
        #arcpy.MakeFeatureLayer_management(service_url, temp_layer)
        #print("Feature layer created from the REST service.")

        # Optionally, save the layer to the geodatabase
        arcpy.CopyFeatures_management(temp_layer, grit_layer)
        print("Feature layer copied to the geodatabase.")
    return

def delete_sj():
    """
    Cleanup Files that are no longer needed
    Delete everything but in use lockfiles.
    
    """
    output_folder = os.path.join(project_folder, "SJO")
    shape_folder = os.path.join(project_folder, "ConvertedSHP")
    holder_folder = os.path.join(project_folder, "HolderFolder")
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)

        if filename.endswith(".lock"):
            continue 
        else:
            os.remove(file_path) 
            
    for filename in os.listdir(holder_folder):
        file_path = os.path.join(holder_folder, filename)

        if filename.endswith(".lock"):
            continue  
        else:
            os.remove(file_path) 
            
    for filename in os.listdir(shape_folder):
        file_path = os.path.join(shape_folder, filename)

        if filename.endswith(".lock"):
            continue  
        else:
            os.remove(file_path)
        


@time_it
def spatial_join():
    """
    Function Name: perform_spatial_join_and_export
    
    Purpose:
    This function performs a spatial join between a target feature layer and a join feature layer, and exports the results to a specified output location. 
    If the spatial join results in an empty feature class, the function exports the relevant attributes from the target feature layer to a SQLite database.
    
    Outputs:
    - A feature class resulting from the spatial join, saved in the SJO folder within the project folder.
    - If the spatial join results in an empty feature class, a SQLite database (output_file.db) 
      is created or updated with the relevant attributes from the target feature layer.
    """
    
    db = os.path.join(project_folder, "Default.gdb")
    target_features = os.path.join(db, "PointLayer")
    
    join_features = project_folder + "\\" + "ND_Index.shp"
    
    output_folder = os.path.join(project_folder, "SJO")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_features = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
       
    field_mappings = arcpy.FieldMappings()

    # Add the FID field from Extract_temp
    fid_field_map = arcpy.FieldMap()
    fid_field_map.addInputField(target_features, "FID_ID")
    field_mappings.addFieldMap(fid_field_map)
    
    # Add the File ID field from Extract_temp
    fid_field_map = arcpy.FieldMap()
    fid_field_map.addInputField(target_features, "MainFileID") 
    field_mappings.addFieldMap(fid_field_map)
    
    fid_field_map = arcpy.FieldMap()
    fid_field_map.addInputField(target_features, "PL_Process")
    field_mappings.addFieldMap(fid_field_map)

    # Add the Tilename field from IndexLidar
    tilename_field_map = arcpy.FieldMap()
    tilename_field_map.addInputField(join_features, "Tilename")
    field_mappings.addFieldMap(tilename_field_map)

    # Add the ASCII_Path field from IndexLidar
    ascii_path_field_map = arcpy.FieldMap()
    ascii_path_field_map.addInputField(join_features, "ASCII_Path")
    field_mappings.addFieldMap(ascii_path_field_map)
    
    # Make a layer from target_features
    arcpy.management.MakeFeatureLayer(target_features, "target_layer")

    # Select points in target_features that are within join_features
    arcpy.management.SelectLayerByLocation(target_features, "WITHIN", join_features)

    # Perform the spatial join
    arcpy.analysis.SpatialJoin(
        target_features="target_layer",
        join_features=join_features,
        out_feature_class=output_features,
        join_type="KEEP_COMMON",  
        field_mapping=field_mappings
    )
    
    # Check if the output_features is empty
    count_result = arcpy.management.GetCount(output_features)
    count = int(count_result.getOutput(0))
    
    if count == 0:
        rows_list = []
        with arcpy.da.SearchCursor("target_layer", ["FID_ID", "MainFileID", "PL_Process"]) as cursor:
            for row in cursor:
                rows_list.append(row)
        df = pd.DataFrame(rows_list, columns=["FID_ID", "MainFileID", "PL_Process"])
        
        output_file = os.path.join(project_folder, "output_file.db")
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()
        
        print(f"printing df: {df}")
        for index, row in df.iterrows():
            #print(f"printing row: {row}")
            cursor.execute("""
                INSERT INTO recordManager (FID, File_name, File_path, FTP, Processed, MainFileID, DateTimeAdded)
                VALUES (?, ?, ?, ?, ?,?, ?)
            """, (int(row[0]) - 0, "Out of Bounds", "Out of Bounds", "Out of Bounds", int(row[2]), int(row[1]), current_timestamp))
        conn.commit()
        cursor.close()
        sys.exit()

    
    


def get_df(text):
    """
    Creates a Pandas dataframe 
    
    """
    project = arcpy.mp.ArcGISProject("CURRENT")
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")

    rows = []
    switch = 1
    if switch == 1:
        cols1 = ['FID_ID', 'TARGET_FID', 'TILENAME', 'ASCII_PATH', 'PL_Process', 'MainFileID', 'SHAPE@']
        with arcpy.da.SearchCursor(fp, cols1) as cursor:
            for row in cursor:
                shape = row[-1]
                x, y = shape.centroid.X, shape.centroid.Y
                rows.append(list(row[:-1]) + [x, y])
        cols = ['FID_ID', 'TARGET_FID', 'TILENAME', 'ASCII_PATH', 'PL_Process', 'MainFileID', 'X_LOCATION', 'Y_LOCATION']
        df = pd.DataFrame(rows, columns=['FID_ID', 'TARGET_FID', 'TILENAME', 'ASCII_PATH', 'PL_Process', 'MainFileID', 'X_LOCATION', 'Y_LOCATION'])
        return df



@time_it
def populate_folder(df_ret): 
    """
    Processes and manages hillshade raster files associated with spatial datasets. 
    It interacts with an SQLite database to track processed files and downloads missing ones from an FTP server for the Department of Water Resources.

    Key Features
    - File Tracking: Ensures missing files are identified and downloaded.
    - Database Logging: Records file details and status in a SQLite database.
    - FTP Download Automation: Fetches files from a remote server if missing.
    - Error Handling: Catches exceptions and logs invalid records.
    """
    
    destination_folder = check_gdb()
    hs_folder = os.path.join(project_folder, 'Hillshade_Converted')
    check_file = project_folder + "\\" + "output_file.db"
    if os.path.exists(check_file):
        output_file = os.path.join(project_folder, "output_file.db")
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()

        arr = df_ret['FID_ID']
        results_arr = []
        count = 0
        cursor.execute("SELECT File_name FROM recordManager")
        results = cursor.fetchall()
        for row in results:
            results_arr.append(row[0])
        missing_fids = [fid for fid in arr if fid not in results_arr]
    
        output_file = os.path.join(project_folder, "output_file.db")
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()

        for x in missing_fids:
            # You need the -1 for the IDs from the original file.
            try:
                count += 1
                lookup = df_ret[df_ret['FID_ID'] == x]
                check = 'dem' + lookup['TILENAME'].iloc[0] + '_Hillshade.tif'
                check_fp = os.path.join(hs_folder, check)
                fn = hs_folder + '\\' + lookup['ASCII_PATH'].iloc[0].split('/')[-1]
                fn_final = fn.replace('.tif', '_Hillshade.tif')
                res = 'ftp://swc:water@lidarftp.swc.nd.gov/' + lookup['ASCII_PATH'].iloc[0]
                file_name = os.path.basename(res)
                destination_path = os.path.join(destination_folder, file_name)
                if os.path.exists(check_fp):
                    print(f"{file_name} already exists. Skipping...")
                else:
                    urllib.request.urlretrieve(res, destination_path)
                    if count % 5 == 0:
                        print(f"{file_name} processed...")
                cursor.execute("""
                    INSERT INTO recordManager (FID, File_name, File_path, FTP, Processed, MainFileID, X_LOCATION, Y_LOCATION, DateTimeAdded)
                    VALUES (?, ?, ?, ?, ?,?, ?, ?, ?)
                """, (int(lookup['FID_ID'].iloc[0]) - 0, lookup['TILENAME'].iloc[0], fn_final, res, int(lookup['PL_Process'].iloc[0]), int(lookup['MainFileID'].iloc[0]), lookup['X_LOCATION'].iloc[0], lookup['Y_LOCATION'].iloc[0], current_timestamp))
            except Exception as e:
                cursor.execute("""
                    INSERT INTO recordManager (FID, File_name, File_path, FTP, Processed, MainFileID, X_LOCATION, Y_LOCATION, DateTimeAdded)
                    VALUES (?, ?, ?, ?, ?,?, ?, ?, ?)
                """, (int(lookup['FID_ID'].iloc[0]) - 0, lookup['TILENAME'].iloc[0], "Invalid", "Invalid", int(lookup['PL_Process'].iloc[0]), int(lookup['MainFileID'].iloc[0]), lookup['X_LOCATION'].iloc[0], lookup['Y_LOCATION'].iloc[0], current_timestamp))
                print(f"An unexpected error occurred for {file_name}: {e}. Skipping...")
                continue
        conn.commit()
        conn.close()

@time_it
def convert_crs(name):
    """
    Description:
    This function processes raster files in the 'Hold_DEMs' directory, applies a hillshade transformation, and saves the results in the 'Hillshade_Converted' folder. The function utilizes ArcPy to check the spatial reference and process only raster datasets.
    
    Parameters:
    - name (str): The input name used to retrieve spatial reference details.
    
    Process:
    1. Determines the geodatabase folder using `check_gdb()`.
    2. Defines the output directory for converted hillshade files (`Hillshade_Converted`).
    3. Iterates through all raster files in the 'Hold_DEMs' directory.
    4. Verifies that the file is a raster dataset using `arcpy.Describe()`.
    5. Applies a hillshade transformation using `arcpy.Raster()` and `Hillshade()`.
    6. Saves the output with a `_Hillshade.tif` suffix.
    7. Skips processing if the output file already exists.
    8. Prints progress updates for every 5 processed files.
    9. Handles and logs errors encountered during processing.
    """
    
    fp = check_gdb()
    fp_converted = Path(fp).parent / 'Hillshade_Converted'

    spatial_ref = arcpy.Describe(name).spatialReference
    count = 0

    new_fp = Path(project_folder) / 'Hold_DEMs'

    for file in new_fp.iterdir():
        if file.suffix.lower() in ['.tif', '.img', '.jpg', '.png']:  
            try:
                if arcpy.Describe(str(file)).datatype == 'RasterDataset':
                    output_name = f"{file.stem}_Hillshade.tif"
                    output_path = fp_converted / output_name
                    if output_path.exists():
                        continue
                        
                    in_raster = arcpy.Raster(str(file))
                    out_hillshade = Hillshade(in_raster)
                    
                    out_hillshade.save(str(output_path))
                    count += 1
                    
                    if count % 5 == 0:
                        print(f"Processing complete for {output_name}")

            except Exception as e:
                print(f"Error processing {file.name}: {e}")
    return

@time_it
def process_geospatial_data(buffer_distance="35 Meters"):
    """
    Function: process_geospatial_data
    
    Description:
    This function processes geospatial data by filtering features based on an SQLite database, applying spatial operations, and converting multipart features to single-part.
    
    Steps:
    1. Sets the workspace and connects to an SQLite database.
    2. Retrieves valid FID values from the 'recordManager' table.
    3. Creates a where clause to filter features in 'SpatialJoin_Output'.
    4. Buffers the selected features by a specified distance.
    5. Performs an erase operation using the 'Joined_Roads_Buffer' layer.
    6. Converts multipart features to single-part with retries on failure.
    
    Parameters:
    - buffer_distance (str, optional): The buffer distance (default: "35 Meters").
    """
    
    # Access the current project and workspace
    attempt, max_tries = 0, 3
    path = os.path.join(project_folder, "HolderFolder")
    gdb_path = os.path.join(project_folder, "default.gdb")

    # Ensure workspace is set
    arcpy.env.workspace = project_folder

    # SQLite database connection
    output_file = os.path.join(project_folder, "output_file.db")
    conn = sqlite3.connect(output_file)
    cursor = conn.cursor()

    # Retrieve file paths from SQLite
    arr_FID = []
    query = "SELECT FID FROM recordManager WHERE FTP != 'Invalid';"
    for row in cursor.execute(query):
        arr_FID.append(row[0])
    conn.close()

    # Generate where clause for FID filtering
    where_clause = "TARGET_FID IN ("
    for x in arr_FID:
        where_clause = where_clause + str(x) + ','
    where_clause = where_clause[:-1] + ')'

    # Locate the "SpatialJoin_Output" layer
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")

    if not fp:
        print("Error: 'SpatialJoin_Output' layer not found.")
        return

    # Create a feature layer with the where clause
    arcpy.MakeFeatureLayer_management(fp, "Extract_temp", where_clause)

    # Buffer the selected features
    buffered_output = os.path.join(path, "buffered_output")
    arcpy.analysis.Buffer("Extract_temp", buffered_output, buffer_distance)

    # Locate the "Roads_Line_Buffer" feature class
    roads_fc = os.path.join(gdb_path, "Joined_Roads_Buffer")

    # Perform erase operation
    erased_output = os.path.join(path, "erased_output")
    arcpy.analysis.Erase(buffered_output, roads_fc, erased_output)
    
    time.sleep(1)
    
    # Convert multipart features to single part
    singlepart_output = os.path.join(path, f"singlepart_output{singlepart_rand}.shp")
    while attempt < max_tries:
        try:
            # Attempt the operation
            arcpy.management.MultipartToSinglepart(erased_output, singlepart_output)
            break  
        except arcpy.ExecuteError as e:
            print(f"Error: {e}")
            attempt += 1
            if attempt < max_retries:
                print(f"Retrying... (Attempt {attempt} of {max_tries})")
                time.sleep(2)  # Wait for 2 seconds before retrying
            else:
                print(f"Failed after {max_tries} attempts.")

@time_it
def correct_output():
    """
    Function: process_spatial_ranking
    
    Description:
    This function analyzes a shapefile, calculates polygon areas, determines centroid coordinates in WGS 84, ranks features by area, and removes lower-ranked polygons.
    
    Steps:
    1. Reads polygon features from 'singlepart_output' and calculates their area.
    2. Projects centroids to WGS 84 and extracts coordinates.
    3. Stores results in a DataFrame and ranks polygons within each 'TARGET_FID'.
    4. Adds a 'RANK' field to the feature class if missing.
    5. Updates rank values and deletes polygons ranked lower than 2.
        path = os.path.join(project_folder, "HolderFolder")
        feature_class = os.path.join(path, f"singlepart_output{singlepart_rand}.shp")
    """

    geo_sr = arcpy.SpatialReference(4326)

    results = []

    with arcpy.da.SearchCursor(feature_class, ['FID', 'TARGET_FID', 'TILENAME', 'SHAPE@']) as cursor:
        for row in cursor:
            polygon_geometry = row[3]  # Access the geometry
            area = polygon_geometry.area  # Get the area of the polygon

            centroid_geom = arcpy.PointGeometry(polygon_geometry.centroid, polygon_geometry.spatialReference)
            # Project the centroid to WGS 84
            projected_centroid = centroid_geom.projectAs(geo_sr)

            x = projected_centroid.centroid.X  # Longitude
            y = projected_centroid.centroid.Y  # Latitude
            results.append({
                "FID": row[0],       
                "TARGET_FID": row[1],   
                "TILENAME": row[2],      
                "Polygon_Area": area,
                "X": x,
                "Y": y,
                "MATCH": ""
            })

    df = pd.DataFrame(results)

    df['RANK'] = df.groupby('TARGET_FID')['Polygon_Area'].rank(ascending=False)

    rank_field = "RANK"

    if rank_field not in [f.name for f in arcpy.ListFields(feature_class)]:
        arcpy.AddField_management(feature_class, rank_field, "DOUBLE")  # Add the Rank field


    with arcpy.da.UpdateCursor(feature_class, ['FID', rank_field]) as cursor:
        for row in cursor:
            fid = row[0]
            # Look up the Rank value from the DataFrame
            rank_value = df.loc[df['FID'] == fid, 'RANK'].values[0]
            row[1] = rank_value
            if rank_value > 2:
                cursor.deleteRow()
            else:
                cursor.updateRow(row)


def group_tracker(total_groups):
    """
    Method of creating unique ID's on Raster Masks.
    
    """
    starting_number = 0
    result = []

    # Loop through each group and calculate the group numbers
    for group in range(1, total_groups + 1):
        # Store the numbers for each group in the result list
        result.append([starting_number, starting_number + 1])
        # Increment the starting number for the next group
        starting_number += 2
        
    return result

def get_ids():
    """
    Gets last processed batch of ID points.
    
    """
    
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


def get_FIDS(shapefile_path, number):
    max_FID = float('-inf')
    min_FID = float('inf')
    
    
    with arcpy.da.SearchCursor(shapefile_path, ['FID_ID', 'FID']) as cursor:
        for row in cursor:
            if row[0] == number:
                if row[1] > max_FID:
                    max_FID = row[1]
                if row[1] < min_FID:
                    min_FID = row[1]

    
    return max_FID, min_FID

@time_it
def mask_extraction():
    """
    Function: mask_extraction
    
    Description:
    Extracts raster data using polygon masks based on the latest processed records in an SQLite database.
    
    Steps:
    1. Connects to the SQLite database and retrieves the latest processed file paths.
    2. Sorts and iterates through valid raster files.
    3. Creates a temporary raster layer for processing.
    4. Uses a grouping function to determine polygon masks.
    5. Generates feature layers for selected polygons.
    6. Performs mask extraction on the raster using these layers.
    7. Saves extracted rasters to a geodatabase.
    
    Returns:
    - List of processed FIDs.
    """
    
    folder = "HolderFolder"
    dest = os.path.join(project_folder, folder)
    so = os.path.join(dest, f'singlepart_output{singlepart_rand}.shp')
    gdb_path = os.path.join(project_folder, "lidar.gdb")

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)
    cursor = conn.cursor()

    # Retrieve FID and file paths from the database
    query = """
    SELECT FID, File_path 
    FROM recordManager 
    WHERE Processed = (SELECT MAX(Processed) FROM recordManager)
    AND FTP != 'Invalid';
    """
    filepaths = {}
    for row in cursor.execute(query):
        filepaths[row[0]] = row[1]
    conn.close()

    sorted_filepaths = dict(sorted(filepaths.items()))
    count = 0
    for key, val in sorted_filepaths.items():
        if arcpy.Exists(val):

            count = int(key)
            filepath = val

            # Create a raster layer
            raster_layer_result = arcpy.MakeRasterLayer_management(filepath, f"temp_raster_layer_{count}")
            raster_layer = raster_layer_result.getOutput(0)
            data_source = raster_layer

            # Group tracker logic (dummy example, replace with actual implementation)
            arr = group_tracker(count)
            partA = arr[-1][0]
            partB = arr[-1][1]

            filtered_layerA = f"Group_{count}_FID_{partA}"
            filtered_layerB = f"Group_{count}_FID_{partB}"

            ids = get_FIDS(so, count)

            # Create feature layers
            arcpy.MakeFeatureLayer_management(so, filtered_layerA, f"FID = {ids[0]}")
            arcpy.MakeFeatureLayer_management(so, filtered_layerB, f"FID = {ids[1]}")

            # Perform mask extraction
            mask_1 = ExtractByMask(data_source, filtered_layerA, "INSIDE")
            mask_2 = ExtractByMask(data_source, filtered_layerB, "INSIDE")
            subfolder_path = os.path.join(project_folder, "lidar.gdb")

            mask_1.save(os.path.join(subfolder_path, filtered_layerA))
            mask_2.save(os.path.join(subfolder_path, filtered_layerB))
            count += 1
            if count % 5 == 0:
                print(f"Mask Extraction: {filtered_layerA} and {filtered_layerB}.")
            
    keys = sorted_filepaths.keys()
    keys_list = list(keys)
    return keys_list

@time_it
def convert_shp(keys_list):
    """
    Function: convert_shp
    
    Description:
    Converts raster datasets from a geodatabase to shapefiles if their names match provided keys.
    
    Steps:
    1. Ensures the output folder exists.
    2. Sets the geodatabase as the workspace and retrieves raster names.
    3. Iterates through rasters, checking for key matches.
    4. Converts matching rasters to polygons and saves them as shapefiles.
    5. Handles errors if conversion fails.
    
    Returns:
    - None (shapefiles are saved in the output folder).
    """

    gdb_path = os.path.join(project_folder, 'lidar.gdb')
    output_folder = os.path.join(project_folder, 'ConvertedSHP')

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Set the workspace to the geodatabase
    arcpy.env.workspace = gdb_path

    # List all rasters in the geodatabase
    rasters = arcpy.ListRasters()

    if rasters:
        for raster in rasters:
            if any(str(key) in raster for key in keys_list):
                try:
                    valid_raster_name = arcpy.ValidateTableName(raster, output_folder)
                    in_raster = os.path.join(gdb_path, raster)
                    output_shp = os.path.join(output_folder, f"{valid_raster_name}.shp")

                    arcpy.RasterToPolygon_conversion(
                        in_raster=in_raster,
                        out_polygon_features=output_shp,
                        simplify="NO_SIMPLIFY",
                        raster_field="Value"
                    )
                except Exception as e:
                    print(f"Error converting {raster}: {e}")
            else:
                continue
    else:
        print("No rasters found in the geodatabase.")


def get_id_sjo(key):
    """
    Gets PointGeometry object for input ID point.
    
    """
    
    sr = arcpy.SpatialReference(4326)
    id = int(key)
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
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

@time_it
def cleanup(keys_list) -> None:
    """
    Cleans up files that are no longer needed
    """
    ConvertedSHP = os.path.join(project_folder, 'ConvertedSHP')
    Hold_DEMs = os.path.join(project_folder, 'Hold_DEMs')
    HolderFolder = os.path.join(project_folder, 'HolderFolder')
    #shutil.rmtree(ConvertedSHP)
    shutil.rmtree(Hold_DEMs)
    delete_sj()
    
    lidar_gdb = os.path.join(project_folder, 'lidar.gdb')
    arcpy.env.workspace = lidar_gdb
    datasets = arcpy.ListRasters()
    itter = 0
    for data in datasets:
        if any(str(key) in data for key in keys_list):
            continue
        else:
            fp = os.path.join(lidar_gdb, data)
            arcpy.Delete_management(fp)
            itter += 1
            
            
def clean_hillshades():
    """
    Creates zipped file of Hillshades for later use.
    
    """
    hillshade_Con = os.path.join(project_folder, "Hillshade_Converted")
    archive_folder = os.path.join(project_folder, "Archive")

    # Get list of .tif files in hillshade_Con
    tif_files = [f for f in os.listdir(hillshade_Con) if f.endswith('.tif')]

    if len(tif_files) > 50:
        # Get existing folders in archive_folder
        existing_folders = [f for f in os.listdir(archive_folder) if os.path.isdir(os.path.join(archive_folder, f))]

        # Extract indices from folders with the format hillshade_holder_X
        indices = [int(f.split('_')[-1]) for f in existing_folders if f.startswith('hillshade_holder_') and f.split('_')[-1].isdigit()]

        # Determine the next index
        next_index = max(indices) + 1 if indices else 1

        # Create new folder with the next index
        hillshade_holder = os.path.join(archive_folder, f"hillshade_holder_{next_index}")
        os.makedirs(hillshade_holder, exist_ok=True)

        # Move all .tif files to hillshade_holder
        for tif_file in tif_files:
            shutil.move(os.path.join(hillshade_Con, tif_file), hillshade_holder)

        # Zip the new folder
        zip_file_path = f"{hillshade_holder}.zip"
        with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(hillshade_holder):
                for file in files:
                    if not file.endswith('.lock'): 
                        zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), hillshade_holder))

        print(f"Moved {len(tif_files)} .tif files to {hillshade_holder} and zipped the folder to {zip_file_path}")
    else:
        pass

def get_point_sjo(key):
    sr = arcpy.SpatialReference(4326)
    id = int(key)
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
    cols = ['TARGET_FID', 'SHAPE@', 'FID_ID']
    numpy_array = arcpy.da.TableToNumPyArray(fp, cols)
    df = pd.DataFrame(numpy_array)
    df_filter = df[df['FID_ID'] == id]

    cell = df_filter.iloc[0, 1]
    return cell

def working_copy():
    """
    Description:
    Creates a point feature class from selected records in an Excel file and spatial data.
    
    Steps:
    1. Reads the 'Id' column from an Excel file.
    2. Retrieves the data source of the "GRIT_Minor_Structures" layer.
    3. Extracts geometries and filters records based on matching 'Id' values.
    4. Prompts user input for the number of records to process.
    5. Creates a new point feature class in a geodatabase.
    6. Adds required fields: 'MainFileID', 'FID_ID', and 'PL_Process'.
    7. Inserts selected records as point geometries into the feature class.
        sr = arcpy.SpatialReference(4326)
    """
    
    db = os.path.join(project_folder, "Default.gdb")
    output_excel_path_ip = os.path.join(project_folder, "IP_PointLayer.xlsx")
    df = pd.read_excel(output_excel_path_ip)

    # Extract only the 'Id' column
    df_id_only = df[['Id']]

    for map in project.listMaps():
            for layer in map.listLayers():
                if layer.name == "GRIT_Minor_Structures":
                    fp2 = layer.dataSource

    cols = ['Id', 'SHAPE@']
    numpy_array = arcpy.da.TableToNumPyArray(fp2, cols)
    df = pd.DataFrame(numpy_array)
    filtered_df = df[df['Id'].isin(df_id_only['Id'])]
    geo_arr = []
    id_arr = []
    FID = []

    sp = input("Please input number to process: ")
    count, max_processed = get_sp() 
    end = count + int(sp)
    max_processed += 1

    itter = count 
    print(itter, count, end)
    for index, val in filtered_df.iloc[count-1:end-1].iterrows():
        geo = val['SHAPE@']
        id_inst = val['Id']
        x = geo.centroid.X
        y = geo.centroid.Y
        point = arcpy.Point(x,y)
        geo_arr.append(point)
        id_arr.append(id_inst)
        FID.append(itter)
        itter += 1

    arcpy.env.workspace = db 

    feature_class_name = "PointLayer"
    arcpy.CreateFeatureclass_management(
        out_path=arcpy.env.workspace,
        out_name=feature_class_name,
        geometry_type="POINT",
        spatial_reference=sr
    )

    # Add the MainFileID field to the feature class
    arcpy.AddField_management(
        in_table=feature_class_name,
        field_name="MainFileID",
        field_type="TEXT"
    )

    arcpy.AddField_management(
        in_table=feature_class_name,
        field_name="FID_ID",
        field_type="LONG"
    )

    arcpy.AddField_management(
        in_table=feature_class_name,
        field_name="PL_Process",
        field_type="LONG"
    )

    with arcpy.da.InsertCursor(feature_class_name, ["SHAPE@", "MainFileID", "FID_ID", "PL_Process"]) as cursor:
        for point, main_file_id, idf in zip(geo_arr, id_arr, FID):
            poly_point = arcpy.PointGeometry(point, sr)
            cursor.insertRow([poly_point, main_file_id, idf, max_processed])

@time_it
def insert_loop():
    """
    Processes shapefiles, extracts the closest gridcode point, and creates polylines.
    """

    # Process shapefiles
    for filename in os.listdir(directory):
        if filename.endswith('.shp'):
            key = filename.split("_")[1]
            if int(key) in num_arr:
                shapefile_path = os.path.join(directory, filename)

                # Load shapefile into a DataFrame
                cols = ['FID', 'Id', 'gridcode', 'SHAPE@']
                numpy_array = arcpy.da.TableToNumPyArray(shapefile_path, cols)
                df = pd.DataFrame(numpy_array)

                # Get the smallest gridcode values
                df_sorted = df.nsmallest(n=30, columns=['gridcode'])

                # Process each point to find the closest one
                start_point = get_id_sjo(key)
                min_dist = float('inf')
                master_val = None

                for _, row in df_sorted.iterrows():
                    cell = row['SHAPE@']
                    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
                    projected_centroid = centroid_geom.projectAs(spatial_reference)

                    point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
                    end_point = arcpy.PointGeometry(point, spatial_reference)

                    distance = start_point.distanceTo(end_point)
                    if distance < min_dist:
                        min_dist = distance
                        master_val = end_point

                # Add to the dictionary
                try:
                    if key not in data_dict:
                        data_dict[key] = arcpy.Array()
                    data_dict[key].add(master_val.centroid)
                except Exception as e:
                    arcpy.AddMessage(f"Centroid error, Hillshade cutoff for {key}: {e}")


    # Insert the polylines into the feature class
    line_shp = os.path.join(gdb, 'PolylineLayer')
    arcpy.env.workspace = gdb
    count = 0

    with arcpy.da.InsertCursor(line_shp, ["SHAPE@", "FID"]) as cursor:
        for key, val in data_dict.items():
            try:
                # Try creating the polyline
                polyline = arcpy.Polyline(val, spatial_reference)
                cursor.insertRow([polyline, mainfileFID_arr[count]])
            except Exception as e:
                # If it fails, print the error and the value of val
                print(f"Error creating polyline for key {key} with value {val}: {e}")
            finally:
                count += 1

def extract_coordinates(shape):
    sr = arcpy.SpatialReference(4326)
    projected_shape = shape.projectAs(sr)
    point = arcpy.Point(projected_shape.centroid.X, projected_shape.centroid.Y)
    return point


def test_line_correction(point1, start_key, point2):
    """
    Adjusts line geometry to align with the most significant axis (horizontal or vertical).
    """
    
    start_point  = {}
    id_point = {}
    end_point = {}
    id_coords = extract_coordinates(start_key)
    start_point['Start'] = [point1.X, point1.Y]
    id_point['ID'] = [id_coords.X, id_coords.Y]
    end_point ['End'] = [point2.X, point2.Y]
    
    start_x, start_y = start_point['Start']
    id_x, id_y = id_point['ID']
    end_x, end_y = end_point['End']

    diff_x = abs(end_x - start_x)
    diff_y = abs(end_y - start_y)

    # Determine if the line is vertical or horizontal and adjust Start and End points
    if diff_x > diff_y:
        # More horizontal; adjust Y of Start and End to match ID
        start_y = id_y
        end_y = id_y
    else:
        # More vertical; adjust X of Start and End to match ID
        start_x = id_x
        end_x = id_x


    adjusted_start = arcpy.Point(start_x, start_y)
    adjusted_id = arcpy.Point(id_x, id_y)
    adjusted_end = arcpy.Point(end_x, end_y)

    array = arcpy.Array([adjusted_start, adjusted_id, adjusted_end])
    return array

if __name__ == "__main__":
    # Environment setup
    current_timestamp = datetime.now()
    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False

    # Get the current ArcGIS Pro project and its folder
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    gdb = os.path.join(project_folder, 'Default.gdb')
    
    # Create necessary folders and geodatabases
    #holderFolder()

    # Generate random number and get spatial information
    singlepart_rand = random.randint(1, 99999)
    x, sjo_Number = get_sp() 
    sjo_Number += 1
    spatial_reference = arcpy.SpatialReference(4326)  

    # Perform preliminary data processing
    working_copy()
    holderFolder()
    check_sqlite()
    spatial_join()

    # Get DataFrame from the spatial join result and populate folders
    df_ret = get_df("SpatialJoin_Output")
    populate_folder(df_ret)

    # CRS conversion
    name = "PointLayer"
    convert_crs(name)

    # Geospatial data processing
    process_geospatial_data()
    correct_output()
    
    # Mask extraction and shapefile conversion
    keys_list = mask_extraction()
    convert_shp(keys_list)

    # Prepare to process the shapefiles
    directory = os.path.join(project_folder, 'ConvertedSHP')
    file_desc = os.path.join(project_folder, f'SJO/SpatialJoin_Output{sjo_Number}.shp')
    
    desc = arcpy.Describe(file_desc)
    spatial_reference = desc.spatialReference

    # Initialize data storage
    data_dict = {}
    num_arr, mainfileFID_arr = get_ids()
    mainfileFID_arr = sorted(mainfileFID_arr)

    insert_loop()

    # Cleanup and finalize
    cleanup(keys_list)
    #clean_hillshades()
    print('done_main_file')
