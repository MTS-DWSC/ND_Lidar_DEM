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
    check_file = os.path.join(project_folder, "Hold_DEMs")
    make_file = os.path.join(project_folder, "Hillshade_Converted")

    if not os.path.exists(check_file):
        os.makedirs(check_file)
    if not os.path.exists(make_file):
        os.makedirs(make_file)
    return check_file

def check_sqlite():
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
        Y_LOCATION REAL
    )
    """)
    #print(f"SQLite database and table are ready at: {output_file}")
    conn.close()
    return    

def get_sp():
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
    # Create the join folder
    folder = "HolderFolder"
    dest = os.path.join(project_folder, folder)
    if not os.path.exists(dest):
        os.makedirs(dest)
        print(f"Folder created at {dest}")
    
    # Create folder to house the geodatabase
    lidar_gdb_path = os.path.join(project_folder, 'lidar.gdb')
    if not arcpy.Exists(lidar_gdb_path):
        arcpy.CreateFileGDB_management(project_folder, 'lidar.gdb')
        print("File geodatabase created.")

    # URL to the raw ZIP file
    url = "https://github.com/MTS-DWSC/ND_Lidar_DEM/raw/main/ND_Index.zip"
    check_file = os.path.join(project_folder, 'ND_Index.shp')
    
    if not os.path.exists(check_file):
        response = requests.get(url, verify=False)
        if response.status_code == 200:
            # Extract the ZIP file
            with zipfile.ZipFile(io.BytesIO(response.content)) as z:
                z.extractall(project_folder)  
            print("ZIP file downloaded and extracted successfully.")
        else:
            print(f"Failed to download the file. Status code: {response.status_code}")
    else:
        print(f"{check_file} already exists.")

    # Grab data source from API
    grit_layer = os.path.join(lidar_gdb_path, "GRIT_Minor_Structures")
    if arcpy.Exists(grit_layer):
        print("GRIT_Minor_Structures already exists in the geodatabase.")
        return
    else:
        # URL of the REST service
        service_url = "https://dotsc.ugpti.ndsu.nodak.edu:6443/arcgis/rest/services/GRIT_all/grit20_bridges_all_feature/MapServer/0"
        
        # Create a temporary layer from the REST service
        temp_layer = "in_memory/GRIT_Minor_Structures"
        arcpy.MakeFeatureLayer_management(service_url, temp_layer)
        print("Feature layer created from the REST service.")

        # Optionally, save the layer to the geodatabase
        arcpy.CopyFeatures_management(temp_layer, grit_layer)
        print("Feature layer copied to the geodatabase.")
    return

def delete_sj():
    output_folder = os.path.join(project_folder, "SJO")
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)

        if filename.endswith(".lock"):
            continue  # Skip this file if it ends with .lock
        else:
            os.remove(file_path) 
            
    output_holder = os.path.join(project_folder, "HolderFolder")
    for filename in os.listdir(output_holder):
        file_path = os.path.join(output_holder, filename)

        if filename.endswith(".lock"):
            continue  # Skip this file if it ends with .lock
        else:
            os.remove(file_path) 
        


@time_it
def spatial_join():
    db = os.path.join(project_folder, "Default.gdb")
    target_features = os.path.join(db, "PointLayer")
    
    join_features = project_folder + "\\" + "ND_Index.shp"
    
    # Specify a valid output path, ensuring no invalid characters
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
    fid_field_map.addInputField(target_features, "MainFileID") # Change for ID
    field_mappings.addFieldMap(fid_field_map)
    
    fid_field_map = arcpy.FieldMap()
    fid_field_map.addInputField(target_features, "PL_Process") # Change for ID
    field_mappings.addFieldMap(fid_field_map)

    # Add the Tilename field from IndexLidar
    tilename_field_map = arcpy.FieldMap()
    tilename_field_map.addInputField(join_features, "Tilename")
    field_mappings.addFieldMap(tilename_field_map)

    # Add the ASCII_Path field from IndexLidar
    ascii_path_field_map = arcpy.FieldMap()
    ascii_path_field_map.addInputField(join_features, "ASCII_Path")
    field_mappings.addFieldMap(ascii_path_field_map)

    # Perform the spatial join
    arcpy.analysis.SpatialJoin(
        target_features=target_features,
        join_features=join_features,
        out_feature_class=output_features,
        join_type="KEEP_COMMON",  
        field_mapping=field_mappings
    )
    arcpy.env.addOutputsToMap = False
    print("done_spatial_join")



def get_df(text):
    ## Get DataFrame for Master File
    project = arcpy.mp.ArcGISProject("CURRENT")
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
    
    #switch = 0
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
    destination_folder = check_gdb()
    hs_folder = os.path.join(project_folder, 'Hillshade_Converted')
    check_file = project_folder + "\\" + "output_file.db"
    if os.path.exists(check_file):
        output_file = os.path.join(project_folder, "output_file.db")
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()
        #df_ret = get_df("SpatialJoin_Output")

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
                    INSERT INTO recordManager (FID, File_name, File_path, FTP, Processed, MainFileID, X_LOCATION, Y_LOCATION)
                    VALUES (?, ?, ?, ?, ?,?, ?, ?)
                """, (int(lookup['FID_ID'].iloc[0]) - 0, lookup['TILENAME'].iloc[0], fn_final, res, int(lookup['PL_Process'].iloc[0]), int(lookup['MainFileID'].iloc[0]), lookup['X_LOCATION'].iloc[0], lookup['Y_LOCATION'].iloc[0]))
            except Exception as e:
                cursor.execute("""
                    INSERT INTO recordManager (FID, File_name, File_path, FTP, Processed, MainFileID, X_LOCATION, Y_LOCATION)
                    VALUES (?, ?, ?, ?, ?,?, ?, ?)
                """, (int(lookup['FID_ID'].iloc[0]) - 0, lookup['TILENAME'].iloc[0], "Invalid", "Invalid", int(lookup['PL_Process'].iloc[0]), int(lookup['MainFileID'].iloc[0]), lookup['X_LOCATION'].iloc[0], lookup['Y_LOCATION'].iloc[0]))
                print(f"An unexpected error occurred for {file_name}: {e}. Skipping...")
                continue
        conn.commit()
        conn.close()
        

## Not looping through it
@time_it
def raster_crs_conversion():
    arcpy.env.addOutputsToMap = False
    fp = os.path.join(project_folder, "Hold_DEMs")
    count = 0
    # -- Specify the reference dataset for spatial reference --
    spatial_ref = arcpy.SpatialReference(4326)
    # -- Define the output folder --
    output_folder = os.path.join(project_folder, "Projected")
    os.makedirs(output_folder, exist_ok=True)  # Ensure the output folder exists
    for file in os.listdir(fp):
        full_path = os.path.join(fp, file)
        # -- Check if the file is a valid raster --
        if arcpy.Describe(full_path).datatype == 'RasterDataset':
            count += 1
            fn = f"{os.path.splitext(file)[0]}"
            out_raster = os.path.join(output_folder, f"{fn}.tif")
        try:
            # -- Project the raster --
            #arcpy.management.ProjectRaster(full_path, out_raster, spatial_ref)
            temp_projected_raster = arcpy.management.ProjectRaster(full_path, "in_memory/temp_projected", spatial_ref)
            exaggerated_dem = Times(temp_projected_raster, 1.5)

            # -- Save the exaggerated DEM --
            exaggerated_dem.save(out_raster)

            #print(f"Projected: {fn} to {spatial_ref.name}")
        except Exception as e:
            print(f"Error processing {full_path}: {e}")
    print("Finished CRS conversion.")
    return

@time_it
def convert_crs(name):
    fp = check_gdb()
    fp_split = fp.split('\\')[:-1]
    fp_converted = '\\'.join(fp_split) + '\\' + 'Hillshade_Converted'
    #print(os.listdir(fp))
    
    dataset = name  # Fully qualify or adjust as needed
    spatial_ref = arcpy.Describe(dataset).spatialReference
    count = 0
    
    new_fp = os.path.join(project_folder, 'Hold_DEMs')
    
    
    for file in os.listdir(new_fp):
        file_name, file_extension = os.path.splitext(file)
        full_path = os.path.join(new_fp, file)
        # -- Check if the file is a valid raster --
        if arcpy.Describe(full_path).datatype == 'RasterDataset':
            in_raster = full_path
            out_hillshade = Hillshade(in_raster)
            #out_hillshade = Hillshade(in_raster, model_shadows = "SHADOWS", altitude = 38)
            
            output_name = f"{file_name}_Hillshade.tif"
            output_path = os.path.join(fp_converted, output_name)
            if os.path.exists(output_path):
                continue
            out_hillshade.save(output_path)
            count += 1
            if count % 5 == 0:
                print(f"Processing complete for {output_name}")
    return

@time_it
def process_geospatial_data(buffer_distance="50 Meters"):
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
            print("Operation succeeded")
            break  # Exit the loop if successful
        except arcpy.ExecuteError as e:
            # If arcpy raises an error, catch it
            print(f"Error: {e}")
            attempt += 1
            if attempt < max_retries:
                print(f"Retrying... (Attempt {attempt} of {max_tries})")
                time.sleep(2)  # Wait for 2 seconds before retrying
            else:
                print(f"Failed after {max_tries} attempts.")

@time_it
def correct_output():
    path = os.path.join(project_folder, "HolderFolder")
    feature_class = os.path.join(path, f"singlepart_output{singlepart_rand}.shp")

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
    output_file = os.path.join(project_folder, "output_file.db")
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
    #print(sorted_filepaths)
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
            #os.makedirs(subfolder_path, exist_ok=True)

            mask_1.save(os.path.join(subfolder_path, filtered_layerA))
            mask_2.save(os.path.join(subfolder_path, filtered_layerB))
            print(f"Mask Extraction: {filtered_layerA} and {filtered_layerB}.")


@time_it
def convert_shp():
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
        print("No rasters found in the geodatabase.")


def get_id_sjo(key):
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


def cleanup(count: int) -> None:
    ConvertedSHP = os.path.join(project_folder, 'ConvertedSHP')
    Hold_DEMs = os.path.join(project_folder, 'Hold_DEMs')
    HolderFolder = os.path.join(project_folder, 'HolderFolder')
    lidar_gdb = os.path.join(project_folder, 'lidar.gdb')
    arcpy.env.workspace = lidar_gdb
    datasets = arcpy.ListRasters()
    itter = 0
    if len(datasets) > count:
        for data in datasets:
            fp = os.path.join(lidar_gdb, data)
            arcpy.Delete_management(fp)
            itter += 1
            if itter == count:
                print(f"Deleted {count} number of files.")
                break

        shutil.rmtree(ConvertedSHP)
        shutil.rmtree(Hold_DEMs)
        delete_sj()
        print("End clean.")

def working_copy():
    sr = arcpy.SpatialReference(4326)
    db = os.path.join(project_folder, "Default.gdb")

    for map in project.listMaps():
            for layer in map.listLayers():
                if layer.name == "GRIT_Minor_Structures":
                    fp2 = layer.dataSource

    cols = ['Id', 'SHAPE@']
    numpy_array = arcpy.da.TableToNumPyArray(fp2, cols)
    df = pd.DataFrame(numpy_array)


    geo_arr = []
    id_arr = []
    FID = []


    sp = input("Please input number to process: ")
    count, max_processed = get_sp() 
    end = count + int(sp)
    max_processed += 1

    itter = count 
    print(itter, count, end)
    #print(df.iloc[0])
    for index, val in df.iloc[count-1:end-1].iterrows():
        #if itter == end:
        #    break

        geo = val['SHAPE@']
        id_inst = val['Id']
        #print(id_inst)
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

    '''
    with arcpy.da.InsertCursor(feature_class_name, ["SHAPE@"]) as cursor:
        for x in geo_arr:
            poly_point = arcpy.PointGeometry(x, sr)
            cursor.insertRow(poly_point)
    '''
    # Add the feature class to the current map
    aprx = arcpy.mp.ArcGISProject("CURRENT")  # Use

if __name__ == "__main__":
    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False
    
    holderFolder()
    singlepart_rand = random.randint(1, 99999)
    working_copy()
    
    # Get the current ArcGIS Pro project
    project = arcpy.mp.ArcGISProject("CURRENT")

    # Get the folder containing the ArcGIS Pro project
    project_folder = os.path.dirname(project.filePath)
    
    gdb = os.path.join(project_folder, 'Default.gdb')
    
    x, sjo_Number = get_sp() 
    sjo_Number += 1
    spatial_reference = arcpy.SpatialReference(4326)  
    # Begin Sequence
    name = "PointLayer"
    create_lidargdb()
    check_sqlite()
    spatial_join()
    df_ret = get_df("SpatialJoin_Output")
    populate_folder(df_ret)

    ##raster_crs_conversion() Need for ML, breaks current itteration
    convert_crs(name)
    
    # ------------------------
    process_geospatial_data()
    correct_output()
    
    mask_extraction()
    convert_shp()
    
    # ------------------------

    # Set the directory and spatial reference
    directory = os.path.join(project_folder, 'ConvertedSHP')
    file_desc = os.path.join(project_folder, f'SJO\SpatialJoin_Output{sjo_Number}.shp')
    
    desc = arcpy.Describe(file_desc)
    spatial_reference = desc.spatialReference
    data_dict = {}
    num_arr, mainfileFID_arr = get_ids()[0], get_ids()[1]
    
    # Loop through the files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.shp'):
            key = filename.split("_")[1]
            if int(key) in num_arr: 
                shapefile_path = os.path.join(directory, filename)

                # Read shapefile into a pandas DataFrame
                cols = ['FID', 'Id', 'gridcode', 'SHAPE@']
                numpy_array = arcpy.da.TableToNumPyArray(shapefile_path, cols)
                df = pd.DataFrame(numpy_array)

                # Sort DataFrame by 'gridcode' column to get the smallest values
                df_sorted = df.nsmallest(n=30, columns=['gridcode'])

                # Extract the key from the filename
                key = filename.split("_")[1]
                min_dist = float('inf')  # Initialize min_dist as infinity
                start_point = get_id_sjo(key)  

                # Iterate over the sorted DataFrame
                for index, row in df_sorted.iterrows():
                    cell = row['SHAPE@']
                    #id_FID = row['MainFileID']
                    geo_sr = spatial_reference
                    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
                    # Project the centroid to WGS 84
                    projected_centroid = centroid_geom.projectAs(geo_sr)

                    x = projected_centroid.centroid.X  # Longitude
                    y = projected_centroid.centroid.Y  # Latitude
                    point = arcpy.Point(x, y)
                    # Create a PointGeometry object for the current point
                    end_point = arcpy.PointGeometry(point, spatial_reference)

                    # Calculate the distance from the start point
                    distance = start_point.distanceTo(end_point)

                    # Update min_dist if the current distance is smaller
                    if distance < min_dist:
                        min_dist = distance
                        master_val = end_point

                # Add the master_val to the dictionary under the key
                if key not in data_dict:
                    data_dict[key] = arcpy.Array()
                data_dict[key].add(master_val.centroid)  

    # Fix this to work with existing
    line_gdb = os.path.join(project_folder, 'Default.gdb')
    line_shp = os.path.join(line_gdb, 'PolylineLayer')

    # Set the workspace to the geodatabase containing the feature class
    arcpy.env.workspace = line_gdb
    
    mainfileFID_arr = sorted(mainfileFID_arr) 
    count = 0
    with arcpy.da.InsertCursor(line_shp, ["SHAPE@", "FID"]) as cursor:
        for key, val in data_dict.items():
            # Create a polyline geometry
            polyline = arcpy.Polyline(val, spatial_reference)
            # Insert the new polyline into the feature class
            cursor.insertRow([polyline, mainfileFID_arr[count]])
            count += 1
    cleanup(10)
    print('done_main_file')



