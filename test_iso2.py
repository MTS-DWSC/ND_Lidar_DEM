def extract_coordinates(shape):
    sr = arcpy.SpatialReference(4326)
    projected_shape = shape.projectAs(sr)
    point = arcpy.Point(projected_shape.centroid.X, projected_shape.centroid.Y)
    return point

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

"""
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
"""

def true_value():
    
    # -----------------------
    cols = ['OID@', 'diss', 'SHAPE@']
    numpy_array = arcpy.da.TableToNumPyArray('dissolved_output', cols)
    df = pd.DataFrame(numpy_array)
    dis_df = df[df['diss'] == 'True']
    shape_value = dis_df['SHAPE@'].iloc[0]
    centroid_geom = arcpy.PointGeometry(shape_value.centroid, shape_value.spatialReference)
    projected_centroid = centroid_geom.projectAs(sr)
    point_focus = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
    to_start_focus = arcpy.PointGeometry(point_focus, sr)
    # -----------------------
    
    
    
    
    # -----------------------
    output_folder = os.path.join(project_folder, "SJO")
    fp = os.path.join(output_folder, f"SpatialJoin_Output{sjo_Number}.shp")
    cols = ['TARGET_FID', 'SHAPE@', 'FID_ID']
    numpy_array = arcpy.da.TableToNumPyArray(fp, cols)
    df1 = pd.DataFrame(numpy_array)
    df_filter = df1[df1['FID_ID'] == 4]

    cell = df_filter.iloc[0, 1]
    x = cell.centroid.X
    y = cell.centroid.Y
    point = arcpy.Point(x, y)
    newT = arcpy.PointGeometry(point, sr)
    # -----------------------
    
    
    
    arr = []
    df['Eliminate'] = True
    df['Distance_lookup'] = 0
    df['Distance_central'] = 0
    for index, row in df.iterrows():
        cell = row['SHAPE@']
        cg = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
        pc = cg.projectAs(sr)

        p = arcpy.Point(pc.centroid.X, pc.centroid.Y)
        ep = arcpy.PointGeometry(p, sr)
        
        arr.append(ep)
        
        dist_to_lookup = ep.distanceTo(to_start_focus)
        dist_to_central = ep.distanceTo(newT)
        
        if dist_to_central < dist_to_lookup:
            df.at[index, 'Eliminate'] = False
        else:
            df.at[index, 'Eliminate'] = True
            
        df.at[index, 'Distance_lookup'] = dist_to_lookup
        df.at[index, 'Distance_central'] = dist_to_central
        
    filtered_df = df[df['Eliminate'] == False]
        
    min_distance_row = filtered_df.loc[filtered_df['Distance_central'].idxmin()]

    # Get the 'OID@' of that row
    min_oid = min_distance_row['SHAPE@']
    
    #spoint = extract_coordinates(shape_start)
    epoint = extract_coordinates(min_oid)

    #print(df[['OID@', 'diss', 'Eliminate', 'Distance_lookup', 'Distance_central']])
    print(epoint)
    return epoint

def isolated_points():
    count = 1
    sr = arcpy.SpatialReference(4326)
    for filename in os.listdir(directory):
        if filename.endswith('.shp'):
            key = filename.split("_")[1]
            if int(key) in num_arr:
                shapefile_path = os.path.join(directory, filename)
                arcpy.AddField_management(shapefile_path, 'diss', 'TEXT', field_length=5)
                
                with arcpy.da.UpdateCursor(shapefile_path, ['diss']) as cursor:
                    for row in cursor:
                        row[0] = 'False' 
                        cursor.updateRow(row)
                


                # Load shapefile into a DataFrame
                cols = ['FID', 'Id', 'gridcode', 'SHAPE@', 'diss']
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
                filtered_df['diss'] = 'False'
                for index, row in filtered_df.iterrows():
                    cell = row['SHAPE@']
                    centroid_geom = arcpy.PointGeometry(cell.centroid, cell.spatialReference)
                    projected_centroid = centroid_geom.projectAs(sr)

                    point = arcpy.Point(projected_centroid.centroid.X, projected_centroid.centroid.Y)
                    end_point = arcpy.PointGeometry(point, sr)

                    distance = start_point.distanceTo(end_point)
                    filtered_df.at[index, 'distance'] = distance

                    # ----

                # First, find the row with the minimum distance
                lowest_gridcode_df = filtered_df.nsmallest(5, 'gridcode')
                fid_of_min_value = lowest_gridcode_df.at[lowest_gridcode_df['distance'].idxmin(), 'Id']
                shape_start = filtered_df.loc[filtered_df['Id'] == fid_of_min_value, 'SHAPE@'].values[0]
                #spoint = extract_coordinates(shape_start)
                
                #filtered_df.loc[filtered_df['Id'] == fid_of_min_value, 'diss'] = 'True'
                with arcpy.da.UpdateCursor(shapefile_path, ['Id', 'diss']) as cursor:
                    for row in cursor:
                        if row[0] == fid_of_min_value:
                            row[1] = 'True' 
                            cursor.updateRow(row)
                
                # Dissolve operation
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
                arcpy.analysis.PairwiseDissolve("temp_layer", output_dissolve_path, dissolve_field="diss", multi_part = "SINGLE_PART")
                
                oid = true_value()
                
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
                
                if count == 1:
                    break

if __name__ == "__main__":
    # Environment setup
    current_timestamp = datetime.now()
    arcpy.env.overwriteOutput = True
    #arcpy.env.addOutputsToMap = False

    # Get the current ArcGIS Pro project and its folder
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    gdb = os.path.join(project_folder, 'Default.gdb')
    
    # Create necessary folders and geodatabases
    #holderFolder()

    # Generate random number and get spatial information
    singlepart_rand = random.randint(1, 99999)
    x, sjo_Number = get_sp() 
    sjo_Number += 0#1
    spatial_reference = arcpy.SpatialReference(4326)  

    # Perform preliminary data processing
    '''
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
    '''
    # Prepare to process the shapefiles
    directory = os.path.join(project_folder, 'ConvertedSHP')
    file_desc = os.path.join(project_folder, f'SJO/SpatialJoin_Output{sjo_Number}.shp')
    
    desc = arcpy.Describe(file_desc)
    spatial_reference = desc.spatialReference

    # Initialize data storage
    data_dict = {}
    num_arr, mainfileFID_arr = get_ids()
    mainfileFID_arr = sorted(mainfileFID_arr)

    # Process shapefiles
    isolated_points()

    # Cleanup and finalize
    #cleanup(keys_list)
    #clean_hillshades()
    print('done_main_file')
            
