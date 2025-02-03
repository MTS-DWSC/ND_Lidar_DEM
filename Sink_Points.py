import arcpy
from arcpy.sa import *
import os
import os
import sys
import traceback
import arcpy
from arcpy import env
from arcpy.sa import *

def JoinOneField(tbl, key1, jtbl, key2, jfield, ftype="", fprec="", flen="",
                 fname="", msg=True):
    """Joins one field from another table (like Join Field, but fast,
    and field values are not modified for no-matches).

    tbl - source table
    key1 - source table key field
    jtbl - join table
    key2 - join field (must be unique in jtbl)
    jfield - field to copy over
    ftype - field type
    fprec - field precision
    flen - field length (for ftype="TEXT")
    fname - field name (if different from jfield)
    msg - print messages to gp (including warnings)

    The function returns a two-integer tuple of number of matches and no-matches.
    ** The field data values are not changed for no-matches **

    JoinOneField(srctbl, key1, jtbl, key2, NHDPlusID", "DOUBLE")"""

    def AddFieldType(tbl, field_name):
        """Determines keyword to create a field with the AddField tool of
        the same type as the specified field.

        Integer > "LONG", Float -> "SINGLE" etc

        example

        AddFieldType("test.dbf", "Shape_Leng")
        "DOUBLE"
        arcpy.AddField_management("test.dbf", "Shape_Len2",
                                  AddFieldType("test.dbf", "Shape_Leng"))
        """
        dtypes = {'OID': 'LONG', 'SmallInteger': 'SHORT', 'Integer': 'LONG',
                  'Double': 'DOUBLE', 'Single': 'FLOAT', 'String': 'TEXT',
                  'Date': 'DATE'}
        fld_list = arcpy.Describe(tbl).Fields
        ftype = [f.type for f in fld_list
                 if f.name.upper() == field_name.upper()]
        if not ftype:
            raise Exception("Field {} not found".format(field_name))
        ftype = ftype[0]
        try:
            return dtypes[ftype]
        except KeyError:
            msg = "{}: Field type {} not supported".format(field_name, ftype)
            raise Exception(msg)

    # Set default field type to be added
    # (Type can be converted by specifying this argument)
    if not ftype:
        ftype = AddFieldType(jtbl, jfield)

    # By default the join field and destination field are the same name
    # You can use a different name for the new field by specifying it
    if not fname:
        fname = jfield

    # Make sure (for the case of shapefile/dbf tables)
    # we have enough precision for double.
    # Default fprec can truncate large #s
    if ftype.upper() == "DOUBLE" and not fprec:
        fprec = 15

    if not arcpy.ListFields(tbl, key1):
        raise Exception("Key Field {} not found in {}".format(key1, tbl))
    if not arcpy.ListFields(jtbl, key2):
        raise Exception("Key Field {} not found in {}".format(jfield, jtbl))
    if not arcpy.ListFields(jtbl, jfield):
        raise Exception("Join Field {} not found in {}".format(jfield, jtbl))
    if (arcpy.ListFields(jtbl, jfield)[0].type == "OID" and
        fname[:3].upper() in ["OBJ", "FID"]):
        # we can't name the field OBJECTID or FID, so use JOIN_OID
        fname = "JOIN_OID"

    ff = arcpy.ListFields(tbl, fname)
    if ff:
        print("Copying matching data to existing field {}".format(fname))
    else:
        print("Adding field {} as {}".format(fname, ftype))
        arcpy.AddField_management(tbl, fname, ftype, fprec, "", flen)

    # load data into dictionary
    jdict = {}
    k = 0
    first_dupe = True
    with arcpy.da.SearchCursor(jtbl, [key2, jfield]) as rows:
        for row in rows:
            key, val = row
            if first_dupe and key in jdict:
                if msg:
                    # display a single warning message for one to many joins
                    print("w", "Key field {} not unique in {}".format(key2, jtbl))
                first_dupe = False
            else:
                jdict[key] = val

    # copy data values to output table
    num_match = 0
    num_nomatch = 0

    with arcpy.da.UpdateCursor(tbl, [key1, fname]) as rows:
        for row in rows:
            try:
                key = row[0]
                row[1] = jdict[key]
                num_match += 1
                rows.updateRow(row)
            except KeyError:
                # no match in dictionary (join no-match)
                # do nothing -- do not change existing value
                num_nomatch += 1
                ## print ("{}:{}".format(key, "no match"))
            except RuntimeError as msg:
                if str(msg).find("The value type is incompatible") != -1:
                    raise Exception(
                        "Invalid value for field {}: {!r}".format(
                            jfield, jdict[key]))
                else:
                    # some other error
                    raise

    # return number of rows matched and not matched

    return num_match, num_nomatch

def FieldNameMap(tbl, maps):
    """Create a field mappings object to rename, drop, or re-order fields.

    Arguments

      tbl - input feature class, table, or table view
      maps - field names map list (';'-delimited string)
            If a field is left out, it will not be included in the field map

    Example

      import arcpy
      Maps = "Shape_Area AREA;BID BID2;AREASQMI AREAMI2"
      Mapper = FieldNameMap("temp.dbf", Maps)
      # to debug: print out the field mapping
      print Mapper.exportToString().replace(";","\n")
      # copy table, keeping only renamed fields: AREA, BID2, AREAMI2
      arcpy.Merge_management("temp.dbf","temp2.dbf",Mapper)
    """

    field_mappings = arcpy.FieldMappings()
    mapList = maps.split(';')

    for rec in mapList:
        fromName,toName = rec.split()
        # create a new field map
        field_map = arcpy.FieldMap()
        # populate it and add to field_mappings
        try:
            field_map.addInputField(tbl, fromName)
            field = field_map.outputField
            field.name = toName
            field_map.outputField = field
            field_mappings.addFieldMap(field_map)
        except:
            raise Exception("FieldNameMap: Cannot not map fields ({}) in {}".format(rec,tbl))

    return field_mappings

def process_lidar_data():
    # Set environment settings
    #name = "Extract_dem12"
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    hold = os.path.join(project_folder, "resfolder")
    gdb = os.path.join(project_folder, 'lidar_1.gdb')

    # Get the raster and parameters
    elev_raster = Raster(os.path.join(gdb, name))
    min_area_process_lidar = 1000
    min_area_units = "CELLS"
    
    #cell_size = elev_raster.meanCellHeight
    #extent = elev_raster.extent
    #thresh = float(min_area_process_lidar)
    
    #print(f"Cell size: {cell_size}")
    #print(f"Raster extent: {extent}")
    
    arcpy.env.extent = "MAXOF"  # Ensures maximum possible extent
    arcpy.env.cellSize = elev_raster.meanCellHeight
    arcpy.env.snapRaster = elev_raster
    arcpy.env.outputCoordinateSystem = elev_raster.spatialReference

    # Fill the raster
    out_fill = Fill(elev_raster)
    out_fill.save(os.path.join(hold, f"output_fill_{name}.tif"))

    # Flow Direction
    fdr = FlowDirection(out_fill, "FORCE")
    fdr.save(os.path.join(hold, f"output_fdr_{name}.tif"))

    # Difference Raster
    diff = out_fill - elev_raster
    diff.save(os.path.join(hold, f"output_diff_{name}.tif"))

    # Flow Accumulation
    fac = FlowAccumulation(fdr)
    fac.save(os.path.join(hold, f"output_acc_{name}.tif"))
    print('Finished processing DEM.')

def calculate_pour_sink():   
    #min_area = 1000
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    hold = os.path.join(project_folder, "resfolder")
    gdb = os.path.join(project_folder, 'lidar_1.gdb')
    min_cells = int(min_area)
    #name = "Extract_dem12"
    elev = Raster(os.path.join(gdb, name))
    fill_raster = os.path.join(hold, ("output_fill_" + name + ".tif"))
    # Filters for only values that matter
    fillcells = Con(Diff(elev, Raster(fill_raster)), 1)
    fillreg = RegionGroup(fillcells, "EIGHT")
    if fillreg.minimum == None:
        print("No fill zones found")

    arcpy.BuildRasterAttributeTable_management(fillreg)
    count = arcpy.management.GetCount(fillreg)[0]
    print("{} fill zones found".format(count))

    # ----

    where = '"COUNT" > {0}'.format(min_cells)
    fillzones = Int(ExtractByAttributes(fillreg, where)) # force to integer, 10.6 gave me a float
    fillzone_raster = os.path.join(hold, f"raster_{name}.tif")
    fillzones.save(fillzone_raster)
    arcpy.BuildRasterAttributeTable_management(fillzone_raster)

    print("Section one finished.")

    # ----

    # Find pour points (cell with maximum COUNT) in each fill zone
    flowaccum_raster = os.path.join(hold, ("output_acc_" + name + ".tif"))
    fac = Raster(flowaccum_raster)
    # Notice here it uses flowaccumulation instead of Elevation
    maxfac = ZonalStatistics(fillzones, "VALUE", fac, "MAXIMUM")
    pourcells = SetNull(Diff(maxfac, fac), fillzones)

    # Find minimum elevation cells in each fill zone
    minele = ZonalStatistics(fillzones, "VALUE", elev, "MINIMUM")
    mincells = SetNull(Diff(elev, minele), fillzones)
    #print(maxfac)

    # ----

    xxpoly = arcpy.CreateScratchName("xxpoly", "", "featureclass", gdb)
    xxpoly1 = arcpy.CreateScratchName("xxpoly1", "", "featureclass", gdb)
    arcpy.RasterToPolygon_conversion(fillzones, xxpoly, "NO_SIMPLIFY")
    GRIDCODE = arcpy.ListFields(xxpoly, "GRID*")[0].name
    arcpy.Dissolve_management(xxpoly, xxpoly1, GRIDCODE)

    # ----

    scr = gdb
    # Calculate fill zone elevation statistics
    xxstats = arcpy.CreateScratchName("xxstats", "", "table", scr)
    ZonalStatisticsAsTable(fillzones, "VALUE", elev, xxstats)

    # ----

    lyr = "tmplyr"
    arcpy.MakeFeatureLayer_management(xxpoly1, lyr)
    GRIDCODE = arcpy.ListFields(lyr, "GRID*")[0].name
    arcpy.AddJoin_management(lyr, GRIDCODE, xxstats, "VALUE")

    # ----

    # want to do a field mapping here to make a nicer output table
    fillzone_poly = os.path.join(hold, f"fillzone_poly_{name}.shp")
    pp = arcpy.Describe(xxpoly1).name
    ss = arcpy.Describe(xxstats).name
    maps = ("{0}.{2} {2};"
            "{1}.COUNT COUNT;"
            "{1}.AREA  AREA;"
            "{1}.MIN MIN_ELE;"
            "{1}.MAX FILL_ELE;"
            "{1}.RANGE FILL_DEP").format(
                arcpy.Describe(xxpoly1).name,
                arcpy.Describe(xxstats).name,
                GRIDCODE)

    fms = FieldNameMap(lyr, maps)
    print('---------')

    arcpy.Merge_management(lyr, fillzone_poly, fms)
    arcpy.Delete_management(lyr)
    print('Section two finished.')

    # ----

    pour_points = os.path.join(hold, f"pour_points_{name}.shp")
    # Check tablehere
    min_points = os.path.join(hold, f"min_points_{name}.shp")
    arcpy.RasterToPoint_conversion(pourcells, pour_points)
    arcpy.RasterToPoint_conversion(mincells, min_points)

    # ----

    # add elevation values to point feature classes
    stats = arcpy.Describe(xxstats).baseName
    arcpy.MakeFeatureLayer_management(pour_points, lyr)
    GRIDCODE = arcpy.ListFields(lyr, "GRID*")[0].name

    # ----

    arcpy.AddField_management(lyr, "FZONE", "LONG")
    expr = "!{}!".format(GRIDCODE)
    arcpy.CalculateField_management(lyr, "FZONE", expr)

    # ----

    arcpy.DeleteField_management(lyr, ["POINTID", GRIDCODE])
    arcpy.AddField_management(lyr, "FILLELEV", "DOUBLE")
    arcpy.AddJoin_management(lyr, "FZONE", xxstats, "VALUE")
    expr = "!{0}.{1}!".format(stats, "MAX")
    arcpy.CalculateField_management(lyr, "FILLELEV", expr)
    arcpy.Delete_management(lyr)
    # Comment this out
    arcpy.MakeFeatureLayer_management(min_points, lyr)

    # ---- 

    # convert grid codes to FZONE field (LONG)
    GRIDCODE = arcpy.ListFields(lyr, "GRID*")[0].name
    arcpy.AddField_management(lyr, "FZONE", "LONG")
    expr = "!{}!".format(GRIDCODE)
    arcpy.CalculateField_management(lyr, "FZONE", expr)
    arcpy.DeleteField_management(lyr, GRIDCODE)
    arcpy.AddField_management(lyr, "ELEV", "DOUBLE")
    arcpy.AddJoin_management(lyr, "FZONE", xxstats, "VALUE")
    expr = "!{0}.{1}!".format(stats, "MIN")
    arcpy.CalculateField_management(lyr, "ELEV", expr)
    arcpy.RemoveJoin_management(lyr, arcpy.Describe(xxstats).name)

    # ----

    tmp_near = "in_memory\\xxnear"

    #print(arcpy.Exists(lyr))
    #print(arcpy.Exists(pour_points))

    arcpy.env.workspace = project_folder
    arcpy.env.scratchWorkspace = project_folder

    # ----

    # Generate in-memory name explicitly
    tmp_near = os.path.join(gdb, 'xxnear')

    # Validate inputs
    if arcpy.Exists(lyr) and arcpy.Exists(pour_points):
        arcpy.GenerateNearTable_analysis(lyr, pour_points, tmp_near, closest_count=5)
    else:
        raise ValueError("Input layers do not exist!")

    # Remove rows that are not min point -> pour point links
    print(arcpy.Exists(tmp_near))
    if arcpy.Exists(tmp_near):
        lyr = arcpy.MakeTableView_management(tmp_near, "tmp_near")
        tmp_near_nm = arcpy.Describe(tmp_near).name
        d_mp = arcpy.Describe(min_points)
        d_pp = arcpy.Describe(pour_points)
    else:
        raise ValueError("Near table was not created successfully!")

    print('Generate Near Table.')

    # ----

    arcpy.AddJoin_management(lyr, "IN_FID", min_points,
                                         d_mp.OIDFieldName)
    arcpy.AddJoin_management(lyr, "{}.NEAR_FID".format(tmp_near_nm),
                             pour_points,
                             d_pp.OIDFieldName)

    # ----

    where = "{0}.{1} <> {2}.{1}".format(d_mp.baseName, "FZONE",
                                        d_pp.baseName, "FZONE")
    arcpy.SelectLayerByAttribute_management(lyr, "", where)
    arcpy.RemoveJoin_management(lyr, d_mp.baseName)
    arcpy.DeleteRows_management(lyr)
    # Copy distances over to min_points
    arcpy.JoinField_management(min_points, d_mp.OIDFieldName,
                               tmp_near, "IN_FID", "NEAR_DIST")
    print('Join on nearest.')

    # ----

    # count number of min points per zone
    tmpcount = os.path.join(gdb, 'xmpcount')
    arcpy.Statistics_analysis(min_points, tmpcount, "FZONE COUNT", "FZONE")
    JoinOneField(min_points, "FZONE", tmpcount, "FZONE", "COUNT_FZONE", fname="MPCOUNT")
    # count zones with dupes
    ndup = 0
    with arcpy.da.SearchCursor(tmpcount, ["FREQUENCY", "FZONE"], "FREQUENCY > 1") as rows:
        for row in rows:
            ndup += 1

    # ----

    # If any found, report count and remove duplicates
    if ndup > 1:
        msg =  ("{} fill zones have more than one minimum point.\n"
                "Number of minimum points found in zone assigned to field MPCOUNT.\n"
                "Removing duplicate minimum points...")
        print(msg.format(ndup))
        # create list of dupe min points, sorted by distance from pp
        tmp1 = os.path.join(gdb, 'tmp1')
        lyrMP = arcpy.MakeFeatureLayer_management(min_points, "lyrMP")
        arcpy.SelectLayerByAttribute_management(lyrMP, "", "MPCOUNT > 1")
        arcpy.Sort_management(lyrMP, tmp1, "FZONE;NEAR_DIST")
        # delete all except the min point closest to the pp
        fz = -1
        with arcpy.da.UpdateCursor(tmp1, ["FZONE", "NEAR_DIST","MPCOUNT"]) as rows:
            for row in rows:
                fz1 = row[0]
                if fz != fz1:
                    fz = fz1 # this is the first one, skip
                else:
                    row[2] = -1 # dupe, tag for delete
                    rows.updateRow(row)
        # select all tagged min points and delete
        arcpy.AddJoin_management(lyrMP, "pointid", tmp1, "pointid")
        arcpy.SelectLayerByAttribute_management(lyrMP, "", '"tmp1.MPCOUNT" = -1')
        arcpy.RemoveJoin_management(lyrMP, "tmp1")
        arcpy.DeleteFeatures_management(lyrMP)
        for k in [lyrMP, tmpcount, tmp1]:
            arcpy.Delete_management(k)

    print('Final.')

def draw_flowpath():
    ptlist = []
    buf = 200
    with arcpy.da.SearchCursor(in_point,("OID@", "SHAPE@XY")) as Rows:
        for Row in Rows:
            ptlist.append(list(Row))
    numpts = len(ptlist)
    print(numpts)

    # make raster objects
    elev = Raster(elev_raster)
    fz = Raster(fill_raster)

    # set raster environment
    env.snapRaster = elev.catalogPath
    env.cellSize = elev.meanCellHeight
    OID = arcpy.Describe(in_point).OIDFieldName

    k = 0

    env.workspace =  env.scratchFolder
    env.scratchWorkspace =  env.scratchFolder
    tback = arcpy.CreateScratchName("xxback", "", "RasterDataset ")
    tline = arcpy.CreateScratchName("xxflowline", "", "FeatureClass", "in_memory")
    env.overwriteOutput = True
    
    # ----
    
    # Initialize loop counter and flag
    first = True
    count = 0
    # Define the geodatabase and feature class path for storing flowlines
    dgdb = os.path.join(project_folder, "Default.gdb")
    tline = f"{dgdb}\\xxflowline"

    # Create clusters of highest densed minpoints via buffer zones on self-join
    # Areas outside of buffer zones will not be considered
    # Isolated loners ignored
    # Loop through each point in the point list (ptlist)
    for kk, pt in enumerate(ptlist):

        # Get point id (oid) and coordinates (x, y) for the current point
        oid = pt[0]
        x, y = pt[1]

        # Set up the raster environment (buffer size and extent around the point)
        buf1 = buf * 2  # Buffer is doubled for extent calculation
        env.extent = arcpy.Extent(x - buf1, y - buf1, x + buf1, y + buf1)

        # Select point from feature class using its ID
        where = "%s = %s" % (OID, oid)
        pnt = arcpy.Point(x, y)  # Create a point object for the location
        ptgeo = arcpy.PointGeometry(pnt)  # Create a PointGeometry object

        # Extract elevation data within the buffer zone for the point
        tbuf = ExtractByCircle(elev, pnt, buf)

        # Try to get the zone value at the point location from the zone raster (fz)
        try:
            # Get the value of the raster at the point (x, y) and convert it to integer
            zoneval = int(arcpy.GetCellValue_management(fz, "%s %s" % (x, y)).getOutput(0))
            print(f"OID: {oid}")
            print(f'-- zone min ele {zoneval} --')
        except:
            # If zone value cannot be determined, skip to the next point
            print("w", "Could not determine zone value")
            continue

        # Find the flow-to cell(s) by masking the elevation buffer with the zone raster
        # Only cells outside the current zone value are considered
        tmsk = Con(IsNull(SetNull(fz, 1, "VALUE <> %s" % zoneval)), 1)
        telev = ExtractByMask(tbuf, tmsk)  # Mask the buffer raster using the flow mask

        # Try to calculate the minimum elevation difference within the buffer
        try:
            # Set cells where the elevation difference is close to the minimum as flow cells
            tmincells = Con(Abs(telev - telev.minimum) < .001, 1)
        except:
            # If flow-to cells cannot be determined, skip to the next point
            print("w", "Could not determine flow-to point")
            continue

        # Calculate cost distance from the starting point to flow-to cells (using the elevation raster)
        tdist = CostDistance(ptgeo, elev, "", tback)  # Compute cost distance
        # Alternatively, distance accumulation can be used (commented out)
        # tdist = DistanceAccumulation(ptgeo, elev, "", "", "", "", "", tback)

        # Calculate the optimal flow path using the cost distance raster
        tpath = CostPath(tmincells, tdist, Raster(tback), "EACH_ZONE")  # Compute the flow path

        # Assign the zone value to all cells in the flow path raster
        tpath1 = SetNull(IsNull(tpath), zoneval)  # Set cells that are null to the zone value

        # Convert the raster-based flow path to polylines
        # RasterToPolyline_conversion takes the raster and outputs polyline geometries
        arcpy.RasterToPolyline_conversion(tpath1, tline, "NODATA", "0", "NO_SIMPLIFY")  # Convert to polylines

        # Append the generated polyline to the existing feature class (output_lines)
        arcpy.Append_management(tline, output_lines, "NO_TEST")  # Append the line to output feature class

        # Clean up by deleting the temporary polyline feature class (tline)
        arcpy.Delete_management(tline)

        # Stop after processing 10 points
        #if count == 10:
        #    break
        count += 1
    print('Final.') #zoneval is min ele

if __name__ == "__main__":
    arcpy.env.addOutputsToMap = False
    name = "Extract_dem12" # input
    
    # ---- Define file paths ----
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    hold = os.path.join(project_folder, "resfolder")
    gdb = os.path.join(project_folder, 'lidar_1.gdb')
    name_min = f"min_points_{name}.shp"
    in_point = os.path.join(hold, name_min)
    fill_raster = os.path.join(hold, f'raster_{name}.tif')
    elev_raster = os.path.join(gdb, 'Extract_dem12')
    dgdb = os.path.join(project_folder, "Default.gdb")
    output_lines = os.path.join(dgdb, 'upload_path')
    """
    project = arcpy.mp.ArcGISProject("CURRENT")
    project_folder = os.path.dirname(project.filePath)
    hold = os.path.join(project_folder, "resfolder")
    gdb = os.path.join(project_folder, 'lidar_1.gdb')
    name_min = f"min_points_{name}.shp"
    in_point = os.path.join(hold, name_min)
    fz = os.path.join(hold, f'raster_{name}.tif')
    elev_raster = os.path.join(gdb, 'Extract_dem12')
    dgdb = os.path.join(project_folder, "Default.gdb")
    output_lines = os.path.join(dgdb, 'upload_path')
    fillzone_raster = os.path.join(hold, f"raster_{name}.tif")
    """
    # ---- Call functions ----
    min_area = 350 # Count of neighboring cells to select by
    min_area_process_lidar = 10000 # Not needed
    buf = 25 # Buffer zone on each min cell
    
    process_lidar_data()
    calculate_pour_sink()
    draw_flowpath()
    
    print('Done - Final.')
    
    ## Further process points that are on both halves of the road.
    ## Buffer to get viscinity
    ## Check why not processing for all

