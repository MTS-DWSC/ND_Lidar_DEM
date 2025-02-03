~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
README: ArcPy Workflow for Geospatial Data Processing

Overview

This workflow is designed to process geospatial data using ArcPy within an ArcGIS Pro environment. It involves creating and managing spatial data, performing spatial joins, converting coordinate systems, and generating output shapefiles and feature classes. The following sections provide a step-by-step guide to understanding and executing the script.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Script Breakdown

This workflow provides a comprehensive approach to geospatial data processing in ArcGIS Pro using ArcPy. 
Proper setup and execution of the script will result in processed spatial data ready for further analysis or visualization.

The script uses the Department of Water Resource's LiDAR data for DEM rasters that are then converted into hillshades.
* A buffer zone around the points for each ID in the database is used as a mask.
* An Polyline is drawn from each half of the mask to areas of lowest elevation.
  * This processes each point then moves on to the next one.
* The Lat, Long of each half of the mask is stored then compared to the central point where the closest position on each half to the central point is recorded
* Final pass to the Polyline actually draws the line, then it is recorded in a database for future reference before prcoessing the next one.
