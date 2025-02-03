~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
README: ArcPy Workflow for Geospatial Data Processing

Overview

This workflow is designed to process geospatial data using ArcPy within an ArcGIS Pro environment.
It involves creating and managing spatial data, performing spatial joins, converting coordinate systems,
and generating output shapefiles and feature classes. The following sections provide a step-by-step guide to understanding and executing the script.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Script Breakdown

This workflow provides a comprehensive approach to geospatial data processing in ArcGIS Pro using ArcPy. 
Proper setup and execution of the script will result in processed spatial data ready for further analysis or visualization.

#### Data Sources

The script utilizes LiDAR data from the Department of Water Resources (DWR) to generate Digital Elevation Models (DEMs). 
These DEMs are then converted into hillshades for further analysis.

#### Workflow Steps
Buffer Zone Creation
* A buffer zone is created around each point in the database. This buffer zone acts as a mask for subsequent processing steps.

Polyline Generation
* From each half of the buffer mask, a polyline is drawn to the areas of lowest elevation. This process is repeated for each point in the database.

Latitude and Longitude Storage
* The latitude and longitude of each half of the buffer mask are stored. These coordinates are then compared to the central point to determine the closest position on each half to the central point.

Final Polyline Drawing
* A final pass is made to draw the polyline based on the closest positions identified in the previous step. The polyline is then recorded in a database for future reference before moving on to the next point.
