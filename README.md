# HimalayasFuture

A project where we are predicting the future of the himalayas. The majority of scripts are intended to run within blender.

# Creating NPZ Files for Blender

Using the scripts *GeneratingFlatEarthsNPZ.py*, we can generate NPZ files that are suitable to be loaded into blender by our *LoadNPZToBlender.py* scripts. To do so, download and extract *ETOPO1_Ice_g_geotiff.tif* from [this website](https://www.ngdc.noaa.gov/mgg/global/#:~:text=ETOPO1%20is%20a%201%20arc,base%20of%20the%20ice%20sheets\).), change the directories within the code appropriately and run it.

# Load Mesh into Blender

Once you have created an NPZ file, copy and paste the code from *LoadNPZToBlender.py* into blender (change the directories appropriately). A mesh representing the Himalayas should appear in blender.

# Move Tectonic Plates

You can run the *WarpFlatMesh.py* to move tectonic plates. You can change the parameter *time* parameter to select what time you want to read GPlates velocity data from.

# Save Data

Runnning the code *ExportData.py* in blender will export the data as an NPZ file containing lon/lat coordinates and the elevation of the topography at those coordinates.