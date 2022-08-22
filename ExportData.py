
import bpy
import bmesh
import numpy as np

# Change this to your desired output directory
saveDir = 'D:/UniSydney/Research/20220726Tibet/BlenderFlatEarth/Output/testOutput.npz'

# Scaling of coordinates (needs to be the same as those used to import mesh into blender)
lonLatScaling = 1e-1
elevsScaling = 1e-5

# Get XYZ coordinates of vertices of currently selected mesh
bm = bmesh.new()
topoMesh = bpy.context.active_object
bm.from_mesh(topoMesh.data)
XYZ = np.array([[v.co.x, v.co.y, v.co.z] for v in bm.verts])
bm.free()

# Undo the scaling done to our coordinates
lon = XYZ[:, 0] / lonLatScaling
lat = XYZ[:, 1] / lonLatScaling
elevations = XYZ[:, 2] / elevsScaling

# Save data as NPZ file
np.savez(saveDir, lon=lon, lat=lat, elevs=elevations)

