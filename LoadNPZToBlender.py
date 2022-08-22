
import bpy
import bmesh
import numpy as np

# Directory of the NPZ file to load into blender
#npzDir = 'D:/UniSydney/Research/20220726Tibet/BlenderFlatEarth/NPZData/HimaHighRes.npz'
#npzDir = 'D:/UniSydney/Research/20220726Tibet/BlenderFlatEarth/NPZData/HimaMedRes.npz'
npzDir = 'D:/UniSydney/Research/20220726Tibet/BlenderFlatEarth/NPZData/HimaLowRes.npz'

# After loading the mesh into blender, this will be the name of the imported mesh
meshName = 'Landscape'

# Scaling of coordinates (so that our mesh is not rediculously large in blender)
lonLatScaling = 1e-1
elevsScaling = 1e-5
includeElevations = True

# Load data
npzFile = np.load(npzDir)
lon = npzFile["lon"] * lonLatScaling
lat = npzFile["lat"] * lonLatScaling
elevs = npzFile["elevs"] * elevsScaling
faces = npzFile["faces"]
edges = []

if includeElevations:
    XYZ = np.stack((lon, lat, elevs), axis=1)
else:
    zero = np.zeros(lon.shape)
    XYZ = np.stack((lon, lat, zero), axis=1)

# Create a new mesh object, load data into it, and link the mesh to the current scene
mesh = bpy.data.meshes.new(meshName)
mesh.from_pydata(list(XYZ), edges, list(faces))
obj = bpy.data.objects.new(meshName, mesh) 
bpy.context.scene.collection.objects.link(obj)
bpy.context.view_layer.objects.active = obj

# Get the active mesh
bm = bmesh.new()
selectedMesh = bpy.context.object.data
bm.from_mesh(selectedMesh)

# Create a color map
color = elevs - np.min(elevs)
color = color / (np.max(elevs) - np.min(elevs))
colorMap = {v : [color[i], color[i], color[i], 1] for i, v in enumerate(bm.verts)}

# Create a new color layer and apply colormap to it
color_layer = bm.loops.layers.color.new("color")
for face in bm.faces:
    for loop in face.loops:
        loop[color_layer] = colorMap[loop.vert]

# Apply changes made to selected mesh
bm.to_mesh(selectedMesh) 


