
import numpy as np
import rasterio
import pyvista as pv

class Topography:
    def __init__(self, 
                    # Specifies which parts of earth to include
                    latLimits = np.array([10, 60]),
                    lonLimits = np.array([50, 130]),

                    # Data file directories and other parameters
                    topoDir = 'D:/UniSydney/Research/20220726Tibet/ETOPO1_Ice_g_geotiff.tif',
                    rawDataDims = np.array([10801, 21601]),
                    earthRadius = 6.371, # Thousand km
                    resolutionReductionFactor = 15 # Not all numbers will work ()
                    ):
        
        # Store input parameters
        self.latLimits = latLimits 
        self.lonLimits = lonLimits 
        self.rawDataDims = (((rawDataDims-1)/resolutionReductionFactor) + 1).astype(int)
        self.topoDir = topoDir
        self.earthRadius = earthRadius
        self.resolutionReductionFactor = resolutionReductionFactor

        # Get indices for where to slice dataset
        self.latIdx, self.lonIdx, self.sliceDims = self.getSliceIndices()

        # Load data, all arrays are 1D of same length
        self.elevations = self.readTopoTIF()
        self.lon, self.lat = self.getLonLatCoords()
        self.flatSphereXYZ = polarToCartesian(earthRadius, self.lon, self.lat)
    
    # Use lon/lat limits to calculate slice indices
    def getSliceIndices(self):
        datDims = self.rawDataDims
        latLim = self.latLimits
        lonLim = self.lonLimits
        lat = np.linspace(90, -90, datDims[0])
        lon = np.linspace(-180, 180, datDims[1])
        latIdx = np.array([
            np.argwhere(lat==latLim[1])[0, 0],
            np.argwhere(lat==latLim[0])[0, 0]+1])
        lonIdx = np.array([
            np.argwhere(lon==lonLim[0])[0, 0],
            np.argwhere(lon==lonLim[1])[0, 0]+1])
        sliceDims = np.array([
            latIdx[1] - latIdx[0],
            lonIdx[1] - lonIdx[0]])
        return latIdx, lonIdx, sliceDims
    
    # Returns np.array of sliced topography data
    def readTopoTIF(self):
        topotTIF = rasterio.open(self.topoDir)
        topo = topotTIF.read()
        topotTIF.close()
        topo = topo[0, ::self.resolutionReductionFactor, ::self.resolutionReductionFactor]
        topo = topo[self.latIdx[0]:self.latIdx[1], self.lonIdx[0]:self.lonIdx[1]]
        return topo.flatten()
    
    # Create lonLat coordinates for elevations in topo
    def getLonLatCoords(self):
        dims = self.sliceDims
        lat = np.linspace(self.latLimits[1], self.latLimits[0], dims[0])
        lon = np.linspace(self.lonLimits[0], self.lonLimits[1], dims[1])
        lonLatGrid = np.meshgrid(lat, lon)
        lonLatGrid = np.array(lonLatGrid)
        lonLatGrid = lonLatGrid.T.reshape(dims[0] * dims[1], 2)
        lon = lonLatGrid[:, 1]
        lat = lonLatGrid[:, 0]
        return lon, lat
    
    # Create XYZ coordinates of flat earth suitable for plotting
    def getFlatEarthXYZ(self, amplificatinFactor=1e-4):
        elevs = amplificatinFactor * self.elevations
        XYZ = np.stack((self.lon, self.lat, elevs), axis=1)
        return XYZ
    
    # Create XYZ coordinates of spherical earth suitable for plotting
    def getSphericalXYZ(self, amplificatinFactor=1e-5):
        elevs = self.elevations
        r = self.earthRadius + elevs * amplificatinFactor
        return polarToCartesian(r, self.lon, self.lat)
    
    # Get mesh of earth (either flat or round)
    def getEarthMesh(self, sphericalMesh=True, amplificatinFactor=1e-5):
        mesh = pv.Plane(i_resolution=self.sliceDims[1]-1, j_resolution=self.sliceDims[0]-1)
        if sphericalMesh:
            mesh.points = self.getSphericalXYZ(amplificatinFactor=amplificatinFactor)
        else:
            mesh.points = self.getFlatEarthXYZ(amplificatinFactor=amplificatinFactor)
        mesh['elevations'] = self.elevations
        return mesh
    
# ================================================ Spherical Coordinate Transformations =============================
# Coordinate transformation function between polar and cartesian
def polarToCartesian(radius, theta, phi, useLonLat=True):
    if useLonLat == True:
        theta, phi = np.radians(theta+180.), np.radians(90. - phi)
    X = radius * np.cos(theta) * np.sin(phi)
    Y = radius * np.sin(theta) * np.sin(phi)
    Z = radius * np.cos(phi)
    
    #Return data either as a list of XYZ coordinates or as a single XYZ coordinate
    if (type(X) == np.ndarray):
        return np.stack((X, Y, Z), axis=1)
    else:
        return np.array([X, Y, Z])

# Coordinate transformation function between polar and cartesian
def cartesianToPolarCoords(XYZ, useLonLat=True):
    X, Y, Z = XYZ[:, 0], XYZ[:, 1], XYZ[:, 2]
    R = (X**2 + Y**2 + Z**2)**0.5
    theta = np.arctan2(Y, X)
    phi = np.arccos(Z / R)
    
    #Return results either in spherical polar or leave it in radians
    if useLonLat == True:
        theta, phi = np.degrees(theta), np.degrees(phi)
        lon, lat = theta - 180, 90 - phi
        lon[lon < -180] = lon[lon < -180] + 360
        return R, lon, lat
    else:
        return R, theta, phi

# Given an axis of rotation and angle, create a quaternion
def quaternion(axis, angle):
    return [np.sin(angle/2) * axis[0], 
            np.sin(angle/2) * axis[1], 
            np.sin(angle/2) * axis[2], 
            np.cos(angle/2)]

# Takes pyvista mesh and saves NPZ for blender
def saveForBlender(pyvistaMesh, saveDir='Mesh.npz'):
    XYZ = pyvistaMesh.points
    faces = pyvistaFacesToBlenderFaces(pyvistaMesh.faces)
    np.savez(saveDir, XYZ=XYZ, faces=faces)

# Converts pyvista faces to blender faces
def pyvistaFacesToBlenderFaces(pvFaces):
    faces = []
    i = 0
    while (i < len(pvFaces)):
        a = pvFaces[i]
        i += 1
        face = []
        for j in range(a):
            face.append(pvFaces[i])
            i += 1
        faces.append(face)
    return np.array(faces)


# Parameters you might want to change
resolutionReductionFactor = 15 # Change this to 1 (High res) 4 (Medium res) 15 (Low res)

# Change this to where you have saved the ETOPO1 file
topoDir = 'D:/UniSydney/Research/20220726Tibet/ETOPO1_Ice_g_geotiff.tif'

# Change this to where you would like to save the data to
saveDir = "D:/UniSydney/Research/20220726Tibet/Github/HimalayasFuture/HimalayasFuture/NPZData/HimaLowRes.npz"

# Create NPZ file of topography and save it
topo = Topography(
    resolutionReductionFactor=resolutionReductionFactor,
    topoDir = topoDir
)

# Get faces for blender mesh
mesh = topo.getEarthMesh(sphericalMesh=False)
faces = pyvistaFacesToBlenderFaces(mesh.faces)

# Save data as NPZ
np.savez(saveDir, lon=topo.lon, lat=topo.lat, elevs=topo.elevations, faces=faces)
