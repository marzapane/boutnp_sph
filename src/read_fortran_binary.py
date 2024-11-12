# Function to read Fortran unformatted binary files
def readfbin(file_path):
    with open(file_path, 'rb') as f:
        # Skip the record marker (4 bytes), read the first scalar (1) and skip the marker
        f.read(4)
        scalar = np.frombuffer(f.read(4), dtype=np.int32)
        f.read(4)
        # Read input dimension
        f.read(4)
        n1 = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n2 = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n3 = np.frombuffer(f.read(4), dtype=np.int32)[0]
        f.read(4)
        # Read epsil, re and time variables
        f.read(4)
        metadata = np.frombuffer(f.read(32), dtype=np.float64)
        f.read(4)
        # Now read the actual data block (all fields in order: E, vx, vy, vz, pr)
        f.read(4)
        vector_size = n1 * n2 * n3
        fields = np.frombuffer(f.read(vector_size *5 *8), dtype=np.float64)
        f.read(4)
        # Now reshape the data array into 5 separate arrays
        fields = fields.reshape(5, n1, n2, n3)

    with open(file_path, 'rb') as f:
        # Skip the record marker (4 bytes), read the first scalar (1) and skip the marker
        f.read(4)
        scalar = np.frombuffer(f.read(4), dtype=np.int32)
        f.read(4)
        # Read input dimension
        f.read(4)
        n1g = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n2g = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n3g = np.frombuffer(f.read(4), dtype=np.int32)[0]
        f.read(4)
        if (n1g==n1 or n2g==n2 or n3g==n3):
            # Read coordinates values
            f.read(4)
            coordinates = np.frombuffer(f.read(vector_size *3 *8), dtype=np.float64)
            f.read(4)
            coordinates = coordinates.reshape(3, n1, n2, n3)
        else:
            raise Exception("Field size incompatible with grid size")
        
    return {
        'n1':   n1,
        'n2':   n2,
        'n3':   n3,
        'x':    coordinates[0],
        'y':    coordinates[1],
        'z':    coordinates[2],
        'E':    fields[0],
        'vx':   fields[1],
        'vy':   fields[2],
        'vz':   fields[3],
        'pr':   fields[4],
    }

(n1,n2,n3,xc,yc,zc,E,vx,vy,vz,pr) = readfbin().values()


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Example data: replace with actual arrays
n1, n2, n3 = vx.shape  # Assume vx, xc, yc, zc have the same shape

# Create a 3D scatter plot of vx at given (x, y, z) coordinates
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Flatten the arrays for plotting
xc_flat = xc.flatten()
yc_flat = yc.flatten()
zc_flat = zc.flatten()
vx_flat = vx.flatten()

# Scatter plot with color representing vx values
sc = ax.scatter(xc.flatten(), yc.flatten(), zc.flatten(), c=vx.flatten(), cmap='viridis', marker='o')

# Add color bar to show the magnitude of vx
plt.colorbar(sc, ax=ax, label='v')

# Set labels
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Set title
# plt.title('3D Scatter plot of vx')

# Show plot
plt.show()

