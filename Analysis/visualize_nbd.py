import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from nbd_utils_code import nbd_read  # Importing your provided code

# Define the directory containing the snapshots
snapshot_directory = r'/mnt/c/Users/shyam/Desktop/dissertation/code/bin'

# List all snapshot files
snapshot_files = sorted([f for f in os.listdir(snapshot_directory) if f.startswith('evolv_binary')])

# Read and plot each snapshot
for snapshot_file in snapshot_files:
    filepath = os.path.join(snapshot_directory, snapshot_file)
    data = nbd_read(filepath)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data.x, data.y, data.z, s=1)
    
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title(f"Snapshot at time {data.time:.2f}")
    
    plt.show()
