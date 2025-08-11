# Import necessary libraries and modules
import os 
import tempfile 
from pylab import * 
from CSXCAD import ContinuousStructure
from openEMS import openEMS  
from openEMS.physical_constants import * 

# Add a custom path for easyMesher
sys.path.append('/home/opt/easyMesh4openEMS')
# Import functions for automatic mesh generation and optimization
from easyMesher import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh

# Define the simulation path and global parameters
Sim_Path = os.path.realpath(os.path.join('.', 'Patch_Antenna'))  # Path to save simulation results
post_proc_only = False  # If True, only post-processing will be performed
unit = 1e-6  # drawing unit (um)

# Patch antenna parameters
patch_length = 20000  
patch_thickness = 10  
substrate_thickness = 635  
substrate_epr = 2.2  
substrate_length = 40000 

# Frequency parameters
f_min = 1.25e9 
f_max = 2.6e9  

# Initialize the FDTD simulator
FDTD = openEMS(EndCriteria=1e-4)           # Set the simulation end criteria at -40 dB
FDTD.SetGaussExcite(f_max / 2, f_max / 2)  # Gaussian excitation with center frequency
# Set boundary conditions: PML (Perfectly Matched Layer) and PEC (Perfect Electric Conductor)
FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PEC', 'PML_8'])

# Create the ContinuousStructure (CSX) for geometry definition
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)  
mesh = CSX.GetGrid()  
mesh.SetDeltaUnit(unit)  

# Mesh setup parameters
primitives_mesh_setup = {}  
properties_mesh_setup = {}  
global_mesh_setup = {
    'drawing_unit': unit, 
    'start_frequency': 0,                                                 
    'stop_frequency': f_max,                                              
    'smooth_metal_edge': 'one_third_two_thirds',                          # useful for thin metal layers, Options: False, 'one_third_two_thirds', 'extra_lines'
    'mesh_resolution': 'medium',                                          # Options: 'low', 'medium', 'high', 'very_high'
    'boundary_distance': ['auto', 'auto', 'auto', 'auto', None, 'auto'],  # Options: value, 'auto', or None
}

# Enhance the CSX and FDTD objects for automatic mesh optimization
CSX = enhance_csx_for_auto_mesh(CSX, primitives_mesh_setup)
FDTD = enhance_FDTD_for_auto_mesh(FDTD, primitives_mesh_setup)

# Add the substrate to the geometry
substrate = CSX.AddMaterial('RO5880', epsilon=substrate_epr) 
start = [-substrate_length, -substrate_length, 0]  
stop = [+substrate_length, +substrate_length, substrate_thickness]  
substrate.AddBox(start, stop, priority=10)  

# Add the patch to the geometry
pec = CSX.AddMetal('patch')  
pec.SetAttributeValue('setup', '5') 
start = [-patch_length, -patch_length, substrate_thickness]  
stop = [patch_length, patch_length, substrate_thickness]  
pec.AddBox(start, stop, priority=100)

# Add a lumped port to the FDTD simulation
portstart = [-7500, 0, 0]  
portstop = [-7500, 0, substrate_thickness]  
port = FDTD.AddLumpedPort(1, 50, portstart, portstop, 'z', 1, priority=100)

# Generate the mesh
GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup, properties_mesh_setup)

if 1: 
    CSX_file = os.path.join(Sim_Path, 'Patch_Antenna.xml')  
    if not os.path.exists(Sim_Path):  #
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)  
    from CSXCAD import AppCSXCAD_BIN  
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file)) 

# Run the simulation if post-processing is not enabled
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)  
    f = linspace(1e6, f_max, 1601)  

port.CalcPort(Sim_Path, f, ref_impedance=50)

s11 = port.uf_ref / port.uf_inc 

plot(f / 1e9, 20 * log10(abs(s11)), 'k-', linewidth=2, label='$S_{11}$')
grid() 
legend()  
ylabel('S-Parameter (dB)')  
xlabel('frequency (GHz)')  
show()  