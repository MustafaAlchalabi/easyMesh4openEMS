import numpy as np
from CSXCAD.Utilities import GetMultiDirs
from .decorateOriginalMethods import CSXWrapper, FDTDWrapper
import easyMesh.geometryUtils as geometryUtils
import easyMesh.meshingUtils as meshingUtils
import easyMesh.meshProcessing as meshProcessing

def GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup=None, properties_mesh_setup=None, **kw):
    automesher = Automesher()
    return automesher.GenMesh(CSX, global_mesh_setup, primitives_mesh_setup, properties_mesh_setup, **kw)

def enhance_csx_for_auto_mesh(original_csx, primitives_mesh_setup):
    return CSXWrapper(original_csx, primitives_mesh_setup)

def enhance_FDTD_for_auto_mesh(original_FDTD, primitives_mesh_setup):
    return FDTDWrapper(original_FDTD, primitives_mesh_setup)

class Automesher:
    def __init__(self):
        self.properties_mesh_setup = {}
        self.primitives_mesh_setup = {}
        self.global_mesh_setup = {} 
        self.mesh_data = {}
        self.mesh_res = None
        self.min_cellsize = None
        self.max_res = None
        self.max_cellsize = None
        self.num_lines = None
        self.min_cellsize_changed = False
        self.found_circles = []
        self.mesh_with_max_cell_size = [[], [], []]  # x, y, z
        self.mesh_with_different_mesh_res = []
        self.metal_edge_res = None

    def GenMesh(self, CSX, global_mesh_setup, primitives_mesh_setup=None, properties_mesh_setup=None, **kw):

        self.properties_mesh_setup = properties_mesh_setup if properties_mesh_setup is not None else {}
        self.primitives_mesh_setup = primitives_mesh_setup if primitives_mesh_setup is not None else {}
        self.global_mesh_setup = global_mesh_setup

        self.csx = CSX
        grid = self.csx.GetGrid()
        # unique_primitives = CSX.GetAllPrimitives()
        unique_primitives = list(self.primitives_mesh_setup.keys())
        # unique_primitives.extend(list(self.primitives_mesh_setup.keys()))

        only_edges = []
        for primitive, mesh_setup in self.primitives_mesh_setup.items():
            if mesh_setup.get('edges_only', False):
                unique_primitives.remove(primitive)
                only_edges.append(primitive)
        if only_edges:
            if len(only_edges) == 1:
                self.collect_mesh_data(only_edges[0], grid, **kw)
            else:
                self.collect_mesh_data_for_multiple_primitives(only_edges, grid, **kw)

        if len(unique_primitives) == 1:
            self.collect_mesh_data(unique_primitives[0], grid, **kw)
        else:
            self.collect_mesh_data_for_multiple_primitives(unique_primitives, grid, **kw)

        self.add_mesh_lines(grid)

    def collect_mesh_data_for_multiple_primitives(self, primitives, grid, **kw):
        (mesh_data,dirs,metal_edge_res) = self.get_mesh_data(primitives, grid, **kw)
        if mesh_data is not None:
            self.mesh_data[tuple(primitives)] = (mesh_data,dirs,metal_edge_res)

    def collect_mesh_data(self, primitive, grid, **kw):
        (mesh_data,dirs,metal_edge_res) = self.get_mesh_data(primitive, grid, **kw)
        if mesh_data is not None:
            self.mesh_data[primitive] = (mesh_data,dirs,metal_edge_res)

    def add_mesh_lines(self, grid):
        for primitive in self.mesh_data:
            for n in range(3):
                if self.mesh_data[primitive][0][n]:
                    grid.AddLine(n, self.mesh_data[primitive][0][n])

    def get_mesh_data(self, polygon, grid, **kw):

        # Initialize the mesh_data list for x, y, and z directions
        mesh_data = [[], [], []]

        # Initialize lists to store edges and coordinates
        x_edges, y_edges, diagonal_edges = [], [], []
        x_coords, y_coords, z_coords = [], [], []

        # Check if there are manually added mesh lines in grid 
        meshingUtils.check_grid(self, grid, mesh_data)

        mesh_map = meshingUtils.get_mesh_map(self)

        geometryUtils.process_polygon(self, polygon, x_coords, y_coords, z_coords, x_edges, y_edges, diagonal_edges, mesh_data)

        # Get the mesh parameters
        meshingUtils.get_mesh_parameters(self)

        # Get unique vertical and horizontal edges
        unique_xedges, unique_yedges = geometryUtils.get_unique_edges(x_edges), geometryUtils.get_unique_edges(y_edges)

        meshingUtils.adjust_mesh_parameters(self, unique_xedges, unique_yedges, z_coords, diagonal_edges, mesh_data)
        
        # Sort x and y coordinates
        sorted_x, sorted_y = np.sort(x_coords), np.sort(y_coords)        
        
        # Sort edges by their starting coordinates
        x_edges.sort(key=lambda edge: edge[0])
        y_edges.sort(key=lambda edge: edge[0]) 
        z_coords = [(z[0], None, None, z[1], False) for z in z_coords]
        z_coords.sort(key=lambda edge: edge[0])

        geometryUtils.metal_edge(self, x_edges, x_coords, y_coords, mesh_data[0], 'x')
        geometryUtils.metal_edge(self, y_edges, x_coords, y_coords, mesh_data[1], 'y')

        # Add graded meshlines at material transitions
        meshingUtils.add_graded_mesh_lines_at_material_transitions(self, x_edges, mesh_data[0], mesh_map[0], 'x')
        meshingUtils.add_graded_mesh_lines_at_material_transitions(self, y_edges, mesh_data[1], mesh_map[1], 'y')
        meshingUtils.add_graded_mesh_lines_at_material_transitions(self, z_coords, mesh_data[2], mesh_map[2], 'z')

        # Handle diagonal edges and add mesh lines for x and y directions
        meshingUtils.handle_diagonal_edges(self, diagonal_edges, x_edges, y_edges, mesh_data[0], 'x')
        meshingUtils.handle_diagonal_edges(self, diagonal_edges, x_edges, y_edges, mesh_data[1], 'y')

        # Refine mesh in small gaps for x and y edges
        meshingUtils.mesh_small_gaps(self, x_edges, mesh_data, 'x')
        meshingUtils.mesh_small_gaps(self, y_edges, mesh_data, 'y')
        meshingUtils.mesh_small_gaps(self, z_coords, mesh_data, 'z')
        
        # Recalculate unique edges after refinement 
        unique_xedges, unique_yedges = geometryUtils.get_unique_edges(x_edges), geometryUtils.get_unique_edges(y_edges)

        # Add missing mesh lines between points and edges for x and y directions
        meshingUtils.add_missing_mesh_lines(self, unique_xedges, sorted_x, diagonal_edges, mesh_data[0], 'x')
        meshingUtils.add_missing_mesh_lines(self, unique_yedges, sorted_y, diagonal_edges, mesh_data[1], 'y')
        
        # Sort edges by their starting coordinates
        x_edges.sort(key=lambda edge: edge[0])
        y_edges.sort(key=lambda edge: edge[0])   

        meshingUtils.add_edges_to_mesh_mesh_data(self, mesh_data[0], x_edges, 'x')
        meshingUtils.add_edges_to_mesh_mesh_data(self, mesh_data[1], y_edges, 'y')
        meshingUtils.add_edges_to_mesh_mesh_data(self, mesh_data[2], z_coords, 'z')

        # Handle circular segments in the polygon
        geometryUtils.detect_all_circles_in_polygon(self, polygon)
        meshingUtils.handle_circular_segments(self, polygon, mesh_data)

        # Add Ports to the mesh_data
        meshingUtils.add_ports_to_mesh_data(self, mesh_data[0], x_edges, 'x')
        meshingUtils.add_ports_to_mesh_data(self, mesh_data[1], y_edges, 'y')
        meshingUtils.add_ports_to_mesh_data(self, mesh_data[2], z_coords, 'z')

        # Smooth and process the mesh lines
        meshProcessing.smooth_and_process_mesh_lines(self, mesh_data, polygon, grid, x_edges, y_edges, z_coords, unique_xedges, unique_yedges, z_coords, mesh_map) 

        # If no z-direction mesh_datas exist, set it to None
        if mesh_data[2] == []:
            mesh_data[2] = None

        # Initialize the final mesh_data with None for all directions
        final_mesh_data = [None, None, None]

       # Assign mesh_datas to the appropriate directions based on global setup 
        dirs = self.global_mesh_setup.get('dirs','xyz')
        if dirs is not None:
            for ny in GetMultiDirs(dirs):
                final_mesh_data[ny] = mesh_data[ny] 

        metal_edge_res = None

        return (final_mesh_data, dirs, metal_edge_res)
    