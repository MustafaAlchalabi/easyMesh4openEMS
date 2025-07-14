from CSXCAD.SmoothMeshLines import SmoothMeshLines
from openEMS.physical_constants import C0
import numpy as np
import geometryUtils

def check_grid(automesher, grid, mesh_data):
    x, y, z = grid.GetLines(0), grid.GetLines(1), grid.GetLines(2)
    grid.ClearLines('x')
    grid.ClearLines('y')
    grid.ClearLines('z')

    if x.size > 0:
        mesh_data[0].extend(x)
    if y.size > 0:
        mesh_data[1].extend(y)
    if z.size > 0:
        mesh_data[2].extend(z)

def get_mesh_parameters(automesher):

        def get_mesh_res():

            fstart = automesher.global_mesh_setup.get('start_frequency', None)
            fstop = automesher.global_mesh_setup.get('stop_frequency', None)
            f0 = automesher.global_mesh_setup.get('f0', None)
            fc = automesher.global_mesh_setup.get('fc', None)
            unit = automesher.global_mesh_setup.get('drawing_unit', 1e-6)
            if fstart is not None and fstop is not None:
                automesher.wave_length = (C0/unit) / fstop
            elif f0 is not None and fc is not None:
                automesher.wave_length = (C0/unit) / (f0+fc)
            else:
                raise ValueError('Please provide start and stop frequency or f0 and fc in the global mesh setup')
            epsilon = 1
            mesh_res = automesher.global_mesh_setup.get('mesh_resolution', 'medium')
            if mesh_res == 'low':
                mesh_res = automesher.wave_length / (15 * epsilon**0.5)
                automesher.max_cellsize_air = automesher.wave_length / 15 
                num_lines = 4
            elif mesh_res == 'medium':
                mesh_res = automesher.wave_length / (20 * epsilon**0.5) 
                automesher.max_cellsize_air = automesher.wave_length / 20
                num_lines = 5
            elif mesh_res == 'high':
                mesh_res = automesher.wave_length / (25 * epsilon**0.5)
                automesher.max_cellsize_air = automesher.wave_length / 25
                num_lines = 6
            elif mesh_res == 'very_high':
                mesh_res = automesher.wave_length / (30 * epsilon**0.5)
                automesher.max_cellsize_air = automesher.wave_length / 30
                num_lines = 7
            else:
                mesh_res = automesher.wave_length / (20 * epsilon**0.5)   
                automesher.max_cellsize_air = automesher.wave_length / 20
                num_lines = 5
            # print('primitives_with_epsilon:', automesher.primitives_with_epsilon)
            return mesh_res, num_lines
        
        automesher.mesh_res = automesher.global_mesh_setup.get('refined_cellsize', None)
        # if automesher.mesh_res is not None:
        #     automesher.mesh_res = automesher.global_mesh_setup.get('mesh_resolution', None)
        if automesher.mesh_res is not None:
            automesher.num_lines = get_mesh_res()[1]
        else:
            automesher.num_lines = 5
        if automesher.mesh_res is None:
            automesher.mesh_res, automesher.num_lines = get_mesh_res()

        automesher.max_cellsize = automesher.global_mesh_setup.get('max_cellsize', get_mesh_res()[0])    
        print('automesher.mesh_res:', automesher.mesh_res)
        automesher.min_cellsize = automesher.global_mesh_setup.get('min_cellsize', automesher.mesh_res / 4)
        automesher.max_res = automesher.min_cellsize + 0.25 * automesher.min_cellsize

def adjust_mesh_parameters(automesher, unique_xedges, unique_yedges, diagonal_edges, mesh_data):
            print('min_cellsize:', automesher.min_cellsize)
            print('mesh_res:', automesher.mesh_res)
            print('num_lines:', automesher.num_lines)
            print('max_res:', automesher.max_res)
            print('max_cellsize:', automesher.max_cellsize)

            distance_smaller_than_min_cellsize = []
            automesher.min_cellsize_changed = False

            if automesher.global_mesh_setup.get('min_cellsize', None) is None:
                for i in range(len(unique_xedges) - 1):
                    condition1 = abs(unique_xedges[i + 1][0] - unique_xedges[i][0]) <= automesher.min_cellsize and abs(unique_xedges[i + 1][0] - unique_xedges[i][0]) >= 1.5
                    condition2 = abs(unique_xedges[i + 1][0] - unique_xedges[i][0]) > automesher.min_cellsize and abs(unique_xedges[i + 1][0] - unique_xedges[i][0]) < automesher.max_res
                    if condition1 or condition2:
                        distance_smaller_than_min_cellsize.append([abs(unique_xedges[i + 1][0] - unique_xedges[i][0]), unique_xedges[i][0], unique_xedges[i + 1][0]])

                if distance_smaller_than_min_cellsize:
                    n = min(distance_smaller_than_min_cellsize)
                    lines = np.linspace(n[1], n[2], automesher.num_lines)
                    ds = abs(np.min(np.diff(lines)))
                    if ds < automesher.min_cellsize:
                        automesher.min_cellsize = ds
                        automesher.min_cellsize_changed = True
                        automesher.mesh_res = round(ds * (automesher.num_lines - 1))
                        automesher.max_res = automesher.min_cellsize + 0.25 * automesher.min_cellsize
                        epsilon = None
                        for primitive in automesher.primitives_mesh_setup.keys():
                            if hasattr(primitive, 'GetProperty') and hasattr(primitive.GetProperty(), 'GetMaterialProperty') and primitive.GetProperty().GetMaterialProperty('epsilon') > 1:
                                current_epsilon = primitive.GetProperty().GetMaterialProperty('epsilon')
                                if epsilon is None or current_epsilon > epsilon:
                                    epsilon = current_epsilon
                        if epsilon:
                            automesher.max_cellsize = automesher.max_cellsize / (epsilon ** 0.5)
                # if automesher.max_cellsize == automesher.mesh_res:
                #     automesher.max_cellsize = automesher.mesh_res * 2

            if not automesher.min_cellsize_changed:
                check_max_resolution(automesher, diagonal_edges, unique_xedges, unique_yedges, automesher.mesh_res, automesher.max_res, mesh_data)
            if automesher.min_cellsize_changed:
                print('min_cellsize changed')
            print('min_cellsize:', automesher.min_cellsize)
            print('mesh_res:', automesher.mesh_res)
            print('num_lines:', automesher.num_lines)
            print('max_res:', automesher.max_res)
            print('max_cellsize:', automesher.max_cellsize)

def get_mesh_map(automesher):
    '''
    Returns a list containing x, y, and z boundaries with their respective properties.
    '''
    properties= automesher.csx.GetAllProperties()
    mesh_map = [[], [], []]  # x, y, z boundaries
    for prop in properties:
        if (hasattr(prop, 'GetMaterialProperty') and prop.GetMaterialProperty('epsilon') is not None) or prop.GetTypeString() == 'Metal':
            epsilon = prop.GetMaterialProperty('epsilon') if not prop.GetTypeString() == 'Metal' else 1.0
            primitives = prop.GetAllPrimitives()
            tmp_x,tmp_y,tmp_z,tmp_x_edges,tmp_y_edges, tmp_diagonal_edges, tmp_mesh_data = [], [], [], [], [], [], []
            if primitives:
                if isinstance(primitives, list):
                    # Process each primitive
                    for prim in primitives:
                        geometryUtils.process_primitive(prim, tmp_x, tmp_y, tmp_x_edges, tmp_y_edges, tmp_diagonal_edges)

                    # Collect z-coordinates from the polygon
                    tmp_z.extend(geometryUtils.collect_z_coordinates(primitives))

                else:
                    # If the polygon is a single primitive, get its x and y coordinates                            
                    prim_x_coords, prim_y_coords = primitives.GetCoords()[0], primitives.GetCoords()[1]
                    tmp_x.extend(prim_x_coords)
                    tmp_y.extend(prim_y_coords)
                    # Process the single polygon to extract edges and coordinates
                    geometryUtils.process_single_polygon(primitives, tmp_x, tmp_y, tmp_x_edges, tmp_y_edges, tmp_diagonal_edges)

                    # Collect z-coordinates from the single polygon
                    tmp_z.extend(geometryUtils.collect_z_coordinates([primitives]))
                tmp_z = [(z[0]) for z in tmp_z]

                z_boundaries = [min(tmp_z), max(tmp_z), epsilon, prop, None, None, None, None, None]
                x_boundaries = [min(tmp_x), max(tmp_x), epsilon, prop, tmp_x, tmp_x_edges, tmp_y, tmp_y_edges, z_boundaries]
                y_boundaries = [min(tmp_y), max(tmp_y), epsilon, prop, tmp_y, tmp_y_edges, tmp_x, tmp_x_edges, z_boundaries]
                mesh_map[0].extend([x_boundaries])
                mesh_map[1].extend([y_boundaries])
                mesh_map[2].extend([z_boundaries])

    return mesh_map

def handle_otheredges(automesher, otheredges, unique_xedges, unique_yedges, mesh_res, max_res, mesh_data, direction, check_max_resolution=False):
    diagonal_edges = otheredges
    if not check_max_resolution:
        other_edges_in_range = []
        if direction == 'x':
            unique_edges = remove_close_unique_edges(automesher, unique_xedges, max_res)
            unique_edges = [unique_xedges[0] for unique_xedges in unique_edges]

        if direction == 'y':
            unique_edges = remove_close_unique_edges(automesher, unique_yedges, max_res)
            unique_edges = [unique_yedges[0] for unique_yedges in unique_edges]
        for edge in diagonal_edges:            
            if direction == 'x':
                start , end = edge[0], edge[1]
                other_edges_in_range = [other_edge for other_edge in diagonal_edges if (start <= other_edge[0] <= end or start <= other_edge[1] <= end or start >= other_edge[0] >= end or start >= other_edge[1] >= end)]
            if direction == 'y':
                start , end = edge[2], edge[3]
                other_edges_in_range = [other_edge for other_edge in diagonal_edges if (start <= other_edge[2] <= end or start <= other_edge[3] <= end or start >= other_edge[2] >= end or start >= other_edge[3] >= end)]
            x_start, x_end, y_start, y_end, prim = edge
            lines_in_range = [direction for direction in unique_edges if start <= direction <= end or start >= direction >= end]

            if not other_edges_in_range:
                alpha = np.atan(abs((y_end-y_start))/abs((x_end-x_start)))
                resolution = mesh_res * np.cos(alpha)
                resolution = automesher.max_res if resolution < automesher.max_res else resolution
                if not lines_in_range:
                    lines=SmoothMeshLines([start, end], resolution)    
                if lines_in_range:
                    lines_in_range.extend([start, end])
                    lines=SmoothMeshLines(lines_in_range, resolution)
                mesh_data.extend(lines)
            if other_edges_in_range:
                alpha = np.round(np.rad2deg(np.atan(abs((y_end-y_start))/abs((x_end-x_start)))),2)
                if direction == 'x':
                    other_edges_in_range_here = [other_edge[0:2] for other_edge in other_edges_in_range]
                    resolution = mesh_res * np.cos(np.deg2rad(alpha))
                    other_edges_in_range_here_all_coords = [other_edge for other_edge in other_edges_in_range if (start <= other_edge[0] <= end or start <= other_edge[1] <= end or start >= other_edge[0] >= end or start >= other_edge[1] >= end)]
                if direction == 'y':
                    other_edges_in_range_here = [other_edge[2:4] for other_edge in other_edges_in_range]
                    resolution = mesh_res * np.sin(np.deg2rad(alpha))
                    other_edges_in_range_here_all_coords = [other_edge for other_edge in other_edges_in_range if (start <= other_edge[2] <= end or start <= other_edge[3] <= end or start >= other_edge[2] >= end or start >= other_edge[3] >= end)]
                min_line = np.min([start,end, np.min(np.min(other_edges_in_range_here))])
                max_line = np.max([start,end, np.max(np.max(other_edges_in_range_here))])
                alphas_in_range = [(np.round(np.rad2deg(np.atan(abs((edge[3] - edge[2])) / abs((edge[1] - edge[0])))),2), edge) for edge in other_edges_in_range_here_all_coords]
                lines_in_mesh_data_range = [direction for direction in mesh_data if min_line <= direction <= max_line] 
                for line in lines_in_mesh_data_range:
                    mesh_data.remove(line)
                if resolution < automesher.max_res:
                    resolution = automesher.max_res
                lines_in_range = [direction for direction in unique_edges if min_line < direction < max_line]
                if lines_in_range:
                    lines_in_range.extend([min_line, max_line])  
                    lines=SmoothMeshLines(lines_in_range, resolution)
                if not lines_in_range:
                    lines=SmoothMeshLines([min_line,max_line], resolution)
                mesh_data.extend(lines)
                if direction == 'x':
                    for alpha_val, edge in alphas_in_range:
                        if alpha_val > alpha:
                            lines_in_mesh_data_range = [line for line in mesh_data if edge[0] <= line <= edge[1] or edge[0] >= line >= edge[1]]
                            for line in lines_in_mesh_data_range:
                                mesh_data.remove(line)
                            resolution = mesh_res * np.cos(np.deg2rad(alpha_val))
                            resolution = automesher.max_res if resolution < automesher.max_res else resolution
                            xlines = SmoothMeshLines([edge[0], edge[1]], resolution)
                            mesh_data.extend(xlines)
                if direction == 'y':
                    for alpha_val, edge in alphas_in_range:
                        if alpha_val < alpha:
                            lines_in_mesh_data_range = [line for line in mesh_data if edge[2] <= line <= edge[3] or edge[2] >= line >= edge[3]]
                            for line in lines_in_mesh_data_range:
                                mesh_data.remove(line)
                            resolution = mesh_res * np.sin(np.deg2rad(alpha_val))
                            resolution = automesher.max_res if resolution < automesher.max_res else resolution
                            ylines = SmoothMeshLines([edge[2], edge[3]], abs(resolution))
                            mesh_data.extend(ylines)
    circ_segments = []
    do_sort = True
    for edge in diagonal_edges:
        circ_segments = geometryUtils.detect_all_circles_in_polygon(automesher, edge[4])
        if circ_segments:
            do_sort = False
    if do_sort:
        if direction == 'x':
            diagonal_edges.sort(key=lambda edge: edge[0])
        if direction == 'y':
            diagonal_edges.sort(key=lambda edge: edge[2])
    dist = []
    for i in range(len(diagonal_edges)):
        for j in range(i + 1, len(diagonal_edges)):
            line1 = diagonal_edges[i]
            line2 = diagonal_edges[j]
            dot_product = np.dot([line1[1] - line1[0], line1[3] - line1[2]], [line2[1] - line2[0], line2[3] - line2[2]])
            norm_product = np.linalg.norm([line1[1] - line1[0], line1[3] - line1[2]]) * np.linalg.norm([line2[1] - line2[0], line2[3] - line2[2]])
            cos_angle = np.clip(dot_product / norm_product, -1.0, 1.0)
            angle = np.round(np.rad2deg(np.arccos(cos_angle)), 2)
            if diagonal_edges[i][4].GetElevation() != diagonal_edges[j][4].GetElevation():
                continue
            if diagonal_edges[i][4] == diagonal_edges[j][4] and (angle not in [0, 180]):
                if not (np.isclose(angle, 0, atol=10) or np.isclose(angle, 180, atol=10)):
                    continue
                else:
                    if direction == 'x':
                        if not (line2[0] < line1[0] < line2[1] or  line2[0] > line1[0] > line2[1]
                        or line2[0] < line1[1] < line2[1] or  line2[0] < line1[1] < line2[1]
                        or line1[0] < line2[0] < line1[1] or  line1[0] > line2[0] > line1[1]
                        or line1[0] < line2[1] < line1[1] or  line1[0] < line2[1] < line1[1]):
                            continue
                    if direction == 'y':
                        if not (line2[2] < line1[2] < line2[3] or  line2[2] > line1[2] > line2[3]
                        or line2[2] < line1[3] < line2[3] or  line2[2] < line1[3] < line2[3]
                        or line1[2] < line2[2] < line1[3] or  line1[2] > line2[2] > line1[3]
                        or line1[2] < line2[3] < line1[3] or  line1[2] < line2[3] < line1[3]):
                            continue
                        
            if np.isclose(angle, 90, atol=1e-2):
                continue
            p1 = np.array([line1[0], line1[2]])  # (x1, y1)
            p2 = np.array([line1[1], line1[3]])  # (x2, y2)
            q1 = np.array([line2[0], line2[2]])  # (x1, y1)
            q2 = np.array([line2[1], line2[3]])  # (x2, y2)
            distance = geometryUtils.distance_between_segments(p1, p2, q1, q2)
            distance = [small_dist for small_dist in distance if small_dist[0] <= mesh_res]
            if distance:
                dist.append(distance)
            alpha = np.round(np.rad2deg(np.atan(abs((q2[1] - q1[1]) / abs((q2[0] - q1[0]))))), 2)
            if np.min(np.diff([line1[0:2], line2[0:2]])) > mesh_res or np.min(np.diff([line1[2:4], line2[2:4]])) > mesh_res:
                continue
    if dist and check_max_resolution:
        for dist0 in dist:
            # discretize the distance with the number of lines and calculate the minimal distance between the lines for the maximum resolution
            x = np.round(np.diff(np.linspace(0, np.min([item[0] for item in dist0]), automesher.num_lines)),1)
            if np.min(x) < automesher.max_res and np.min(x) > 0:
                automesher.max_res = np.min(x)

    if dist and not check_max_resolution:    
        for dist1 in dist:
            if direction == 'x':
                coords_of_p = [item[1][0] for item in dist1]
                resolution = mesh_res * np.cos(np.deg2rad(alpha))
                # start_and_end_points = [line1[0], line1[1], line2[0], line2[1]]
                start_and_end_points = [dist1[0][2][0], dist1[0][3][0], dist1[0][4][0], dist1[0][5][0]]
            if direction == 'y':
                coords_of_p = [item[1][1] for item in dist1]
                resolution = mesh_res * np.sin(np.deg2rad(alpha))
                # start_and_end_points = [line1[2], line1[3], line2[2], line2[3]]
                start_and_end_points = [dist1[0][2][1], dist1[0][3][1], dist1[0][4][1], dist1[0][5][1]]
            lines_in_range = [lines for lines in mesh_data if np.min(coords_of_p) <= lines <= np.max(coords_of_p)]
            for line in lines_in_range:
                mesh_data.remove(line)
            lines_before_min = [line for line in mesh_data if line < np.min(coords_of_p) and abs(line - np.min(coords_of_p)) < max_res]
            lines_after_max = [line for line in mesh_data if line > np.max(coords_of_p) and abs(line - np.max(coords_of_p)) < max_res]

            if lines_before_min and lines_after_max:
                min_line = min(min(lines_before_min), min(lines_after_max))
                max_line = max(max(lines_before_min), max(lines_after_max))
                lines = check_edges_in_range(unique_edges, min_line, max_line, start_and_end_points, mesh_data, max_res)
                mesh_data.extend(lines)                        
            elif lines_before_min:
                min_line = min(min(lines_before_min), np.min(coords_of_p))
                max_line = max(max(lines_before_min), np.max(coords_of_p))
                lines = check_edges_in_range(unique_edges, min_line, max_line, start_and_end_points, mesh_data, max_res)
                mesh_data.extend(lines)
            elif lines_after_max:
                min_line = min(min(lines_after_max), np.min(coords_of_p))
                max_line = max(max(lines_after_max), np.max(coords_of_p))
                lines = check_edges_in_range(unique_edges, min_line, max_line, start_and_end_points, mesh_data, max_res)
                mesh_data.extend(lines)
            else:
                min_line = np.min(coords_of_p)
                max_line = np.max(coords_of_p)
                lines = check_edges_in_range(unique_edges, min_line, max_line, start_and_end_points, mesh_data, max_res)
                mesh_data.extend(lines)

def check_edges_in_range(unique_edges, min_line, max_line, start_and_end_points, mesh_data, resolution):
    edges_in_range = [edge for edge in unique_edges if min_line < edge < max_line]
    if edges_in_range:
        lines_in_range = [line for line in mesh_data if min(min_line, max_line, min(start_and_end_points)) <= line <= max(min_line, max_line, max(start_and_end_points))]
        if lines_in_range:
            for line in lines_in_range:
                mesh_data.remove(line)
        edges_in_range.extend([min_line, max_line])
        edges_in_range.extend(start_and_end_points)  
        lines=SmoothMeshLines(edges_in_range, resolution/1)
    if not edges_in_range:
        lines_in_range = [line for line in mesh_data if min(min_line, max_line) <= line <= max(min_line, max_line)]
        if lines_in_range:
            for line in lines_in_range:
                mesh_data.remove(line)
        # edges_in_range.extend(start_and_end_points)
        edges_in_range.extend([min_line, max_line])  
        lines=SmoothMeshLines(edges_in_range, resolution/1)
    return lines

def add_ports_to_mesh_data(automesher, mesh_data, edges, direction):
    x, y = [], []
    zedges = []
    x_edges, y_edges = [], []
    diagonal_edges = []
    for prim in automesher.primitives_mesh_setup:
        if hasattr(prim, 'priority'):
            port_coords_x, port_coords_y, port_coords_z = geometryUtils.transfer_port_to_polygon(prim.start, prim.stop)
            # if prim.measplane_shift:
            #     print('prim.meas_plane_shift:', prim, prim.measplane_shift)
            x.extend(port_coords_x)
            y.extend(port_coords_y)
            geometryUtils.collect_edges(port_coords_x, port_coords_y, prim, x_edges, y_edges, diagonal_edges)
            z = [(z, None, None, prim) for z in port_coords_z]
            zedges.extend(z)
    if direction == 'x':
        if x_edges:
            edges.extend(x_edges)
    if direction == 'y':
        if y_edges:
            edges.extend(y_edges)
    if direction == 'z': 
        if zedges:
            edges.extend(zedges)
    for edge in edges:
        if hasattr(edge[3], 'priority'):
            mesh_data.append(edge[0])
            for line in mesh_data:
                if abs(line - edge[0]) < automesher.min_cellsize/2:
                    if any([e[0] == line and hasattr(e[3], 'priority') for e in edges]):
                        continue
                    elif not any([e[0] == line and hasattr(e[3], 'priority') for e in edges]):
                        if line in mesh_data:
                            mesh_data.remove(line)
    return mesh_data
    
def add_edges_to_mesh_mesh_data(automesher, mesh_data, edges, mesh_res, min_cellsize, direction):

    remove_close_edges(automesher, edges, direction)
    for edge in edges:
        dirs = automesher.primitives_mesh_setup.get(edge[3], {}).get('dirs') or \
                automesher.properties_mesh_setup.get(edge[3].GetProperty() if hasattr(edge[3], 'GetProperty') else None, {}).get('dirs') or \
                automesher.global_mesh_setup.get('dirs', 'xyz')
        if direction == 'x':
            if 'x' in dirs:              
                mesh_data.append(edge[0])
        if direction == 'y':
            if 'y' in dirs:
                mesh_data.append(edge[0])
        if direction == 'z':
            if 'z' in dirs:
                mesh_data.append(edge[0])
    # mesh_data.extend(edge[0] for edge in edges)

def remove_close_edges(automesher, edges, direction):

    edges_to_remove = []
    for i in range(len(edges) - 1):
        if abs(edges[i+1][0] - edges[i][0]) < automesher.min_cellsize and abs(edges[i+1][0] - edges[i][0]) > 0:
            if hasattr(edges[i][3], 'priority') and hasattr(edges[i + 1][3], 'priority'):
                if edges[i][3].priority == edges[i + 1][3].priority:
                    continue
            elif (getattr(edges[i][3], 'GetPriority', lambda: None)() or edges[i][3].priority) > (getattr(edges[i + 1][3], 'GetPriority', lambda: None)() or edges[i + 1][3].priority):
                edges_to_remove.append(edges[i + 1])
            elif (getattr(edges[i][3], 'GetPriority', lambda: None)() or edges[i][3].priority) < (getattr(edges[i + 1][3], 'GetPriority', lambda: None)() or edges[i + 1][3].priority):
                edges_to_remove.append(edges[i])
            else:
                print(f"\033[91mWarning: Detected closely spaced edges at ({direction} = {edges[i][0]} and {direction} = {edges[i+1][0]}, d{direction} = {abs(edges[i][0]-edges[i+1][0])}). This configuration may lead to prolonged simulation times. Consider assigning different priorities or modifying the structure to optimize performance.\033[0m")
                continue
                # if abs(edges[i][2]-edges[i][1]) < abs(edges[i+1][2]-edges[i+1][1]):
                #     edges_to_remove.append(edges[i])
                # else:
                #     edges_to_remove.append(edges[i + 1])
    edges_to_remove = list({edge[0]: edge for edge in edges_to_remove}.values())
    if edges_to_remove:
        edges_to_remove_first_elements = {edge[0] for edge in edges_to_remove}
        edges[:] = [edge for edge in edges if edge[0] not in edges_to_remove_first_elements]

def remove_close_unique_edges(automesher, unique_edges, max_res):

    unique_edges_to_remove = []
    for i in range(len(unique_edges) - 1):
        if abs(unique_edges[i+1][0] - unique_edges[i][0]) < automesher.min_cellsize:
            if hasattr(unique_edges[i][1], 'priority') and hasattr(unique_edges[i + 1][1], 'priority'):
                if unique_edges[i][1].priority == unique_edges[i + 1][1].priority:
                    continue
            elif (getattr(unique_edges[i][1], 'GetPriority', lambda: None)() or unique_edges[i][1].priority) > (getattr(unique_edges[i + 1][1], 'GetPriority', lambda: None)() or unique_edges[i + 1][1].priority):
                unique_edges_to_remove.append(unique_edges[i + 1])
            elif (getattr(unique_edges[i][1], 'GetPriority', lambda: None)() or unique_edges[i][1].priority) < (getattr(unique_edges[i + 1][1], 'GetPriority', lambda: None)() or unique_edges[i + 1][1].priority):
                unique_edges_to_remove.append(unique_edges[i])
            else:
                if unique_edges[i][0] < unique_edges[i + 1][0]:
                    unique_edges_to_remove.append(unique_edges[i])
                else:
                    unique_edges_to_remove.append(unique_edges[i+1])
    unique_edges_to_remove = list({edge[0]: edge for edge in unique_edges_to_remove}.values())
    if unique_edges_to_remove:
        unique_edges_to_remove_first_elements = {edge[0] for edge in unique_edges_to_remove}
        unique_edges[:] = [edge for edge in unique_edges if edge[0] not in unique_edges_to_remove_first_elements]
                            
    return unique_edges


def mesh_small_gaps(automesher, unique_edges, mesh_res, max_res, num_lines, mesh_data, direction):
    unique_edges.sort(key = lambda x: x[0])
    min_cellsize = automesher.global_mesh_setup.get('min_cellsize', None)
    use_num_lines = True if min_cellsize is None else False
    max_res_list= []
    if direction == 'z':
        for i in range(len(unique_edges) - 1):
            if abs(np.diff([unique_edges[i][0], unique_edges[i + 1][0]])) <= mesh_res and abs(np.diff([unique_edges[i][0], unique_edges[i + 1][0]])) >= max_res:
                z_in_range = [z for z in mesh_data[2] if unique_edges[i][0] <= z <= unique_edges[i + 1][0]]
                for z in z_in_range:
                    mesh_data[2] = list(mesh_data[2])  
                    mesh_data[2].remove(z)
                # if use_num_lines:
                #     new_max_res = np.diff(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], num_lines))
                #     new_max_res = np.max(new_max_res)
                #     # max_res_list.append(new_max_res)
                # else:
                #     new_max_res = max_res
                zlines = SmoothMeshLines([unique_edges[i][0], unique_edges[i + 1][0]], automesher.max_res)
                # if len(zlines) <= 4:
                #     mesh_data[2] = list(mesh_data[2]) 
                #     mesh_data[2].extend(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], 4))
                # else:
                mesh_data[2].extend(zlines)
    else:
        for i in range(len(unique_edges) - 1):
            if abs(np.diff([unique_edges[i][0], unique_edges[i + 1][0]])) <= mesh_res and abs(np.diff([unique_edges[i][0], unique_edges[i + 1][0]])) >= max_res and abs(np.diff([unique_edges[i][0], unique_edges[i + 1][0]])) >= 1.5:
                y1, y2 = unique_edges[i][1], unique_edges[i][2]
                y1_next, y2_next = unique_edges[i + 1][1], unique_edges[i + 1][2]
                if (y1 <= y1_next <= y2 or y1 >= y1_next >= y2 or
                    y1 <= y2_next <= y2 or y1 >= y2_next >= y2 or
                    y1_next <= y1 <= y2_next or y1_next >= y1 >= y2_next or
                    y1_next <= y2 <= y2_next or y1_next >= y2 >= y2_next):
                        # print('unique_edges[i], unique_edges[i + 1]:', unique_edges[i], unique_edges[i + 1])
                        if direction == 'x':
                            x_in_range = [x for x in mesh_data[0] if unique_edges[i][0] <= x <= unique_edges[i + 1][0]]
                            for x in x_in_range:
                                mesh_data[0].remove(x)
                            if use_num_lines:
                                new_max_res = np.diff(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], num_lines))
                                new_max_res = np.max(new_max_res)
                                max_res_list.append(new_max_res)
                            else:
                                new_max_res = max_res
                            xlines = SmoothMeshLines([unique_edges[i][0], unique_edges[i + 1][0]], new_max_res)
                            if len(xlines) <= 4:
                                mesh_data[0].extend(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], 4))
                            else:
                                mesh_data[0].extend(xlines)
                        elif direction == 'y':
                            y_in_range = [y for y in mesh_data[1] if unique_edges[i][0] <= y <= unique_edges[i + 1][0]]
                            for y in y_in_range:
                                mesh_data[1].remove(y)
                            if use_num_lines:
                                new_max_res = np.diff(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], num_lines))
                                new_max_res = np.max(new_max_res)
                                max_res_list.append(new_max_res)
                            else:
                                new_max_res = max_res
                            ylines = SmoothMeshLines([unique_edges[i][0], unique_edges[i + 1][0]], new_max_res)
                            if len(ylines) <= 4:
                                mesh_data[1].extend(np.linspace(unique_edges[i][0], unique_edges[i + 1][0], 4))
                            else:
                                mesh_data[1].extend(ylines)
    if max_res_list:
        new_min_cellsize = min(max_res_list) - 0.25*min(max_res_list)
        if new_min_cellsize < automesher.min_cellsize:
            automesher.min_cellsize = new_min_cellsize
            # automesher.min_cellsize_changed = True
            automesher.max_res = min(max_res_list)


def check_max_resolution(automesher, diagonal_edges, unique_xedges, unique_yedges, mesh_res, max_res, mesh_data):
    # x-direction
    handle_otheredges(automesher, diagonal_edges, unique_xedges, unique_yedges, mesh_res, max_res, mesh_data[0], 'x', check_max_resolution=True)
    # y-direction
    handle_otheredges(automesher, diagonal_edges, unique_xedges, unique_yedges, mesh_res, max_res, mesh_data[1], 'y', check_max_resolution=True)

def handle_circular_segments(automesher, polygon, mesh_data):

    circ_segments = automesher.found_circles

    if circ_segments:
        print(f"Gefundene Kreisabschnitte: {len(circ_segments)}")
        for i, (x_seg, y_seg) in enumerate(circ_segments):
            print(f"  Kreis #{i+1}:")
            print("    X:", x_seg)
            print("    Y:", y_seg)
    else:
        print("Kein kreisförmiger Abschnitt gefunden.")

    # delete mesh lines that are inside the circle
    if circ_segments:
        for x_seg, y_seg in circ_segments:
            # Remove x lines inside the circle
            mesh_data[0] = [line for line in mesh_data[0] if not (min(x_seg) < line < max(x_seg))]
            # Add x lines inside the circle
            mesh_data[0].extend(SmoothMeshLines([min(x_seg), max(x_seg)], automesher.max_res))
            # Remove y lines inside the circle
            mesh_data[1] = [line for line in mesh_data[1] if not (min(y_seg) < line < max(y_seg))]
            # Add y lines inside the circle
            mesh_data[1].extend(SmoothMeshLines([min(y_seg), max(y_seg)], automesher.max_res))

def add_missing_mesh_lines(automesher, unique_edges, sorted_points, diagonal_edges, mesh_res, mesh_data, direction):
    'Check if the first and last point are x or y edges, if not it adds the missing mesh lines between the point and the edge'
    # if unique_edges.size > 0:
    if unique_edges:
        if unique_edges[-1][0] < sorted_points[-1]:
            for other_edge in diagonal_edges:
                if direction == 'x':
                    start, end = other_edge[0], other_edge[1]
                if direction == 'y':
                    start, end = other_edge[2], other_edge[3]
                if start <= sorted_points[-1] <= end or end <= sorted_points[-1] <= start:
                    if abs(np.diff([unique_edges[-1][0], min(start, end)])) < mesh_res:
                        lines = np.linspace(unique_edges[-1][0], min(start, end), 5)[1:]
                    else:
                        lines = SmoothMeshLines([unique_edges[-1][0], min(start, end)], mesh_res)[1:]
                    mesh_data.extend(lines)
        if unique_edges[0][0] > sorted_points[0]:
            for other_edge in diagonal_edges:
                if direction == 'x':
                    start, end = other_edge[0], other_edge[1]
                if direction == 'y':
                    start, end = other_edge[2], other_edge[3]
                if start <= sorted_points[0] <= end or end <= sorted_points[0] <= start:
                    if abs(np.diff([unique_edges[0][0], max(start, end)])) < mesh_res:
                        lines = np.linspace(unique_edges[0][0], max(start, end), 5)[1:]
                    else:
                        lines = SmoothMeshLines([max(start, end), unique_edges[0][0]], mesh_res)
                    mesh_data.extend(lines)    

def add_graded_mesh_lines(automesher, start, end, start_res, target_cellsize, growth):
    # print(f"Adding graded mesh lines from {start} to {end} with start resolution {start_res}, target cell size {target_cellsize}, and growth factor {growth}")
    if start_res <= 0:
        start_res = automesher.min_cellsize
    if target_cellsize <= 0:
        target_cellsize = automesher.mesh_res
    lines = []
    if start < end:
        current = start
        lines.append(current)
        temp_cellsize = start_res
        while current + temp_cellsize < end:
            current += temp_cellsize
            lines.append(current)
            if temp_cellsize < target_cellsize:
                temp_cellsize *= growth
            else:
                temp_cellsize = target_cellsize

    else:
        current = start
        lines.append(current)
        temp_cellsize = start_res
        while current - temp_cellsize > end:
            current -= temp_cellsize
            lines.append(current)
            if temp_cellsize < target_cellsize:
                temp_cellsize *= growth
            else:
                temp_cellsize = target_cellsize

    return lines

def add_graded_mesh_lines_at_material_transitions(automesher, edges, mesh_data, mesh_res, mesh_map, direction):
    for i in range(len(edges) - 1):
        if abs(np.diff([edges[i][0], edges[i + 1][0]])) == 0 and edges[i][0] > np.min(edges[0][0]) and edges[i][0] < np.max(edges[-1][0]):
            if hasattr(edges[i][3], 'GetProperty') and hasattr(edges[i + 1][3], 'GetProperty'):
                if edges[i][3].GetProperty()!= edges[i + 1][3].GetProperty():
                    if hasattr(edges[i][3].GetProperty(),'GetMaterialProperty') and not hasattr(edges[i + 1][3].GetProperty(), 'GetMaterialProperty') or \
                        hasattr(edges[i + 1][3].GetProperty(),'GetMaterialProperty') and not hasattr(edges[i][3].GetProperty(), 'GetMaterialProperty'):
                        if hasattr(edges[i][3].GetProperty(),'GetMaterialProperty'):
                            epsilon = edges[i][3].GetProperty().GetMaterialProperty('epsilon')
                        elif hasattr(edges[i + 1][3].GetProperty(), 'GetMaterialProperty'):
                            epsilon = edges[i + 1][3].GetProperty().GetMaterialProperty('epsilon')
                        if not automesher.min_cellsize_changed:
                            target_size = automesher.max_cellsize_air / epsilon**0.5
                        else:
                            target_size = mesh_res
                        lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]-target_size, automesher.min_cellsize, target_size, 1.3)
                        mesh_data.extend(lines)
                        lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]+target_size, automesher.min_cellsize, target_size, 1.3)
                        mesh_data.extend(lines)

                    if hasattr(edges[i][3].GetProperty(),'GetMaterialProperty') and hasattr(edges[i + 1][3].GetProperty(), 'GetMaterialProperty'):
                        # print('edges[i], epsilon, edges[i + 1], epsilon:', edges[i], edges[i][3].GetProperty().GetMaterialProperty('epsilon'), edges[i + 1], edges[i + 1][3].GetProperty().GetMaterialProperty('epsilon'))
                        if edges[i][3].GetProperty().GetMaterialProperty('epsilon') != edges[i + 1][3].GetProperty().GetMaterialProperty('epsilon'):
                            epsilon = max(edges[i][3].GetProperty().GetMaterialProperty('epsilon'), edges[i + 1][3].GetProperty().GetMaterialProperty('epsilon'))
                        else:
                            continue
                            epsilon = edges[i][3].GetProperty().GetMaterialProperty('epsilon')
                        if not automesher.min_cellsize_changed:
                            target_size = automesher.max_cellsize_air / epsilon**0.5
                        else:
                            target_size = mesh_res
                        lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]-target_size, automesher.min_cellsize, target_size, 1.3)
                        mesh_data.extend(lines)
                        lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]+target_size, automesher.min_cellsize, target_size, 1.3)
                        mesh_data.extend(lines)
                        
        # if abs(np.diff([edges[i][0], edges[i + 1][0]])) != 0 and edges[i][0] > np.min(edges[0][0]) and edges[i][0] < np.max(edges[-1][0]):
        #     if hasattr(edges[i][3], 'GetProperty') and hasattr(edges[i + 1][3], 'GetProperty'):
        #         if edges[i][3].GetProperty() != edges[i + 1][3].GetProperty():
        #             if hasattr(edges[i][3].GetProperty(),'GetMaterialProperty') and not hasattr(edges[i + 1][3].GetProperty(), 'GetMaterialProperty') or \
        #                 hasattr(edges[i + 1][3].GetProperty(),'GetMaterialProperty') and not hasattr(edges[i][3].GetProperty(), 'GetMaterialProperty'):
        #                 if hasattr(edges[i][3].GetProperty(),'GetMaterialProperty'):
        #                     epsilon = edges[i][3].GetProperty().GetMaterialProperty('epsilon')
        #                 elif hasattr(edges[i + 1][3].GetProperty(), 'GetMaterialProperty'):
        #                     epsilon = edges[i + 1][3].GetProperty().GetMaterialProperty('epsilon')
        #                 if not automesher.min_cellsize_changed:
        #                     target_size = automesher.max_cellsize_air / epsilon**0.5
        #                 else:
        #                     target_size = mesh_res
        #                 lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]-target_size, automesher.min_cellsize, target_size, 1.3)
        #                 if 2*abs(max(lines)-min(lines)) + target_size < abs(edges[i][0] - edges[i-1][0]):
        #                     if edges[i][1] < edges[i + 1][1] < edges[i][2] or edges[i][1] > edges[i + 1][1] > edges[i][2] or \
        #                         edges[i][1] < edges[i + 1][2] < edges[i][2] or edges[i][1] > edges[i + 1][2] > edges[i][2] or\
        #                         edges[i + 1][1] < edges[i][1] < edges[i + 1][2] or edges[i + 1][1] > edges[i][1] > edges[i + 1][2] or \
        #                         edges[i + 1][1] < edges[i][2] < edges[i + 1][2] or edges[i + 1][1] > edges[i][2] > edges[i + 1][2]:
        #                         mesh_data.extend(lines)
        #                 lines = add_graded_mesh_lines(automesher, edges[i][0], edges[i][0]+target_size, automesher.min_cellsize, target_size, 1.3)
        #                 if 2*abs(max(lines)-min(lines)) + target_size < abs(edges[i+1][0] - edges[i][0]):
        #                     if edges[i][1] < edges[i + 1][1] < edges[i][2] or edges[i][1] > edges[i + 1][1] > edges[i][2] or \
        #                         edges[i][1] < edges[i + 1][2] < edges[i][2] or edges[i][1] > edges[i + 1][2] > edges[i][2] or\
        #                         edges[i + 1][1] < edges[i][1] < edges[i + 1][2] or edges[i + 1][1] > edges[i][1] > edges[i + 1][2] or \
        #                         edges[i + 1][1] < edges[i][2] < edges[i + 1][2] or edges[i + 1][1] > edges[i][2] > edges[i + 1][2]:
        #                         mesh_data.extend(lines)
        #                     # mesh_data.extend(lines)
                    

    mesh_map.sort(key=lambda x: x[0])
    for i in range(len(mesh_map)-1):
        if direction == 'z':
            condition  = True
        else:
            condition  = (mesh_map[i][8][0] <= mesh_map[i+1][8][0] <= mesh_map[i][8][1] or mesh_map[i][8][0] <= mesh_map[i+1][8][1] <= mesh_map[i][8][1] or mesh_map[i+1] [8][0] <= mesh_map[i][8][0] <= mesh_map[i+1][8][1] or mesh_map[i+1][8][0] <= mesh_map[i][8][1] <= mesh_map[i+1][8][1])
        if mesh_data and mesh_map[i][0] <= min(mesh_data):
            continue
        if mesh_map[i][0] < mesh_map[i+1][0] and mesh_map[i][2] != mesh_map[i+1][2] and condition:
            if not automesher.min_cellsize_changed:
                target_size = automesher.max_cellsize_air / max(mesh_map[i][2], mesh_map[i+1][2])**0.5
            else:
                target_size = mesh_res
            # target_size = automesher.max_cellsize_air / max(mesh_map[i][2], mesh_map[i+1][2])**0.5
            lines = add_graded_mesh_lines(automesher, mesh_map[i+1][0], mesh_map[i+1][0]-target_size, automesher.min_cellsize, target_size, 1.3)
            mesh_data.extend(lines)
            lines = add_graded_mesh_lines(automesher, mesh_map[i+1][0], mesh_map[i+1][0]+target_size, automesher.min_cellsize, target_size, 1.3)
            mesh_data.extend(lines)
        if mesh_map[i+1][1] < mesh_map[i][1] and mesh_map[i][2] != mesh_map[i+1][2] and condition:
            if not automesher.min_cellsize_changed:
                target_size = automesher.max_cellsize_air / max(mesh_map[i][2], mesh_map[i+1][2])**0.5
            else:
                target_size = mesh_res            
            # target_size = automesher.max_cellsize_air / max(mesh_map[i][2], mesh_map[i+1][2])**0.5
            lines = add_graded_mesh_lines(automesher, mesh_map[i+1][1], mesh_map[i+1][1]-target_size, automesher.min_cellsize, target_size, 1.3)
            mesh_data.extend(lines)
            lines = add_graded_mesh_lines(automesher, mesh_map[i+1][1], mesh_map[i+1][1]+target_size, automesher.min_cellsize, target_size, 1.3)
            mesh_data.extend(lines)

    # for mesh_maps in mesh_map:
    #     start, end, epsilon, prop, x_point, xedges, ypoint, yedges, z_boundaries = mesh_maps
    #     if not automesher.min_cellsize_changed:
    #         target_size = automesher.max_cellsize_air / mesh_maps[2]**0.5
    #     else:
    #         target_size = mesh_res
    #     if xedges:
    #         xedges.sort(key=lambda x: x[0])
    #         for i in range(len(xedges)-1):
    #             if not mesh_data or xedges[i][0] <= min(mesh_data):
    #                 continue
    #             else:
    #                 lines = add_graded_mesh_lines(automesher, xedges[i][0], xedges[i][0]-target_size, automesher.min_cellsize, target_size, 1.3)
    #                 if i > 0 and (abs(max(lines)-min(lines)) + target_size >= xedges[i][0] - xedges[i-1][0] or xedges[i][0] - xedges[i-1][0] == 0):
    #                     lines_to_remove = [line for line in mesh_data if ((xedges[i-1][0] <= line <= xedges[i][0]) or (xedges[i][0] <= line <= xedges[i-1][0]))]
    #                     if lines_to_remove:
    #                         for line in lines_to_remove:
    #                             if line in mesh_data:
    #                                 mesh_data.remove(line)
    #                     # equal_lines = np.linspace(xedges[i-1][0], xedges[i][0], automesher.num_lines)
    #                     # mesh_data.extend(equal_lines)                    
    #                 else:   
    #                     mesh_data.extend(lines)
    #                 lines = add_graded_mesh_lines(automesher, xedges[i][0], xedges[i][0]+target_size, automesher.min_cellsize, target_size, 1.3)
    #                 if abs(max(lines)-min(lines)) + target_size >= xedges[i+1][0] - xedges[i][0] or xedges[i+1][0] - xedges[i][0] == 0:
    #                     lines_to_remove = [line for line in mesh_data if ((xedges[i][0] <= line <= xedges[i+1][0]) or (xedges[i+1][0] <= line <= xedges[i][0]))]
    #                     if lines_to_remove:
    #                         for line in lines_to_remove:
    #                             mesh_data.remove(line)    
    #                     # equal_lines = np.linspace(xedges[i][0], xedges[i+1][0], automesher.num_lines)
    #                     # mesh_data.extend(equal_lines)               
    #                 else:
    #                     mesh_data.extend(lines)




    



                    

        
        