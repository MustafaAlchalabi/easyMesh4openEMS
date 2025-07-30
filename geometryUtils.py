from CSXCAD import CSPrimitives, CSProperties
import numpy as np

def process_polygon(automesher, polygon, x, y, z, x_edges, y_edges, diagonal_edges, mesh_data):
    # Check if the input is a list of primitives
    if isinstance(polygon, list):
        # Process each primitive
        for prim in polygon:
            process_primitive(prim, x, y, x_edges, y_edges, diagonal_edges)

        # Collect z-coordinates from the polygon
        z.extend(collect_z_coordinates(polygon))

        # Process z-coordinates and add them to the mesh_data list
        process_z_coordinates(automesher, z, mesh_data)

    else:
        # If the polygon is a single primitive, get its x and y coordinates
        prim_x_coords, prim_y_coords = polygon.GetCoords()[0], polygon.GetCoords()[1]
        x.extend(prim_x_coords)
        y.extend(prim_y_coords)
        # Process the single polygon to extract edges and coordinates
        process_single_polygon(polygon, x, y, x_edges, y_edges, diagonal_edges)

        # Collect z-coordinates from the single polygon
        z.extend(collect_z_coordinates([polygon]))

        # Process z-coordinates and add them to the mesh_data
        process_z_coordinates(automesher, z, mesh_data)

def get_unique_edges(edges):
        unique_edges = [(edge[0], edge[3]) for edge in edges]
        unique_edges = list(set(unique_edges))
        unique_edges = list({edge[0]: edge for edge in unique_edges}.values())
        unique_edges.sort(key=lambda x: x[0])
        return unique_edges 

def tranfer_box_to_polygon(box):
    start = np.fmin(box.GetStart(), box.GetStop())
    stop = np.fmax(box.GetStart(), box.GetStop())
    x_coords = [start[0], stop[0], stop[0], start[0], start[0]]
    y_coords = [start[1], start[1], stop[1], stop[1], start[1]]
    z_coords = [float(start[2]), float(stop[2])]
    return x_coords, y_coords, z_coords

def transfer_port_to_polygon(start, stop):
    port_coords_x = [start[0], stop[0], stop[0], start[0], start[0]]
    port_coords_y = [start[1], start[1], stop[1], stop[1], start[1]]
    port_coords_z = [start[2], stop[2]]
    return port_coords_x, port_coords_y, port_coords_z

def process_primitive(prim, x, y, x_edges, y_edges, diagonal_edges):
    if not hasattr(prim, 'GetType'):
        port_coords_x, port_coords_y, port_coords_z = transfer_port_to_polygon(prim.start, prim.stop)
        x.extend(port_coords_x)
        y.extend(port_coords_y)
        collect_edges(port_coords_x, port_coords_y, prim, x_edges, y_edges, diagonal_edges)
    elif prim.GetType() == CSPrimitives.BOX:
        box_coords_x, box_coords_y, box_coords_z = tranfer_box_to_polygon(prim)
        x.extend(box_coords_x)
        y.extend(box_coords_y)
        collect_edges(box_coords_x, box_coords_y, prim, x_edges, y_edges, diagonal_edges)
    else:
        xx, yy = prim.GetCoords()[0], prim.GetCoords()[1]
        x.extend(xx)
        y.extend(yy)
        if xx[-1] != xx[0] or yy[-1] != yy[0]:
            xx = np.append(xx, xx[0])
            yy = np.append(yy, yy[0])
        collect_edges(xx, yy, prim, x_edges, y_edges, diagonal_edges)

def collect_edges(x_coords, y_coords, prim, x_edges, y_edges, diagonal_edges):
    for i in range(len(x_coords) - 1):
        if x_coords[i] != x_coords[i + 1] and y_coords[i] != y_coords[i + 1]:
            diagonal_edges.append([x_coords[i], x_coords[i + 1], y_coords[i], y_coords[i + 1], prim])
        if x_coords[i] == x_coords[i + 1]:
            x_edges.append([x_coords[i], y_coords[i], y_coords[i + 1], prim, False])
        if y_coords[i] == y_coords[i + 1]:
            y_edges.append([y_coords[i], x_coords[i], x_coords[i + 1], prim, False])

def collect_z_coordinates(polygon):
    z = [(prim.GetElevation(), prim) for prim in polygon if hasattr(prim, 'GetType') and prim.GetType() != CSPrimitives.BOX]
    z.extend((prim.GetElevation() + prim.GetLength(), prim) for prim in polygon if hasattr(prim, 'GetType') and prim.GetType() == CSPrimitives.LINPOLY)
    box_coords_z = [(tranfer_box_to_polygon(prim)[2][0], prim) for prim in polygon if hasattr(prim, 'GetType') and prim.GetType() == CSPrimitives.BOX]
    box_coords_z.extend((tranfer_box_to_polygon(prim)[2][1], prim) for prim in polygon if hasattr(prim, 'GetType') and prim.GetType() == CSPrimitives.BOX)
    z = list(set(z))
    z.sort(key=lambda x: x[0])
    z.extend(box_coords_z)
    return z

def process_z_coordinates(automesher, z, mesh_data):
    for z_val, prim in z:
        dirs = automesher.primitives_mesh_setup.get(prim, {}).get('dirs') or \
        automesher.properties_mesh_setup.get(prim.GetProperty(), {}).get('dirs') or \
        automesher.global_mesh_setup.get('dirs')
        if dirs is not None and 'z' in dirs:
            mesh_data[2].append(z_val)            

def process_single_polygon(polygon, x, y, x_edges, y_edges, diagonal_edges):
    xx, yy = polygon.GetCoords()[0], polygon.GetCoords()[1]

    x = np.append(x, xx)
    y = np.append(y, yy)
    for i in range(len(xx) - 1):
        if xx[i] != xx[i + 1] and yy[i] != yy[i + 1]:
            diagonal_edges.append([xx[i], xx[i + 1], yy[i], yy[i + 1], polygon])
        if xx[i] == xx[i + 1]:
            x_edges.append([xx[i], yy[i], yy[i + 1], polygon])
        if yy[i] == yy[i + 1]:
            y_edges.append([yy[i], xx[i], xx[i + 1], polygon])       

def distance_between_segments(p1, p2, q1, q2):
    p = np.linspace(p1, p2, 10)
    def point_to_line_distance(p, a, b):
        # Projektion des Punktes p auf die Linie a-b
        ap = p - a
        ab = b - a
        t = np.dot(ap, ab) / np.dot(ab, ab)
        # t = np.clip(t, 0, 1)  # Projektion auf das Segment beschränken
        closest_point = a + t * ab
        return np.linalg.norm(p - closest_point)
    distances = []
    for p_point in p:
        distances.append((point_to_line_distance(p_point, q1, q2), p_point, p1, p2, q1, q2))
    return distances

def detect_all_circles_in_polygon(automesher, polygon, min_points=20, tolerance=0.01):
    
    if automesher.global_mesh_setup.get('use_circle_detection', False) is False:
        return []
    if isinstance(polygon, list):
        coords = [prim.GetCoords() for prim in polygon if hasattr(prim, 'GetCoords')]
        x_coords = np.concatenate([coord[0] for coord in coords])
        y_coords = np.concatenate([coord[1] for coord in coords])
    else:                
        coords = polygon.GetCoords()
        x_coords = np.array(coords[0])
        y_coords = np.array(coords[1])
    N = len(x_coords)

    found_segments = []
    used_indices = set()

    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            indices = tuple(range(start, start + length))

            # Falls diese Punkte schon Teil eines gefundenen Kreises sind, überspringen
            if any(i in used_indices for i in indices):
                continue

            sub_x = x_coords[list(indices)]
            sub_y = y_coords[list(indices)]

            # Mittelpunkt-Näherung
            x0, y0 = np.mean(sub_x), np.mean(sub_y)
            radii = np.sqrt((sub_x - x0)**2 + (sub_y - y0)**2)
            mean_radius = np.mean(radii)
            deviation = np.abs(radii - mean_radius) / mean_radius

            if np.all(deviation < tolerance):
                found_segments.append((sub_x, sub_y))
                used_indices.update(indices)
                break  # Nicht überlappend: nächster Startpunkt
            
    return found_segments

def fit_circle_least_squares(x, y):
    """
    Least-squares Kreis-Fit nach algebraischer Methode.

        np.c_[np.array([1,2,3]), np.array([4,5,6])]:
            array([[1, 4],
                    [2, 5],
                    [3, 6]])
    """
    A = np.c_[2*x, 2*y, np.ones(len(x))]
    b = x**2 + y**2
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    xc, yc = sol[0], sol[1]
    r = np.sqrt(sol[2] + xc**2 + yc**2)
    return xc, yc, r

def detect_all_arcs_in_polygon(automesher, polygon, min_points=10, tolerance=0.03, max_radius_dev=0.03):
    
    if automesher.global_mesh_setup.get('use_arc_detection', False) is False:
        return []
    if isinstance(polygon, list):
        coords = [prim.GetCoords() for prim in polygon if hasattr(prim, 'GetCoords')]
        x_coords = np.concatenate([coord[0] for coord in coords])
        y_coords = np.concatenate([coord[1] for coord in coords])
    else:
        coords = polygon.GetCoords()
        x_coords, y_coords = np.array(coords[0]), np.array(coords[1])
    N = len(x_coords)

    found_arcs = []
    used_indices = set()

    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            indices = tuple(range(start, start + length))

            if any(i in used_indices for i in indices):
                continue

            sub_x = x_coords[list(indices)]
            sub_y = y_coords[list(indices)]

            try:
                xc, yc, r = fit_circle_least_squares(sub_x, sub_y)
            except Exception:
                continue

            distances = np.sqrt((sub_x - xc)**2 + (sub_y - yc)**2)
            mean_radius = np.mean(distances)
            deviation = np.abs(distances - mean_radius) / mean_radius

            if np.max(deviation) < max_radius_dev:
                found_arcs.append((sub_x, sub_y))
                used_indices.update(indices)
                break  # Nicht überlappen, nächster Start

    return found_arcs


def is_circle_in_polygon(polygon, min_points=8, tolerance=0.01, return_segment=False):
    coords = polygon.GetCoords()
    x_coords, y_coords = np.array(coords[0]), np.array(coords[1])
    N = len(x_coords)
    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            sub_x = x_coords[start:start+length]
            sub_y = y_coords[start:start+length]
            if len(sub_x) < 3:
                continue
            x0, y0 = np.mean(sub_x), np.mean(sub_y)
            radii = np.sqrt((sub_x - x0) ** 2 + (sub_y - y0) ** 2)
            mean_radius = np.mean(radii)
            deviation = np.abs(radii - mean_radius) / mean_radius
            if np.all(deviation < tolerance):
                if return_segment:
                    return True, sub_x, sub_y
                return True
    return (False, None, None) if return_segment else False


def calc_min_distance( x):
    min_distance = float('inf')
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            distance = abs(x[i] - x[j])
            if distance > 0 and distance < min_distance:
                min_distance = distance
    return min_distance

def point_in_polygon(polygon, point):
    """
    Raycasting Algorithm to find out whether a point is in a given polygon.
    Performs the even-odd-rule Algorithm to find out whether a point is in a given polygon.
    This runs in O(n) where n is the number of edges of the polygon.
    *
    :param polygon: an array representation of the polygon where polygon[i][0] is the x Value of the i-th point and polygon[i][1] is the y Value.
    :param point:   an array representation of the point where point[0] is its x Value and point[1] is its y Value
    :return: whether the point is in the polygon (not on the edge, just turn < into <= and > into >= for that)
    """

    # A point is in a polygon if a line from the point to infinity crosses the polygon an odd number of times
    odd = False
    # For each edge (In this case for each point of the polygon and the previous one)
    i = 0
    j = len(polygon[0]) - 1
    while i < len(polygon[0]) - 1:
        i = i + 1
        # If a line from the point into infinity crosses this edge
        # One point needs to be above, one below our y coordinate
        # ...and the edge doesn't cross our Y corrdinate before our x coordinate (but between our x coordinate and infinity)

        if (((polygon[1][i] > point[1]) != (polygon[1][j] > point[1])) and (point[0] < ((polygon[0][j] - polygon[0][i]) * (point[1] - polygon[1][i]) / (polygon[1][j] - polygon[1][i])) +polygon[0][i])):                # Invert odd
            odd = not odd
        j = i
    # If the number of crossings was odd, the point is in the polygon
    return odd

def yline_in_polygon(polygon, x_point, y_start, y_end):

    loop = np.linspace(y_start, y_end, 10)
    loop = loop[1:-1] 

    for y_val in loop:
        point = [x_point, y_val]
        if not point_in_polygon(polygon, point):
            return False

    return True

def xline_in_polygon(polygon, x_start, x_end, y_point):

    loop = np.linspace(x_start, x_end, 10)
    loop = loop[1:-1] 

    for x_val in loop:
        point = [x_val, y_point]
        if not point_in_polygon(polygon, point):
            return False

    return True

def metal_edge(automesher, edges, x_coords, y_coords, mesh_data, direction):

    # if metal_edge_res is not None:
    #     if unique_xedges[0] <= sorted_x[0]:
    #         mesh_data_in_range =  [mesh_data for mesh_data in mesh_data[0] if unique_xedges[0]-mer[1] <= mesh_data <= unique_xedges[0]-mer[0]]
    #         if not mesh_data_in_range:
    #             mesh_data[0].append(unique_xedges[0]-mer[1])
    #             mesh_data[0].append(unique_xedges[0]-mer[0])
    #         else:
    #             mesh_data[0] = [h for h in mesh_data[0] if h not in mesh_data_in_range]
    #             mesh_data[0].append(unique_xedges[0]-mer[1])
    #             mesh_data[0].append(unique_xedges[0]-mer[0])
    #     if unique_xedges[-1] >= sorted_x[-1]:
    #         mesh_data_in_range =  [mesh_data for mesh_data in mesh_data[0] if unique_xedges[-1]+mer[0] <= mesh_data <= unique_xedges[-1]+mer[1]]
    #         if not mesh_data_in_range:
    #             mesh_data[0].append(unique_xedges[-1]+mer[0])
    #             mesh_data[0].append(unique_xedges[-1]+mer[1])
    #         else:
    #             mesh_data[0] = [h for h in mesh_data[0] if h not in mesh_data_in_range]
    #             mesh_data[0].append(unique_xedges[-1]+mer[0])
    #             mesh_data[0].append(unique_xedges[-1]+mer[1])
    #     if unique_yedges[0] <= sorted_y[0]:
    #         mesh_data_in_range =  [mesh_data for mesh_data in mesh_data[1] if unique_yedges[0]-mer[1] <= mesh_data <= unique_yedges[0]-mer[0]]
    #         if not mesh_data_in_range:
    #             mesh_data[1].append(unique_yedges[0]-mer[1])
    #             mesh_data[1].append(unique_yedges[0]-mer[0])
    #         else:
    #             mesh_data[1] = [h for h in mesh_data[1] if h not in mesh_data_in_range]
    #             mesh_data[1].append(unique_yedges[0]-mer[1])
    #             mesh_data[1].append(unique_yedges[0]-mer[0])
    #     if unique_yedges[-1] >= sorted_y[-1]:
    #         mesh_data_in_range =  [mesh_data for mesh_data in mesh_data[1] if unique_yedges[-1]+mer[0] <= mesh_data <= unique_yedges[-1]+mer[1]]
    #         if not mesh_data_in_range:
    #             mesh_data[1].append(unique_yedges[-1]+mer[0])
    #             mesh_data[1].append(unique_yedges[-1]+mer[1])
    #         else:
    #             mesh_data[1] = [h for h in mesh_data[1] if h not in mesh_data_in_range]
    #             mesh_data[1].append(unique_yedges[-1]+mer[0])
    #             mesh_data[1].append(unique_yedges[-1]+mer[1])
    min_distance_x = calc_min_distance(x_coords)
    min_distance_y = calc_min_distance(y_coords)
    mer = np.array([-1.0, 2.0]) / 3 * automesher.min_cellsize
    edges_to_add = []
    edges_to_remove = []
    for edge in list(edges):  # Iterate over a copy of the list
        if hasattr(edge[3], 'GetProperty') and isinstance(edge[3].GetProperty(), CSProperties.CSPropMetal):
            x, y, x_edges, y_edges, diagonal_edges = [], [], [], [], []
            process_primitive(edge[3], x, y, x_edges, y_edges, diagonal_edges)
            coords = [x, y]

            if direction == 'x':
                condition1 = yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2]) and not yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2])
                condition2 = yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2]) and yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2])
            if direction == 'y':
                condition1 = xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2) and not xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2)
                condition2 = xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2) and xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2)
            if condition1:
                mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]-mer[1] <= mesh_data <= edge[0]-mer[0]]
                if not mesh_data_in_range:
                    edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                    edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                    edges_to_remove.append(edge)
                else:
                    mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                    mesh_data.append(edge[0]-mer[1])
                    mesh_data.append(edge[0]-mer[0])
                    edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                    edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                    edges_to_remove.append(edge)
            elif condition2:
                continue
            else:
                mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]+mer[0] <= mesh_data <= edge[0]+mer[1]]
                if not mesh_data_in_range:
                    edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                    edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                    edges_to_remove.append(edge)
                else:
                    mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                    mesh_data.append(edge[0]+mer[0])
                    mesh_data.append(edge[0]+mer[1])
                    edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                    edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                    edges_to_remove.append(edge)
    edges.extend(edges_to_add)
    for edge in edges_to_remove:
        if edge in edges:  
            edges.remove(edge)
    # if direction == 'x':
    #     min_distance_x = calc_min_distance(x)
    # if direction == 'y':
    #     min_distance_y = calc_min_distance(y)
    # for i in range(len(edges) - 1):
    #     if direction == 'x':
    #         condition1 = yline_in_polygon(coords, edges[i][0]+min_distance_x/2, edges[i][1], edges[i][2]) and not yline_in_polygon(coords, edges[i][0]-min_distance_x/2, edges[i][1], edges[i][2])
    #         condition2 = yline_in_polygon(coords, edges[i][0]-min_distance_x/2, edges[i][1], edges[i][2]) and yline_in_polygon(coords, edges[i][0]+min_distance_x/2, edges[i][1], edges[i][2])
    #     if direction == 'y':
    #         condition1 = xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]+min_distance_y/2) and not xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]-min_distance_y/2)
    #         condition2 = xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]-min_distance_y/2) and xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]+min_distance_y/2)
    #         if i > 0 and abs(edges[i][0] - edges[i + 1][0]) > mesh_res and abs(edges[i][0] - edges[i - 1][0]) > mesh_res:
    #             if condition1:
    #                 mesh_data_in_range =  [mesh_data for mesh_data in mesh_data if edges[i][0]-mer[1] <= mesh_data <= edges[i][0]-mer[0]]
    #                 if not mesh_data_in_range:
    #                     mesh_data.append(edges[i][0]-mer[1])
    #                     mesh_data.append(edges[i][0]-mer[0])
    #                 else:
    #                     mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
    #                     mesh_data.append(edges[i][0]-mer[1])
    #                     mesh_data.append(edges[i][0]-mer[0])
    #             elif condition2:
    #                 continue
    #             else:
    #                 mesh_data_in_range =  [mesh_data for mesh_data in mesh_data if edges[i][0]+mer[0] <= mesh_data <= edges[i][0]+mer[1]]
    #                 if not mesh_data_in_range:
    #                     mesh_data.append(edges[i][0]+mer[0])
    #                     mesh_data.append(edges[i][0]+mer[1])
    #                 else:
    #                     mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
    #                     mesh_data.append(edges[i][0]+mer[0])
    #                     mesh_data.append(edges[i][0]+mer[1])