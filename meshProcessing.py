from CSXCAD.SmoothMeshLines import SmoothMeshLines
import numpy as np
import geometryUtils, meshingUtils

def smooth_and_process_mesh_lines(automesher, mesh_data, polygon, grid, x_edges, y_edges, z, unique_xedges, unique_yedges, z_coords, mesh_map):

    mesh_data[0] = sorted(mesh_data[0])
    mesh_data[1] = sorted(mesh_data[1])
    mesh_data[2] = sorted(mesh_data[2])
    automesher.mesh_with_max_cell_size = [[], [], []]
    xmax, xmin, ymax, ymin, zmax, zmin = max(mesh_data[0]), min(mesh_data[0]), max(mesh_data[1]), min(mesh_data[1]), max(mesh_data[2]), min(mesh_data[2])

    if isinstance(polygon, list):
        if automesher.min_cellsize_changed:
            if not any(automesher.primitives_mesh_setup.get(prim, {}).get('edges_only', False) for prim in polygon):
                for i in range(len(mesh_data[0]) - 1):
                    if mesh_data[0][i + 1] - mesh_data[0][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[0].append((mesh_data[0][i], mesh_data[0][i + 1]))
                for i in range(len(mesh_data[1]) - 1):
                    if mesh_data[1][i + 1] - mesh_data[1][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[1].append((mesh_data[1][i], mesh_data[1][i + 1]))
                for i in range(len(mesh_data[2]) - 1):
                    if mesh_data[2][i + 1] - mesh_data[2][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[2].append((mesh_data[2][i], mesh_data[2][i + 1]))

                for idx in range(3):
                    mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
                    for start, end in automesher.mesh_with_max_cell_size[idx]:
                        mesh_data[idx] = [line for line in mesh_data[idx] if not (start < line < end)]
        else:
            for idx in range(3):
                mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
    else:
        if automesher.min_cellsize_changed:
            if not automesher.primitives_mesh_setup.get(polygon, {}).get('edges_only', False):
                for i in range(len(mesh_data[0]) - 1):
                    if mesh_data[0][i + 1] - mesh_data[0][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[0].append((mesh_data[0][i], mesh_data[0][i + 1]))
                for i in range(len(mesh_data[1]) - 1):
                    if mesh_data[1][i + 1] - mesh_data[1][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[1].append((mesh_data[1][i], mesh_data[1][i + 1]))
                for i in range(len(mesh_data[2]) - 1):
                    if mesh_data[2][i + 1] - mesh_data[2][i] > automesher.max_cellsize / 2:
                        automesher.mesh_with_max_cell_size[2].append((mesh_data[2][i], mesh_data[2][i + 1]))

                for idx in range(3):
                    mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
                    for start, end in automesher.mesh_with_max_cell_size[idx]:
                        mesh_data[idx] = [line for line in mesh_data[idx] if not (start < line < end)]
        else:
            for idx in range(3):
                mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()

    if mesh_map[0]:
        mesh_map[0].sort(key=lambda epsilon: epsilon[2], reverse=True)
    if mesh_map[1]:
        mesh_map[1].sort(key=lambda epsilon: epsilon[2], reverse=True)
    if mesh_map[2]:
        mesh_map[2].sort(key=lambda epsilon: epsilon[2], reverse=True)
    
    for map in mesh_map[0]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = [line for line in mesh_data[0] if map[0] <= line <= map[1]]
        if lines_to_be_smoothed:
            mesh_data[0].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize).tolist())
    for map in mesh_map[1]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = [line for line in mesh_data[1] if map[0] <= line <= map[1]]
        if lines_to_be_smoothed:
            mesh_data[1].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize).tolist())
    for map in mesh_map[2]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = [line for line in mesh_data[2] if map[0] <= line <= map[1]]
        if lines_to_be_smoothed:
            mesh_data[2].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize).tolist())

    if automesher.min_cellsize_changed:
        if automesher.mesh_with_max_cell_size[0]:
            lines = SmoothMeshLines(mesh_data[0], automesher.max_cellsize/2, 1.3).tolist()
            lines = sorted(set(lines))
            lines_in_range = []
            mean_resolution = []
            skipping_list = []
            for start, end in automesher.mesh_with_max_cell_size[0]:
                lines_to_add_in_lines_in_range = [line for line in lines if start < line < end]
                if lines_to_add_in_lines_in_range:
                    lines_to_add_in_lines_in_range = sorted(set(lines_to_add_in_lines_in_range))
                    lines_in_range_in_mesh_data = [line for line in mesh_data[0] if start < line < end]
                    mean_resolution_in_mesh_data = np.mean(np.diff(lines_in_range_in_mesh_data)) if lines_in_range_in_mesh_data else automesher.max_cellsize/2
                    if lines_in_range_in_mesh_data and mean_resolution_in_mesh_data > automesher.mesh_res:
                        skipping_list.append(True)
                    else:
                        skipping_list.append(False)
                else:
                    lines_to_add_in_lines_in_range = []
                lines_in_range.append(lines_to_add_in_lines_in_range)
                mean_resolution.append(np.mean(np.diff(lines_to_add_in_lines_in_range)) if lines_to_add_in_lines_in_range else automesher.max_cellsize/2)
            if lines_in_range and mean_resolution:
                for i in range(len(lines_in_range)):
                    if lines_in_range[i]:
                        if mean_resolution[i] <= automesher.max_cellsize/2 and skipping_list[i]:
                            lines = [line for line in lines if not( min(lines_in_range[i]) < line < max(lines_in_range[i]))]
            
            # check if there are lines in range between  min lines and max lines alredy in mesh_data with the mean resolution <= automesher.max_cellsize/2. if so than delete these lines from the lines list
            for start, end in automesher.mesh_with_max_cell_size[0]:
                lines_to_add = [line for line in lines if (start < line < end)]
                mesh_data[0].extend(lines_to_add)
        if automesher.mesh_with_max_cell_size[1]:
            lines = SmoothMeshLines(mesh_data[1], automesher.max_cellsize/2, 1.3).tolist()
            lines = sorted(set(lines))
            lines_in_range = []
            mean_resolution = []
            skipping_list = []
            for start, end in automesher.mesh_with_max_cell_size[1]:
                lines_to_add_in_lines_in_range = [line for line in lines if start < line < end]
                if lines_to_add_in_lines_in_range:
                    lines_to_add_in_lines_in_range = sorted(set(lines_to_add_in_lines_in_range))
                    lines_in_range_in_mesh_data = [line for line in mesh_data[1] if start < line < end]
                    mean_resolution_in_mesh_data = np.mean(np.diff(lines_in_range_in_mesh_data)) if lines_in_range_in_mesh_data else automesher.max_cellsize/2
                    if lines_in_range_in_mesh_data and mean_resolution_in_mesh_data > automesher.mesh_res:
                        skipping_list.append(True)
                    else:
                        skipping_list.append(False)
                else:
                    lines_to_add_in_lines_in_range = []
                lines_in_range.append(lines_to_add_in_lines_in_range)
                mean_resolution.append(np.mean(np.diff(lines_to_add_in_lines_in_range)) if lines_to_add_in_lines_in_range else automesher.max_cellsize/2)
            if lines_in_range and mean_resolution:
                for i in range(len(lines_in_range)):
                    if lines_in_range[i]:
                        if mean_resolution[i] <= automesher.max_cellsize/2 and skipping_list[i]:
                            lines = [line for line in lines if not(min(lines_in_range[i]) < line < max(lines_in_range[i]))]
            for start, end in automesher.mesh_with_max_cell_size[1]:
                lines_to_add = [line for line in lines if (start < line < end)]
                mesh_data[1].extend(lines_to_add)
        if automesher.mesh_with_max_cell_size[2]:
            lines = SmoothMeshLines(mesh_data[2], automesher.max_cellsize/2, 1.3).tolist()
            lines = sorted(set(lines))
            lines_in_range = []
            mean_resolution = []
            skipping_list = []
            for start, end in automesher.mesh_with_max_cell_size[2]:
                lines_to_add_in_lines_in_range = [line for line in lines if start < line < end]
                if lines_to_add_in_lines_in_range:
                    lines_to_add_in_lines_in_range = sorted(set(lines_to_add_in_lines_in_range))
                    lines_in_range_in_mesh_data = [line for line in mesh_data[2] if start < line < end]
                    mean_resolution_in_mesh_data = np.mean(np.diff(lines_in_range_in_mesh_data)) if lines_in_range_in_mesh_data else automesher.max_cellsize/2
                    if lines_in_range_in_mesh_data and mean_resolution_in_mesh_data > automesher.mesh_res:
                        skipping_list.append(True)
                    else:
                        skipping_list.append(False)
                else:
                    lines_to_add_in_lines_in_range = []
                lines_in_range.append(lines_to_add_in_lines_in_range)
                mean_resolution.append(np.mean(np.diff(lines_to_add_in_lines_in_range)) if lines_to_add_in_lines_in_range else automesher.max_cellsize/2)
            if lines_in_range and mean_resolution:
                for i in range(len(lines_in_range)):
                    if lines_in_range[i]:
                        if mean_resolution[i] <= automesher.max_cellsize/2 and skipping_list[i]:
                            lines = [line for line in lines if not(min(lines_in_range[i]) < line < max(lines_in_range[i]))]
            for start, end in automesher.mesh_with_max_cell_size[2]:
                lines_to_add = [line for line in lines if (start < line < end)]
                mesh_data[2].extend(lines_to_add)

    # for i in range(1, len(np.diff(lines[2])) - 1):
    #     # check if the difference between two consecutive z values is greater than 2 times the difference between the next two consecutive z values
    #     if i + 1 < len(lines[2][0]) and np.round(np.diff(lines[2][0])[i] / np.diff(lines[2][0])[i + 1], 1) > 2 and np.diff(lines[2][0])[i] > automesher.min_cellsize:
    #         lines[2][0] = list(lines[2][0])  # Convert to list
    #         lines[2][0].extend(SmoothMeshLines([lines[2][0][i], lines[2][0][i + 1]], automesher.mesh_res/2, 1.3))

    # # Check lines between x edges
    # for i in range(len(x_edges) - 1):
    #     if abs(x_edges[i][0] - x_edges[i + 1][0]) > automesher.mesh_res:
    #         lines_in_range = [line for line in mesh_data[0] if x_edges[i][0] < line < x_edges[i + 1][0]]
    #         if not lines_in_range:
    #             mesh_data[0] = np.append(mesh_data[0], np.linspace(x_edges[i][0], x_edges[i + 1][0], automesher.num_lines))
    # # Check lines between y edges
    # for i in range(len(y_edges) - 1):
    #     if abs(y_edges[i][0] - y_edges[i + 1][0]) > automesher.mesh_res:
    #         lines_in_range = [line for line in mesh_data[1] if y_edges[i][0] < line < y_edges[i + 1][0]]
    #         if not lines_in_range:
    #             mesh_data[1] = np.append(mesh_data[1], np.linspace(y_edges[i][0], y_edges[i + 1][0], automesher.num_lines))

#     automesher.global_mesh_setup:'boundary_distance': [ 1000, 1000, 1000, 1000, 1000, 1000 ], # value, auto or None
    graded_lines_y = []
    graded_lines_x = []
    graded_lines_z = []
    # Ensure all numbers in mesh_data are converted to the same format (float)
    mesh_data[0] = [float(num) if isinstance(num, np.float64) else num for num in mesh_data[0]]
    mesh_data[1] = [float(num) if isinstance(num, np.float64) else num for num in mesh_data[1]]
    mesh_data[2] = [float(num) if isinstance(num, np.float64) else num for num in mesh_data[2]]
    mesh_data[0] = sorted(set(mesh_data[0]))
    mesh_data[1] = sorted(set(mesh_data[1]))
    mesh_data[2] = sorted(set(mesh_data[2]))

    distance = automesher.global_mesh_setup.get('boundary_distance', [0, 0, 0, 0, 0, 0])
    for i in range(len(distance)):
        if distance[i] == 'auto':
            distance[i] = automesher.wave_length
        elif distance[i] is None:
            distance[i] = 0
    if xmax in mesh_data[0]:
        graded_lines_x.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[0]), xmax+distance[0], abs(np.max(mesh_data[0])- mesh_data[0][np.argmax(mesh_data[0]) - 1]), automesher.max_cellsize_air, 1.3))
    if xmin in mesh_data[0]:
        graded_lines_x.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[0]), xmin-distance[1], abs(np.min(mesh_data[0]) - mesh_data[0][np.argmin(mesh_data[0]) + 1]) , automesher.max_cellsize_air, 1.3))
    if ymax in mesh_data[1]:
        graded_lines_y.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[1]), ymax+distance[2], abs(np.max(mesh_data[1])- mesh_data[1][np.argmax(mesh_data[1]) - 1]), automesher.max_cellsize_air, 1.3))
    if ymin in mesh_data[1]:
        graded_lines_y.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[1]), ymin-distance[3], abs(np.min(mesh_data[1]) - mesh_data[1][np.argmin(mesh_data[1]) + 1]), automesher.max_cellsize_air, 1.3))
    if mesh_data[2] and zmax in mesh_data[2]:
        graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[2]), zmax+distance[4], abs(np.max(mesh_data[2])- mesh_data[2][np.argmax(mesh_data[2]) - 1]), automesher.max_cellsize_air, 1.3))
    if mesh_data[2] and zmin in mesh_data[2]:
        if np.argmin(mesh_data[2]) + 1 >= len(mesh_data[2]):
            graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[2]), zmin-distance[5], 0, automesher.max_cellsize_air, 1.3))
        else:
            graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[2]), zmin-distance[5], abs(np.min(mesh_data[2]) - mesh_data[2][np.argmin(mesh_data[2]) + 1]), automesher.max_cellsize_air, 1.3))

    # add  graded lines to lines list
    mesh_data[0] = np.append(mesh_data[0], graded_lines_x)
    mesh_data[1] = np.append(mesh_data[1], graded_lines_y)
    mesh_data[2] = np.append(mesh_data[2], graded_lines_z)

    mesh_data[0] = mesh_data[0].tolist() if isinstance(mesh_data[0], np.ndarray) else mesh_data[0]
    mesh_data[1] = mesh_data[1].tolist() if isinstance(mesh_data[1], np.ndarray) else mesh_data[1]
    mesh_data[2] = mesh_data[2].tolist() if isinstance(mesh_data[2], np.ndarray) else mesh_data[2]

    # if automesher.global_mesh_setup.get('min_cellsize', None) is not None or automesher.min_cellsize_changed:
    mesh_data[0] = process_mesh_data(mesh_data[0], automesher.min_cellsize, unique_xedges)
    mesh_data[1] = process_mesh_data(mesh_data[1], automesher.min_cellsize, unique_yedges)
    mesh_data[2] = process_mesh_data(mesh_data[2], automesher.min_cellsize, z_coords)

def process_mesh_data(mesh_data, min_cellsize, unique_edges):
    mesh_data = sorted(mesh_data)
    while True:  
        new_mesh_data = []
        skip_next = False
        changed = False 

        for i in range(len(mesh_data) - 1):
            if skip_next:
                skip_next = False
                continue

            if abs(mesh_data[i+1] - mesh_data[i]) < min_cellsize / 2:
                changed = True  
                if any(mesh_data[i] == edge[0] for edge in unique_edges) and not any(mesh_data[i+1] == edge[0] for edge in unique_edges):
                    new_mesh_data.append(mesh_data[i])
                    skip_next = True
                elif any(mesh_data[i+1] == edge[0] for edge in unique_edges) and not any(mesh_data[i] == edge[0] for edge in unique_edges):
                    new_mesh_data.append(mesh_data[i+1])
                    skip_next = True
                elif any(mesh_data[i] == edge[0] for edge in unique_edges) and any(mesh_data[i+1] == edge[0] for edge in unique_edges):
                    # new_mesh_data.append((mesh_data[i] + mesh_data[i+1]) / 2)
                    skip_next = False
                    continue
                else:
                    new_mesh_data.append((mesh_data[i] + mesh_data[i+1]) / 2)
                    skip_next = True
            else:
                new_mesh_data.append(mesh_data[i])

        if not skip_next and mesh_data:
            new_mesh_data.append(mesh_data[-1])

        if not changed:
            break  

        mesh_data = new_mesh_data  

    return mesh_data
