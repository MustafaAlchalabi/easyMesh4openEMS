from CSXCAD.SmoothMeshLines import SmoothMeshLines
import numpy as np
import geometryUtils, meshingUtils

def smooth_and_process_mesh_lines(automesher, mesh_data, polygon, grid, x_edges, y_edges, z, unique_xedges, unique_yedges, z_coords, mesh_map):

    mesh_data[0] = sorted(mesh_data[0])
    mesh_data[1] = sorted(mesh_data[1])
    mesh_data[2] = sorted(mesh_data[2])
    z_coords.sort(key=lambda edge: edge[0])
    x_edges.sort(key=lambda edge: edge[0])
    y_edges.sort(key=lambda edge: edge[0])  
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
                    if len(mesh_data[idx]) > 1 and automesher.mesh_res > 0:
                        mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
                    else:
                        continue
                        # raise ValueError(f"Invalid input for SmoothMeshLines: len(mesh_data[{idx}])={len(mesh_data[idx])}, automesher.mesh_res={automesher.mesh_res}")
                    for start, end in automesher.mesh_with_max_cell_size[idx]:
                        mesh_data[idx] = [line for line in mesh_data[idx] if not (start < line < end)]
        # else:
        #     for idx in range(3):
        #         mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
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
                    if len(mesh_data[idx]) > 1 and automesher.mesh_res > 0:
                        mesh_data[idx] = SmoothMeshLines(mesh_data[idx], automesher.mesh_res).tolist()
                    else:
                        continue
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
    already_smoothed = [[], [], []]
    tmp_mesh_with_max_cell_size = automesher.mesh_with_max_cell_size

    for i, (start, end) in enumerate(tmp_mesh_with_max_cell_size[0]):
        if not any(start == x_edges[0] for x_edges in unique_xedges):
            # search for the closest edge in xedges and replace it
            closest_edge = min(unique_xedges, key=lambda edge: abs(edge[0] - start))
            if closest_edge and abs(start-closest_edge[0]) > automesher.mesh_res:
                continue
            start = closest_edge[0]
        if not any(end == x_edges[0] for x_edges in unique_xedges):
            closest_edge = min(unique_xedges, key=lambda edge: abs(edge[0] - end))
            if closest_edge and abs(end- closest_edge[0]) > automesher.mesh_res:
                continue
            end = closest_edge[0]
        tmp_mesh_with_max_cell_size[0][i] = (start, end)

    for i, (start, end) in enumerate(tmp_mesh_with_max_cell_size[1]):
        if not any(start == y_edges[0] for y_edges in unique_yedges):
            # search for the closest edge in yedges and replace it
            closest_edge = min(unique_yedges, key=lambda edge: abs(edge[0] - start))
            if closest_edge and abs(start-closest_edge[0]) > automesher.mesh_res:
                continue
            start = closest_edge[0]
        if not any(end == y_edges[0] for y_edges in unique_yedges):
            closest_edge = min(unique_yedges, key=lambda edge: abs(edge[0] - end))
            if closest_edge and abs(end-closest_edge[0]) > automesher.mesh_res:
                continue
            end = closest_edge[0]
        tmp_mesh_with_max_cell_size[1][i] = (start, end)

    for map in mesh_map[0]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = sorted(set(float(line) for line in mesh_data[0] if map[0] <= line <= map[1]))
        if any(map[0] == start or map[0] == end for start, end in tmp_mesh_with_max_cell_size[0]) or \
           any(map[1] == start or map[1] == end for start, end in tmp_mesh_with_max_cell_size[0]):
            continue
        if lines_to_be_smoothed:
            already_smoothed[0].append((map[0], map[1]))
            mesh_data[0].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize, 1.3).tolist())

    for map in mesh_map[1]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = sorted(set(float(line) for line in mesh_data[1] if map[0] <= line <= map[1]))
        if any(map[0] == start or map[0] == end for start, end in tmp_mesh_with_max_cell_size[1]) or \
              any(map[1] == start or map[1] == end for start, end in tmp_mesh_with_max_cell_size[1]):    
            continue
        if lines_to_be_smoothed:
            already_smoothed[1].append((map[0], map[1]))
            mesh_data[1].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize, 1.3).tolist())
            
    for map in mesh_map[2]:
        max_cellsize = automesher.max_cellsize_air / map[2]**0.5
        lines_to_be_smoothed = sorted(set(float(line) for line in mesh_data[2] if map[0] <= line <= map[1]))
        if any(map[0] == start or map[0] == end for start, end in tmp_mesh_with_max_cell_size[2]) or \
              any(map[1] == start or map[1] == end for start, end in tmp_mesh_with_max_cell_size[2]):
            continue
        if lines_to_be_smoothed:
            mesh_data[2].extend(SmoothMeshLines(lines_to_be_smoothed, max_cellsize, 1.3).tolist())

    if automesher.min_cellsize_changed:
        if automesher.mesh_with_max_cell_size[0]:
            lines = SmoothMeshLines(mesh_data[0], automesher.max_cellsize/2, 1.3).tolist()
            lines = sorted(set(lines))
            lines_in_range = []
            mean_resolution = []
            skipping_list = [False] * len(automesher.mesh_with_max_cell_size[0])
            for idx, (start, end) in enumerate(automesher.mesh_with_max_cell_size[0]):
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
            skipping_list = [False] * len(automesher.mesh_with_max_cell_size[1])
            for idx, (start, end) in enumerate(automesher.mesh_with_max_cell_size[1]):
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
            skipping_list = [False] * len(automesher.mesh_with_max_cell_size[2])
            for idx, (start, end) in enumerate(automesher.mesh_with_max_cell_size[2]):
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

    mesh_data[0] = sorted(mesh_data[0])
    mesh_data[1] = sorted(mesh_data[1])
    mesh_data[2] = sorted(mesh_data[2])
    z_coords.sort(key=lambda edge: edge[0])
    x_edges.sort(key=lambda edge: edge[0])
    y_edges.sort(key=lambda edge: edge[0])  
    # Check lines between x edges
    for i in range(len(x_edges) - 1):
        if abs(x_edges[i][0] - x_edges[i + 1][0]) > automesher.mesh_res:
            lines_in_range = [line for line in mesh_data[0] if x_edges[i][0] < line < x_edges[i + 1][0]]
            lines_in_range = sorted(set(lines_in_range))  # Ensure unique lines 
            if not lines_in_range:
                lines_to_add = np.linspace(x_edges[i][0], x_edges[i + 1][0], automesher.num_lines)
                lines_to_add = list(lines_to_add)   
                mesh_data[0].extend(lines_to_add)

    # Check lines between y edges
    for i in range(len(y_edges) - 1):
        if abs(y_edges[i][0] - y_edges[i + 1][0]) > automesher.mesh_res:
            lines_in_range = [line for line in mesh_data[1] if y_edges[i][0] < line < y_edges[i + 1][0]]
            lines_in_range = sorted(set(lines_in_range))  # Ensure unique lines
            if not lines_in_range:
                lines_to_add = np.linspace(y_edges[i][0], y_edges[i + 1][0], automesher.num_lines)
                lines_to_add = list(lines_to_add)   
                mesh_data[1].extend(lines_to_add)

    # Check lines between z edges
    for i in range(len(z_coords) - 1):
        if abs(z_coords[i][0] - z_coords[i + 1][0]) > automesher.mesh_res_z:
            lines_in_range = [line for line in mesh_data[2] if z_coords[i][0] < line < z_coords[i + 1][0]]
            lines_in_range = sorted(set(lines_in_range))  # Ensure unique lines
            if not lines_in_range:
                lines_to_add = np.linspace(z_coords[i][0], z_coords[i + 1][0], automesher.num_lines)
                mesh_data[2].extend(lines_to_add.tolist())
            if len(lines_in_range) < automesher.num_lines-2:
                lines_to_add = np.linspace(z_coords[i][0], z_coords[i + 1][0], automesher.num_lines)
                if lines_in_range:
                    lines_to_add = list(lines_to_add)
                mesh_data[2].extend(lines_to_add)

    # automesher.global_mesh_setup:'boundary_distance': [ 1000, 1000, 1000, 1000, 1000, 1000 ], # value, auto or None
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
    if xmin in mesh_data[0]:
        graded_lines_x.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[0]), xmin-distance[0], abs(np.min(mesh_data[0]) - mesh_data[0][np.argmin(mesh_data[0]) + 1]) , automesher.max_cellsize_air, 1.3))
    if xmax in mesh_data[0]:
        graded_lines_x.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[0]), xmax+distance[1], abs(np.max(mesh_data[0])- mesh_data[0][np.argmax(mesh_data[0]) - 1]), automesher.max_cellsize_air, 1.3))
    if ymin in mesh_data[1]:
        graded_lines_y.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[1]), ymin-distance[2], abs(np.min(mesh_data[1]) - mesh_data[1][np.argmin(mesh_data[1]) + 1]), automesher.max_cellsize_air, 1.3))
    if ymax in mesh_data[1]:
        graded_lines_y.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[1]), ymax+distance[3], abs(np.max(mesh_data[1])- mesh_data[1][np.argmax(mesh_data[1]) - 1]), automesher.max_cellsize_air, 1.3))
    if mesh_data[2] and zmin in mesh_data[2]:
        if np.argmin(mesh_data[2]) + 1 >= len(mesh_data[2]):
            graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[2]), zmin-distance[4], 0, automesher.max_cellsize_air, 1.3))
        else:
            graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.min(mesh_data[2]), zmin-distance[4], abs(np.min(mesh_data[2]) - mesh_data[2][np.argmin(mesh_data[2]) + 1]), automesher.max_cellsize_air, 1.3))
    if mesh_data[2] and zmax in mesh_data[2]:
        graded_lines_z.extend(meshingUtils.add_graded_mesh_lines(automesher, np.max(mesh_data[2]), zmax+distance[5], abs(np.max(mesh_data[2])- mesh_data[2][np.argmax(mesh_data[2]) - 1]), automesher.max_cellsize_air, 1.3))

    # add  graded lines to lines list
    mesh_data[0] = np.append(mesh_data[0], graded_lines_x)
    mesh_data[1] = np.append(mesh_data[1], graded_lines_y)
    mesh_data[2] = np.append(mesh_data[2], graded_lines_z)

    mesh_data[0] = mesh_data[0].tolist() if isinstance(mesh_data[0], np.ndarray) else mesh_data[0]
    mesh_data[1] = mesh_data[1].tolist() if isinstance(mesh_data[1], np.ndarray) else mesh_data[1]
    mesh_data[2] = mesh_data[2].tolist() if isinstance(mesh_data[2], np.ndarray) else mesh_data[2]

    # if automesher.global_mesh_setup.get('min_cellsize', None) is not None or automesher.min_cellsize_changed:
    mesh_data[0] = process_mesh_data(mesh_data[0], automesher.min_cellsize, x_edges)
    mesh_data[1] = process_mesh_data(mesh_data[1], automesher.min_cellsize, y_edges)
    mesh_data[2] = process_mesh_data(mesh_data[2], automesher.min_cellsize_z, z_coords)

    lines = [sorted(set(mesh_data[0])), sorted(set(mesh_data[1])), sorted(set(mesh_data[2]))]
    changed = True
    iteration_count = 0  # Initialize iteration counter
    max_iterations = 10  # Set maximum iterations to avoid infinite loop
    while changed and iteration_count < max_iterations:
        differences = np.diff(lines[0])  # Precompute differences
        changed = False
        for i in range(1, len(differences) - 1):
            # Check if the difference between two consecutive x values is greater than 2 times the difference between the next two consecutive z values
            if abs(np.round(differences[i + 1] / differences[i], 1)) == 2:
                # Add extra lines between these two x values
                lines_to_add = lines[0][i+1] + abs(lines[0][i+1]-lines[0][i])
                lines[0].append(lines_to_add)
                mesh_data[0].append(lines_to_add)
                lines[0] = sorted(set(lines[0]))  # Re-sort and remove duplicates
                differences = np.diff(lines[0])  # Recompute differences
                changed = True
            if abs(np.round(differences[i + 1] / differences[i], 1)) > 2:
                # Add extra lines between these two x values
                lines_to_add = lines[0][i+1] + 1.3*abs(lines[0][i+1]-lines[0][i])
                lines[0].append(lines_to_add)
                mesh_data[0].append(lines_to_add)
                lines[0] = sorted(set(lines[0]))  # Re-sort and remove duplicates
                differences = np.diff(lines[0])  # Recompute differences
                changed = True
        iteration_count += 1  # Increment iteration counter

    lines = [sorted(set(mesh_data[0]), reverse=True), sorted(set(mesh_data[1]), reverse=True), sorted(set(mesh_data[2]), reverse=True)]
    changed = True
    iteration_count = 0  # Initialize iteration counter
    max_iterations = 10  # Set maximum iterations to avoid infinite loop
    while changed and iteration_count < max_iterations:
        # lines = [sorted(set(mesh_data[0]), reverse=True), sorted(set(mesh_data[1]), reverse=True), sorted(set(mesh_data[2]), reverse=True)]
        differences = np.diff(lines[0])  # Precompute differences (reverse order for descending)
        changed = False
        for i in range(1, len(differences) - 1):
            # Check if the difference between two consecutive x values is greater than 2 times the difference between the next two consecutive z values
            if abs(np.round(differences[i + 1] / differences[i], 1)) == 2:
                # Add extra lines between these two x values
                lines_to_add = lines[0][i+1] - abs(lines[0][i+1]-lines[0][i])
                lines[0].append(lines_to_add)
                mesh_data[0].append(lines_to_add)
                lines[0] = sorted(set(lines[0]), reverse=True)  # Re-sort in descending order and remove duplicates
                differences = np.diff(lines[0])  # Recompute differences (reverse order for descending)
                changed = True
            if abs(np.round(differences[i + 1] / differences[i], 1)) > 2:
                # Add extra lines between these two x values
                lines_to_add = lines[0][i+1] - 1.3*abs(lines[0][i+1]-lines[0][i])
                lines[0].append(lines_to_add)
                mesh_data[0].append(lines_to_add)
                lines[0] = sorted(set(lines[0]), reverse=True)  # Re-sort in descending order and remove duplicates
                differences = np.diff(lines[0])
                changed = True
        iteration_count += 1  # Increment iteration counter

    lines = [sorted(set(mesh_data[0])), sorted(set(mesh_data[1])), sorted(set(mesh_data[2]))]
    changed = True
    iteration_count = 0  # Initialize iteration counter
    max_iterations = 10  # Set maximum iterations to avoid infinite loop
    while changed and iteration_count < max_iterations:
        differences = np.diff(lines[1])
        changed = False
        for i in range(1, len(differences) - 1):
            if abs(np.round(differences[i + 1] / differences[i], 1)) == 2:
                lines_to_add = lines[1][i+1] + abs(lines[1][i+1]-lines[1][i])
                lines[1].append(lines_to_add)
                mesh_data[1].append(lines_to_add)
                lines[1] = sorted(set(lines[1]))
                differences = np.diff(lines[1])
                changed = True
            if abs(np.round(differences[i + 1] / differences[i], 1)) > 2:
                lines_to_add = lines[1][i+1] + 1.3*abs(lines[1][i+1]-lines[1][i])
                lines[1].append(lines_to_add)
                mesh_data[1].append(lines_to_add)
                lines[1] = sorted(set(lines[1]))
                differences = np.diff(lines[1])
                changed = True
        iteration_count += 1  # Increment iteration counter
    lines = [sorted(set(mesh_data[0]), reverse=True), sorted(set(mesh_data[1]), reverse=True), sorted(set(mesh_data[2]), reverse=True)]
    changed = True
    iteration_count = 0  # Initialize iteration counter
    max_iterations = 10  # Set maximum iterations to avoid infinite loop
    while changed and iteration_count < max_iterations:
        differences = np.diff(lines[1])
        changed = False
        for i in range(1, len(differences) - 1):
            if abs(np.round(differences[i + 1] / differences[i], 1)) == 2:
                lines_to_add = lines[1][i+1] - abs(lines[1][i+1]-lines[1][i])
                lines[1].append(lines_to_add)
                mesh_data[1].append(lines_to_add)
                lines[1] = sorted(set(lines[1]), reverse=True)
                differences = np.diff(lines[1])
                changed = True
            if abs(np.round(differences[i + 1] / differences[i], 1)) > 2:
                lines_to_add = lines[1][i+1] - 1.3*abs(lines[1][i+1]-lines[1][i])
                lines[1].append(lines_to_add)
                mesh_data[1].append(lines_to_add)
                lines[1] = sorted(set(lines[1]), reverse=True)
                differences = np.diff(lines[1])
                changed = True
        iteration_count += 1  # Increment iteration counter
        
    mesh_data[0] = process_mesh_data(mesh_data[0], automesher.min_cellsize, x_edges)
    mesh_data[1] = process_mesh_data(mesh_data[1], automesher.min_cellsize, y_edges)
    mesh_data[2] = process_mesh_data(mesh_data[2], automesher.min_cellsize_z, z_coords)

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
            # if any(mesh_data[i] == edge[0] for edge in unique_edges) and any(mesh_data[i+1] == edge[0] for edge in unique_edges):
            #     matching_edge_i = next(edge for edge in unique_edges if mesh_data[i] == edge[0])
            #     matching_edge_i_plus_1 = next(edge for edge in unique_edges if mesh_data[i+1] == edge[0])
            #     if abs(mesh_data[i+1] - mesh_data[i]) > 0:
            #         if hasattr(matching_edge_i[3], 'priority') and hasattr(matching_edge_i_plus_1[3], 'priority'):
            #             print('matching edges with same value found:', mesh_data[i+1]- mesh_data[i], matching_edge_i_plus_1, matching_edge_i)
            #             new_mesh_data.append(mesh_data[i])
            #             new_mesh_data.append(mesh_data[i+1])
            #         # if hasattr(matching_edge_i_plus_1, 'priority'):
            #         #     new_mesh_data.append(mesh_data[i+1])
            #         print('new mesh data:', new_mesh_data)
            #         continue
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
    # print(f"Mesh data changed: {new_mesh_data}")

    return mesh_data
