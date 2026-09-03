# AutoMesher for openEMS (easyMesh4openEMS)

A small, pragmatic **automatic mesh generator** for the Python bindings of **openEMS/CSXCAD**. It scans your geometry, infers where resolution is needed, and builds **smooth, graded** mesh lines in **x/y/z** — so you can focus on modeling, not manual meshing.

---

## Features

* **One‑call meshing:** `GenerateMesh(CSX, global_mesh_setup, ...)` builds mesh lines and writes them to the CSX grid.
* **Material‑aware smoothing:** mesh density adapts by region (air vs. dielectric/metal).
* **Edge, gap & diagonal handling:** special heuristics for metal edges, small gaps, and diagonal/circular segments.
* **Ports supported:** waveguide/MSL/lumped/rect ports are added to the meshing hints automatically.
* **PML-aware boundary meshing:** boundary conditions set through the wrapped `FDTD` object are tracked, and the extra cells required by `PML_N` boundaries are added automatically outside the requested air margin.

---

## Installation (local, not published on PyPI)

At the moment the mesher is **not published on PyPI**.  
You install it **locally from the project folder**:

1. Clone or download this repository.
2. Open a terminal in the `easyMesh4openEMS` folder (the folder that contains `pyproject.toml`).
3. Run:

```bash
pip install .
```

After that, you can import the package like this:

```python
from easyMesh import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh
```

> **Note:**  
> - The **distribution** is installed from the local folder with `pip install .`.  
> - The **Python package name** (for `import`) is `easyMesh`.

### Requirements

You need:

- Python ≥ 3.7  
- `numpy` (installed automatically via `pip install .`)  
- Working Python bindings for:
  - `openEMS`
  - `CSXCAD`

Make sure you can do e.g.:

```python
from openEMS import openEMS
from CSXCAD import ContinuousStructure
```

before using `easyMesh4openEMS`.

---

## Quick Start

```python
from openEMS import openEMS
from CSXCAD import ContinuousStructure
from easyMesh import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh

# 1) Create your FDTD/CSX as usual
FDTD = openEMS()
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)

# 2) Wrap CSX/FDTD before adding geometry, ports, or boundary conditions
primitives_mesh_setup = {}
properties_mesh_setup = {}
CSX = enhance_csx_for_auto_mesh(CSX, primitives_mesh_setup)
FDTD = enhance_FDTD_for_auto_mesh(FDTD, primitives_mesh_setup)

# 3) Set boundary conditions AFTER wrapping FDTD.
#    This lets easyMesh detect PML_N boundaries and add their cells automatically.
FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'])

# 4) Describe your global meshing intent
global_mesh_setup = {
    # Either provide start/stop OR f0/fc (frequencies in Hz)
    'start_frequency': 1e9,
    'stop_frequency': 3e9,
    # alternative: 'f0': 2e9, 'fc': 1e9,

    'drawing_unit': 1e-6,        # geometry unit (meters per drawing unit); 1e-6 => um units
    'mesh_resolution': 'medium', # one of: 'low'|'medium'|'high'|'very_high'

    # Optional wavelength reference for automatic wavelength-based mesh settings.
    # For antennas, use the antenna design frequency rather than necessarily f_stop.
    # 'target_frequency': 2e9,  # optional; for antennas, typically use the design frequency

    # Boundary order: [xmin, xmax, ymin, ymax, zmin, zmax]
    # boundary_distance is the AIR MARGIN only. PML_N cells are appended automatically.
    # 'auto' currently means lambda/3 using the selected target wavelength.
    'boundary_distance': ['auto', 'auto', 'auto', 'auto', 'auto', 'auto'],

    # Optional knobs
    'min_cellsize': None,        # computed from geometry if None
    'max_cellsize': None,        # defaults derived from resolution & epsilon
    'refined_cellsize': None,    # override nominal cellsize if desired

    # Heuristics/toggles
    'smooth_metal_edge': 'one_third_two_thirds', # useful for thin metal layers, Options: False, 'one_third_two_thirds', 'extra_lines'
    'use_circle_detection': False,               # detect circles for better angular resolution
    'handle_closely_placed_edges': True,  # if True, then mesher will try to handle close placed edges by merging them
}

# 5) Create your structure
# substrate = CSX.AddMaterial('RO5880', epsilon=substrate_epr)
# substrate.AddBox(start, stop, priority=10) etc....

# 6) (Optional) Provide per‑primitive/property hints
# Example: later, when you add geometry (if CSX is wrapped), hints are auto‑collected.
# You can also add entries manually, e.g. to restrict directions:
# primitives_mesh_setup[my_prim] = { 'dirs': 'xy', 'edges_only': False, 'metal_edge_res': None }

# 7) Generate and write mesh lines to CSX
GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup, properties_mesh_setup)

# 8) Continue with your usual openEMS workflow (run, post-processing ...)
```

> **Important for PML boundaries:** `enhance_FDTD_for_auto_mesh(...)` must be called **before** `FDTD.SetBoundaryCond(...)`. The wrapper records the boundary conditions for the mesher; boundary conditions that were set on the original `FDTD` object before wrapping cannot be detected retroactively.

---

## API Overview

### `GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup=None, properties_mesh_setup=None, **kw)`

Runs the full automeshing pipeline. It inspects geometry, computes meshlines per direction, smooths them, and writes them to `CSX.GetGrid()`.

**Parameters**

* `CSX`: your CSXCAD object.
* `global_mesh_setup` *(dict, required)* — see **Global parameters** below.
* `primitives_mesh_setup` *(dict, optional)* — per‑primitive options (see **Per‑primitive/property options**).
* `properties_mesh_setup` *(dict, optional)* — per‑property options (material/metal groups).


---

### `enhance_csx_for_auto_mesh(original_csx, primitives_mesh_setup)`

Wraps your `CSX` so that **any new primitive you add later** is automatically registered with default mesh hints. Useful when you don’t want to manually manage `primitives_mesh_setup`.

---

### `enhance_FDTD_for_auto_mesh(original_FDTD, primitives_mesh_setup)`

Wraps your `FDTD` so **ports** (`AddLumpedPort`, `AddWaveGuidePort`, `AddRectWaveGuidePort`, `AddMSLPort`) are auto‑registered with default port hints. The wrapper also records calls to `FDTD.SetBoundaryCond(...)`, which allows the mesher to detect `PML_N` boundaries and append the corresponding PML cells. Therefore, wrap the `FDTD` object **before** calling `SetBoundaryCond(...)`.

---

## Configuration reference

### Global parameters (keys for `global_mesh_setup`)

* **Frequencies** *(choose one pair)*

  * `start_frequency` + `stop_frequency`
  * `f0` + `fc`

  These provide the default wavelength reference used to derive wavelength-based mesh sizes and limits.

  * `target_frequency` *(optional)*: overrides that wavelength reference for the automatic meshing calculations. For antenna simulations, this can be the antenna design frequency (for example 868 MHz) even when `stop_frequency` is higher. It does not replace the required start/stop or `f0`/`fc` pair.

* **Units**

  * `drawing_unit` (default `1e-6`): meters per drawing unit.

* **Resolution preset**

  * `mesh_resolution`: `'low' | 'medium' | 'high' | 'very_high'`
    Controls the nominal cell size and number of intermediate lines used by the smoother. Rough intuition:

    * `low`   → coarser mesh
    * `medium` (default)
    * `high` / `very_high` → finer

* **Boundary distance and PML** *(optional)*

  * `boundary_distance`: six entries in the order `[xmin, xmax, ymin, ymax, zmin, zmax]`. Numeric values are in drawing units.
  * A numeric `boundary_distance` is interpreted as the **air margin only** between the structure/antenna and the start of the absorbing boundary region. Do **not** add the PML thickness manually.
  * `'auto'` currently uses a conservative air margin of `lambda/3`, based on the mesher wavelength. If `target_frequency` is set, that wavelength is based on `target_frequency`.
  * `None` means no added air margin on that side. With the current implementation, automatic PML-cell extension is only applied on sides whose resulting boundary distance is greater than zero; for a `PML_N` side, use `'auto'` or a positive numeric air margin if the PML cells should be appended automatically.
  * If a captured boundary condition is `PML_N` (for example `PML_8`), easyMesh appends the additional PML mesh cells **after** the air margin. The user should therefore pass only the desired air distance, such as an explicit `lambda/4`, not `lambda/4 + PML thickness`.

  > To make this work, call `enhance_FDTD_for_auto_mesh(...)` before `FDTD.SetBoundaryCond(...)`.

* **Direct overrides** *(optional)*

  * `refined_cellsize`: override nominal cell size if you need a specific resolution.
  * `min_cellsize`: minimal spacing allowed between lines (auto‑tightened if small gaps are found).
  * `max_cellsize`: maximum spacing allowed (scaled by dielectric `epsilon`).
  * `dirs`: a string subset of `'x'`, `'y'`, `'z'` to allow meshing in only certain directions.

* **Heuristics/toggles** *(optional)*

  * `smooth_metal_edge`: `'one_third_two_thirds'` or `'extra_lines'` or `False` enables extra lines near metal/port edges using a 1/3–2/3 rule.
  * `use_circle_detection`: `True/False` to detect full circles for refined angular meshing.

> The mesher automatically tightens `min_cellsize` and related limits when it finds close edges or small gaps.

---

## How it works 

Take a look for the examples in the Tutorials folder

<!-- 1. **Collect** current grid lines (if any) and existing geometry; clear the grid temporarily.
2. **Parse** primitives/ports into edge sets (vertical, horizontal, diagonal) and z‑boundaries.
3. **Derive** nominal `mesh_res`, `min_cellsize`, `max_cellsize` from frequency & materials.
4. **Refine** where needed: metal edges, small gaps, diagonals, circles/arcs, material transitions.
5. **Smooth** lines with `SmoothMeshLines` and **respect** per‑primitive/property `dirs` filters.
6. **Write back** final lines to the CSX grid. -->

---

## Troubleshooting

- **`ImportError: No module named 'openEMS'` or `'CSXCAD'`**  
  Make sure the Python bindings for openEMS and CSXCAD are installed and on your `PYTHONPATH`.

- **Mesh is too coarse / too fine**  
  - Adjust `mesh_resolution` (`low` ↔ `very_high`) or set `refined_cellsize`/`min_cellsize` directly.

- **PML cells are missing at a boundary**  
  - Make sure `FDTD = enhance_FDTD_for_auto_mesh(...)` is called before `FDTD.SetBoundaryCond(...)`.
  - Use a `PML_N` boundary string such as `PML_8` and use `'auto'` or a positive numeric `boundary_distance` on that side.
  - Treat `boundary_distance` as the air margin only; do not add the PML cells manually.

- **Very dense mesh near tiny gaps**  
  This is usually intentional to resolve small features.  
  If it is too dense:
  - Increase any minimum cell size parameters you use.
  - Slightly enlarge gaps in your geometry if physically acceptable.

---
## Acknowledgments

This project is funded by the Federal Ministry of Research, Technology and Space (BMFTR) under the DI-DEMICO project as part of the DE:Sign – Design Initiative for Microelectronics.

## Contributing

The commit masseges have to follow the rules defined in https://www.conventionalcommits.org/en/v1.0.0/. Tests and minimal examples are welcome.

---

