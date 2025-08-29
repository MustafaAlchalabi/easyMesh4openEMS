# AutoMesher for openEMS (easyMesh4openEMS)

A small, pragmatic **automatic mesh generator** for the Python bindings of **openEMS/CSXCAD**. It scans your geometry, infers where resolution is needed, and builds **smooth, graded** mesh lines in **x/y/z** — so you can focus on modeling, not manual meshing.

---

## Features

* **One‑call meshing:** `GenerateMesh(CSX, global_mesh_setup, ...)` builds mesh lines and writes them to the CSX grid.
* **Material‑aware smoothing:** mesh density adapts by region (air vs. dielectric/metal).
* **Edge, gap & diagonal handling:** special heuristics for metal edges, small gaps, and diagonal/circular segments.
* **Ports supported:** waveguide/MSL/lumped/rect ports are added to the meshing hints automatically.

---

## Installation

Clone the repository with

```bash
git clone https://github.com/MustafaAlchalabi/easyMesh4openEMS.git
```

to `<path of cloned repository>`. In order to import the mesher you have to tell python whese to search for it. For this purpose (as there is no module installation so far) either add 

```bash
export PYTHONPATH=/<path with cloned repository>/easyMesh4openEMS:$PYTHONPATH
```
to your `.bashrc` or `.profile`, etc. or append 

```python
import sys
sys.path.append('/<path with cloned repository>/easyMesh4openEMS')
```

in your Python script before importing. Then import:

```python
import sys
sys.path.append('/<path-to>/easyMesh4openEMS')
```

Then import:

```python
from easyMesher import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh
```

---

## Quick Start

```python
import openEMS
from CSXCAD import CSXCAD
from easyMesher import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh

# 1) Create your FDTD/CSX as usual
FDTD = openEMS()
CSX = ContinuousStructure()  # or however you initialize your project

# 2) Wrap CSX/FDTD so newly added primitives/ports are auto‑tracked
primitives_mesh_setup = {}
properties_mesh_setup = {}
CSX = enhance_csx_for_auto_mesh(CSX, primitives_mesh_setup={})
FDTD = enhance_FDTD_for_auto_mesh(FDTD, primitives_mesh_setup={})

# 3) Describe your global meshing intent
global_mesh_setup = {
    # Either provide start/stop OR f0/fc (unit = drawing units)
    'start_frequency': 1e9,
    'stop_frequency': 3e9,
    # alternative: 'f0': 2e9, 'fc': 1e9,

    'drawing_unit': 1e-6,        # geometry unit (meters per drawing unit); 1e-6 => um units
    'mesh_resolution': 'medium', # one of: 'low'|'medium'|'high'|'very_high'

    # Optional knobs
    'min_cellsize': None,        # computed from geometry if None
    'max_cellsize': None,        # defaults derived from resolution & epsilon
    'refined_cellsize': None,    # override nominal cellsize if desired

    # Heuristics/toggles
    'smooth_metal_edge': 'one_third_two_thirds', # useful for thin metal layers, Options: False, 'one_third_two_thirds', 'extra_lines'
    'use_circle_detection': False,               # detect circles for better angular resolution

}

# 4) create your structure

# substrate = CSX.AddMaterial('RO5880', epsilon=substrate_epr) 
# substrate.AddBox(start, stop, priority=10) etc.... 

# 5) (Optional) Provide per‑primitive/property hints

# Example: later, when you add geometry (if CSX is wrapped), hints are auto‑collected.
# You can also add entries manually, e.g. to restrict directions:
# primitives_mesh_setup[my_prim] = { 'dirs': 'xy', 'edges_only': False, 'metal_edge_res': None }

# 6) Generate and write mesh lines to CSX
GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup, properties_mesh_setup)

# 7) Continue with your usual openEMS workflow (run, post-processing ...)
```

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

Wraps your `FDTD` so **ports** (`AddLumpedPort`, `AddWaveGuidePort`, `AddRectWaveGuidePort`, `AddMSLPort`) are also auto‑registered with default port hints.

---

## Configuration reference

### Global parameters (keys for `global_mesh_setup`)

* **Frequencies** *(choose one pair)*

  * `start_frequency` + `stop_frequency`
  * `f0` + `fc`

  These determine a wavelength used to derive a nominal cell size and limits.

* **Units**

  * `drawing_unit` (default `1e-6`): meters per drawing unit.

* **Resolution preset**

  * `mesh_resolution`: `'low' | 'medium' | 'high' | 'very_high'`
    Controls the nominal cell size and number of intermediate lines used by the smoother. Rough intuition:

    * `low`   → coarser mesh
    * `medium` (default)
    * `high` / `very_high` → finer

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

* **ImportError for openEMS/CSXCAD**: ensure their Python modules are installed and on `PYTHONPATH`.
* **Mesh looks too coarse/fine**: adjust `mesh_resolution` or set `refined_cellsize`/`min_cellsize` directly.
* **No lines in a direction**: check global/primitive `dirs` filters.
* **Very dense mesh near tiny gaps**: expected; increase `min_cellsize` or simplify geometry.

--- -->

## Contributing

The commit masseges have to follow the rules defined in https://www.conventionalcommits.org/en/v1.0.0/. Tests and minimal examples are welcome.

---

