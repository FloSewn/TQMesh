# TQMesh: Simplifying Two-Dimensional Mesh Generation
-----------------------
<img src="doc/Example_6.png" alt="TQMesh-Example-6" width="750"/>

**TQMesh** is a C++ library and standalone application for generating two-dimensional meshes 
composed of triangular and quadrilateral elements. It is intended for use in scientific and 
engineering applications such as finite element and finite volume methods, where control over 
element size, topology, and boundary representation is required.

The mesh generation is based on an advancing-front approach. Users describe the computational 
domain by boundary segments and optional interior constraints, while element sizes can be controlled 
through boundary attributes or user-defined size functions.

---

## Features

- **Advancing‑front mesh generation**  
  Meshes are generated from user‑defined boundary segments, providing explicit control over the domain geometry.

- **Triangular and quadrilateral elements**  
  The mesher supports pure triangular meshes as well as mixed meshes with quadrilateral elements, including quad layers near boundaries.

- **Local mesh size control**  
  Element sizes can be influenced by boundary‑based sizing factors or by user‑defined size functions.

- **Boundary and interior constraints**  
  Boundary markers, fixed interior vertices, and fixed interior edges can be used to guide the meshing process.

- **Multiple meshes and merging**  
  Several meshes can be generated and merged while preserving conformity along shared boundaries.

---

## Used in Research

This library has been used in the following academic work:

- **Quadrilateral Mesh-Based Reactive Transport Modeling in Non-Orthogonal Random Fracture-Matrix Systems**  
  Su, Danyang, et al.  
  SSRN 5251117

---

## Installation

**TQMesh** is a header‑only library. To use the library components, include the headers located in:

- `src/algorithm`
- `src/utils`

The directory `src/app` contains the source code for the standalone application.

### Build and install the application

```sh
git clone https://github.com/FloSewn/TQMesh
cd TQMesh
mkdir build
cd build
cmake ..
make install
```

To select a specific compiler or build type, pass the appropriate options to CMake, for example:

```sh
cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug ..
```
---

## Library Usage

Example programs demonstrating how to use the **TQMesh** library are provided in `src/examples`.

After building the project, the examples can be executed using:

```sh
./bin/run_examples <example-id>
```

where `<example-id>` refers to one of the available example configurations.

---

## Application Usage

The **TQMesh** application reads simple parameter files describing the mesh geometry and options. Example input files are provided in the `input/` directory.

Run the application from the project root using:

```sh
./bin/TQMesh <input-file>
```

---

## Examples and Capabilities

Mesh generation in **TQMesh** is based on a small set of global parameters, 
boundary definitions (nodes and edges), and optional interior constraints. Boundary edges can be assigned 
markers that are propagated to the resulting mesh.

The following examples illustrate typical use cases.

### Local mesh refinement

Element sizes can be controlled through boundary‑based sizing factors or user‑defined size functions.

<img src="doc/simple_triangular_mesh.png" alt="Simple triangular mesh" width="300"/>
<img src="doc/thin_fracture.png" alt="Thin fracture mesh" width="338"/>

### Quadrilateral layers and refinement

Quadrilateral elements and quad layers near boundaries can be generated using a paving‑style approach, followed by optional quad refinement.

<img src="doc/square_in_channel.png" alt="Quad-dominated mesh" width="650"/>
<img src="doc/MeshRefinement_1.png" alt="Quad refinement" width="250"/>

### Boundary shapes

Meshes can be created from predefined boundary shapes such as rectangles, circles, and triangles.

<img src="doc/boundary_shapes.png" alt="Boundary shapes" width="250"/>
<img src="doc/MeshRefinement_2.png" alt="Mesh refinement example" width="250"/>
<img src="doc/MeshRefinement_3.png" alt="Mesh refinement example" width="250"/>

### Element coloring and interior control

Each element carries an associated color value that can be modified during meshing. Fixed interior vertices may be defined to locally influence mesh resolution.

<img src="doc/fixed_vertices.png" alt="Fixed interior vertices" width="250"/>

### Merging meshes

Multiple meshes can be generated and merged while preserving boundary conformity, enabling different regions to be assigned different element attributes.

<img src="doc/merge_meshes.png" alt="Merged meshes" width="250"/>
<img src="doc/multiple_meshes.png" alt="Multiple meshes" width="250"/>

### Importing boundaries from CSV

Boundary definitions can be imported from CSV files to support automated meshing workflows.

<img src="doc/airfoil.png" alt="Airfoil mesh" width="350"/>

### Fixed interior edges

Interior edges can be specified explicitly to guide the advancing‑front process.

<img src="doc/fixed_edges.png" alt="Fixed interior edges" width="500"/>

---

## Output Formats

Meshes can be exported in:

- **VTU** format (readable by ParaView)
- A **plain text** format

The text format includes vertex coordinates, element connectivity, boundary information, adjacency data, and size function values.

---

## OpenFOAM Conversion

A helper script is provided to convert **TQMesh** output into a format suitable for OpenFOAM by extruding the 2D mesh in the z‑direction:

```sh
python scripts/convert2foam.py [-e EXTRUSION] Mesh.txt export-prefix
```

---

## Plotting Meshes

A Python script for visualizing meshes is available at `scripts/plot_mesh.py`:

```sh
python scripts/plot_mesh.py Mesh.txt (-s -c -v -e -b -f)
```

Optional flags enable visualization of size functions, element colors, indices, and boundaries.

<img src="doc/ExampleMesh_SizeFunction.png" alt="Mesh visualization" width="350"/>

---

## Tests and Benchmarks

**TQMesh** uses a customized QuadTree data structure to store and access mesh entities. Benchmark plots illustrating performance characteristics are provided below.

<img src="doc/BenchmarkPlot_QTree.png" alt="QuadTree benchmark" width="400"/>
<img src="doc/BenchmarkPlot_Mesh.png" alt="Mesh benchmark" width="400"/>

---

## Roadmap

- Improved triangle‑to‑quad morphing
- Boundary definitions via splines
- Extended documentation and tests

Contributions are welcome.

---

## Third‑Party Libraries

- exprtk — C++ Mathematical Expression Toolkit Library

---

## References

- O'Rourke, J. *Computational Geometry in C*. Cambridge University Press, 1998.
- Shewchuk, J. R. *Lecture Notes on Delaunay Mesh Generation*, 2012.
- Lo, D. S. H. *Finite Element Mesh Generation*. CRC Press, 2014.
- Blazek, J. *Computational Fluid Dynamics: Principles and Applications*. Butterworth‑Heinemann, 2015.
- Zhou, T., & Shimada, K. *An Angle‑Based Approach to Two‑Dimensional Mesh Smoothing*, IMR 2000.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

