# Changelog

## [1.3.2] - 2024-06-10

### Fixed

- Add missing call to update of face-connectivity prior to mesh smoothing strategies, 
  which let to invalid meshing attempts in example 2 and example 6 (based on issue [#24](https://github.com/FloSewn/TQMesh/issues/24)) 
  ([`0ce9393`](https://github.com/FloSewn/TQMesh/commit/0ce9393))

### Changed

- Remove `ve_intersection_` attribute in `src/algorithm/FrontUpdate.h`, which prevented quad layer heights 
  to be removed adequately in the vicinity of small element sizes and which had no impact on the actual mesh generation
  ([#26](https://github.com/FloSewn/TQMesh/issues/26)) ([`1762bdd`](https://github.com/FloSewn/TQMesh/commit/1762bdd))
- Change namespace `TQMesh::TQAlgorithm` simply to `TQMesh`
- Add single include header files `src/algorithm/TQMesh.h` and `src/utils/CppUtils.h` to reduce include statements 
  within exampes / TQMesh application

### Added

- **convert2foam.py** to convert the TQMesh text-format to the OpenFOAM compatible mesh format ([#19](https://github.com/FloSewn/TQMesh/issues/19)) ([`df2a078`](https://github.com/FloSewn/TQMesh/commit/df2a078))
- New example **thin fracture**, which involves merging two meshes that model a thin fraction interface (based on [#26](https://github.com/FloSewn/TQMesh/issues/26)).
- **MeshChecker** class in `src/algorithm/MeshChecker.h`, which acts as interface for general mesh validity checks ([`7ca7320`](https://github.com/FloSewn/TQMesh/commit/7ca7320))
  (might be extened in future for mesh quality checks)

