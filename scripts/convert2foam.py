import sys, os, argparse
import numpy as np

def io_create_directories(export_dir):
    constant_dir = os.path.join(export_dir, 'constant')
    mesh_dir = os.path.join(constant_dir, 'polyMesh')
    if not os.path.exists(export_dir):
        os.mkdir(export_dir)
    if not os.path.exists(constant_dir):
        os.mkdir(constant_dir)
    if not os.path.exists(mesh_dir):
        os.mkdir(mesh_dir)
    open(os.path.join(export_dir, 'tqmesh.foam'), 'a').close()

    return mesh_dir

def io_write_foamfile_header(foamfile, filetype, note=None):

    class_type = {
            'owner'     : 'labelList',
            'neighbour' : 'labelList',
            'points'    : 'vectorField',
            'faces'     : 'faceList',
            'boundary'  : 'polyBoundaryMesh',
            }

    output_str  = 'FoamFile\n{\n'
    output_str += '\tversion     2.0;\n'.expandtabs(4)
    output_str += '\tformat      ascii;\n'.expandtabs(4)
    output_str += '\tclass       {:};\n'.format(class_type[filetype]).expandtabs(4)
    if note is not None:
        output_str += '\tnote        \"{:}\";\n'.format(note).expandtabs(4)
    output_str += '\tlocation    \"constant/polyMesh\";\n'.expandtabs(4)
    output_str += '\tobject      {:};\n'.format(filetype).expandtabs(4)
    output_str += '}\n'

    foamfile.write(output_str)




def io_clear_comments(lines):
    ''' Clears all comments from the input string
    '''
    for i in range(len(lines)-1, 0, -1):
        line = lines[i].replace(' ','')
        if line[0] == '#':
            lines.pop( i )

def io_read_vertices(mesh_id, lines):
    ''' Read vertices from the input string
    '''
    vertices = []
    start, n_vert = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the vertices
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "VERTICES":
            start  = i + 1
            n_vert = int(line[1])
            break

    # Read the vertices
    for i in range(start, start+n_vert):
        line = lines[i]
        (x,y) = ( float(s) for s in line.split(',') )
        vertices.append( (x,y) )

    return vertices

def io_read_interior_edges(mesh_id, lines):
    ''' Read interior edges from the input string
    '''
    edges, nbr_elements = [], []
    start, n_edge = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the boundary edges
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "INTERIOREDGES":
            start = i + 1
            n_edge = int(line[1])
            break

    # Read the boundary edges and their markers
    for i in range(start, start+n_edge):
        line = lines[i]
        (v1,v2,e1,e2) = ( int(s) for s in line.split(',') )
        edges.append( (v1,v2) )
        nbr_elements.append( (e1,e2) )

    return edges, nbr_elements


def io_read_boundary_edges(mesh_id, lines):
    ''' Read boundary edges from the input string
    '''
    edges = []
    bdry_elements = []
    markers = []
    start, n_edge = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the boundary edges
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "BOUNDARYEDGES":
            start = i + 1
            n_edge = int(line[1])
            break

    # Read the boundary edges and their markers
    for i in range(start, start+n_edge):
        line = lines[i]
        (v1,v2,e,m) = ( int(s) for s in line.split(',') )
        edges.append( (v1,v2) )
        bdry_elements.append( e )
        markers.append( m )

    return edges, bdry_elements, markers


def io_read_triangles(mesh_id, lines):
    ''' Read triangles from the input string
    '''
    tris, tri_colors = [], []
    start, n_tris = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the triangles
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "TRIANGLES":
            start = i + 1
            n_tris = int(line[1])
            break

    # Read triangles and their corresponding mesh ID
    for i in range(start, start+n_tris):
        line = lines[i]
        (v1,v2,v3,color) = ( int(s) for s in line.split(',') )
        tris.append( (v1,v2,v3) )
        tri_colors.append( color )

    return tris, tri_colors

def io_read_quads(mesh_id, lines):
    ''' Read quadrilaterals from the input string
    '''
    quads, quad_colors = [], []
    start, n_quads = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the quads
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "QUADS":
            start = i + 1
            n_quads = int(line[1])
            break

    # Read quads and their corresponding mesh ID
    for i in range(start, start+n_quads):
        line = lines[i]
        (v1,v2,v3,v4,color) = ( int(s) for s in line.split(',') )
        quads.append( (v1,v2,v3,v4) )
        quad_colors.append( color )

    return quads, quad_colors





class TQMesh:

    def __init__(self, mesh_id, lines):
        ''' Initialize a new TQMesh from a file, passed as an array
            of lines. The mesh_id corresponds to the index number
            of the mesh in the file.
        '''
        self.mesh_id = mesh_id
        self.vertices = np.array( io_read_vertices( mesh_id, lines ) )
        self.interior_edges, self.interior_elements = io_read_interior_edges( mesh_id, lines )
        self.boundary_edges, self.boundary_elements, self.boundary_markers = io_read_boundary_edges( mesh_id, lines )
        self.triangles, self.tri_colors = io_read_triangles( mesh_id, lines )
        self.quads, self.quad_colors = io_read_quads( mesh_id, lines )


class OpenFOAMMesh:

    def __init__(self, tqmesh, z_offset=1.0):

        points_2d = np.array(tqmesh.vertices)
        tris_2d = np.array(tqmesh.triangles)
        quads_2d = np.array(tqmesh.quads)
        intr_edges_2d = np.array(tqmesh.interior_edges)
        intr_elems_2d = np.array(tqmesh.interior_elements)
        bdry_edges_2d = np.array(tqmesh.boundary_edges)
        bdry_elems_2d = np.array(tqmesh.boundary_elements)
        bdry_marker_2d = np.array(tqmesh.boundary_markers)

        self._swap_intr_edges(intr_edges_2d, intr_elems_2d)

        bdry_ids = np.argsort(bdry_marker_2d)
        bdry_marker_2d = bdry_marker_2d[bdry_ids]
        bdry_edges_2d = bdry_edges_2d[bdry_ids]
        bdry_elems_2d = bdry_elems_2d[bdry_ids]

        n_intr_edges = intr_edges_2d.shape[0]
        n_tris, n_quads = tris_2d.shape[0], quads_2d.shape[0]

        self.points = self._init_points(points_2d, z_offset)
        self.faces  = self._init_faces(tris_2d, quads_2d, intr_edges_2d, bdry_edges_2d)
        self.owner = self._init_owner(intr_elems_2d, bdry_elems_2d, n_tris, n_quads)
        self.neighbour = self._init_neighbour(intr_elems_2d)
        self.boundaries = self._init_boundaries(bdry_marker_2d, n_intr_edges, n_tris, n_quads)


    def export(self, mesh_dir):
        self._write_points(mesh_dir)
        self._write_faces(mesh_dir)
        self._write_owner(mesh_dir)
        self._write_neighbour(mesh_dir)
        self._write_boundaries(mesh_dir)

    def _swap_intr_edges(self, intr_edges, intr_elements):
        ''' Swap direction of edges (and their neighbours),
        where the left neighboring element index is larger than the
        right neighboring element index '''
        swap = np.argwhere(intr_elements[:,0] < intr_elements[:,1])

        tmp = intr_elements[swap,0]
        intr_elements[swap,0] = intr_elements[swap,1]
        intr_elements[swap,1] = tmp

        tmp = intr_edges[swap,0]
        intr_edges[swap,0] = intr_edges[swap,1]
        intr_edges[swap,1] = tmp


    def _init_points(self, points_2d, z_offset=1.0):
        n_points = points_2d.shape[0]
        z0, z1 = np.zeros((n_points,1)), z_offset*np.ones((n_points,1))
        base_points = np.hstack((points_2d, z0))
        top_points = np.hstack((points_2d, z1))
        return np.vstack((base_points, top_points))

    def _init_faces(self, tris, quads, intr_edges, bdry_edges):
        ''' All faces are oriented in CCW direction towards the cell owner,
        such that the resulting normal vectors are pointing outwards from it
        '''
        n_offset = self.points.shape[0] // 2

        intr_faces_base = intr_edges
        intr_faces_top = intr_edges + n_offset
        intr_faces = np.hstack( (intr_faces_base[:,::-1], intr_faces_top) )

        bdry_faces_base = bdry_edges
        bdry_faces_top = bdry_edges + n_offset
        bdry_faces = np.hstack( (bdry_faces_base[:,::-1], bdry_faces_top) )

        base_faces = [ quad[::-1].tolist() for quad in quads ] \
                   + [ tri[::-1].tolist() for tri in tris ]

        top_faces = [ quad.tolist() for quad in (quads + n_offset) ] \
                  + [ tri.tolist() for tri in (tris + n_offset) ]

        faces = intr_faces.tolist() + bdry_faces.tolist() \
              + base_faces + top_faces

        return faces

    def _init_owner(self, intr_elements, bdry_elements, n_tris, n_quads):
        intr_owner = intr_elements[:,0].tolist()
        bdry_owner = bdry_elements.tolist()
        base_owner = np.arange(n_tris + n_quads).tolist()
        top_owner  = base_owner
        return intr_owner + bdry_owner + base_owner + top_owner

    def _init_neighbour(self, intr_elements):
        neighbours = intr_elements[:,1].tolist()
        return neighbours

    def _init_boundaries(self, bdry_marker, n_intr_edges, n_tris, n_quads):
        boundaries = {}

        # Add boundaries from 2D mesh
        offset = n_intr_edges
        count = 1
        for marker in np.unique(bdry_marker):
            name = 'boundary_{:}'.format(count)
            n_faces = np.sum(bdry_marker == marker)
            boundaries.update({name : (offset, n_faces)})
            offset += n_faces
            count += 1

        # Add additional base and top boundaries for 3D mesh
        n_faces = n_tris + n_quads
        boundaries.update({
            'boundary_{:}'.format(count) : (offset, n_faces)
        })
        offset += n_faces
        count += 1
        boundaries.update({
            'boundary_{:}'.format(count) : (offset, n_faces)
        })

        return boundaries


    def _write_points(self, mesh_dir):
        point_str = '\n{:}\n(\n'.format(self.points.shape[0])
        for (x,y,z) in self.points:
            point_str += '({:g} {:g} {:g})\n'.format(x,y,z)
        point_str += ')\n'

        with open(os.path.join(mesh_dir, 'points'), 'w') as file:
            io_write_foamfile_header(file, 'points')
            file.write(point_str)

    def _write_faces(self, mesh_dir):
        face_str = '\n{:}\n(\n'.format(len(self.faces))
        for face in self.faces:
            n = len(face)
            if n == 3:
                (v1,v2,v3) = face
                face_str += '{:d}({:d} {:d} {:d})\n'.format(n,v1,v2,v3)
            if n == 4:
                (v1,v2,v3,v4) = face
                face_str += '{:d}({:d} {:d} {:d} {:d})\n'.format(n,v1,v2,v3,v4)
        face_str += ')\n'

        with open(os.path.join(mesh_dir, 'faces'), 'w') as file:
            io_write_foamfile_header(file, 'faces')
            file.write(face_str)


    def _write_owner(self, mesh_dir):
        owner_str = '\n{:}\n(\n'.format(len(self.owner))
        for e in self.owner:
            owner_str += '{:d}\n'.format(e)
        owner_str += ')\n'

        with open(os.path.join(mesh_dir, 'owner'), 'w') as file:
            io_write_foamfile_header(file, 'owner')
            file.write(owner_str)

    def _write_neighbour(self, mesh_dir):
        neighbour_str = '\n{:}\n(\n'.format(len(self.neighbour))
        for e in self.neighbour:
            neighbour_str += '{:d}\n'.format(e)
        neighbour_str += ')\n'

        with open(os.path.join(mesh_dir, 'neighbour'), 'w') as file:
            io_write_foamfile_header(file, 'neighbour')
            file.write(neighbour_str)

    def _write_boundaries(self, mesh_dir):
        boundary_str = '\n{:}\n(\n'.format(len(self.boundaries))

        # Ensure that boundaries are written in sorted order
        boundary_data = []
        for name, (start_face, n_faces) in self.boundaries.items():
            boundary_data.append([name, (start_face, n_faces)])
        boundary_data.sort(key = lambda entry: entry[0])

        for name, (start_face, n_faces) in boundary_data:
            boundary_str += '    {:}\n    {{\n'.format(name)
            boundary_str += '        type            patch;\n'
            boundary_str += '        nFaces          {:};\n'.format(n_faces)
            boundary_str += '        startFace       {:};\n'.format(start_face)
            boundary_str += '    }\n'
        boundary_str += ')\n'

        with open(os.path.join(mesh_dir, 'boundary'), 'w') as file:
            io_write_foamfile_header(file, 'boundary')
            file.write(boundary_str)












def main(argv):

    parser = argparse.ArgumentParser(
        prog='convert2foam',
        description='Convert TQMesh.txt files to the openFOAM mesh format',
    )

    general = parser.add_argument_group('general arguments')
    general.add_argument(dest='input-meshfile', help='The TQMesh grid file to convert')
    general.add_argument(dest='export-prefix', help='The export prefix for the output openFOAM mesh')

    additional = parser.add_argument_group('additional arguments')
    additional.add_argument('-e', '--extrusion', dest='extrusion', default=1.0, type=float, required=False,
                            help='The applied extrusion height in the z-direction')

    args = vars(parser.parse_args())

    input_meshfile = args['input-meshfile']
    export_dir = args['export-prefix']
    z_offset = args['extrusion']

    mesh_dir = io_create_directories(export_dir)

    with open(input_meshfile, 'r') as f:
        lines = f.readlines()

    io_clear_comments( lines )

    tqmesh = TQMesh(0, lines)
    foammesh = OpenFOAMMesh(tqmesh, z_offset)

    foammesh.export(mesh_dir)


if __name__ == '__main__': main(sys.argv)
