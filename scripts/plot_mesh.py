import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import collections as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, os


class TQMesh:

    def __init__(self, mesh_id, lines):
        ''' Initialize a new TQMesh from a file, passed as an array
            of lines. The mesh_id corresponds to the index number
            of the mesh in the file.
        '''
        self.mesh_id = mesh_id
        self.vertices = np.array( io_read_vertices( mesh_id, lines ) )
        self.front, self.front_markers = io_read_front( mesh_id, lines, self.vertices )
        self.boundary, self.boundary_elements, self.boundary_markers = io_read_boundary( mesh_id, lines, self.vertices )
        self.triangles, self.tri_colors = io_read_triangles( mesh_id, lines, self.vertices )
        self.quads, self.quad_colors = io_read_quads( mesh_id, lines, self.vertices )
        self.size_func = np.array( io_read_sizefunction( mesh_id, lines ) )

    def get_extent(self):
        ''' Return the mesh extents
        '''
        xy_max = np.max( self.vertices )
        xy_min = np.min( self.vertices )

        return (xy_min, xy_max)

    def get_element_colors(self):
        ''' Returns an array of all unique element color IDs
        '''
        t_colors = np.unique( self.tri_colors ).astype(int)
        q_colors = np.unique( self.quad_colors ).astype(int)
        return np.unique( np.hstack( (t_colors, q_colors) ) ).tolist()

    def plot_vertices(self, ax, indices=False):
        ''' Plot all mesh vertices
        '''
        ax.scatter(self.vertices[:,0], self.vertices[:,1],
                   marker='o', color='k', s=3)

        if indices:
            for i, v in enumerate(self.vertices):
                ax.text(v[0], v[1], str(i))

    def plot_front(self, ax):
        ''' Plot the advancing front of the mesh
        '''
        front_collection = mc.LineCollection( self.front, colors=[(.7,.4,.4)], lw=1.0, ls='--' )
        ax.add_collection( front_collection )

    def plot_boundaries(self, ax):
        ''' Plot the mesh boundary edges
        '''
        bdry_collection = mc.LineCollection( self.boundary, colors=[(.7,.4,.4)], lw=2.0, ls='--' )
        ax.add_collection( bdry_collection )


    def plot_triangles(self, ax, element_colors=[], indices=False):
        ''' Plot all mesh triangles
        '''
        n_colors = len(element_colors)

        if n_colors > 0:
            colors = plt.cm.Set1( np.linspace(0.0,1.0,n_colors) )
            elem_colors = colors[self.tri_colors]
        else:
            elem_colors = ['None']

        tri_collection = mc.PolyCollection(self.triangles, edgecolors=['k'],
                                           facecolors=elem_colors, lw=0.5 )
        ax.add_collection( tri_collection )

        if indices:
            for i, t in enumerate(self.triangles):
                c = np.mean(t, axis=0)
                ax.text(c[0], c[1], str(i+len(self.quads)), c=(.3,.3,.3))


    def plot_quads(self, ax, element_colors=[], indices=False):
        ''' Plot all mesh quads
        '''
        n_colors = len(element_colors)

        if n_colors > 0:
            colors = plt.cm.Set1( np.linspace(0.0,1.0,n_colors) )
            elem_colors = colors[self.quad_colors]
        else:
            elem_colors = ['None']

        quad_collection = mc.PolyCollection(self.quads, edgecolors=['k'],
                                            facecolors=elem_colors, lw=0.5 )
        ax.add_collection( quad_collection )

        if indices:
            for i, q in enumerate(self.quads):
                c = np.mean( q, axis=0 )
                ax.text(c[0], c[1], str(i), c=(.3,.3,.3))

    def plot_sizefunction(self, ax):
        tris = np.array(self.triangles)
        s = self.size_func
        x = self.vertices[:,0]
        y = self.vertices[:,1]
        cont = ax.tricontourf(x, y, s, levels=50, cmap='Spectral')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cont, cax=cax)



def io_gather_mesh_ids(lines):
    ''' Returns a list of mesh-IDs that are defined in the mesh file.
        This is used to estimate the number of meshes in the file.
    '''
    mesh_ids = []

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            mesh_ids.append( int(line[1]) )

    return mesh_ids

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

def io_read_front(mesh_id, lines, vertices):
    ''' Read advancing front edges from the input string
    '''
    edges = []
    markers = []
    start, n_edge = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the front edges
    # that are associated to the given mesh-ID
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "FRONT":
            start = i + 1
            n_edge = int(line[1])
            break

    # Read the front edges and their markers
    for i in range(start, start+n_edge):
        line = lines[i]
        (v1,v2,m) = ( int(s) for s in line.split(',') )
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        edges.append( [(x1,y1),(x2,y2)] )
        markers.append( m )

    return edges, markers

def io_read_boundary(mesh_id, lines, vertices):
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
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        edges.append( [(x1,y1),(x2,y2)] )
        bdry_elements.append( e )
        markers.append( m )

    return edges, bdry_elements, markers


def io_read_triangles(mesh_id, lines, vertices):
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
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        (x3, y3) = vertices[v3]
        tris.append( [(x1,y1),(x2,y2),(x3,y3)] )
        tri_colors.append( color )

    return tris, tri_colors

def io_read_quads(mesh_id, lines, vertices):
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
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        (x3, y3) = vertices[v3]
        (x4, y4) = vertices[v4]
        quads.append( [(x1,y1),(x2,y2),(x3,y3),(x4,y4)] )
        quad_colors.append( color )

    return quads, quad_colors

def io_read_qtree(lines):
    ''' Read qtree from the input string
    '''
    quads = []
    start, n_quads = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "QTREE-LEAFS":
            start = i + 1
            n_quads = int(line[1])
            break

    for i in range(start, start+n_quads):
        line = lines[i]
        (c1,c2,scale,size) = ( float(s) for s in line.split(',') )
        halfscale = 0.5*scale
        (x1, y1) = (c1-halfscale,c2-halfscale)
        (x2, y2) = (c1+halfscale,c2-halfscale)
        (x3, y3) = (c1+halfscale,c2+halfscale)
        (x4, y4) = (c1-halfscale,c2+halfscale)
        quads.append( [(x1,y1),(x2,y2),(x3,y3),(x4,y4)] )

    return quads

def io_read_sizefunction(mesh_id, lines):
    ''' Read the size function from the input string
    '''
    size_fun = []
    start, n_lines = 0, 0
    mesh_passed = False

    # Estimate start and ending line to read the size function
    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "MESH":
            if int(line[1]) == mesh_id:
                mesh_passed = True
        if mesh_passed and line[0] == "SIZEFUNCTION":
            start  = i + 1
            n_lines = int(line[1])
            break

    # Read the size function values
    for i in range(start, start+n_lines):
        line = lines[i]
        size_fun.append( float(line) )

    return size_fun

def io_read_sizefunction_DEPRECATED(lines):
    ''' Read the size function from the input string
    '''
    start, n_lines = 0, 0

    xy_min, xy_max = [0.0,0.0], [0.0,0.0]
    Nx, Ny = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "SIZE-FUNCTION":
            start = i + 1
            xy_min[0] = float(line[1])
            xy_min[1] = float(line[2])
            xy_max[0] = float(line[3])
            xy_max[1] = float(line[4])
            Nx = int(line[5])
            Ny = int(line[6])
            n_lines = int(np.ceil(Nx * Ny / 10))

    sizefun = []
    for i in range(start, start+n_lines):
        line = lines[i]
        sizefun += [ float(s) for s in line.split(',') ]


    x = np.linspace( xy_min[0], xy_max[0], Nx )
    y = np.linspace( xy_min[1], xy_max[1], Ny )

    X,Y = np.meshgrid( x, y )
    Z = np.array( sizefun ).reshape( (Ny,Nx) )

    return X, Y, Z



def main():
    if len(sys.argv) < 2:
        print("plot_mesh.py <mesh.txt> <flags>")
        sys.exit(1)

    file_name = sys.argv[1]

    with open(file_name, 'r') as f:
        lines = f.readlines()

    io_clear_comments( lines )

    mesh_ids = io_gather_mesh_ids( lines )

    meshes = []
    for i in mesh_ids:
        meshes.append( TQMesh(i, lines) )

    # Read the quadtree structure
    qtree = io_read_qtree( lines )

    # Initialize face colors for mesh entities
    element_colors = []
    if '-c' in sys.argv:
        for m in meshes:
            element_colors += m.get_element_colors()

    # Create the mesh
    fig, ax = plt.subplots(1,1,dpi=200)

    if '-q' in sys.argv:
        qtree_collection = mc.PolyCollection( qtree, edgecolors=['b'], facecolors=['None'], lw=1, ls='--' )
        ax.add_collection( qtree_collection )


    xy_min, xy_max = [], []

    for i, mesh in enumerate(meshes):
        xy_min_i, xy_max_i = mesh.get_extent()
        xy_min.append( xy_min_i )
        xy_max.append( xy_max_i )

        if '-v' in sys.argv:
            mesh.plot_vertices(ax, indices=True)

        if '-f' in sys.argv:
            mesh.plot_front(ax)

        if '-b' in sys.argv:
            mesh.plot_boundaries(ax)

        if '-s' in sys.argv:
            mesh.plot_sizefunction(ax)

        mesh.plot_triangles(ax, element_colors, '-e' in sys.argv)
        mesh.plot_quads(ax, element_colors, '-e' in sys.argv)


    xy_min = np.min(np.vstack( xy_min ).T)
    xy_max = np.max(np.vstack( xy_max ).T)

    ax.set_xlim( (xy_min,xy_max))
    ax.set_ylim( (xy_min,xy_max))

    ax.set_axis_off()
    ax.set_aspect('equal')



    fig.tight_layout()
    export_fig_path = "MeshPlot.png"
    print("Writing {:}.".format(export_fig_path))
    fig.savefig(export_fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)






if __name__ == '__main__': main()
