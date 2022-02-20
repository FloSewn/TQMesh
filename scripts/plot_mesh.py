import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import collections as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, os


def io_clear_comments(lines):
    ''' Clears all comments from the input string
    '''
    for i in range(len(lines)-1, 0, -1):
        line = lines[i].replace(' ','')
        if line[0] == '#':
            lines.pop( i )

def io_read_vertices(lines):
    ''' Read vertices from the input string
    '''
    vertices = []
    start, n_vert = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "VERTICES":
            start  = i + 1
            n_vert = int(line[1])
            break

    for i in range(start, start+n_vert):
        line = lines[i]
        (x,y) = ( float(s) for s in line.split(',') )
        vertices.append( (x,y) )

    return vertices

def io_read_front(lines, vertices):
    ''' Read advancing front edges from the input string
    '''
    edges = []
    markers = []
    start, n_edge = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "FRONT":
            start = i + 1
            n_edge = int(line[1])
            break

    for i in range(start, start+n_edge):
        line = lines[i]
        (v1,v2,m) = ( int(s) for s in line.split(',') )
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        edges.append( [(x1,y1),(x2,y2)] )
        markers.append( m )

    return edges, markers

def io_read_boundary(lines, vertices):
    ''' Read boundary edges from the input string
    '''
    edges = []
    markers = []
    start, n_edge = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "BOUNDARYEDGES":
            start = i + 1
            n_edge = int(line[1])
            break

    for i in range(start, start+n_edge):
        line = lines[i]
        (v1,v2,m) = ( int(s) for s in line.split(',') )
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        edges.append( [(x1,y1),(x2,y2)] )
        markers.append( m )

    return edges, markers


def io_read_triangles(lines, vertices):
    ''' Read triangles from the input string
    '''
    tris = []
    start, n_tris = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "TRIANGLES":
            start = i + 1
            n_tris = int(line[1])
            break

    for i in range(start, start+n_tris):
        line = lines[i]
        (v1,v2,v3) = ( int(s) for s in line.split(',') )
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        (x3, y3) = vertices[v3]
        tris.append( [(x1,y1),(x2,y2),(x3,y3)] )

    return tris

def io_read_quads(lines, vertices):
    ''' Read quadrilaterals from the input string
    '''
    quads = []
    start, n_quads = 0, 0

    for i, line in enumerate(lines):
        line = line.split(" ")
        if line[0] == "QUADS":
            start = i + 1
            n_quads = int(line[1])
            break

    for i in range(start, start+n_quads):
        line = lines[i]
        (v1,v2,v3,v4) = ( int(s) for s in line.split(',') )
        (x1, y1) = vertices[v1]
        (x2, y2) = vertices[v2]
        (x3, y3) = vertices[v3]
        (x4, y4) = vertices[v4]
        quads.append( [(x1,y1),(x2,y2),(x3,y3),(x4,y4)] )

    return quads

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

def io_read_sizefunction(lines):
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
    vertices = np.array( io_read_vertices( lines ) )
    front, front_markers = io_read_front( lines, vertices )
    boundary, boundary_markers = io_read_boundary( lines, vertices )
    triangles = io_read_triangles( lines, vertices )
    quads = io_read_quads( lines, vertices )
    qtree = io_read_qtree( lines )
    X,Y,Z = io_read_sizefunction( lines )




    fig, ax = plt.subplots(1,1,dpi=200)

    if '-q' in sys.argv:
        qtree_collection = mc.PolyCollection( qtree, edgecolors=['b'], facecolors=['None'], lw=1, ls='--' )
        ax.add_collection( qtree_collection )

    if '-s' in sys.argv:
        size_fun = ax.contourf(X,Y,Z, levels=50, cmap='Spectral')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(size_fun, cax=cax)


    quad_collection = mc.PolyCollection( quads, edgecolors=['k'], facecolors=['None'], lw=0.5 )
    ax.add_collection( quad_collection )

    tri_collection = mc.PolyCollection( triangles, edgecolors=['k'], facecolors=['None'], lw=0.5 )
    ax.add_collection( tri_collection )

    if '-i' in sys.argv:

        for i, q in enumerate(quads):
            c = np.mean( q, axis=0 )
            ax.text(c[0], c[1], str(i), c=(.3,.3,.3))

        for i, t in enumerate(triangles):
            c = np.mean( t, axis=0 )
            ax.text(c[0], c[1], str(i+len(quads)), c=(.3,.3,.3))

    # Plot the mesh boundaries

    # Plot remaining advancing front edges
    if '-f' in sys.argv:
        front_collection = mc.LineCollection( front, colors=[(.7,.4,.4)], lw=1.0, ls='--' )
        ax.add_collection( front_collection )

    # Plot the boundary edges
    if '-b' in sys.argv:
        bdry_collection = mc.LineCollection( boundary, colors=[(.7,.4,.4)], lw=2.0, ls='--' )
        ax.add_collection( bdry_collection )

    if len(vertices) > 0:
        #ax.scatter( vertices[:,0], vertices[:,1], marker='o', color='k', s=3 )

        if '-v' in sys.argv:
            for i, v in enumerate(vertices):
                ax.text(v[0], v[1], str(i))

    ax.set_aspect('equal')

    xy_max = np.max( vertices )
    xy_min = np.min( vertices )
    ax.set_xlim( (xy_min,xy_max))
    ax.set_ylim( (xy_min,xy_max))

    ax.set_ylim( (np.min(vertices[:,1]), np.max(vertices[:,1]) ) )

    ax.set_axis_off()

    #plt.show()

    fig.tight_layout()
    export_fig_path = "MeshPlot.png"
    print("Writing {:}.".format(export_fig_path))
    fig.savefig(export_fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)






if __name__ == '__main__': main()
