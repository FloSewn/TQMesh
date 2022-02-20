import sys
import numpy as np
from matplotlib import pyplot as plt

def main():

    if len(sys.argv) < 2:
        print("plot_mesh_benchmark.py <data>")
        sys.exit(1)

    data = np.loadtxt( sys.argv[1], delimiter=',', skiprows=2 )

    h        = data[:,0]
    L        = data[:,1]
    n_v      = data[:,2]
    n_t      = data[:,3]
    n_q      = data[:,4]
    n_ie     = data[:,5]
    n_be     = data[:,6]
    t_layer  = data[:,7]
    t_mesh   = data[:,8]
    t_smooth = data[:,9]

    n = n_t + n_q

    t_ref1 = n * 5.0E-5
    t_ref2 = n**2 * 2.0E-7

    fig, ax = plt.subplots(1,1,figsize=(4,4))
    ax.plot(n, t_mesh, marker='s', mfc='w', c='r', label='Mesh generation')
    ax.plot(n, t_smooth, marker='o', mfc='w', c='b', label='Mesh smoothing')

    ax.plot(n, t_ref1, marker='', c='k', ls='--', label=r'$\mathcal{O}(n)$')
    #ax.plot(n, t_ref2, marker='', c='k', ls='-.', label=r'$\mathcal{O}(n^2)$')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylim((5.0E-03,2.0E+02))

    ax.set_ylabel('Time in seconds')
    ax.set_xlabel('Number of elements')

    ax.legend()

    #plt.show()

    fig.tight_layout()
    export_fig_path = "BenchmarkPlot_Mesh.png"
    print("Writing {:}.".format(export_fig_path))
    fig.savefig(export_fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__': main()
