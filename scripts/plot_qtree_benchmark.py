import sys
import numpy as np
from matplotlib import pyplot as plt


def main():

    if len(sys.argv) < 2:
        print("plot_qtree_benchmark.py <data>")
        sys.exit(1)

    data = np.loadtxt( sys.argv[1], delimiter=',' )

    n  = data[:,0]
    t_qtree = data[:,1]
    t_brute = data[:,2]

    t_ref1 = n * 5.0E-6
    t_ref2 = n**2 * 2.0E-7

    fig, ax = plt.subplots(1,1,figsize=(4,4))
    ax.plot(n, t_qtree, marker='s', mfc='w', c='r', label='QTree')
    ax.plot(n, t_brute, marker='o', mfc='w', c='b', label='Brute force')

    ax.plot(n, t_ref1, marker='', c='k', ls='--', label=r'$\mathcal{O}(n)$')
    ax.plot(n, t_ref2, marker='', c='k', ls='-.', label=r'$\mathcal{O}(n^2)$')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylabel('Search time in seconds')
    ax.set_xlabel('Number of samples')

    ax.legend()

    #plt.show()

    fig.tight_layout()
    export_fig_path = "BenchmarkPlot_QTree.png"
    print("Writing {:}.".format(export_fig_path))
    fig.savefig(export_fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__': main()
