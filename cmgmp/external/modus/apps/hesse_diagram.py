#!/usr/bin/env python3
import os
import sys
from time import time

import glfw
import numpy as np

import modus as cm
import modus._modus as _cm


def hesse_diagram(graph, dot_file):
    with open(dot_file, 'w') as dot:
        G = []
        G.append('graph {')
        G.append('node [shape=none, label="", style=none, color="0 0 0", margin=0, width=0.2, height=0]')
        G.append('splines=false')
        G.append('layout=neato')
        # G.append('P')
        # G.append('E')

        # Create ranks manually
        rank_sep = 0.9
        node_sep = 0.8
        names = dict()
        for k in range(0, graph.dim() + 2):
            R_k = graph.rank(k)
            n_k = len(R_k)
            len_k = (n_k - 1) * node_sep
            for i in range(n_k):
                f = R_k[i]
                sv = '"' + f.sign_vector() + '"'
                if k == -1:
                    sv = '0'
                elif k == graph.dim() + 1:
                    sv = '1'
                x = -len_k/2 + i * node_sep
                y = -rank_sep * k

                svf = sv.replace('0', '00')
                svf = svf.replace('-', '−0')
                svf = svf.replace('+', '0−')
                svf = svf.replace('−', 't')
                svf = svf.replace('0', '−')
                svf = svf.replace('t', '0')

                if sv == '1':
                    svf = '0'
                if k == 0:
                    svf = '1'

                G.append(sv + ' [label=%s,pos="%f,%f!"]' % (svf,x,y))
                # G.append(sv + ' [pos="%f,%f!"]' % (x,y))

        # Create lattice
        for k in range(0, graph.dim() + 2):
            R_k = graph.rank(k)
            n_k = len(R_k)
            for i in range(n_k):
                f = R_k[i]
                f_sv = '"' + f.sign_vector() + '"'
                if k == -1:
                    f_sv = '0'
                for g in f.superfaces():
                    g_sv = '"' + g.sign_vector() + '"'
                    if g.rank() ==  graph.dim() + 1:
                        g_sv = '1'
                    G.append(f_sv + ' -- ' + g_sv)
        
        G.append('}')

        dot.write('\n'.join(G) + '\n')

def main():
    if sys.argv[1] == 'cube':
        # Create cube hyperplanes.
        np.zeros((6, 3))
        # Create graph.
        # Create Hesse diagram.
    elif sys.argv[1] == 'S2':
        # Create hyperplanes.
        H = np.identity(3)
        d = np.zeros((3,1))
        # Create arrangement.
        A = _cm.build_arrangement(H, d, 1e-8)
        A.update_sign_vectors(1e-8)
        # Print sign vectors by rank.
        # for k in range(0, A.dim() + 1):
        #     n = len(A.sign_vectors(k))
        #     R = A.rank(k)
        #     for i in range(n):
        #         print(R[i].sign_vector())
        # Create Hesse diagram.
        hesse_diagram(A, sys.argv[2])

if __name__=='__main__':
    main()