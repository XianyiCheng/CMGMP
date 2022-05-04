#!/usr/bin/env python3
import os
from time import time

import glfw
import numpy as np

import modus as cm
import modus._modus as _cm
from modus.geometry import *

np.set_printoptions(suppress=True, precision=5, linewidth=300, sign=' ')

def gen_cs_modes():
    # Case: Box wall.
    for i in range(1, 6):
        system = cm.box_case(i)
        A, b = cm.build_normal_velocity_constraints(system.collider.manifolds)
        T, c = cm.build_tangential_velocity_constraints(system.collider.manifolds, 2)

        t_start = time()
        cs_modes, lattice, info = cm.enumerate_contact_separating_3d(A, b)
        print('case box', i)
        print('time python', time() - t_start)

        cs_modes = np.array(cs_modes).view(np.uint32)

        test_case_path = os.path.join('test', 'cs_modes', 'box-%d.npz' % i)
        print(test_case_path)
        print(cs_modes)

        np.savez_compressed(test_case_path, A=A, b=b, T=T, cs_modes=cs_modes)

def gen_boxbox_cases():
    # Case: Box-Box
    system = cm.box_box_case(1)
    A, b = cm.build_normal_velocity_constraints(system.collider.manifolds)
    T, c = cm.build_tangential_velocity_constraints(system.collider.manifolds, 2)

    t_start = time()
    cs_modes, lattice, info = cm.enumerate_contact_separating_3d(A, b)

    print('box-box-1')
    print('time python')

    print(cs_modes)

    cs_modes = np.array(cs_modes).view(np.uint32)

    print(cs_modes)

    test_case_path = os.path.join('test', 'cs_modes', 'box-box-%d.npz' % 1)

    np.savez_compressed(test_case_path, A=A, b=b, T=T, cs_modes=cs_modes)

def gen_ss_modes():
    # Case: Box wall.
    for i in range(1, 3):
        system = cm.box_case(i)
        A, b = cm.build_normal_velocity_constraints(system.collider.manifolds)
        T, c = cm.build_tangential_velocity_constraints(system.collider.manifolds, 2)

        t_start = time()
        cs_modes, lattice, info = cm.enumerate_contact_separating_3d(A, b)
        print('case box', i)
        print('time python', time() - t_start)
        # print(cs_modes)

        ss_modes = cm.enumerate_sliding_sticking_3d_exponential(A, b, T)

        for cs_mode, ss_mode in ss_modes.items():
            # print(cs_mode)
            # print(np.array(ss_mode))
            ss_modes[cs_mode] = np.array(ss_mode).view(np.uint32)
            # print(np.array(ss_mode).view(np.uint32))
            # print(np.array(ss_mode).view(cha))

        ss_modes['A'] = A
        ss_modes['b'] = b
        ss_modes['T'] = T

        test_case_path = os.path.join('test', 'box-%d.npz' % i)
        print(test_case_path)

        # np.savez('test/box-1.npz', 'A', A, 'b', b, 'T', T, **ss_modes)
        np.savez_compressed(test_case_path, **ss_modes)

        npz = np.load(test_case_path)

        # print(list(np.load(test_case_path).keys()))

        # print(ss_modes)

    # Case: 

if not glfw.init():
    raise RuntimeError('Error: Failed to initialize GLFW.')
win = glfw.create_window(400, 400, 'hello', None, None)
glfw.make_context_current(win)
glfw.swap_interval(0)

# gen_ss_modes()
# gen_cs_modes()
gen_boxbox_cases()

glfw.terminate()
