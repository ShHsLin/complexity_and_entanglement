import numpy as np
import pickle
import os, sys
sys.path.append('..')

from quantumCircuit.methods.aux_methods import *
from quantumCircuit.methods.states import *
from quantumCircuit.methods.gates import *

from tensor_network_functions.mps_func import MPS_2_state, state_2_MPS, lpr_2_plr, plr_2_lpr
from tensor_network_functions.mps_func import init_mps, right_canonicalize, MPS_compression_variational

import tensor_network_functions.mps_func as mps_func
import tensor_network_functions.circuit_func as circuit_func

from circuit.ed import get_H_Ising

from scipy.sparse.linalg import eigsh, expm, expm_multiply

import time

from circuit.parse_args import parse_args


if __name__ == "__main__":
    
    create_path("HPC_data")
    create_path("HPC_data/ED_evolution")
    create_path("HPC_data/circuits_evolution")

    args = parse_args()

    option = args.option
    filename = args.filename

    with open(filename) as f:
        line = f.readlines()[option-1]

    args = parse_args(line)
    L = args.L
    g = args.g
    h = args.h
    T = args.T
    depth = args.depth
    N_iter = args.N_iter  # should be a multiple of 100

    save_filename = f"Ising_evolution_L{L}_g{g}_h{h}_T{T}".replace(".","-")

    print(f"Starting L={L}, g={g}, h={h}, T={T}, depth={depth}")

    N = L  # I should probably just go with L

    t1 = time.time()
    #H = get_H_Ising(g, h, 1., L)

    H = sp.csr_matrix((2**N, 2**N))  # using sparse is far faster!!!
    for ii in range(N-1):
        H = H + Gate.xx([ii,ii+1]).toSparseArray(N)
    for ii in range(N):
        H = H + g*Gate.z(ii).toSparseArray(N) + h*Gate.x(ii).toSparseArray(N)

    t2 = time.time()


    print(f"Hamiltonian constructed in {t2-t1} seconds")

    t1 = time.time()
    V = np.zeros((2**N,1), dtype=complex)
    V[0,0] = 1.
    V = expm_multiply(-1j*T*H, V)
    t2 = time.time()

    save_obj(V, "HPC_data/ED_evolution/"+save_filename+"_evolution")

    print(f"Evolution solved in {t2-t1} seconds")
