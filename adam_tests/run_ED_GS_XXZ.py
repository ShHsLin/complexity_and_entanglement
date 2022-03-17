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
    create_path("HPC_data/ED")
    create_path("HPC_data/circuits")

    args = parse_args()

    option = args.option
    filename = args.filename

    with open(filename) as f:
        line = f.readlines()[option-1]

    args = parse_args(line)
    L = args.L
    delta = args.delta
    N_iter = args.N_iter  # should be a multiple of 100

    save_filename = f"XXZ_GS_L{L}_delta{delta}".replace(".","-")

    print(f"Starting L={L}, delta={delta}")

    N = L  # I should probably just go with L

    t1 = time.time()
    #H = get_H_Ising(g, h, 1., L)

    H = sp.csr_matrix((2**N, 2**N))  # using sparse is far faster!!!
    for ii in range(N-1):
        H = H + Gate.xx([ii,ii+1]).toSparseArray(N) 
        H = H + Gate.yy([ii,ii+1]).toSparseArray(N)
        H = H + delta*Gate.zz([ii,ii+1]).toSparseArray(N)

    t2 = time.time()


    print(f"Hamiltonian constructed in {t2-t1} seconds")

    t1 = time.time()
    # for ground state
    _,V = eigsh(H, 1, which='SA')
    V = V.reshape((-1,1)) * (1 + 0*1j)  # make vector complex for MPS!
    t2 = time.time()

    save_obj(V, "HPC_data/ED/"+save_filename+"_ED")

    print(f"Eigenproblem solved in {t2-t1} seconds")
