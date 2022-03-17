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


def initial_circuit(L, depth):
    #np.random.seed(1)
    np.set_printoptions(linewidth=2000, precision=5, threshold=4000)

    ############### PARAMETERS INITIALIZATION #################
    ## We should test identity initialization and
    ## trotterization initialization

    idx = 0
    my_circuit = []

    ################# CIRCUIT INITIALIZATION  ######################
    # product_state = [np.array([1., 0.]).reshape([2, 1, 1]) for i in range(L)]
    product_state = np.zeros([2**L], dtype=np.complex128)
    product_state[0] = 1.
    for dep_idx in range(depth):

        random_layer = []
        for idx in range(L-1):
            if (idx + dep_idx) % 2 == 0:
                #random_layer.append(np.eye(4).reshape([2,2,2,2]).reshape([4, 4]))
                random_layer.append(circuit_func.random_2site_U(2).reshape([4, 4]))  # random close to identity
            else:
                random_layer.append(np.eye(4).reshape([2,2,2,2]).reshape([4, 4]))

        my_circuit.append(random_layer)
        current_depth = dep_idx + 1

    return my_circuit


def run_optimization(target_state, my_circuit, L, depth, N_iter):
    np.set_printoptions(linewidth=2000, precision=5, threshold=4000)

    tol = 1e-12
    cov_crit = tol * 0.1
    max_N_iter = N_iter

    product_state = np.zeros([2**L], dtype=np.complex128)
    product_state[0] = 1.

    iter_state = circuit_func.circuit_2_state(my_circuit, product_state)

    fidelity_reached = np.abs(circuit_func.overlap_exact(target_state, iter_state))**2
    #print("fidelity reached : ", fidelity_reached)
    error_list = []
    error_list.append(1. - fidelity_reached)

    assert np.isclose(circuit_func.overlap_exact(target_state, target_state), 1.)
    for idx in range(0, N_iter):
        #################################
        #### variational optimzation ####
        #################################
        iter_state, my_circuit = circuit_func.var_circuit_exact(target_state, iter_state,
                                                         my_circuit, product_state, brickwall=True)
        #################
        #### Measure ####
        #################
        assert np.isclose(circuit_func.overlap_exact(iter_state, iter_state), 1.)

        fidelity_reached = np.abs(circuit_func.overlap_exact(target_state, iter_state))**2

        #print("fidelity reached : ", fidelity_reached)
        error_list.append(1. - fidelity_reached)

        ################
        ## Forcing to stop if already converge
        ################
        if (fidelity_reached > 1 - 1e-12):
            print("breaking!")
            break

    print("fidelity reached : ", fidelity_reached)

    iter_state = circuit_func.circuit_2_state(my_circuit, product_state)

    return iter_state, my_circuit, error_list


if __name__ == "__main__":
    
    create_path("HPC_data")
    create_path("HPC_data/ED")
    create_path("HPC_data/circuits")

    args = parse_args()

    L = args.L
    option = args.option
    filename = args.filename

    with open(filename) as f:
        line = f.readlines()[option-1]

    args = parse_args(line)
    delta = args.delta
    depth = args.depth
    N_iter = args.N_iter  # should be a multiple of 100

    save_filename = f"XXZ_GS_L{L}_delta{delta}".replace(".","-")

    print(f"Starting L={L}, delta={delta}, depth={depth}")

    N = L  # I should probably just go with L

    if os.path.exists("HPC_data/ED/"+save_filename+"_ED.pkl"):
        print("Loading ED data")
        V = load_obj("HPC_data/ED/"+save_filename+"_ED")
    else:
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

    circuit = initial_circuit(L, depth)

    error_list = []
    for ii in range(int(N_iter/100)):
        t1 = time.time()
        state, circuit, temp_error_list  = run_optimization(V.flatten(), circuit, L, depth, 100)  # save every 100
        error_list = error_list + temp_error_list[1:]
        t2 = time.time()
        print(f"100 iterations in {t2-t1} seconds")

        data = {"state": state, "circuit": circuit, "errors": error_list}
        save_obj(data, "HPC_data/circuits/"+save_filename+f"_depth{depth}_N_iter{N_iter}_circuit")