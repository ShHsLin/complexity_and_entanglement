import numpy as np
import pickle
import os, sys
sys.path.append('..')

from quantumCircuit.methods.states import *
from quantumCircuit.methods.gates import *

from tensor_network_functions.mps_func import MPS_2_state, state_2_MPS, lpr_2_plr, plr_2_lpr
from tensor_network_functions.mps_func import init_mps, right_canonicalize, MPS_compression_variational

import tensor_network_functions.mps_func as mps_func
import tensor_network_functions.circuit_func as circuit_func

from scipy.sparse.linalg import eigsh, expm, expm_multiply

import time

from matplotlib import pyplot as plt


def run_optimization(target_state, L, depth, N_iter):
    #np.random.seed(1)
    np.set_printoptions(linewidth=2000, precision=5, threshold=4000)

    """
    args = parse_args.parse_args()

    L = args.L
    J = 1.
    depth = args.depth
    N_iter = args.N_iter
    T = args.T  # the target state is corresponding to time T.
    """

    save_each = 100
    tol = 1e-12
    cov_crit = tol * 0.1
    max_N_iter = N_iter

    ############### PARAMETERS INITIALIZATION #################
    ## We should test identity initialization and
    ## trotterization initialization

    idx = 0
    my_circuit = []
    t_list = [0]
    error_list = []
    Sz_array = np.zeros([N_iter, L], dtype=np.complex)
    ent_array = np.zeros([N_iter, L-1], dtype=np.double)
    num_iter_array = np.zeros([N_iter], dtype=np.int)

    ################# CIRCUIT INITIALIZATION  ######################
    # product_state = [np.array([1., 0.]).reshape([2, 1, 1]) for i in range(L)]
    product_state = np.zeros([2**L], dtype=np.complex128)
    product_state[0] = 1.
    for dep_idx in range(depth):
        # identity_layer = [np.eye(4, dtype=np.complex).reshape([2, 2, 2, 2]) for i in range(L-1)]
        # my_circuit.append(identity_layer)

        random_layer = []
        for idx in range(L-1):
            if (idx + dep_idx) % 2 == 0:
                #random_layer.append(np.eye(4).reshape([2,2,2,2]).reshape([4, 4]))
                random_layer.append(circuit_func.random_2site_U(2).reshape([4, 4]))
            else:
                random_layer.append(np.eye(4).reshape([2,2,2,2]).reshape([4, 4]))

        my_circuit.append(random_layer)
        current_depth = dep_idx + 1

    iter_state = circuit_func.circuit_2_state(my_circuit, product_state)
    '''
    Sz_array[0, :] = mps_func.expectation_values_1_site(mps_of_layer[-1], Sz_list)
    ent_array[0, :] = mps_func.get_entanglement(mps_of_last_layer)
    '''
    fidelity_reached = np.abs(circuit_func.overlap_exact(target_state, iter_state))**2
    #print("fidelity reached : ", fidelity_reached)
    error_list.append(1. - fidelity_reached)

    stop_crit = 1e-1
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
        num_iter_array[idx] = idx
        t_list.append(idx)

        ################
        ## Forcing to stop if already converge
        ################
        if (fidelity_reached > 1 - 1e-12 or np.abs((error_list[-1] - error_list[-2])/error_list[-1]) < 1e-4) and idx > save_each:
            break

    print("fidelity reached : ", fidelity_reached)

    num_data = len(error_list)
    num_iter_array = num_iter_array[:num_data]

    iter_state = circuit_func.circuit_2_state(my_circuit, product_state)

    return iter_state


N = 16

H = sp.csr_matrix((2**N, 2**N))  # using sparse is far faster!!!
for ii in range(N-1):
    H = H + Gate.xx([ii,ii+1]).toSparseArray(N)
for ii in range(N):
    H = H + 1.4*Gate.z(ii).toSparseArray(N) + 0.9045*Gate.x(ii).toSparseArray(N)
#H = H + 0.0001*Gate.x(0).toSparseArray(N)
#H = H.todense()

print("Hamiltonian constructed")

# for ground state
#E,V = eigsh(H, 1, which='SA')
#V = V.reshape((-1,1)) * (1 + 0*1j)  # make vector complex for MPS!

# for excited state
#E,V = eigsh(H, int(0.125*2**N), which='SA')
#V = V[:,-1].reshape((-1,1)) * (1 + 0*1j)  # make vector complex for MPS!

# for time evolution
V = np.zeros((2**N,1), dtype=complex)
V[0,0] = 1.
t = 4
V = expm_multiply(-1j*t*H, V)

psi = State.vector(V)

print("Eigenproblem solved")


entropy_exact = [psi.bipartiteEE([x for x in range(y)]) for y in range(N+1) ]
entropy_list = [entropy_exact]

depth_list = [1,2,3,4,5,6]
for depth in depth_list:
    psi_circuit = run_optimization(V.flatten(), N, depth, 5000)
    psi_circuit = State.vector(psi_circuit.reshape((-1,1)))
    entropy_list.append([psi_circuit.bipartiteEE([x for x in range(y)]) for y in range(N+1) ])

for entropy in entropy_list:
    plt.plot(entropy)

plt.show()