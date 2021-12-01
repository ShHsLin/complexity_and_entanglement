import numpy as np
import pickle
import os, sys
sys.path.append('..')

from quantumCircuit.methods.states import *
from quantumCircuit.methods.gates import *

from tensor_network_functions.mps_func import MPS_2_state, state_2_MPS, lpr_2_plr, plr_2_lpr
from tensor_network_functions.mps_func import init_mps, right_canonicalize, MPS_compression_variational

from scipy.sparse.linalg import eigsh, expm

import time

from matplotlib import pyplot as plt

def MPS_compressed_state(psi_mps_exact, N, chi):
    psi_chi = init_mps(N,chi,2)
    psi_chi = lpr_2_plr(psi_chi)  # convert to plr for right canonicalization
    psi_chi, _ = right_canonicalize(psi_chi, chi=chi)  # need right canonical form for trial function in MPS compression
    psi_chi = plr_2_lpr(psi_chi)  # convert to lpr for MPS compression
    trunc_error = MPS_compression_variational(psi_chi,psi_mps_exact)
    print(trunc_error)
    psi_chi = lpr_2_plr(psi_chi)  # convert to plr for MPS_2_state
    psi_chi = MPS_2_state(psi_chi)

    print(f"fidelity = {np.abs(np.vdot(psi.toVector().flatten(),psi_chi.flatten()))**2}")

    psi_chi = State.vector(psi_chi)

    return psi_chi


N = 8

H = sp.csr_matrix((2**N, 2**N))  # using sparse is far faster!!!
for ii in range(N-1):
    H = H + Gate.xx([ii,ii+1]).toSparseArray(N)
for ii in range(N):
    H = H + 0.8*Gate.z(ii).toSparseArray(N) + 0.1*Gate.x(ii).toSparseArray(N)
#H = H + 0.0001*Gate.x(0).toSparseArray(N)
#H = H.todense()

print("Hamiltonian constructed")


E,V = eigsh(H, int(0.125*2**N), which='SA')
V = V[:,-1].reshape((-1,1)) * (1 + 0*1j)  # make vector complex for MPS!
psi = State.vector(V)

print("Eigenproblem solved")


psi_mps_exact = state_2_MPS(V, N, 2**(N//2))
psi_mps_exact = plr_2_lpr(psi_mps_exact)  # convert to lpr for MPS compression


entropy_exact = [psi.bipartiteEE([x for x in range(y)]) for y in range(N+1) ]
entropy_list = [entropy_exact]
chi_list = [16,8,4,2]
for chi in chi_list:
    psi_chi = MPS_compressed_state(psi_mps_exact, N, chi)
    entropy_list.append([psi_chi.bipartiteEE([x for x in range(y)]) for y in range(N+1) ])

for entropy in entropy_list:
    plt.plot(entropy)

entropy_list = []
chi_list = [16,8,4,2]
for chi in chi_list:
    psi_chi = state_2_MPS(V, N, chi)
    psi_chi = MPS_2_state(psi_chi)
    psi_chi = State.vector(psi_chi)
    entropy_list.append([psi_chi.bipartiteEE([x for x in range(y)]) for y in range(N+1) ])


for entropy in entropy_list:
    plt.plot(entropy, '--')
plt.legend(['exact'] + ['chi='+str(chi) for chi in chi_list]  + ['chi='+str(chi) for chi in chi_list])
plt.show()
