import numpy as np
import pickle
import os, sys
sys.path.append('..')

from quantumCircuit.methods.aux_methods import *
from quantumCircuit.methods.states import *
from quantumCircuit.methods.gates import *

from matplotlib import pyplot as plt

L = 12
g = 1.2
h = 0.1
depth = 4
N_iter = 2000

save_filename = f"Ising_GS_L{L}_g{g}_h{h}".replace(".","-")
V = load_obj("HPC_data/ED/"+save_filename+"_ED")

psi = State.vector(V)

entropy_exact = [psi.bipartiteEE([x for x in range(y)]) for y in range(L+1) ]
entropy_list = [entropy_exact]

data = load_obj("HPC_data/circuits/"+save_filename+f"_depth{depth}_N_iter{N_iter}_circuit")

plt.plot(data["errors"])
plt.show()

psi_circuit = State.vector(data["state"].reshape((-1,1)))
entropy_list.append([psi_circuit.bipartiteEE([x for x in range(y)]) for y in range(L+1) ])

for entropy in entropy_list:
    plt.plot(entropy)

plt.show()