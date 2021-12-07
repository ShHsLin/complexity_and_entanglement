import numpy as np
import pickle
import os, sys
sys.path.append('..')

from quantumCircuit.methods.aux_methods import *
from quantumCircuit.methods.states import *
from quantumCircuit.methods.gates import *

from matplotlib import pyplot as plt

L = 20
g = 1.4
h = 0.9045
depth = 8
N_iter = 5000

save_filename = f"Ising_GS_L{L}_g{g}_h{h}".replace(".","-")
V = load_obj("HPC_data/ED/"+save_filename+"_ED")

psi = State.vector(V)

entropy_exact = [psi.bipartiteEE([x for x in range(y)]) for y in range(L+1) ]
entropy_list = [entropy_exact]

data = load_obj("HPC_data/circuits/"+save_filename+f"_depth{depth}_N_iter{N_iter}_circuit")

#plt.plot(np.log(data["errors"]))
#plt.show()


psi_circuit = State.vector(data["state"].reshape((-1,1)))


print(psi.geometricEntanglement_svd())
print(psi_circuit.geometricEntanglement_svd())
#print(psi_circuit.geometricEntanglement())

entropy_list.append([psi_circuit.bipartiteEE([x for x in range(y)]) for y in range(L+1) ])

for entropy in entropy_list:
    plt.plot(entropy)

plt.show()




for kk in range(1,10):
    qubits = [(kk-1)/2 + x for x in range(L-kk)]
    MI_exact = np.array([psi.mutualInformation([y],[y+kk]) for y in range(L-kk)])
    #plt.plot(MI)
    MI = np.array([psi_circuit.mutualInformation([y],[y+kk]) for y in range(L-kk)])
    plt.plot(qubits,MI/MI_exact)

plt.legend([str(x) for x in range(9)])
plt.show()