# complexity_and_entanglement

Our goal is to study the complexity and entanglement properties of various quantum states.


### Cloning the Project
The project involves submodule and should be cloned by
```
git clone --recursive git@github.com:ShHsLin/complexity_and_entanglement.git
```
After cloning, it might be good to consider attach back to the HEAD of the submodule, i.e.
```
cd complexity_and_entanglement
cd auto_isoTNS; git checkout master; cd ..
cd tensor_network_functions; git checkout main; cd ..
```


It has dependencies on
- [tensor_network_functions](https://github.com/ShHsLin/tensor_network_functions)
- [auto_isoTNS](https://github.com/ShHsLin/auto_isoTNS/)


For more instruction, see [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules)


It is recommended to always update to the newest version before working. This can be done by
```
git fetch  # the standard part
git pull
git submodule foreach git fetch  # the part for submodule
git submodule foreach git pull
```


### Working
- circuit directory contains script which can turn a given exact state into QC.
- tensor_network_functions contains script that can turn a given exact state into MPS.


### TODO
- Add the exact diagonalization code here
- Setup the pipeline (ED -> QC, MPS -> savestate)
- Setup the measurement.py


