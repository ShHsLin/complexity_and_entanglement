
## Goal:
Given an exact state, we learn the circuit representation of such state.


### Requirement
numpy, jax

### Current support circuit type:
brickwall and staircase

### Current support optimization methods
Exact contraction + Polar decomposition
```
python run_optimization_exact.py --N_iter 1000 --depth 4
```

Exact contraction + Riemannian gradient descent
```
python run_optimization_exact_jax.py --N_iter 1000 --depth 4
```

MPS contraction + Polar decomposition
```          
python run_optimization_mps.py --N_iter 1000 --depth 4
```

* MPS and Exact contraction + Polar decomposition should give same result.
Because currently, MPS contraction is set to have large truncation bond dimension.
MPS contraction could deal with larger system size with low entanglement.

* Riemannian gradient descent gives slightly worse result for simple circuit.
Have not yet tested fully for all cases.




### Detail of possible flags:
You can use the --help to see the detail or as following.
```
Quantum Circuit Simulation

optional arguments:
  -h, --help            show this help message and exit
  --L L                 system size. Default: 10
  --depth DEPTH         depth of the circuitDefault: 2
  --N_iter N_ITER       (maximum) number of iteration in the
                        optimizationDefault: 1
  --brickwall BRICKWALL
                        Whether or not using brickwalltype in 1 for true, 0
                        for falseDefault: False
```


## TODO:

###  adding MPS contraction + Riemannian gradient descent
* add mps function to circuit\_func\_jax.py
* overcome the svd problem; truncation no problem? stability?
* fix the convention to U\_ij instead of U\_ji
* after fixing the convention in circuit\_func\_jax.py;
we can then change the circuit\_func.py and then whole project to
use the U\_ij notation.

* Having stable MPS + RGD, we can then apply this on PEPS.

  

