# CRoCoDiL 🐊

Welcome to the *$\,$**C**onvection-**R**eaction **o**f **C**arb**o**n **D**iodixe **i**n **L**UCiFEx $\,$* package!

CRoCoDiL solves partial differential equations modelling flow in a porous medium coupled to the solutal and/or thermal transport of dissolved CO<sub>2</sub> and porosity evolution due to the dissolution of capillary-trapped CO<sub>2</sub>. It is based on the core package [LUCiFEx](https://github.com/george-poole/LUCiFEx), which itself is based on [FEniCSx](https://github.com/FEniCS/dolfinx). For more details on the mathematical formulation, see `doc`. User-defined choices of domain, constitutive relations, boundary conditions and initial conditions may be prescribed to investigate both novel and classic systems.

## Installation (macOS)

See `https://github.com/george-poole/LUCiFEx` to first install the `lucifex` package.

`git clone https://github.com/george-poole/CRoCoDiL.git`

`pip install .` (or `pip install -e .` for editable mode).

## Simulations

### Specifying a simulation

Examples of classic systems (e.g. Rayleigh-Bénard, Rayleigh-Taylor) specified in this framework can be found in the module `crocodil.dns`. 

```python
# file `dns_specific.py` 
from lucifex.sim import configure_simulation, run_from_cli
from crocodil.dns import dns_generic

@configure_simulation(
    # arguments setting the simulation's default I/O and compiler options
)
def dns_specific(
    # arguments for the domain, boundary conditions, etc.
):
    # code creating objects for domain, boundary conditions, etc.
    return dns_generic(
        # arguments specifying domain, boundary conditions, etc.
    )

if __name__ == "__main__":
    run_from_cli(dns_specific)
```

### Running a simulation

#### In an iPython notebook

For small-scale testing and prototyping, it can be beneficial to run the simulation in an interactive environment and have the simulation data available in memory for immediate visualization or further computations.

```python
# file `dns_specific.ipynb` in same directory as `dns_specific.py`
from lucifex.sim import integrate
from dns_specific import dns_specific

simulation = dns_specific(
    # arguments overriding the simulation's default I/O and low-level behaviour
)(
    # arguments for the domain, boundary conditions, etc.
)
run(
    simulation,
    # arguments for the time loop
)
```

In particular, the `store_delta` argument determines the interval size, in either number of integration steps if of type `int` or simulation time if of type `float`, at which to keep simulation data in memory. If `None`, no simulation data is kept in memory. If unspecified in the command line, `store_delta` will assume its value configured in the `configure_simulation` decorator function.

#### From the command line

For long-running simulations and batches of simulations, it is best to run from the command line.

`python dns_specific.py --help` will list all argument names and default values. 

The complete set of command-line arguments consists of those passed to `configure_simulation` to override the simulation's default I/O and compiler options, those passed to `dns_specific` which are required for the domain, constitutive relations, boundary conditions and initial conditions, and those passed to `run` for the time loop. In particular, the `write_delta` argument is much like the `store_delta` except that it determined the interval size at which to write simulation data to file.

To run a simulation with the argument `X` set to $X=X_0$ and `Y` set to $Y=Y_0$
<br>
`python dns_specific.py --X X0 --Y Y0`<br>

To run multiple simulations in parallel, the GNU command line utility `parallel` is recommended.

To run `N_PROC` simultaneous processes with the arguments `X` and `Y` taking values $(X, Y)\in\{(X_0, Y_0), (X_0, Y_1), (X_1, Y_0), (X_1, Y_1)\}$
<br>
`parallel -j N_PROC "python dns_specific.py --X {1} --Y {2}" ::: X0 X1 ::: Y0 Y1` <br>

To run `N_PROC` simultaneous processes with the arguments `X` and `Y` taking values $(X, Y)\in\{(X_0, Y_0), (X_1, Y_1), (X_2, Y_2)\}$
<br>
`parallel -j N_PROC --link "python dns_specific.py --X {1} --Y {2}" ::: X0 X1 X2 ::: Y0 Y1 Y2` <br>

Further command line utilities:
* `caffeinate` e.g. `caffeinate -d -i -s -t <SECONDS> <COMMAND>` to prevent sleeping
* `nohup` e.g. `nohup <COMMAND> & disown` to run without interruption
* `htop` and `kill` for process mangement
* (Ctrl + Z) followed by `bg` to move process to background <br>

### Postprocessing

In addition to the core functionality provided by the `lucifex.io` and `lucifex.viz` modules, CRoCoDiL contains the `crocodil.post` module for postprocessing utilties such as computing contours and plume statistics calculations.