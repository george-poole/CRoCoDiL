# CRoCoDiL 🐊

Welcome to the *&nbsp;**C**onvection-**R**eaction **o**f **C**arb**o**n **D**iodixe **i**n **L**UCiFEx&nbsp;* package!

CRoCoDiL provides a framework to model the geological sequestration of CO<sub>2</sub> in a porous medium. It is based on the core package [LUCiFEx](https://github.com/george-poole/LUCiFEx), which itself is based on [FEniCSx](https://github.com/FEniCS/dolfinx). To get started with LUCiFEx, refer to the [user guide](https://george-poole.github.io/LUCiFEx/notebooks/P01_user_guide/U00_introduction.html).

## What does CRoCoDiL do?

+ solves the equations of thermosolutal convection coupled to Darcy flow and porosity evolution due to the dissolution of capillary-trapped CO<sub>2</sub>
+ lets you customise the domain, initial conditions, boundary conditions and constitutive relations
+ enables easy visualization and postprocessing in iPython notebooks

## Installation (macOS)

Please note that CRoCoDiL is a research code still under active development.

Refer to `README.md` in the [LUCiFEx repository](https://github.com/george-poole/LUCiFEx) to first install the `lucifex` package. Then to install the `crocodil` package

`git clone https://github.com/george-poole/CRoCoDiL.git`

`pip install .` (or `pip install -e .` for editable mode)