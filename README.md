# DLMP under Uncertainty Code Supplement
Repository containing supplementary data and code for the paper "*Distribution Electricity Prices under Uncertainty*" by Robert Mieth and Yury Dvorkin.

**Instructions:**

The optimization model was implemented by using [JuMP](https://github.com/JuliaOpt/JuMP.jl) and auxiliary packages in the [Julia](http://julialang.org/downloads/) programming language.
Additionally, we used Mosek 8.1.0.63 in our numerical experiments. [Mosek](https://www.mosek.com) is a commercial solver which must be installed and licensed. The solver was chosen for its specific features for second-order conic programming. For more information on solvers, see the JuMP documentation.

The experiments require Julia 1.1. All necessary packages and their respective versions are distributed in the `Project.toml` file. The packages can be loaded in a new Julia environment directly via the Julia package manager IDE or running: 
```
julia> import Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

**Running the code:**

All necessary functions are provided by the files `src/model_definition.jl`, `src/input.jl`, `src/output.jl` and `src/tools.jl`. All necessary data is provided in the `data` folder.
An exemplary case study that showcases the model implementation is provided as a Jupyter Notebook `uncertain_dlmp_case_study.ipynb`.
To run the notebook install the `IJulia` package:
```
julia> Pkg.add("IJulia")
```
and then in your terminal:
```
jupyter notebook uncertain_dlmp_case_study.ipynb
```
More explanation on the model is given in the notebook. 
