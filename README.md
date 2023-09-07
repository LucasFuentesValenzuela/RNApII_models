# RNApII_models

This package accompanies the paper `RNA polymerase II dynamics and mRNA stability feedback determine mRNA scaling with cell size`. It contains the code enabling two different parts of the analysis in the paper.

For both parts, one needs to create a `config.toml` file and put it in the root repository of this package. The `config.toml` file should contain the path to a `DATADIR` (containing the data measured in the paper) and a `RESULTSDIR` (containing simulation results run from this package).

## Dynamic equilibrium model
The fitting of the dynamic equilibrium model and its comparison with ChIP-seq data is performed in the `Python` notebook under `DE_model/`.

## TASEP models for RNApII occupancy
The rest of the repository consists in stochastic models of RNApII promoter-binding and elongation. The source code is implemented in `Julia` and is located under `src/`. To run simulations appearing in the paper, `cd` into the root repository of this package. Then add the package and activate the associated environment:
```
]
add .
activate .
```
The package should now be installed and activated.

To run one simulation, run
`julia feasible_analysis.jl --type screen --ntimes 10 --Omega 1`
The meaning of the arguments are the following:
- `type`: dictates the range of values for the `kon` and the initiation rate `alpha` to be considered when running simulations. `screen` is for a large parameter screen to identify feasible points. `narrow` is for exploring scaling in the neighborhood of a feasible point. `wide` is for exploring scaling over a large range of expression levels around a feasible point.
- `ntimes`: the number of times each simulation is run.
- `Omega`: the average residence time of RNApIIs on the promoter, in seconds. 

### Notebooks
A few notebooks accompany the implementation.
- `PromoterBinding_feasibility_dataAnalysis.jl`: this is the main notebook of the paper, comparing measured data with extensive simulations. Note that simulations must have been run and saved in `config["RESULTSDIR"]` for this notebook to work.
- `Theoretical_models.jl` is an exploratory notebook using theoretical models of the system. It was used for debugging and building intuition.
- `Feasibility_theory.jl` is a notebook that implements a few of the results from stochastic simulations on the basis of theoretical models.
