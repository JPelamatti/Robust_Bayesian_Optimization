# Robust Bayesian optimization
The following code was developed during a post-doc in collaboration with Ecole des Mines de Saint-Etienne and Ecole Centrale de Lyon as part of the OQUAIDO chair by Julien Pelamatti, Rodolphe Le Riche, Christophette Blanchet-Scalliet and Céline Helbert. The theory behind the implemented algorithms is detailed in the paper 'Coupling and selecting constraints in Bayesian optimization under uncertainties'.

Part of the code presented in this repository is inspired from the code developed by Reda El Amri et al. for the paper ''A Sampling Criterion for Constrained BayesianOptimization with Uncertainties''

Files :
  - `Auxiliary_functions.R` : contains the auxiliary functions required for the robust Bayesian optimization algorithm
  - `main_MMCS.R` : allows to perform the robust Bayesian optimization relying on Multioutput Modeling and selection of the refined constraint
  - `main_MMCU.R` : allows to perform the robust Bayesian optimization relying on Multioutput Modeling and simultaneous and identical refinement of all problem functions
  - `main_REF.R` : allows to perform the robust Bayesian optimization relying on independent modeling of all outputs and simultaneous and identical refinement of all problem functions
  - `main_SMCS.R` : allows to perform the robust Bayesian optimization relying on independent modeling of all outputs and selection of the refined constraint
  - `simulate_new.R` : contains the function allowing to simulate Gaussian process trajectories
