# epiclomal_real pipeline

epiclomal_real.yaml runs the entire real pipeline:
- runs the non-probabilistic methods
- runs EpiclomalBasic a given number of times. Each run corresponds to one Variational Bayes (VB) iteration.
- runs EpiclomalRegion a given number of times. Each run corresponds to one Variational Bayes (VB) iteration.
- evaluates the EpiclomalBasic results
- evaluates the EpiclomalRegion results

epiclomal_real.yaml runs only the non-probabilistic methods.
