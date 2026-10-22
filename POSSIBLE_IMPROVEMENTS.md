# Possible improvements

Ideas that have been discussed but are not implemented yet. SPOT is meant to stay simple, fast and general; for anything beyond its scope, a flexible forward simulator such as [SLiM](https://messerlab.org/slim/) or simuPOP is the better tool. Items here will be considered depending on user feedback.

## Model

- **Historical DNA of bacteria, viruses and other organisms with very short generations (planned next).** A few centuries of pathogen samples can mean 10⁵ to 10⁶ generations, which takes from half a minute to several minutes per power estimate if every generation is simulated. The planned route is rescaling, as in SLiM: simulate at most about 2,000 generations, with selection made proportionally stronger and the population proportionally smaller, so that Ne·s and the total drift over the period stay the same. The page would choose the scaling factor automatically, keep it small enough for the approximation to hold (scaled s below about 0.1), say when rescaling is used, and accept generation times down to about one hour (0.0001 years). Rescaled results need to be checked against full simulations (for example at 20,000 generations) before release.
- **More flexible demographic histories.** Ne can currently be constant, change exponentially between the first and last sampling points, or be set separately for each interval between sampling points. Possible additions: change points that do not coincide with sampling dates, gradual (rather than abrupt) recovery after a bottleneck, and random fluctuations in population size.
- **Overlapping generations and age structure**, for long-lived species where individuals of different ages breed together.
- **Selection that changes over time**, for example selection that stops, reverses or fluctuates.
- **Population structure**: several connected populations instead of one population with an optional migrant source.
- **Linked selection** (hitchhiking, background selection) and more than one locus.
- **Replicate populations (1 to 100)**: several independent populations starting from the same allele frequency and evolving in parallel (for example several regions or lineages sampled through time), with the evidence combined across replicates, for example by summing their scores. The drift-aware threshold works the same way, and run time grows in proportion to the number of replicates.
- **Sampling a period rather than a date**: spreading the individuals of a time point over a window of years instead of pooling them at one representative year.

## Statistics

- Other test statistics (for example likelihood-based tests of selection from time series) next to the trend test.
- A false discovery rate correction as an alternative to Bonferroni.

## Software

- **Command-line and HPC use**: a small `run.js` wrapper that reads parameters from a JSON file and writes results, so that thousands of settings can be run as an array job on a cluster (for example with SLURM).
- **Python or R package**: a port of the simulation core, if wider use or new features require a full rewrite.
- **Using all processor cores (Web Workers)**: run the replicates, the parameter power curve and the sample-size planner in parallel background threads of the browser, for a speed-up of roughly 4 to 8 times on a typical laptop without changing the model. The number of threads would be chosen automatically from the number of cores, and each block of replicates would get its own fixed random seed, so results stay identical whatever the number of threads.
- Saving and loading parameter sets from the web page, and exporting results as a table.
