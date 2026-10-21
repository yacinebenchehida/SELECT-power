# Possible improvements

Ideas that have been discussed but are not implemented yet. SPOT is meant to stay simple, fast and general; for anything beyond its scope, a flexible forward simulator such as [SLiM](https://messerlab.org/slim/) or simuPOP is the better tool. Items here will be considered depending on user feedback.

## Model

- **Viruses, bacteria and other organisms with very short generations.** This would need sampling dates finer than a year (days or hours), an asexual haploid model with very large population sizes, and a way to handle tens of thousands of generations quickly. The usual approach is rescaling: simulate fewer generations with proportionally stronger selection and a smaller population, keeping Ne·s and the drift per unit time unchanged, as is done in SLiM.
- **More flexible demographic histories.** Ne can currently be constant, change exponentially between the first and last sampling points, or be set separately for each interval between sampling points. Possible additions: change points that do not coincide with sampling dates, gradual (rather than abrupt) recovery after a bottleneck, and random fluctuations in population size.
- **Overlapping generations and age structure**, for long-lived species where individuals of different ages breed together.
- **Selection that changes over time**, for example selection that stops, reverses or fluctuates.
- **Population structure**: several connected populations instead of one population with an optional migrant source.
- **Linked selection** (hitchhiking, background selection) and more than one locus.
- **Sampling a period rather than a date**: spreading the individuals of a time point over a window of years instead of pooling them at one representative year.

## Statistics

- Other test statistics (for example likelihood-based tests of selection from time series) next to the trend test.
- A false discovery rate correction as an alternative to Bonferroni.

## Software

- **Command-line and HPC use**: a small `run.js` wrapper that reads parameters from a JSON file and writes results, so that thousands of settings can be run as an array job on a cluster (for example with SLURM).
- **Python or R package**: a port of the simulation core, if wider use or new features require a full rewrite.
- Saving and loading parameter sets from the web page, and exporting results as a table.
