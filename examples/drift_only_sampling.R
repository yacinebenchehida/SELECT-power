# Drift only (no selection): how the sampled data are generated
# Diploid, 3 sampling points

set.seed(1)

# ---- settings ----
p0          <- 0.10          # allele frequency at the start (generation 0)
N           <- 500           # breeding individuals per generation (diploid: 2N gene copies)
n_gen       <- 90            # number of generations
sample_gens <- c(0, 45, 90)  # generations at which specimens are collected
n_ind       <- c(20, 20, 40) # diploid individuals sampled at each point
n_rep       <- 5             # number of replays

pop_p_at_samples <- matrix(NA, nrow = n_rep, ncol = length(sample_gens))
sample_counts    <- matrix(NA, nrow = n_rep, ncol = length(sample_gens))

for (r in 1:n_rep) {

  # ---- STEP 1: the population (random, re-simulated in EVERY replay) ----
  p <- numeric(n_gen + 1)    # p[g + 1] = population allele frequency in generation g
  p[1] <- p0

  for (g in 1:n_gen) {
    # next generation: draw 2N gene copies from the current generation's frequency
    copies_with_allele <- rbinom(1, size = 2 * N, prob = p[g])
    p[g + 1] <- copies_with_allele / (2 * N)
  }

  # ---- STEP 2: the sampling (from this replay's own population) ----
  for (k in 1:length(sample_gens)) {
    p_now    <- p[sample_gens[k] + 1]   # population frequency at this sampling point
    n_copies <- 2 * n_ind[k]
    pop_p_at_samples[r, k] <- p_now
    sample_counts[r, k]    <- rbinom(1, size = n_copies, prob = p_now)
  }
}

cat("TRUE population frequency at the 3 points (differs between replays = drift):\n")
print(round(pop_p_at_samples, 3))
cat("\nAllele copies found in the sample (sampling noise on top of drift):\n")
print(sample_counts)
