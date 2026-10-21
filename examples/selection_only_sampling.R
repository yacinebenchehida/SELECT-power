# Selection only (no drift): how the sampled data are generated
# Diploid, 3 sampling points

set.seed(1)

# ---- settings ----
p0          <- 0.01          # allele frequency at the start (generation 0)
s           <- 0.05          # selection coefficient
h           <- 0.5           # dominance (0.5 = additive)
n_gen       <- 90            # number of generations (e.g. 1935 -> 2025, one per year)
sample_gens <- c(0, 45, 90)  # generations at which specimens are collected
n_ind       <- c(20, 20, 40) # diploid individuals sampled at each point
n_rep       <- 5             # number of replays

# ---- STEP 1: the population (deterministic, computed ONCE) ----
p <- numeric(n_gen + 1)      # p[g + 1] = allele frequency in generation g
p[1] <- p0

for (g in 1:n_gen) {
  q <- 1 - p[g]

  # genotype frequencies before selection (Hardy-Weinberg)
  AA <- p[g]^2
  Aa <- 2 * p[g] * q
  aa <- q^2

  # fitness of each genotype
  wAA <- 1 + s
  wAa <- 1 + h * s
  waa <- 1
  w_mean <- AA * wAA + Aa * wAa + aa * waa

  # allele frequency after selection = frequency in the next generation
  p[g + 1] <- (AA * wAA + 0.5 * Aa * wAa) / w_mean
}

true_p <- p[sample_gens + 1]   # true population frequency at the 3 sampling points
cat("True population frequency at the 3 points:\n")
print(round(true_p, 4))

# ---- STEP 2: the sampling (random, repeated n_rep times) ----
sample_counts <- matrix(NA, nrow = n_rep, ncol = length(sample_gens))
sample_freq   <- matrix(NA, nrow = n_rep, ncol = length(sample_gens))

for (r in 1:n_rep) {
  for (k in 1:length(sample_gens)) {
    n_copies <- 2 * n_ind[k]   # diploid: 2 gene copies per individual

    # each copy carries the allele with probability true_p[k]
    sample_counts[r, k] <- rbinom(1, size = n_copies, prob = true_p[k])
    sample_freq[r, k]   <- sample_counts[r, k] / n_copies
  }
}

cat("\nAllele copies found in the sample (rows = replays, columns = points):\n")
print(sample_counts)
cat("\nSample allele frequency:\n")
print(round(sample_freq, 3))
