# Selection + drift + sampling, and how a power value is obtained
# Diploid, 3 sampling points. Same logic as the simulator, stripped to the essentials.
#
# Each generation: (1) selection pushes the allele frequency deterministically,
#                  (2) drift redraws 2N gene copies from that pushed frequency.
# At the sampling generations, specimens are drawn from the population of that replay.
# Score = observed change in sample frequency between the first and last points.
# Threshold = score exceeded by only 5% of the no-selection replays.
# Power = share of selection replays whose score is above the threshold.

set.seed(1)

# ---- settings ----
p0          <- 0.30          # allele frequency at the start (generation 0)
N           <- 10000         # breeding individuals per generation (diploid: 2N gene copies)
s           <- 0.01          # selection coefficient
h           <- 0.5           # dominance (0.5 = additive)
n_gen       <- 100           # number of generations
sample_gens <- c(0, 50, 100) # generations at which specimens are collected
n_ind       <- c(20, 20, 40) # diploid individuals sampled at each point
alpha       <- 0.05          # significance level
n_neutral   <- 2000          # replays without selection (set the threshold)
n_selection <- 1000          # replays with selection (estimate power)

# ---- one replay: population history, then sampling; returns the sample frequencies ----
one_replay <- function(s_used) {
  p <- p0
  sample_freq <- numeric(length(sample_gens))

  for (g in 0:n_gen) {

    # sampling, if this generation is a sampling point
    k <- which(sample_gens == g)
    if (length(k) == 1) {
      n_copies <- 2 * n_ind[k]
      sample_freq[k] <- rbinom(1, size = n_copies, prob = p) / n_copies
    }

    if (g == n_gen) break

    # (1) selection: deterministic push
    q <- 1 - p
    AA <- p^2
    Aa <- 2 * p * q
    aa <- q^2
    w_mean <- AA * (1 + s_used) + Aa * (1 + h * s_used) + aa
    p_star <- (AA * (1 + s_used) + 0.5 * Aa * (1 + h * s_used)) / w_mean

    # (2) drift: the next generation's 2N gene copies are drawn from p_star
    p <- rbinom(1, size = 2 * N, prob = p_star) / (2 * N)
  }
  sample_freq
}

# ---- step 1: no-selection replays ----
score_neutral <- numeric(n_neutral)
for (r in 1:n_neutral) {
  f <- one_replay(0)
  score_neutral[r] <- f[length(f)] - f[1]   # observed change, last minus first
}

# ---- step 2: threshold = value exceeded by only alpha of the neutral scores ----
sorted <- sort(score_neutral)
threshold <- sorted[ceiling((1 - alpha) * n_neutral)]

# ---- step 3: selection replays ----
score_selection <- numeric(n_selection)
for (r in 1:n_selection) {
  f <- one_replay(s)
  score_selection[r] <- f[length(f)] - f[1]
}

# ---- step 4: count ----
power <- mean(score_selection > threshold)
false_alarm <- mean(score_neutral > threshold)

cat("Threshold (observed change):", round(threshold, 3), "\n")
cat("Neutral replays above it:   ", round(100 * false_alarm, 1), "%  (should be about", 100 * alpha, "%)\n")
cat("Power:                      ", round(100 * power, 1), "%\n")
