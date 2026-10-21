# Statistical Programming
# Practical 2: Sequence duels and the perils of independence
#
# Group 3:
# Alejandra Campo Archbold (s3023700)
# Weihang Gu (s2932201)
# Weichen Zhao (s2929481)
#
# Contributions:
# - Alejandra Campo Archbold: 33% [designed the question structure, the
#   shared simulation function race and the plots]
# - Weihang Gu: []
# - Weichen Zhao: []
#
# GitHub URL: https://github.com/alejandrarchbold/Statistical_Programming.git


#-----------------------------------------------------------------------
# Overview
#-----------------------------------------------------------------------

# Two players choose different patterns of heads (1) and tails (0).
# Both observe the same coin-toss stream. The first pattern to
# appear wins. We simulate duels, estimate win probabilities and
# calculate exact waiting times. We then compare best replies,
# non-transitive cycles and the effects of coin bias.


#-----------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------


# Question 1: firsthit

firsthit <- function(x, a) {
  # Description: finds the first occurrence of a in x.
  # Inputs: x (coin tosses), a (the pattern to find).
  # Output: the ending position of the first match, or NA.
  # How it works: checks each window of length k in order.

  k <- length(a)
  n <- length(x)

  if (k == 0) {
    stop("Please enter a 0/1 pattern.")
  }

  if (n < k) {
    return(NA_integer_)
  }

  for (m in k:n) {
    w <- x[(m - k + 1):m]

    if (all(w == a)) {
      return(m)
    }
  }

  return(NA_integer_)
}


# Input check shared by race and ewait

check_p <- function(p) {
  # Description: stops unless p is a single probability with 0 < p < 1.
  if (length(p) != 1 || !is.numeric(p) ||
      !is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be between 0 and 1.")
  }
}


# Shared simulation function for Question 2, Question 3 and Question 4

race <- function(l_patterns, p = 0.5, nreps = 1, block = 20) {
  # Description: simulates coin streams until a pattern appears.
  # Inputs: l_patterns (list of patterns), p (chance of heads),
  #         nreps (number of streams), block (batch size).
  # Output: matrix containing winner indices and toss counts.
  # How it works: tosses are generated in blocks, and only unfinished
  #   streams receive a new block. Within a block the tosses are
  #   added one at a time: each stream keeps its last k tosses in
  #   a row of `recent`, which is compared with every pattern.

  # Input checks: equal-length, distinct 0/1 patterns and valid
  # p, nreps and block.
  if (!is.list(l_patterns) || length(l_patterns) == 0) {
    stop("the patterns must be a non-empty list.")
  }

  k <- length(l_patterns[[1]])

  if (k == 0) {
    stop("Patterns cannot be empty. Add 0/1 pattern.")
  }

  for (i in seq_along(l_patterns)) {
    a <- l_patterns[[i]]
    if (!is.numeric(a) || length(a) != k ||
        anyNA(a) || !all(a %in% c(0, 1))) {
      stop("Patterns must have equal lengths and contain 0/1.")
    }
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        if (all(a == l_patterns[[j]])) {
          stop("Patterns must be different.")
        }
      }
    }
  }

  check_p(p)

  if (length(nreps) != 1 || !is.numeric(nreps) ||
      !is.finite(nreps) || nreps < 1 ||
      nreps != floor(nreps)) {
    stop("nreps must be a positive integer.")
  }

  if (length(block) != 1 || !is.numeric(block) ||
      !is.finite(block) || block < 1 ||
      block != floor(block)) {
    stop("block must be a positive integer.")
  }

  recent <- matrix(0, nrow = nreps, ncol = k)
  winner <- rep(NA_integer_, nreps)
  tosses <- integer(nreps)
  active <- seq_len(nreps)

  while (length(active) > 0) {

    block_tosses <- matrix(
      rbinom(length(active) * block, size = 1, prob = p),
      nrow = length(active),
      ncol = block
    )

    for (j in seq_len(block)) {

      # streams that finished earlier in this block are skipped
      unfinished <- which(is.na(winner[active]))

      if (length(unfinished) == 0) {
        break
      }

      rows <- active[unfinished]

      # shift the window left by one and add the new toss at the end
      if (k > 1) {
        recent[rows, 1:(k - 1)] <- recent[rows, 2:k, drop = FALSE]
      }
      recent[rows, k] <- block_tosses[unfinished, j]
      tosses[rows] <- tosses[rows] + 1

      # a full window exists only after k tosses
      ready <- rows[tosses[rows] >= k]

      if (length(ready) == 0) {
        next
      }

      windows <- recent[ready, , drop = FALSE]

      # a window matches pattern q when all k positions are equal
      for (q in seq_along(l_patterns)) {
        target <- matrix(l_patterns[[q]], nrow = length(ready),
                         ncol = k, byrow = TRUE)
        matched <- rowSums(windows == target) == k
        winner[ready[matched]] <- q
      }
    }

    active <- which(is.na(winner))
  }

  return(cbind(winner = winner, tosses = tosses))
}


# Question 2: duel

duel <- function(a, b, p = 0.5) {
  # Description: simulates one duel between a and b.
  # Inputs: a and b (patterns), p (chance of heads).
  # Output: c(winner, tosses), with 1 for a and 0 for b.
  # How it works: calls race for one shared coin stream.

  result <- race(list(a, b), p, nreps = 1)

  return(c(
    winner = as.integer(result[1, "winner"] == 1),
    tosses = unname(result[1, "tosses"])
  ))
}


# Question 3: pwin

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Description: estimates the probability that a beats b.
  # Inputs: a and b (patterns), p (heads), nreps (duels).
  # Output: the proportion of games won by a.
  # How it works: counts the wins for a across simulated duels.

  results <- race(list(a, b), p, nreps = nreps)

  return(mean(results[, "winner"] == 1))
}


# Question 4: ewait

ewait <- function(a, p = 0.5) {
  # Description: calculates the exact mean waiting time for a.
  # Inputs: a (pattern), p (chance of heads).
  # Output: expected tosses until a first appears.
  # How it works: a has an overlap of length j if its first j
  #   outcomes equal its last j outcomes (j = k always counts).
  #   Each overlap adds 1 / P(a[1:j]) to the expected wait.

  k <- length(a)

  if (k == 0 || anyNA(a) || any(!(a %in% c(0, 1)))) {
    stop("Please enter a non-empty 0/1 pattern.")
  }

  check_p(p)

  wait <- 0

  for (j in seq_len(k)) {

    prefix <- a[1:j]
    suffix <- a[(k - j + 1):k]

    if (all(prefix == suffix)) {
      heads <- sum(prefix)
      tails <- j - heads
      prob <- p^heads * (1 - p)^tails
      wait <- wait + 1 / prob
    }
  }

  return(wait)
}


# Question 4: simwait

simwait <- function(a, p = 0.5, nreps = 20000) {
  # Description: simulates waiting times for a single pattern.
  # Inputs: a (pattern), p (heads), nreps (independent runs).
  # Output: the waiting time for each simulated stream.
  # How it works: calls race with one pattern instead of two.

  results <- race(list(a), p, nreps = nreps)

  return(results[, "tosses"])
}


# Question 5: all_patterns

all_patterns <- function(k) {
  # Description: generates all 2^k patterns of length k.
  # Input: k (positive integer).
  # Output: a matrix with one pattern per row. Row i is the number
  #   i - 1 written in binary, with the first toss as the most
  #   significant digit (1 = H, 0 = T). For k = 3 the order is
  #   TTT, TTH, THT, THH, HTT, HTH, HHT, HHH.
  # How it works: expand.grid varies its first column fastest, so
  #   the columns are reversed to put the first toss last to change.

  if (length(k) != 1 || !is.numeric(k) ||
      !is.finite(k) || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }

  patterns <- expand.grid(rep(list(0:1), k))
  patterns <- as.matrix(patterns)
  patterns <- patterns[, k:1, drop = FALSE]

  return(unname(patterns))
}


# Question 5: label

label <- function(a) {
  # Description: converts a binary pattern to H/T notation.
  # Input: a (0/1 pattern).
  # Output: a character label, e.g. c(1, 0, 1) -> "HTH".

  return(paste(c("T", "H")[a + 1], collapse = ""))
}


# Question 5: pmat

pmat <- function(k, p = 0.5, nreps = 10000) {
  # Description: estimates win probabilities for all pattern pairs.
  # Inputs: k (pattern length), p (heads), nreps (games per pair).
  # Output: matrix M, with M[i,j] = P(row i beats column j). Rows
  #   and columns follow the order of all_patterns(k) (for k = 3,
  #   1 = TTT, 2 = TTH, ..., 8 = HHH); the diagonal is NA.
  # How it works: simulates each pair once and fills both entries,
  #   since one of the two patterns always wins.

  patterns <- all_patterns(k)
  n <- nrow(patterns)
  labels <- apply(patterns, 1, label)

  M <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(labels, labels))

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      prob <- pwin(patterns[i, ], patterns[j, ], p = p, nreps = nreps)
      M[i, j] <- prob
      M[j, i] <- 1 - prob
    }
  }

  return(M)
}


# Question 5(a): best_reply

best_reply <- function(M) {
  # Description: finds the best reply to each pattern.
  # Input: M (win probability matrix).
  # Output: each pattern's best reply and its win probability.
  # How it works: column j of M holds the chances of beating
  #   pattern j, so the best reply is the row with the largest
  #   entry in column j.

  n <- ncol(M)

  result <- data.frame(
    pattern = colnames(M),
    reply = character(n),
    probability = numeric(n)
  )

  for (j in seq_len(n)) {
    best <- which.max(M[, j])
    result$reply[j] <- rownames(M)[best]
    result$probability[j] <- M[best, j]
  }

  result$probability <- round(result$probability, 3)

  return(result)
}


# Question 5(b): even_pairs

even_pairs <- function(M, tol = 0.03) {
  # Description: finds pairs with nearly equal win probabilities.
  # Inputs: M (win matrix), tol (allowed distance from 0.5).
  # Output: a table of pairs and their estimated probabilities.
  # How it works: checks the upper triangle of M, so each pair
  #   appears once.

  pos <- which(upper.tri(M) & abs(M - 0.5) < tol, arr.ind = TRUE)

  result <- data.frame(
    pattern_a = rownames(M)[pos[, 1]],
    pattern_b = colnames(M)[pos[, 2]],
    p_a_wins = round(M[pos], 3)
  )

  return(result)
}


# Question 5(b): find_cycle

find_cycle <- function(M, cutoff = 0.55) {
  # Description: finds a shortest clear non-transitive cycle.
  # Inputs: M (win matrix), cutoff (minimum win probability).
  # Output: pattern indices forming the cycle (the first pattern
  #   is repeated at the end), or NULL.
  # How it works: beats[i, j] is TRUE when i beats j by more than
  #   cutoff. From each start pattern, the search visits patterns
  #   in order of how many steps away they are, so the first step
  #   back to the start closes the shortest cycle through it.
  #   parent[v] records the pattern from which v was reached, so
  #   the path can be rebuilt. The shortest cycle over all start
  #   patterns is returned; its length is not fixed in advance.

  n <- nrow(M)

  beats <- M > cutoff
  beats[is.na(beats)] <- FALSE

  best <- NULL

  for (start in seq_len(n)) {

    visited <- rep(FALSE, n)
    visited[start] <- TRUE
    parent <- rep(NA_integer_, n)
    queue <- start
    cycle <- NULL

    while (length(queue) > 0 && is.null(cycle)) {

      u <- queue[1]
      queue <- queue[-1]

      for (v in which(beats[u, ])) {

        if (v == start) {
          path <- u
          while (path[1] != start) {
            path <- c(parent[path[1]], path)
          }
          cycle <- c(path, start)
          break
        }

        if (!visited[v]) {
          visited[v] <- TRUE
          parent[v] <- u
          queue <- c(queue, v)
        }
      }
    }

    if (!is.null(cycle) &&
        (is.null(best) || length(cycle) < length(best))) {
      best <- cycle
    }
  }

  return(best)
}


# Question 5(b): show_cycle

show_cycle <- function(M, cutoff = 0.55) {
  # Description: displays a non-transitive cycle.
  # Inputs: M (win matrix), cutoff (minimum win probability).
  # Output: the cycle and the probability of each step.
  # How it works: prints the results returned by find_cycle.

  cycle <- find_cycle(M, cutoff)

  if (is.null(cycle)) {
    cat("No clear non-transitive cycle found.\n")
    return(invisible(NULL))
  }

  labels <- rownames(M)

  cat("Cycle of length", length(cycle) - 1, ":",
      paste(labels[cycle], collapse = " -> "), "\n")

  for (i in seq_len(length(cycle) - 1)) {
    a <- cycle[i]
    b <- cycle[i + 1]
    cat(labels[a], "beats", labels[b], ":", round(M[a, b], 3), "\n")
  }
}


#-----------------------------------------------------------------------
# Tests and Simulations
#-----------------------------------------------------------------------

start_time <- proc.time()["elapsed"]

set.seed(2026)


# Question 1: Test firsthit on streams that can be checked by hand.
# x is H T T H T H H H T, the stream used as an example in the sheet.

x <- c(1, 0, 0, 1, 0, 1, 1, 1, 0)

cat("Q1: HTH (expected 6):", firsthit(x, c(1, 0, 1)), "\n")
cat("Q1: THH (expected 7):", firsthit(x, c(0, 1, 1)), "\n")
cat("Q1: HHH (expected 8):", firsthit(x, c(1, 1, 1)), "\n")
cat("Q1: TTT absent (expected NA):", firsthit(x, c(0, 0, 0)), "\n")
cat("Q1: stream shorter than pattern (expected NA):",
    firsthit(c(1, 0), c(1, 0, 1)), "\n")
cat("Q1: empty stream (expected NA):", firsthit(c(), c(1, 0, 1)), "\n")


# Question 2: Simulate one duel.

cat("\nQ2: One duel between HTH and THH\n")
print(duel(c(1, 0, 1), c(0, 1, 1)))


# Question 3: Estimate the win probability.
# HHH wins only if the first three tosses are HHH: after any tail,
# THH is completed before HHH. So the exact value is 1/8 = 0.125.

p_fair <- pwin(c(1, 1, 1), c(0, 1, 1))

cat("\nQ3: P(HHH beats THH) (exact 0.125):", round(p_fair, 4), "\n")


# Question 4: Compare exact and simulated waiting times.

cat("\nQ4: Length-2 patterns, fair coin (expected 6, 6, 4, 4):",
    ewait(c(1, 1)), ewait(c(0, 0)), ewait(c(1, 0)), ewait(c(0, 1)), "\n")

patterns4 <- list(c(1, 1, 1), c(1, 1, 0), c(1, 0, 1),
                  c(1, 0, 1), c(1, 1, 0))
probs4 <- c(0.5, 0.5, 0.5, 0.7, 0.3)

exact4 <- numeric(length(patterns4))
mc4 <- numeric(length(patterns4))
se4 <- numeric(length(patterns4))

for (i in seq_along(patterns4)) {
  exact4[i] <- ewait(patterns4[[i]], p = probs4[i])
  waits <- simwait(patterns4[[i]], p = probs4[i])
  mc4[i] <- mean(waits)
  se4[i] <- sd(waits) / sqrt(length(waits))
}

comparison4 <- data.frame(
  pattern = vapply(patterns4, label, character(1)),
  p = probs4,
  exact = round(exact4, 3),
  simulated = round(mc4, 3),
  std_error = round(se4, 3)
)

cat("\nQ4: Exact and simulated waiting times\n")
print(comparison4, row.names = FALSE)


# Question 5: Construct the fair-coin win probability matrix.

M <- pmat(k = 3, p = 0.5, nreps = 10000)

cat("\nQ5: Win probabilities, entry (i, j) = P(row i beats column j)\n")
print(round(M, 3))


# Question 5(a): Best replies.

cat("\nQ5(a): Best replies\n")
print(best_reply(M), row.names = FALSE)


# Question 5(b): Nearly even pairs and a non-transitive cycle.
# With nreps = 10000 each estimate has standard error at most
# sqrt(0.25 / 10000) = 0.005. A pair within 0.03 of 0.5 (6 standard
# errors) is treated as evenly matched, and a win above 0.55
# (10 standard errors) counts as clear.

cat("\nQ5(b): Nearly even pairs\n")
print(even_pairs(M), row.names = FALSE)

cat("\nQ5(b): Shortest non-transitive cycle\n")
show_cycle(M)

# Any three patterns either include an evenly matched pair or are
# transitive (one of them beats both others), so there is no cycle
# of length 3. In the cycle found, the opposite corners TTH/HHT and
# THH/HTT are even pairs, which is why it needs four patterns. Each
# step is won with probability between 2/3 and 3/4.


# Question 5(c): Compare waiting times with duel performance.

patterns3 <- all_patterns(3)
waits3 <- numeric(nrow(patterns3))

for (i in seq_len(nrow(patterns3))) {
  waits3[i] <- ewait(patterns3[i, ])
}

avg_win <- rowMeans(M, na.rm = TRUE)

comparison5 <- data.frame(
  pattern = rownames(M),
  expected_wait = waits3,
  mean_duel_win = round(avg_win, 3)
)

comparison5 <- comparison5[order(comparison5$expected_wait), ]

cat("\nQ5(c): Waiting times and average duel win probabilities\n")
print(comparison5, row.names = FALSE)

cat("\nQ5(c): Same expected wait (8 tosses), very different duels\n")
cat("P(HTT beats TTH):", round(M["HTT", "TTH"], 3), "\n")
cat("P(THH beats HHT):", round(M["THH", "HHT"], 3), "\n")
cat("Q5(c): Expected waits 14 against 8, yet an even duel\n")
cat("P(TTT beats TTH):", round(M["TTT", "TTH"], 3), "\n")

# On average the duel results follow the expected waits (TTT and HHH
# are last), but the average hides the head-to-head results. HTT and
# TTH both wait 8 tosses on average, yet HTT wins 3 duels in 4,
# while TTT, which takes almost twice as long, still draws with TTH.


# Question 6: Duel duration and Monte Carlo variation for THH
# against HHH (THH is the best reply to HHH in Question 5(a)).

a6 <- c(0, 1, 1)
b6 <- c(1, 1, 1)

# Duel lengths for plots 1 and 2.
games <- race(list(a6, b6), nreps = 10000)

duration <- games[, "tosses"]
proportions <- tabulate(duration) / length(duration)
median_tosses <- median(duration)

# P(duel takes <= 5 tosses) is exactly 1/2, so the simulated median
# can come out as 5 or 6 depending on the random draws.

# 50 independent estimates of P(THH beats HHH) for plot 3.
estimates <- numeric(50)

for (i in seq_len(50)) {
  estimates[i] <- pwin(a6, b6, nreps = 2000)
}

oldpar <- par(no.readonly = TRUE)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1))


# Plot 1: Empirical cumulative distribution of the duel length.

freq <- table(duration)
xval <- as.numeric(names(freq))
cumprob <- cumsum(as.numeric(freq)) / length(duration)

plot(xval, cumprob, type = "s", ylim = c(0, 1),
     xlab = "Number of tosses", ylab = "Cumulative probability",
     main = "Cumulative duel length: THH vs HHH")
lines(c(median_tosses, median_tosses), c(0, 1), col = "red", lty = 2)
text(median_tosses, 0.2, paste("Median =", median_tosses),
     pos = 4, col = "red")


# Plot 2: Proportion of duels ending at each toss.

plot(seq_along(proportions), proportions, type = "h",
     xlim = c(0, min(30, length(proportions))),
     xlab = "Number of tosses", ylab = "Proportion of duels",
     main = "Duel finish probabilities: THH vs HHH")
lines(c(median_tosses, median_tosses), c(0, max(proportions)),
      col = "red", lty = 2)
text(median_tosses, max(proportions), paste("Median =", median_tosses),
     pos = 4, col = "red")


# Plot 3: Monte Carlo variation around the exact value 7/8
# (Question 3 shows HHH wins with probability 1/8).

reference <- 7 / 8
limits <- range(c(estimates, reference)) + c(-0.03, 0.03)

plot(estimates, pch = 19, ylim = limits,
     xlab = "Simulation run", ylab = "Estimated win probability",
     main = "P(THH beats HHH) in 50 runs of pwin (nreps = 2000)")
lines(c(1, 50), c(reference, reference), col = "red", lty = 2)
text(25, reference + 0.02, "Exact value = 0.875", col = "red")

layout(1)
par(oldpar)


# Question 7: Investigate the effect of coin bias.

for (bias in c(0.25, 0.7)) {

  M_bias <- pmat(k = 3, p = bias, nreps = 20000)

  cat("\nQ7: Best replies for p =", bias, "\n")
  print(best_reply(M_bias), row.names = FALSE)

  cat("\nQ7: Nearly even pairs for p =", bias, "\n")
  pairs <- even_pairs(M_bias)
  if (nrow(pairs) == 0) {
    cat("None\n")
  } else {
    print(pairs, row.names = FALSE)
  }

  cat("\nQ7: Non-transitive cycle for p =", bias, "\n")
  show_cycle(M_bias)

  cat("\nQ7: Cycle with every step above 0.6 for p =", bias, "\n")
  show_cycle(M_bias, cutoff = 0.6)
}

# With a biased coin the best replies change and become uneven. For
# p = 0.25 the replies to patterns with many heads win almost every
# duel (THH beats HHH with probability about 0.985), while the best
# reply to TTT wins only about 0.57. For p = 0.7 the best reply to
# HHH drops from 0.875 to about 0.66.
# With a fair coin, the symmetry between H and T makes many pairs
# exactly even, which forces the shortest cycle to have length 4.
# The bias breaks most of these ties (only THH/HHT stays close, at
# about 0.51, for p = 0.7), so cycles of length 3 appear. However
# non-transitivity becomes weaker: some steps of these cycles are
# narrow wins (about 0.57). If every step must exceed 0.6, there is
# no cycle for p = 0.25 and only a cycle of length 4 for p = 0.7,
# while the fair-coin cycle has every step above 0.65.


#-----------------------------------------------------------------------
# Conclusion
#-----------------------------------------------------------------------

# Expected waiting times for single patterns do not determine their
# head-to-head results: patterns with the same expected wait can win
# 3 duels in 4 against each other, and the cycles show there is no
# best pattern at all. A duel depends on how the end of one pattern
# overlaps the start of the other, that is, on the joint behaviour
# of the two waiting times on the same stream, not on their separate
# means. Comparing isolated summaries treats the two waiting times
# as independent when they are not, and leads to confident but wrong
# conclusions. The same hazard appears whenever quantities that share
# a common source of randomness are summarised one at a time.


#-----------------------------------------------------------------------
# Running time
#-----------------------------------------------------------------------

cat(
  "\nElapsed time (seconds):",
  round(proc.time()["elapsed"] - start_time, 2),
  "\n"
)
