
# Statistical Programming
# Practical 2: Sequence duels and the perils of independence
#
# Group 3:
# Alejandra Campo Archbold (s3023700)
# Weihang Gu (s2932201)
# Weichen Zhao (s2929481)
#
# Contributions:
# - Alejandra Campo Archbold: []
# - Weihang Gu: []
# - Weichen Zhao: []
#
# GitHub URL:



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


# Q1: firsthit

firsthit <- function(x, a) {
  # Description: finds the first occurrence of a in x.
  # Inputs: x (coin tosses), a (the pattern to find).
  # Output: the ending position of the first match, or NA.
  # How it works: checks each window of length k in order.
  
  k <- length(a)
  n <- length(x)
  
  if (k == 0) {
    stop("Please enter a non-empty pattern.")
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


# Shared simulation function for Q2, Q3 and Q4

race <- function(pats, p = 0.5, nreps = 1, block = 20) {
  # Description: simulates coin streams until a pattern appears.
  # Inputs: pats (list of patterns), p (chance of heads),
  #         nreps (number of streams), block (batch size).
  # Output: matrix containing winner indices and toss counts.
  # How it works: checks patterns using a sliding window.
  #         Only unfinished streams receive new blocks.
  
  # Check that pats is a non-empty list.
  if (!is.list(pats) || length(pats) == 0) {
    stop("pats must be a non-empty list.")
  }
  
  k <- length(pats[[1]])
  
  if (k == 0) {
    stop("Patterns cannot be empty.")
  }
  
  # Check each pattern.
  for (i in seq_along(pats)) {
    
    a <- pats[[i]]
    
    if (!is.numeric(a) || length(a) != k ||
        anyNA(a) || !all(a %in% c(0, 1))) {
      stop("Patterns must have equal lengths and contain 0/1.")
    }
    
    # Check that no two patterns are identical.
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        if (all(a == pats[[j]])) {
          stop("Patterns must be different.")
        }
      }
    }
  }
  
  # Check the simulation parameters.
  if (length(p) != 1 || !is.numeric(p) ||
      !is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be between 0 and 1.")
  }
  
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
  
  # Store the most recent k tosses for each stream.
  recent <- matrix(0, nrow = nreps, ncol = k)
  
  winner <- rep(NA_integer_, nreps)
  tosses <- integer(nreps)
  
  active <- seq_len(nreps)
  
  # Continue until every stream has a winner.
  while (length(active) > 0) {
    
    # Generate new tosses for unfinished streams.
    new <- matrix(
      rbinom(length(active) * block, size = 1, prob = p),
      nrow = length(active),
      ncol = block
    )
    
    for (j in seq_len(block)) {
      
      unfinished <- which(is.na(winner[active]))
      
      if (length(unfinished) == 0) {
        break
      }
      
      rows <- active[unfinished]
      current <- new[unfinished, j]
      
      # Update the sliding window.
      if (k > 1) {
        recent[rows, 1:(k - 1)] <-
          recent[rows, 2:k, drop = FALSE]
      }
      
      recent[rows, k] <- current
      tosses[rows] <- tosses[rows] + 1
      
      # A match is possible after at least k tosses.
      ready <- rows[tosses[rows] >= k]
      
      if (length(ready) == 0) {
        next
      }
      
      windows <- recent[ready, , drop = FALSE]
      
      # Compare the window with each pattern.
      for (q in seq_along(pats)) {
        
        a <- pats[[q]]
        
        target <- matrix(
          a,
          nrow = length(ready),
          ncol = k,
          byrow = TRUE
        )
        
        matched <- rowSums(windows == target) == k
        
        winner[ready[matched]] <- q
      }
    }
    
    # Keep only streams without a winner.
    active <- which(is.na(winner))
  }
  
  return(cbind(winner = winner, tosses = tosses))
}


# Q2: duel

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


# Q3: pwin

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Description: estimates the probability that a beats b.
  # Inputs: a and b (patterns), p (heads), nreps (duels).
  # Output: the proportion of games won by a.
  # How it works: counts the wins for a across simulated duels.
  
  results <- race(list(a, b), p, nreps = nreps)
  
  return(mean(results[, "winner"] == 1))
}


# Q4: ewait

ewait <- function(a, p = 0.5) {
  # Description: calculates the exact mean waiting time for a.
  # Inputs: a (pattern), p (chance of heads).
  # Output: expected tosses until a first appears.
  # How it works: uses matching prefixes and suffixes of a.
  
  k <- length(a)
  
  if (k == 0 || anyNA(a) || any(!(a %in% c(0, 1)))) {
    stop("Please enter a non-empty 0/1 pattern.")
  }
  
  if (length(p) != 1 || !is.numeric(p) ||
      !is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be between 0 and 1.")
  }
  
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


# Q4: simwait

simwait <- function(a, p = 0.5, nreps = 3000) {
  # Description: simulates waiting times for a single pattern.
  # Inputs: a (pattern), p (heads), nreps (independent runs).
  # Output: the waiting time for each simulated stream.
  # How it works: calls race with one pattern instead of two.
  
  results <- race(list(a), p, nreps = nreps)
  
  return(results[, "tosses"])
}


# Q5: all_patterns

all_patterns <- function(k) {
  # Description: generates all patterns of length k.
  # Input: k (positive integer).
  # Output: a matrix of patterns in binary counting order.
  # How it works: creates all 0/1 combinations with expand.grid.
  
  if (length(k) != 1 || !is.numeric(k) ||
      !is.finite(k) || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }
  
  patterns <- expand.grid(rep(list(0:1), k))
  patterns <- as.matrix(patterns)
  patterns <- patterns[, k:1, drop = FALSE]
  
  return(unname(patterns))
}


# Q5: label

label <- function(a) {
  # Description: converts a binary pattern to H/T notation.
  # Input: a (0/1 pattern).
  # Output: a character label.
  # How it works: maps 0 to T and 1 to H.
  
  return(paste(c("T", "H")[a + 1], collapse = ""))
}


# Q5: pmat

pmat <- function(k, p = 0.5, nreps = 10000) {
  # Description: estimates win probabilities for all pattern pairs.
  # Inputs: k (pattern length), p (heads), nreps (games per pair).
  # Output: matrix M, with M[i,j] = P(row i beats column j).
  # How it works: simulates each pair once and fills both entries.
  
  patterns <- all_patterns(k)
  n <- nrow(patterns)
  
  labels <- apply(patterns, 1, label)
  
  M <- matrix(NA_real_, nrow = n, ncol = n)
  
  rownames(M) <- labels
  colnames(M) <- labels
  
  for (i in seq_len(n - 1)) {
    
    for (j in (i + 1):n) {
      
      a <- patterns[i, ]
      b <- patterns[j, ]
      
      prob <- pwin(a, b, p = p, nreps = nreps)
      
      M[i, j] <- prob
      M[j, i] <- 1 - prob
    }
  }
  
  return(M)
}


# Q5(a): best_reply

best_reply <- function(M) {
  # Description: finds the best reply to each pattern.
  # Input: M (win probability matrix).
  # Output: each pattern's best reply and its win probability.
  # How it works: finds the largest entry in each column.
  
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


# Q5(b): even_pairs

even_pairs <- function(M, tol = 0.03) {
  # Description: finds pairs with nearly equal win probabilities.
  # Inputs: M (win matrix), tol (allowed distance from 0.5).
  # Output: a table of pairs and their estimated probabilities.
  # How it works: checks the upper triangle of M.
  
  pos <- which(
    upper.tri(M) & abs(M - 0.5) < tol,
    arr.ind = TRUE
  )
  
  result <- data.frame(
    pattern_a = rownames(M)[pos[, 1]],
    pattern_b = colnames(M)[pos[, 2]],
    p_a_wins = round(M[pos], 3)
  )
  
  return(result)
}


# Q5(b): find_cycle

find_cycle <- function(M, cutoff = 0.55) {
  # Description: finds a shortest clear non-transitive cycle.
  # Inputs: M (win matrix), cutoff (minimum win probability).
  # Output: pattern indices forming the cycle, or NULL.
  # How it works: uses breadth-first search from each pattern.
  
  n <- nrow(M)
  
  if (n < 3) {
    return(NULL)
  }
  
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
    
    if (!is.null(cycle) && length(cycle) >= 4 &&
        (is.null(best) || length(cycle) < length(best))) {
      best <- cycle
    }
  }
  
  return(best)
}


# Q5(b): show_cycle

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
  
  cat(
    "Cycle of length", length(cycle) - 1, ":",
    paste(labels[cycle], collapse = " -> "), "\n"
  )
  
  for (i in seq_len(length(cycle) - 1)) {
    
    a <- cycle[i]
    b <- cycle[i + 1]
    
    cat(
      labels[a], "beats", labels[b], ":",
      round(M[a, b], 3), "\n"
    )
  }
}


#-----------------------------------------------------------------------
# Run code and analysis
#-----------------------------------------------------------------------

start_time <- proc.time()["elapsed"]

set.seed(2026)


# Q1: Test firsthit on simple streams.

x <- c(1, 0, 0, 1, 0, 1, 1)
a <- c(1, 0, 1)

cat("Q1: First hit of HTH:", firsthit(x, a), "\n")
cat("Q1: HHH absent:", firsthit(x, c(1, 1, 1)), "\n")
cat("Q1: Short stream:", firsthit(c(1, 0), a), "\n")


# Q2: Simulate one duel.

cat("\nQ2: One duel between HTH and THH\n")

print(duel(c(1, 0, 1), c(0, 1, 1)))


# Q3: Estimate the win probability.

p_fair <- pwin(c(1, 1, 1), c(0, 1, 1))

cat("\nQ3: Win probability estimate\n")
cat("P(HHH beats THH):", round(p_fair, 3), "\n")


# Q4: Compare exact and simulated waiting times.

patterns4 <- list(
  c(1, 1, 1),
  c(1, 1, 0),
  c(1, 0, 1),
  c(1, 1, 0)
)

probs4 <- c(0.5, 0.5, 0.5, 0.7)

exact4 <- numeric(length(patterns4))
mc4 <- numeric(length(patterns4))

for (i in seq_along(patterns4)) {
  
  exact4[i] <- ewait(patterns4[[i]], p = probs4[i])
  
  mc4[i] <- mean(simwait(patterns4[[i]], p = probs4[i]))
}

comparison4 <- data.frame(
  pattern = vapply(patterns4, label, character(1)),
  p = probs4,
  exact = round(exact4, 3),
  simulated = round(mc4, 3)
)

cat("\nQ4: Exact and simulated waiting times\n")
print(comparison4, row.names = FALSE)


# Q5: Construct the fair-coin win probability matrix.

M <- pmat(k = 3, p = 0.5, nreps = 10000)


# Q5(a): Best replies.

cat("\nQ5(a): Best replies\n")
print(best_reply(M), row.names = FALSE)


# Q5(b): Nearly even pairs and a non-transitive cycle.

cat("\nQ5(b): Nearly even pairs\n")
print(even_pairs(M), row.names = FALSE)

cat("\nQ5(b): Shortest non-transitive cycle\n")
show_cycle(M)


# Q5(c): Compare waiting times with duel performance.

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

comparison5 <- comparison5[
  order(comparison5$expected_wait),
]

cat("\nQ5(c): Waiting times and average duel win probabilities\n")
print(comparison5, row.names = FALSE)

cat(
  "P(TTH beats THH):",
  round(M["TTH", "THH"], 3),
  "\n"
)


# Q6: Duel duration and Monte Carlo variation.

a6 <- c(0, 1, 1)
b6 <- c(1, 1, 1)

# Use the same duel data for the first two plots.
games <- race(list(a6, b6), nreps = 10000)

duration <- games[, "tosses"]
counts <- tabulate(duration)
proportions <- counts / length(duration)
median_tosses <- median(duration)

# Use 50 independent estimates for the third plot.
estimates <- numeric(50)

for (i in seq_len(50)) {
  estimates[i] <- pwin(a6, b6, nreps = 2000)
}

# Arrange the three plots.
oldpar <- par(no.readonly = TRUE)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1))


# Plot 1: Empirical cumulative distribution.

freq <- table(duration)
xval <- as.numeric(names(freq))
cumprob <- cumsum(as.numeric(freq)) / length(duration)

plot(
  xval, cumprob,
  type = "s",
  ylim = c(0, 1),
  xlab = "Number of tosses",
  ylab = "Cumulative probability",
  main = "Cumulative duel length"
)

abline(v = median_tosses, col = "red", lty = 2)

legend(
  "bottomright",
  legend = paste("Median =", median_tosses),
  col = "red",
  lty = 2,
  bty = "n"
)


# Plot 2: Discrete duel finish probabilities.

plot(
  seq_along(proportions), proportions,
  type = "h",
  xlim = c(0, min(30, length(proportions))),
  xlab = "Number of tosses",
  ylab = "Proportion of duels",
  main = "Duel finish probabilities"
)

abline(v = median_tosses, col = "red", lty = 2)

legend(
  "topright",
  legend = paste("Median =", median_tosses),
  col = "red",
  lty = 2,
  bty = "n"
)


# Plot 3: Monte Carlo variation.

reference <- 7 / 8
limits <- range(c(estimates, reference)) + c(-0.03, 0.03)

plot(
  estimates,
  pch = 19,
  ylim = limits,
  xlab = "Simulation run",
  ylab = "Estimated win probability",
  main = "Monte Carlo variation"
)

lines(
  c(1, 50),
  c(reference, reference),
  col = 2,
  lty = 2
)

text(
  25,
  reference + 0.02,
  "Exact value = 0.875",
  col = 2
)

layout(1)
par(oldpar)


# Q7: Investigate the effect of coin bias.

for (bias in c(0.25, 0.7)) {
  
  M_bias <- pmat(k = 3, p = bias, nreps = 20000)
  
  cat("\nQ7: Best replies for p =", bias, "\n")
  print(best_reply(M_bias), row.names = FALSE)
  
  cat("\nQ7: Non-transitive cycle for p =", bias, "\n")
  show_cycle(M_bias)
}

# Coin bias changes best replies and can change the cycles found.


#-----------------------------------------------------------------------
# Conclusion
#-----------------------------------------------------------------------

# Expected waiting times for individual patterns do not determine
# their head-to-head win probabilities. Both patterns share the same
# toss stream, so their waiting times are dependent. Comparing only
# separate mean waiting times ignores this dependence and can give
# misleading rankings.


#-----------------------------------------------------------------------
# Running time
#-----------------------------------------------------------------------

cat(
  "\nElapsed time (seconds):",
  round(proc.time()["elapsed"] - start_time, 2),
  "\n"
)
