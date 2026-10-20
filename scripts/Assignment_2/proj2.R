# Statistical Programming
# Practical 2: Sequence duels and the perils of independence
#
# Group 3:
# Alejandra Campo Archbold (s3023700)
# Weihang Gu (s2932201)
# Weichen Zhao (s2929481)
#
# Contributions (approximate proportions):
# - Alejandra Campo Archbold: []
# - Weihang Gu: []
# - Weichen Zhao: []
#
# GitHub URL:
# https://github.com/alejandrarchbold/Statistical_Programming.git


#-----------------------------------------------------------------------------
# Overview
#-----------------------------------------------------------------------------

# A coin with P(heads) = p is tossed repeatedly (1 = heads, 0 = tails). Two
# players each choose a pattern of k outcomes and the same stream of tosses
# is watched until one of the patterns appears as k consecutive tosses; that
# player wins. We simulate duels, estimate win probabilities between all
# patterns, compute exact expected waiting times of single patterns, look for
# non-transitive cycles and study how the bias p changes the game.


#-----------------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------------


# Q1: firsthit

firsthit <- function(x, a) {
  # Description: finds where pattern a first appears in the toss stream x.
  # Inputs: x (vector of 0/1 tosses), a (0/1 pattern).
  # Output: the toss at which a is first completed, or NA if it never is.
  # How it works: checks each window of k tosses in order, where the window
  #   ending at toss m is x[(m-k+1):m].
  
  k <- length(a)
  n <- length(x)
  if (k == 0) stop("Please add a 0/1 pattern.")
  if (n < k) return(NA)
  
  for (m in k:n) {
    if (all(x[(m - k + 1):m] == a)) return(m)
  }
  NA
}


# Shared simulation engine for Q2, Q3 and Q4

race <- function(pats, p = 0.5, nreps = 1, block = 50) {
  # Description: simulates nreps independent coin streams. Each stream is
  #   watched for every pattern in pats and stops when the first of them
  #   appears.
  # Inputs: pats (list of 0/1 patterns of the same length k), p (P(heads)),
  #   nreps (number of streams), block (tosses added per extension).
  # Output: matrix with one row per stream and columns winner (position in
  #   pats of the pattern that appeared first) and tosses (toss at which it
  #   appeared).
  # How it works: row i of X is stream i. While some streams are unfinished,
  #   only those rows receive a new block of tosses. Each pattern is compared
  #   with every k-window of every row at once; a row finishes when some
  #   window matches, and the winner is the pattern with the earliest match.
  
  k <- length(pats[[1]])
  res <- matrix(NA, nreps, 2, dimnames = list(NULL, c("winner", "tosses")))
  active <- 1:nreps
  X <- matrix(0, nreps, 0)
  
  while (length(active) > 0) {
    new <- matrix(rbinom(length(active) * block, 1, p), length(active))
    X <- cbind(X, new)
    ends <- k:ncol(X)
    
    # hit[i, q] = first toss at which pattern q appears in stream i (Inf if
    # it has not appeared). Window ending at toss e matches pattern a when
    # X[, e - k + j] == a[j] for j = 1, ..., k.
    hit <- sapply(pats, function(a) {
      match <- matrix(TRUE, nrow(X), length(ends))
      for (j in 1:k) match <- match & X[, ends - k + j, drop = FALSE] == a[j]
      ifelse(rowSums(match) > 0, ends[max.col(match, "first")], Inf)
    })
    hit <- matrix(hit, nrow(X))
    
    first <- apply(hit, 1, min)
    done <- first < Inf
    res[active[done], "winner"] <- apply(hit[done, , drop = FALSE], 1,
                                         which.min)
    res[active[done], "tosses"] <- first[done]
    
    X <- X[!done, , drop = FALSE]
    active <- active[!done]
  }
  res
}

check_pair <- function(a, b) {
  # Description: stops with an error unless a and b are distinct 0/1
  #   patterns of the same positive length.
  if (length(a) == 0 || length(a) != length(b)) {
    stop("a and b must have the same positive length.")
  }
  if (!all(c(a, b) %in% c(0, 1))) stop("Patterns must contain only 0 and 1.")
  if (all(a == b)) stop("a and b must be different.")
}


# Q2: duel

duel <- function(a, b, p = 0.5) {
  # Description: simulates one duel between patterns a and b.
  # Inputs: a and b (distinct 0/1 patterns), p (P(heads)).
  # Output: c(winner, tosses), with winner 1 if a wins and 0 if b wins.
  # How it works: runs race for a single stream.
  check_pair(a, b)
  r <- race(list(a, b), p, nreps = 1)
  c(winner = as.numeric(r[1, "winner"] == 1),
    tosses = as.numeric(r[1, "tosses"]))
}


# Q3: pwin

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Description: estimates the probability that a appears before b.
  # Inputs: a and b (distinct 0/1 patterns), p (P(heads)),
  #   nreps (number of simulated duels).
  # Output: the proportion of the nreps duels won by a.
  # How it works: runs race for nreps streams at once.
  check_pair(a, b)
  r <- race(list(a, b), p, nreps)
  mean(r[, "winner"] == 1)
}


# Q4: ewait and simwait

ewait <- function(a, p = 0.5) {
  # Description: exact expected waiting time of pattern a.
  # Inputs: a (0/1 pattern), p (P(heads)).
  # Output: the expected number of tosses until a first appears.
  # How it works: a has an overlap of length j if its first j outcomes equal
  #   its last j outcomes (j = k always counts). Each overlap adds
  #   1 / P(a[1:j]), with P(a[1:j]) = p^h (1-p)^(j-h) and h the number of
  #   heads in a[1:j].
  k <- length(a)
  wait <- 0
  for (j in 1:k) {
    if (all(a[1:j] == a[(k - j + 1):k])) {
      h <- sum(a[1:j])
      wait <- wait + 1 / (p^h * (1 - p)^(j - h))
    }
  }
  wait
}

simwait <- function(a, p = 0.5, nreps = 20000) {
  # Description: simulated waiting times of pattern a, to check ewait.
  # Output: the waiting time in each of nreps independent streams.
  # How it works: runs race with a as the only pattern.
  race(list(a), p, nreps)[, "tosses"]
}


# Q5: pmat

all_patterns <- function(k) {
  # Description: all 2^k patterns of length k, one per row, in binary
  #   counting order with the first toss as the most significant digit
  #   (0 = T, 1 = H). For k = 3: TTT, TTH, THT, THH, HTT, HTH, HHT, HHH.
  unname(as.matrix(expand.grid(rep(list(0:1), k)))[, k:1, drop = FALSE])
}

label <- function(a) {
  # Description: writes a pattern with letters, e.g. c(1, 0, 1) -> "HTH".
  paste(c("T", "H")[a + 1], collapse = "")
}

pmat <- function(k, p = 0.5, nreps = 10000) {
  # Description: estimated win probabilities between all 2^k patterns.
  # Inputs: k (pattern length), p (P(heads)), nreps (duels per pair).
  # Output: matrix M where M[i, j] = P(pattern i beats pattern j), with
  #   patterns ordered as in all_patterns(k) and NA on the diagonal.
  # How it works: each pair is simulated once; since a duel always has one
  #   winner, M[j, i] = 1 - M[i, j].
  pats <- all_patterns(k)
  n <- nrow(pats)
  labs <- apply(pats, 1, label)
  M <- matrix(NA, n, n, dimnames = list(labs, labs))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      M[i, j] <- pwin(pats[i, ], pats[j, ], p, nreps)
      M[j, i] <- 1 - M[i, j]
    }
  }
  M
}


# Q5(a): best_reply

best_reply <- function(M) {
  # Description: best reply to every pattern.
  # Input: M (win probability matrix from pmat).
  # Output: for each pattern, the pattern that beats it with the highest
  #   probability and that probability.
  # How it works: column j of M holds the chances of beating pattern j, so
  #   the best reply is the row of the largest entry in column j.
  best <- apply(M, 2, which.max)
  data.frame(pattern = colnames(M), reply = rownames(M)[best],
             probability = round(M[cbind(best, 1:ncol(M))], 3))
}


# Q5(b): even_pairs, find_cycle and show_cycle

even_pairs <- function(M, tol = 0.015) {
  # Description: pairs of patterns that are (essentially) evenly matched.
  # Inputs: M (win matrix), tol (allowed distance of P from 0.5).
  # Output: prints the pairs with |P(a beats b) - 0.5| < tol.
  pos <- which(upper.tri(M) & abs(M - 0.5) < tol, arr.ind = TRUE)
  if (nrow(pos) == 0) {
    cat("No evenly matched pairs.\n")
  } else {
    print(data.frame(pattern_a = rownames(M)[pos[, 1]],
                     pattern_b = colnames(M)[pos[, 2]],
                     p_a_wins = round(M[pos], 3)), row.names = FALSE)
  }
}

find_cycle <- function(M, cutoff = 0.55) {
  # Description: shortest non-transitive cycle of distinct patterns, each
  #   beating the next with probability above cutoff.
  # Inputs: M (win matrix), cutoff (minimum win probability of each step).
  # Output: pattern indices of the cycle (start repeated at the end), or
  #   NULL if there is none.
  # How it works: beats[i, j] is TRUE if i clearly beats j. A breadth-first
  #   search from each start s visits patterns in order of the number of
  #   steps from s, so the first step back to s closes the shortest cycle
  #   through s. The shortest cycle over all starts is returned.
  beats <- M > cutoff
  beats[is.na(beats)] <- FALSE
  n <- nrow(M)
  best <- NULL
  
  for (s in 1:n) {
    parent <- rep(NA, n)
    queue <- s
    cycle <- NULL
    while (length(queue) > 0 && is.null(cycle)) {
      u <- queue[1]
      queue <- queue[-1]
      for (v in which(beats[u, ])) {
        if (v == s) {
          # follow the parents back from u to s to recover the path
          cycle <- c(u, s)
          while (cycle[1] != s) cycle <- c(parent[cycle[1]], cycle)
          break
        }
        if (is.na(parent[v])) {
          parent[v] <- u
          queue <- c(queue, v)
        }
      }
    }
    if (!is.null(cycle) && (is.null(best) || length(cycle) < length(best))) {
      best <- cycle
    }
  }
  best
}

show_cycle <- function(M, cutoff = 0.55) {
  # Description: prints the cycle from find_cycle and the win probability
  #   of each of its steps.
  cycle <- find_cycle(M, cutoff)
  if (is.null(cycle)) {
    cat("No clear non-transitive cycle found.\n")
    return(invisible(NULL))
  }
  labs <- rownames(M)
  cat("Cycle of length", length(cycle) - 1, ":",
      paste(labs[cycle], collapse = " -> "), "\n")
  for (i in 1:(length(cycle) - 1)) {
    cat(" ", labs[cycle[i]], "beats", labs[cycle[i + 1]], ":",
        round(M[cycle[i], cycle[i + 1]], 3), "\n")
  }
}


#-----------------------------------------------------------------------------
# Run code and analysis
#-----------------------------------------------------------------------------

set.seed(2026)


# Q1: Test firsthit on streams that can be checked by hand

# Stream H T T H T H H H T: HTH is completed at toss 6 and THH at toss 7.
x <- c(1, 0, 0, 1, 0, 1, 1, 1, 0)
cat("Q1: firsthit on H T T H T H H H T\n")
cat("  HTH (expected 6):", firsthit(x, c(1, 0, 1)), "\n")
cat("  THH (expected 7):", firsthit(x, c(0, 1, 1)), "\n")
cat("  HHH (expected 8):", firsthit(x, c(1, 1, 1)), "\n")
cat("  TTT (expected NA):", firsthit(x, c(0, 0, 0)), "\n")
cat("  HTH in H T (expected NA):", firsthit(c(1, 0), c(1, 0, 1)), "\n")


# Q2: Simulate one duel

cat("\nQ2: One duel between HTH and THH\n")
print(duel(c(1, 0, 1), c(0, 1, 1)))


# Q3: Estimate a win probability

cat("\nQ3: P(HHH beats THH), p = 0.5:",
    pwin(c(1, 1, 1), c(0, 1, 1)), "\n")


# Q4: Compare exact and simulated waiting times

cat("\nQ4: Fair coin length-2 patterns (expected 6, 6, 4, 4):",
    ewait(c(1, 1)), ewait(c(0, 0)), ewait(c(1, 0)), ewait(c(0, 1)), "\n")

patterns4 <- list(c(1, 1, 1), c(1, 0, 1), c(1, 1, 0), c(1, 0, 1), c(1, 1, 0))
probs4 <- c(0.5, 0.5, 0.5, 0.3, 0.7)
exact4 <- numeric(5)
sim4 <- numeric(5)
for (i in 1:5) {
  exact4[i] <- ewait(patterns4[[i]], probs4[i])
  sim4[i] <- mean(simwait(patterns4[[i]], probs4[i]))
}
cat("\nQ4: Exact and simulated expected waiting times\n")
print(data.frame(pattern = sapply(patterns4, label), p = probs4,
                 exact = round(exact4, 3), simulated = round(sim4, 3)),
      row.names = FALSE)


# Q5: Fair-coin win matrix, entry (i, j) = P(row i beats column j)

M <- pmat(k = 3, p = 0.5, nreps = 20000)
cat("\nQ5: Win probability matrix (fair coin)\n")
print(round(M, 2))


# Q5(a): Best replies

cat("\nQ5(a): Best replies (fair coin)\n")
print(best_reply(M), row.names = FALSE)


# Q5(b): Evenly matched pairs and the shortest non-transitive cycle

# With nreps = 20000 the standard error of an estimate is at most 0.0035,
# so pairs within 0.015 of 0.5 are treated as evenly matched, and a pattern
# "clearly beats" another when it wins with probability above 0.55.
cat("\nQ5(b): Evenly matched pairs (fair coin)\n")
even_pairs(M)
cat("\nQ5(b): Shortest non-transitive cycle (fair coin)\n")
show_cycle(M)

# The ten evenly matched pairs remove possible steps from a cycle. Every
# cycle of three patterns that never loses passes through one of them, so
# no cycle of clear wins has length 3, and the shortest has length 4:
# TTH -> THH -> HHT -> HTT -> TTH, each step won with probability 2/3-3/4.


# Q5(c): Expected waiting time against head-to-head performance

pats <- all_patterns(3)
ew <- apply(pats, 1, ewait)
avg_win <- rowMeans(M, na.rm = TRUE)
cat("\nQ5(c): Ranking by expected wait and by average duel win probability\n")
print(data.frame(pattern = rownames(M), ewait = ew, ewait_rank = rank(ew),
                 avg_win = round(avg_win, 3), duel_rank = rank(-avg_win)),
      row.names = FALSE)

# Duels between patterns with the same expected waiting time; the most
# one-sided are listed first.
pos <- which(upper.tri(M) & outer(ew, ew, "=="), arr.ind = TRUE)
win <- pmax(M[pos], 1 - M[pos])
winner_first <- M[pos] >= 0.5
same <- data.frame(
  winner = ifelse(winner_first, rownames(M)[pos[, 1]], rownames(M)[pos[, 2]]),
  loser = ifelse(winner_first, rownames(M)[pos[, 2]], rownames(M)[pos[, 1]]),
  ewait = ew[pos[, 1]], p_winner = round(win, 3))
cat("\nQ5(c): Duels between patterns with equal expected waiting time\n")
print(same[order(-win), ], row.names = FALSE)

# The four patterns with the shortest expected wait (8 tosses) are also the
# best in duels, but the expected wait cannot separate them: among them HTT
# beats TTH and THH beats HHT with probability about 0.75, and the cycle
# above shows that none of the four is best.


# Q6: Duel length and Monte Carlo variability, for THH against HHH
# (THH is the best reply to HHH found in Q5(a)).

a <- c(0, 1, 1)
b <- c(1, 1, 1)
par(mfrow = c(1, 2))

# Plot 1: estimated distribution of the number of tosses in a duel
games <- race(list(a, b), nreps = 20000)
freq <- tabulate(games[, "tosses"]) / 20000
med <- median(games[, "tosses"])
plot(seq_along(freq), freq, type = "h", xlim = c(0, 30),
     xlab = "Number of tosses", ylab = "Proportion of duels",
     main = "Duel length: THH vs HHH")
lines(c(med, med), c(0, max(freq)), col = "red", lty = 2, lwd = 2)
text(med, max(freq), paste("Median =", med), pos = 4, col = "red")

# Plot 2: 50 independent estimates of P(THH beats HHH) with nreps = 2000.
# HHH can only win if the first three tosses are HHH: after any tail, THH
# is completed before HHH can be. So the exact value is 1 - 1/8 = 0.875.
# The dashed lines are +- 2 standard errors of a single estimate.
est <- numeric(50)
for (i in 1:50) est[i] <- pwin(a, b, nreps = 2000)
se <- sqrt(0.875 * 0.125 / 2000)
plot(est, pch = 19, ylim = 0.875 + c(-4, 4) * se,
     xlab = "Simulation run", ylab = "Estimated win probability",
     main = "P(THH beats HHH) in 50 runs")
lines(c(1, 50), c(0.875, 0.875), col = "red", lwd = 2)
lines(c(1, 50), 0.875 + c(2, 2) * se, col = "red", lty = 2)
lines(c(1, 50), 0.875 - c(2, 2) * se, col = "red", lty = 2)
text(25, 0.875 + 3.5 * se, "Exact value = 0.875", col = "red")
par(mfrow = c(1, 1))


# Q7: Effect of the coin bias

for (bias in c(0.25, 0.7)) {
  M_bias <- pmat(k = 3, p = bias, nreps = 20000)
  cat("\nQ7: Best replies for p =", bias, "\n")
  print(best_reply(M_bias), row.names = FALSE)
  cat("\nQ7: Evenly matched pairs for p =", bias, "\n")
  even_pairs(M_bias)
  cat("\nQ7: Shortest non-transitive cycle for p =", bias, "\n")
  show_cycle(M_bias)
}

# With p = 0.25 or 0.7 the best replies change (e.g. for p = 0.25 the best
# reply to TTH is TTT and to HTT is THT) and the second player's advantage
# grows for most patterns: at p = 0.25 THH beats HHH about 98% of the time.
# The ties of the fair coin come from its symmetry and disappear (only
# THH vs HHT at p = 0.7 stays close, at about 0.51), so a cycle of length 3
# becomes possible. The game stays non-transitive, but some steps of these
# shorter cycles are only narrow wins (about 0.57).


#-----------------------------------------------------------------------------
# Conclusion
#-----------------------------------------------------------------------------

# Isolated expected waiting times do not determine who wins a duel:
# patterns with the same expected wait win against each other 3 times in 4,
# and the cycles show there is no best pattern to rank at all. The waiting
# times of the two players are dependent, because both watch the same
# stream, and the outcome depends on how the end of one pattern overlaps
# the start of the other. Comparing two separate summaries treats the
# waiting times as if they were independent, which ignores exactly this
# dependence and gives confident but wrong conclusions.
