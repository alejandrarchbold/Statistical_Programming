
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


#-------------------------------------------------------------------------------
# Overview:
#
# Two players choose patterns of heads (1) and tails (0).
# They watch the same coin tosses, and the first pattern to appear wins.
# We simulate duels, find win probabilities and calculate exact waiting
# times. We then compare patterns, look for cycles and change the coin bias.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Function definitions
#-------------------------------------------------------------------------------

# Q1: firsthit

firsthit <- function(x, a){
  # Description: finds the first occurrence of pattern a in x.
  # Inputs: x (coin tosses), a (the pattern to find).
  # Output: the ending position of the first match, or NA.
  # How it works: checks each window of length k in order.
  
  k <- length(a)
  n <- length(x)
  
  if(k == 0){
    stop("Please add a 0/1 pattern.")
  }
  
  if(n < k){
    return(NA)
  }
  
  for(m in k:n){
    w <- x[(m-k+1):m]
    
    if(all(w == a)){
      return(m)
    }
  }
  
  return(NA)
}


# Q2 and Q3: simulate_duels

simulate_duels <- function(a, b, p=0.5, nreps=1, block=20){
  # Description: simulates several duels between a and b.
  # Inputs: a and b (patterns), p (chance of heads),
  #         nreps (number of duels), block (tosses per batch).
  # Output: a matrix of winners and toss counts for all duels.
  # How it works: each row is one duel. Both patterns use the same
  # tosses in that row. Only unfinished rows receive new tosses.
  
  k <- length(a)
  
  if(k == 0 || length(b) != k ||
     anyNA(a) || anyNA(b) ||
     any(!(a %in% c(0,1))) || any(!(b %in% c(0,1))) ||
     all(a == b)){
    stop("a and b must be different 0/1 patterns of equal length.")
  }
  
  if(length(p) != 1 || !is.numeric(p) ||
     !is.finite(p) || p <= 0 || p >= 1){
    stop("p must be between 0 and 1.")
  }
  
  if(length(nreps) != 1 || !is.numeric(nreps) ||
     !is.finite(nreps) || nreps < 1 ||
     nreps != floor(nreps)){
    stop("nreps must be a positive integer.")
  }
  
  if(length(block) != 1 || !is.numeric(block) ||
     !is.finite(block) || block < 1 ||
     block != floor(block)){
    stop("block must be a positive integer.")
  }
  
  # Keep only the most recent k tosses from each duel.
  recent <- matrix(0, nrow=nreps, ncol=k)
  winner <- rep(NA_integer_, nreps)
  tosses <- integer(nreps)
  unfinished <- seq_len(nreps)
  
  while(length(unfinished) > 0){
    
    new <- matrix(
      rbinom(length(unfinished)*block, size=1, prob=p),
      nrow=length(unfinished), ncol=block
    )
    
    for(j in seq_len(block)){
      active_rows <- is.na(winner[unfinished])
      active <- unfinished[active_rows]
      
      if(length(active) == 0){
        break
      }
      
      current <- new[active_rows, j]
      
      if(k > 1){
        recent[active, 1:(k-1)] <- recent[active, 2:k, drop=FALSE]
      }
      
      recent[active, k] <- current
      tosses[active] <- tosses[active] + 1
      ready <- active[tosses[active] >= k]
      
      if(length(ready) > 0){
        windows <- recent[ready, , drop=FALSE]
        
        # Compare the rows with a and b at the same time.
        target_a <- matrix(a, nrow=length(ready),
                           ncol=k, byrow=TRUE)
        target_b <- matrix(b, nrow=length(ready),
                           ncol=k, byrow=TRUE)
        
        match_a <- rowSums(windows == target_a) == k
        match_b <- rowSums(windows == target_b) == k
        
        winner[ready[match_a]] <- 1
        winner[ready[match_b]] <- 0
      }
    }
    
    unfinished <- which(is.na(winner))
  }
  
  return(cbind(winner, tosses))
}


# Q2: duel

duel <- function(a, b, p=0.5){
  # Description: simulates one game between patterns a and b.
  # Inputs: a and b (patterns), p (chance of heads).
  # Output: c(winner, tosses), with 1 for a and 0 for b.
  # How it works: uses simulate_duels for one game.
  
  result <- simulate_duels(a, b, p, nreps=1)
  return(result[1,])
}


# Q3: pwin

pwin <- function(a, b, p=0.5, nreps=10000){
  # Description: estimates the probability that a beats b.
  # Inputs: a and b (patterns), p (chance of heads),
  #         nreps (number of simulated games).
  # Output: the proportion of games won by a.
  
  results <- simulate_duels(a, b, p, nreps=nreps)
  return(mean(results[, "winner"]))
}


# Q4: ewait

ewait <- function(a, p=0.5){
  # Description: calculates the exact mean wait for pattern a.
  # Inputs: a (pattern), p (chance of heads).
  # Output: the expected number of tosses until a first appears.
  # How it works: adds 1/P(prefix) for every self-overlap length.
  
  k <- length(a)
  
  if(k == 0 || anyNA(a) || any(!(a %in% c(0,1)))){
    stop("Please enter a non-empty 0/1 pattern.")
  }
  
  if(length(p) != 1 || !is.numeric(p) ||
     !is.finite(p) || p <= 0 || p >= 1){
    stop("p must be between 0 and 1.")
  }
  
  wait <- 0
  
  for(j in seq_len(k)){
    first <- a[1:j]
    last <- a[(k-j+1):k]
    
    if(all(first == last)){
      h <- sum(first)
      prob <- p^h * (1-p)^(j-h)
      wait <- wait + 1/prob
    }
  }
  
  return(wait)
}


# Q4: simwait

simwait <- function(a, p=0.5, nreps=3000, block=20){
  # Description: estimates the wait for one pattern by simulation.
  # Inputs: a (pattern), p (heads), nreps (runs), block (batch size).
  # Output: the simulated waiting time for each run.
  # How it works: keeps adding tosses until firsthit finds a.
  
  waits <- integer(nreps)
  
  for(i in seq_len(nreps)){
    x <- integer(0)
    hit <- NA
    
    while(is.na(hit)){
      x <- c(x, rbinom(block, size=1, prob=p))
      hit <- firsthit(x, a)
    }
    
    waits[i] <- hit
  }
  
  return(waits)
}


# Q5: pmat

pmat <- function(k, p=0.5, nreps=10000){
  # Description: estimates win probabilities for all pattern pairs.
  # Inputs: k (pattern length), p (heads), nreps (games per pair).
  # Output: matrix M, where M[i,j] is P(pattern i beats pattern j).
  # How it works: generates every 0/1 pattern in binary order.
  # For k=3 the order is 000, 001, ..., 111. Each pair is
  # simulated once. The diagonal is NA and reverse wins add to 1.
  
  if(length(k) != 1 || !is.numeric(k) ||
     !is.finite(k) || k < 1 || k != floor(k)){
    stop("k must be a positive integer.")
  }
  
  # expand.grid creates all 0/1 combinations for k positions.
  patterns <- as.matrix(expand.grid(rep(list(0:1), k)))
  patterns <- patterns[, k:1, drop=FALSE]
  labels <- apply(patterns, 1, paste, collapse="")
  n <- nrow(patterns)
  
  M <- matrix(NA_real_, nrow=n, ncol=n)
  rownames(M) <- labels
  colnames(M) <- labels
  
  for(i in seq_len(n-1)){
    for(j in (i+1):n){
      prob <- pwin(patterns[i,], patterns[j,], p=p, nreps=nreps)
      M[i,j] <- prob
      M[j,i] <- 1-prob
    }
  }
  
  return(M)
}


# Q5(a): best_reply

best_reply <- function(M){
  # Description: finds the best reply to each pattern.
  # Input: M (matrix where each row pattern beats a column pattern).
  # Output: pattern, its best reply and the estimated win chance.
  # How it works: finds the largest value in each column.
  
  n <- ncol(M)
  result <- data.frame(
    pattern=colnames(M),
    reply=character(n),
    probability=numeric(n)
  )
  
  for(j in seq_len(n)){
    i <- which.max(M[,j])
    result$reply[j] <- rownames(M)[i]
    result$probability[j] <- M[i,j]
  }
  
  result$probability <- round(result$probability, 3)
  return(result)
}


# Q5(b): find_cycle

find_cycle <- function(M, cutoff=0.55){
  # Description: finds the shortest clear non-transitive cycle.
  # Inputs: M (win matrix), cutoff (minimum win probability).
  # Output: the positions in a cycle, or NULL if none is found.
  # How it works: builds paths of length 3, then 4, and so on.
  # A path forms a cycle if its last pattern beats its first.
  
  n <- nrow(M)
  
  if(n < 3){
    return(NULL)
  }
  
  beats <- M > cutoff
  beats[is.na(beats)] <- FALSE
  
  for(size in 3:n){
    for(start in seq_len(n)){
      paths <- list(start)
      
      for(step in seq_len(size-1)){
        new_paths <- list()
        
        for(path in paths){
          current <- tail(path, 1)
          options <- which(beats[current,])
          options <- options[!(options %in% path)]
          
          for(j in options){
            new_paths[[length(new_paths)+1]] <- c(path,j)
          }
        }
        
        paths <- new_paths
      }
      
      for(path in paths){
        last <- tail(path, 1)
        
        if(beats[last,start]){
          return(c(path,start))
        }
      }
    }
  }
  
  return(NULL)
}


# Q5(b): show_cycle

show_cycle <- function(M, cutoff=0.55){
  # Description: prints a cycle and its win probabilities.
  # Inputs: M (win matrix), cutoff (minimum win probability).
  # Output: a labelled cycle and the probability of each win.
  
  cycle <- find_cycle(M, cutoff=cutoff)
  
  if(is.null(cycle)){
    cat("No clear non-transitive cycle found.\n")
    return(invisible(NULL))
  }
  
  labels <- rownames(M)
  cat("Cycle:", paste(labels[cycle], collapse=" -> "), "\n")
  
  for(i in seq_len(length(cycle)-1)){
    a <- cycle[i]
    b <- cycle[i+1]
    cat(labels[a], "beats", labels[b], ":",
        round(M[a,b], 3), "\n")
  }
}


#-------------------------------------------------------------------------------
# Run code and analysis
#-------------------------------------------------------------------------------

# A fixed seed makes the simulation results reproducible.
set.seed(2026)


# Q1: Check firsthit on simple streams

x <- c(1,0,0,1,0,1,1)
a <- c(1,0,1)

cat("Q1: First hit of 101:", firsthit(x,a), "(expected 6)\n")
cat("Q1: Pattern 111 absent:", firsthit(x,c(1,1,1)), "\n")
cat("Q1: Short stream:", firsthit(c(1,0),a), "\n")


# Q2: Simulate one duel

a <- c(1,0,1)
b <- c(0,1,1)
result <- duel(a,b)

cat("\nQ2: Duel between 101 and 011\n")
print(result)


# Q3: Estimate the chance that HHH beats THH

p_fair <- pwin(c(1,1,1), c(0,1,1))
p_bias <- pwin(c(1,1,1), c(0,1,1), p=0.7)

cat("\nQ3: Win probability estimates\n")
cat("P(HHH beats THH), p=0.5:", round(p_fair,3), "\n")
cat("P(HHH beats THH), p=0.7:", round(p_bias,3), "\n")

# The exact probabilities for this pair are 0.125 and 0.343.


# Q4: Compare exact and simulated waiting times

patterns4 <- list(
  c(1,1,1), c(1,1,0), c(1,0,1), c(1,1,0)
)
labels4 <- c("HHH", "HHT", "HTH", "HHT")
probs4 <- c(0.5, 0.5, 0.5, 0.7)

exact4 <- numeric(length(patterns4))
mc4 <- numeric(length(patterns4))

for(i in seq_along(patterns4)){
  exact4[i] <- ewait(patterns4[[i]], p=probs4[i])
  mc4[i] <- mean(simwait(patterns4[[i]], p=probs4[i]))
}

comparison4 <- data.frame(
  pattern=labels4,
  p=probs4,
  exact=round(exact4,3),
  simulated=round(mc4,3)
)

cat("\nQ4: Exact and simulated waiting times\n")
print(comparison4, row.names=FALSE)


# Q5: Construct the fair-coin win matrix

M <- pmat(k=3, p=0.5, nreps=10000)


# Q5(a): Find the best reply to every pattern

cat("\nQ5(a): Best replies, p=0.5\n")
print(best_reply(M), row.names=FALSE)


# Q5(b): Find near-ties and a non-transitive cycle

# An estimate within 0.03 of 0.5 is treated as nearly even.
# upper.tri selects each unordered pair only once.
tie_pos <- which(
  upper.tri(M) & abs(M-0.5) < 0.03,
  arr.ind=TRUE
)

ties <- data.frame(
  pattern_a=rownames(M)[tie_pos[,1]],
  pattern_b=colnames(M)[tie_pos[,2]],
  p_a_wins=round(M[tie_pos],3)
)

cat("\nQ5(b): Approximately even pairs\n")
print(ties, row.names=FALSE)

# The near-ties rule out three clear wins in a fair-coin cycle.
# Try different cycle lengths rather than assuming three patterns.
cat("\nQ5(b): Shortest clear cycle, p=0.5\n")
show_cycle(M, cutoff=0.55)


# Q5(c): Compare individual waits and head-to-head results

patterns3 <- as.matrix(expand.grid(rep(list(0:1),3)))
patterns3 <- patterns3[,3:1,drop=FALSE]
waits3 <- numeric(nrow(patterns3))

for(i in seq_len(nrow(patterns3))){
  waits3[i] <- ewait(patterns3[i,])
}

waiting_rank <- data.frame(
  pattern=rownames(M),
  expected_wait=waits3
)
waiting_rank <- waiting_rank[order(waiting_rank$expected_wait),]

cat("\nQ5(c): Ranking by expected waiting time\n")
print(waiting_rank, row.names=FALSE)

# TTH and THH have the same isolated mean wait of 8 tosses,
# but one is much more likely to win when they compete.
cat("\nQ5(c): TTH against THH\n")
cat("Expected wait of TTH:", ewait(c(0,0,1)), "\n")
cat("Expected wait of THH:", ewait(c(0,1,1)), "\n")
cat("P(TTH beats THH):", round(M["001","011"],3), "\n")


# Q6: Plot duel duration and Monte Carlo variation

# Plot 1: empirical distribution of the tosses until a winner.
games <- simulate_duels(c(1,1,1), c(0,1,1), nreps=10000)
duration <- games[,"tosses"]
freq <- table(duration)
xval <- as.numeric(names(freq))
cumprob <- cumsum(as.numeric(freq))/length(duration)
med <- median(duration)

par(mfrow=c(1,2))

plot(xval, cumprob, type="s", ylim=c(0,1),
     xlab="Number of tosses", ylab="Cumulative probability",
     main="Duel length: HHH vs THH")
lines(c(med,med), c(0,1), col=2, lty=2)
text(med, 0.55, paste("Median =",med), pos=4, col=2)

# Plot 2: 50 separate estimates of P(THH beats HHH).
estimates <- numeric(50)

for(i in seq_len(50)){
  estimates[i] <- pwin(c(0,1,1), c(1,1,1), nreps=2000)
}

reference <- 7/8
limits <- range(c(estimates, reference)) + c(-0.03,0.03)

plot(estimates, pch=19, ylim=limits,
     xlab="Simulation run", ylab="Estimated win probability",
     main="Monte Carlo variation")
lines(c(1,50), c(reference,reference), col=2, lty=2)
text(25, reference+0.02, "Exact value = 0.875", col=2)

par(mfrow=c(1,1))


# Q7: Compare the game when the coin is biased

for(bias in c(0.25,0.7)){
  M_bias <- pmat(k=3, p=bias, nreps=20000)
  
  cat("\nQ7: Best replies for p =", bias, "\n")
  print(best_reply(M_bias), row.names=FALSE)
  
  cat("\nQ7: Non-transitive cycle for p =", bias, "\n")
  show_cycle(M_bias, cutoff=0.55)
}

# Bias changes the best replies. For example, 000 becomes the best
# reply to 001 when p=0.25, rather than 100 for a fair coin.
# When p=0.7, 110 becomes the best reply to 010 rather than 001.
# Both biased cases can have three-pattern cycles, while the fair
# coin needs four patterns for a cycle of clear wins.


#-------------------------------------------------------------------------------
# Conclusion
#-------------------------------------------------------------------------------

# TTH and THH have equal mean waiting times, but TTH usually wins
# their direct duel. The patterns use the same toss stream, so the
# waiting times are dependent. Ranking patterns by their separate
# waiting times can therefore give the wrong head-to-head answer.

