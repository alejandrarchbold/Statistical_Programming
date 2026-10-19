# Statistical Programming
# Practical 2: Sequence duels and perils of independence

# Group 3: Alejandra Campo Archbold (s3023700), Weihang Gu (s2932201), Weichen Zhao (s2929481)

# Contributions:


#1.
firsthit <- function(x, a) {
  k <- length(a)
  n <- length(x)
  if (n < k) return(NA)
  
  for (i in k:n) {
    if (identical(x[(i - k + 1):i], a)) {
      return(i)
    }
  }
  return(NA)
}

find_first_match_mat <- function(mat, a) {
  k <- length(a)
  n_cols <- ncol(mat)
  n_rows <- nrow(mat)
  
  if (n_cols < k) return(rep(NA, n_rows))
  
  match_mat <- matrix(TRUE, nrow = n_rows, ncol = n_cols - k + 1)
  for (offset in 1:k) {
    match_mat <- match_mat & (mat[, offset:(n_cols - k + offset), drop = FALSE] == a[offset])
  }
  
  first_pos <- max.col(match_mat, ties.method = "first")
  has_match <- rowSums(match_mat) > 0
  
  res <- rep(NA, n_rows)
  res[has_match] <- first_pos[has_match] + k - 1
  return(res)
}

# 2.
duel <- function(a, b, p = 0.5) {
  k <- length(a)
  block_size <- 50
  tosses_so_far <- sample(c(0, 1), size = 30, replace = TRUE, prob = c(1 - p, p))
  
  repeat {
    pos_a <- firsthit(tosses_so_far, a)
    pos_b <- firsthit(tosses_so_far, b)
    
    if (!is.na(pos_a) || !is.na(pos_b)) {
      if (is.na(pos_b)) return(c(winner = 1, tosses = pos_a))
      if (is.na(pos_a)) return(c(winner = 0, tosses = pos_b))
      if (pos_a < pos_b) return(c(winner = 1, tosses = pos_a))
      return(c(winner = 0, tosses = pos_b))
    }
    tosses_so_far <- c(tosses_so_far, sample(c(0, 1), size = block_size, replace = TRUE, prob = c(1 - p, p)))
  }
}

# 3.
pwin <- function(a, b, p = 0.5, nreps = 10000) {
  k <- length(a)
  block_size <- 50
  
  stream_mat <- matrix(sample(c(0, 1), size = nreps * 30, replace = TRUE, prob = c(1 - p, p)), 
                       nrow = nreps, ncol = 30)
  
  winners <- rep(NA, nreps)
  active <- rep(TRUE, nreps)
  
  repeat {
    pos_a <- find_first_match_mat(stream_mat[active, , drop = FALSE], a)
    pos_b <- find_first_match_mat(stream_mat[active, , drop = FALSE], b)
    
    finished <- !is.na(pos_a) | !is.na(pos_b)
    
    if (any(finished)) {
      active_indices <- which(active)
      finished_indices <- active_indices[finished]
      
      pa <- pos_a[finished]
      pb <- pos_b[finished]
      
      win_a <- ifelse(is.na(pb), TRUE, ifelse(is.na(pa), FALSE, pa < pb))
      winners[finished_indices] <- ifelse(win_a, 1, 0)
      
      active[finished_indices] <- FALSE
    }
    
    if (!any(active)) break
    
    new_blocks <- matrix(sample(c(0, 1), size = nreps * block_size, replace = TRUE, prob = c(1 - p, p)),
                         nrow = nreps, ncol = block_size)
    stream_mat <- cbind(stream_mat, new_blocks)
  }
  
  return(mean(winners))
}

# 4.
ewait <- function(a, p = 0.5) {
  k <- length(a)
  total_ew <- 0
  
  for (j in 1:k) {
    prefix <- a[1:j]
    suffix <- a[(k - j + 1):k]
    if (identical(prefix, suffix)) {
      hj <- sum(prefix)
      prob_j <- (p^hj) * ((1 - p)^(j - hj))
      total_ew <- total_ew + (1 / prob_j)
    }
  }
  return(total_ew)
}

# 5
pmat <- function(k, p = 0.5, nreps = 10000) {
  pat_df <- expand.grid(rep(list(0:1), k))
  pat_mat <- as.matrix(pat_df[, ncol(pat_df):1]) 
  n_pat <- nrow(pat_mat)
  res_mat <- matrix(NA, nrow = n_pat, ncol = n_pat)
  
  for (i in 1:n_pat) {
    for (j in 1:n_pat) {
      if (i < j) {
        a <- pat_mat[i, ]
        b <- pat_mat[j, ]
        p_val <- pwin(a, b, p = p, nreps = nreps)
        res_mat[i, j] <- p_val
        res_mat[j, i] <- 1 - p_val
      }
    }
  }
  
  rownames(res_mat) <- apply(pat_mat, 1, paste0, collapse = "")
  colnames(res_mat) <- apply(pat_mat, 1, paste0, collapse = "")
  return(res_mat)
}

set.seed(42)
k_val <- 3
p_val <- 0.5

W_matrix <- pmat(k = 3, p = 0.5, nreps = 3000)
patterns_names <- rownames(W_matrix)
pat_df <- expand.grid(rep(list(0:1), k_val))
pat_mat <- as.matrix(pat_df[, ncol(pat_df):1])

best_replies <- data.frame(Pattern = character(), BestReply = character(), WinProb = numeric(), stringsAsFactors = FALSE)
for (i in 1:nrow(W_matrix)) {
  row_vals <- W_matrix[i, ]
  best_idx <- which.max(row_vals)
  best_replies <- rbind(best_replies, data.frame(
    Pattern = patterns_names[i],
    BestReply = patterns_names[best_idx],
    WinProb = row_vals[best_idx]
  ))
}
print(best_replies, row.names = FALSE)

evenly_matched <- data.frame(PatternA = character(), PatternB = character(), ProbAB = numeric(), stringsAsFactors = FALSE)
for (i in 1:(nrow(W_matrix)-1)) {
  for (j in (i+1):nrow(W_matrix)) {
    p_ab <- W_matrix[i, j]
    if (!is.na(p_ab) && abs(p_ab - 0.5) < 0.05) { 
      evenly_matched <- rbind(evenly_matched, data.frame(
        PatternA = patterns_names[i],
        PatternB = patterns_names[j],
        ProbAB = p_ab
      ))
    }
  }
}
print(evenly_matched, row.names = FALSE)

find_cycle <- function(mat) {
  n <- nrow(mat)
  names_p <- rownames(mat)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j || mat[i,j] <= 0.5) next
      for (l in 1:n) {
        if (l == i || l == j || mat[j,l] <= 0.5) next
        if (mat[l,i] > 0.5) { 
          return(c(names_p[i], names_p[j], names_p[l]))
        }
      }
    }
  }
  return(NULL)
}

cycle_found <- find_cycle(W_matrix)
if (!is.null(cycle_found)) {
  cat(sprintf("Cycle: %s -> %s -> %s -> %s\n", 
              cycle_found[1], cycle_found[2], cycle_found[3], cycle_found[1]))
}

ew_vals <- sapply(1:nrow(pat_mat), function(i) ewait(pat_mat[i, ], p = p_val))
ew_df <- data.frame(Pattern = patterns_names, EWait = ew_vals)
ew_df <- ew_df[order(ew_df$EWait), ]
print(ew_df, row.names = FALSE)

a_pat <- c(1, 1, 1)
b_pat <- c(0, 1, 1)
duel_tosses <- replicate(2000, duel(a_pat, b_pat, p = 0.5)["tosses"])
dens <- density(duel_tosses)
plot(dens, main = "Estimated Distribution of Duel Length (HHH vs THH)",
     xlab = "Number of Tosses", ylab = "Density", col = "darkblue", lwd = 2)
med_tosses <- median(duel_tosses)
lines(c(med_tosses, med_tosses), c(0, max(dens$y)), col = "red", lty = 2, lwd = 2)
text(x = med_tosses + 3, y = max(dens$y) * 0.8, 
     labels = paste("Median:", med_tosses), col = "red", pos = 4)

sim_runs <- replicate(50, pwin(c(1,1,1), c(0,1,1), p = 0.5, nreps = 2000))
mean_val <- mean(sim_runs)
plot(sim_runs, type = "p", pch = 16, col = "royalblue",
     main = "Variability of P(THH beats HHH) Across 50 Independent Runs",
     xlab = "Simulation Run Index", ylab = "Estimated Win Probability",
     ylim = c(mean_val - 0.05, mean_val + 0.05))
lines(c(1, 50), c(mean_val, mean_val), col = "darkred", lty = 2, lwd = 2)
text(x = 25, y = mean_val + 0.015, 
     labels = paste("Mean Value ~", round(mean_val, 3)), col = "darkred")

W_mat_biased_low <- pmat(k = 3, p = 0.25, nreps = 1000)
W_mat_biased_high <- pmat(k = 3, p = 0.7, nreps = 1000)