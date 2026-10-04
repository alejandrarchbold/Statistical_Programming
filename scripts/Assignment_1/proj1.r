# Practical 1: Music simulator

# Group members: Alejandra Campo Archbold, Weihang Gu, Weichen Zhao s2929481

# Contributions: <who did what, and roughly what % each>

# Overview: a simple Markov model that learns from Palestrina's music
# and generates new pieces. Steps: read and clean the data, turn notes
# into tokens, build a matrix of note sequences, predict notes one by
# one, and compare different weights.


# -------------------------------------------------------------------
# 0. Read the data  (Task 3-4)
# -------------------------------------------------------------------
# Task: read palestrina.txt into R and look at the first lines.
# Load the raw text file (one line per bar or info line).

music_raw <- readLines("palestrina.txt")

head(music_raw, 20)


# -------------------------------------------------------------------
# 1. Clean the data  (Task 4a-4e)
# -------------------------------------------------------------------
# Task: write preprocess() to remove the header and "#" lines, add
# " |" after each bar and " ||" after each piece, and split into notes.
# Remove headers and info lines, mark bar ends ("|") and piece ends
# ("||"), then split everything into single notes.

# Input:  music, raw lines of the file.
# Output: vector of notes, rests, "|" and "||".
preprocess <- function(music) {
  
  # Drop composer, source and blank line.
  music <- music[-(1:3)]
  
  # Drop info lines. Check only the first character: "#" also means sharp.
  music <- music[substr(music, 1, 1) != "#"]
  
  # A blank next line means the piece ends here.
  next_is_blank <- c(music[-1] == "", TRUE)
  is_bar <- music != ""
  
  # "||" for the last bar of a piece, "|" for the rest.
  ending <- ifelse(next_is_blank[is_bar], "||", "|")
  music[is_bar] <- paste(music[is_bar], ending)
  
  # Remove blank lines and split bars into notes.
  unlist(strsplit(music[is_bar], " "))
}

music_clean <- preprocess(music_raw)

# Check: over 1 million notes, ending in "G3b ||".
length(music_clean)
tail(music_clean, 10)


# -------------------------------------------------------------------
# 2. Turn notes into tokens  (Task 5)
# -------------------------------------------------------------------
# Task: use unique() to get the notes and match() to get the tokens.
# Each note becomes a number: its position in the list of unique notes.

notes  <- unique(music_clean)
tokens <- match(music_clean, notes)


# -------------------------------------------------------------------
# 3. Build the sequence matrix  (Task 6)
# -------------------------------------------------------------------
# Task: write make_matrix() to build an (n - mlag) x (mlag + 1) matrix
# of shifted tokens, and remove rows with "||" in the first mlag columns.
# Each row of M holds mlag notes plus the note that followed them.

# Input:  tokens, all tokens; end, token of "||"; mlag, notes of context.
# Output: matrix with mlag + 1 columns.
make_matrix <- function(tokens, end, mlag = 4) {
  
  n <- length(tokens)
  M <- matrix(NA_integer_, nrow = n - mlag, ncol = mlag + 1)
  
  # Column j is the token vector shifted by j - 1.
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # Drop rows that run across two pieces ("||" before the last column).
  crosses_pieces <- rowSums(M[, 1:mlag, drop = FALSE] == end) > 0
  M[!crosses_pieces, , drop = FALSE]
}

M <- make_matrix(tokens = tokens, end = which(notes == "||"), mlag = 4)


# -------------------------------------------------------------------
# 4. Choose the next note  (Task 7)
# -------------------------------------------------------------------
# Task: write next_note() to sample the next token from key, using
# the mixture model with weights w and at most one call to sample().
# Find rows of M matching the last 4, 3, 2 and 1 notes, weight each
# group by w, and sample one of the notes that came next.

# Input:  key, notes so far; M, sequence matrix; w, weight per key length
#         (w[1] = longest key, w[mlag] = last note only).
# Output: one token.
next_note <- function(key, M, w = rep(1, ncol(M) - 1)) {
  
  mlag <- ncol(M) - 1
  
  # Keep only the last mlag notes.
  key <- tail(key, mlag)
  key_len <- length(key)
  
  candidates <- c()
  probs <- c()
  
  # Try the full key first, then shorter ones.
  for (i in 1:key_len) {
    
    sub_key <- key[i:key_len]
    
    # Columns to compare, also the index of the weight.
    mc <- mlag - length(sub_key) + 1
    
    # Skip key lengths with no weight.
    if (w[mc] == 0) next
    
    # Rows with 0 differences are matches.
    n_diff  <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))
    matches <- which(n_diff == 0)
    
    # Save the next notes; share the weight equally between matches.
    if (length(matches) > 0) {
      candidates <- c(candidates, M[matches, mlag + 1])
      probs <- c(probs, rep(w[mc] / length(matches), length(matches)))
    }
  }
  
  # No match: pick any note from the data.
  if (length(candidates) == 0) return(sample(M[, mlag + 1], 1))
  
  # One match: no need to sample.
  if (length(candidates) == 1) return(candidates)
  
  sample(candidates, 1, prob = probs)
}


# -------------------------------------------------------------------
# 5. Generate a piece  (Task 8)
# -------------------------------------------------------------------
# Task: write generate_music() to call next_note() repeatedly and
# return the piece as one string of notes.
# Add notes one at a time until "||", 8 bars, or 200 notes.

# Input:  start, first token; M, sequence matrix; w, weights;
#         notes, note names.
# Output: one string with the piece.
generate_music <- function(start, M, w = rep(1, ncol(M) - 1), notes) {
  
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  
  piece <- start
  n_bars <- 0
  
  repeat {
    new <- next_note(key = piece, M = M, w = w)
    piece <- c(piece, new)
    
    if (new == bar_token) n_bars <- n_bars + 1
    
    # Stopping rules.
    if (new == end_token || n_bars >= 8 || length(piece) >= 200) break
  }
  
  # Tokens back to note names.
  paste(notes[piece], collapse = " ")
}


# -------------------------------------------------------------------
# 6. Compare weights  (Task 8)
# -------------------------------------------------------------------
# Task: test w1 = (1,1,1,1), w2 = (1,0,0,0), w3 = (0,0,0,1) and
# comment on the lag used and how it changes the music.
# Generate one piece with each weight vector and compare them.

set.seed(2026)
start_note <- which(notes == "F#5w")

# w1: mixes all key lengths. Most Palestrina-like: smooth melody,
# bars of mostly sensible length.
w1 <- c(1, 1, 1, 1)
music_w1 <- generate_music(start = start_note, M = M, w = w1, notes = notes)
music_w1

# w2: only 4-note keys. Short early keys get no weight, so notes are
# random, and random keys rarely match later. Result: big jumps and
# empty bars.
w2 <- c(1, 0, 0, 0)
music_w2 <- generate_music(start = start_note, M = M, w = w2, notes = notes)
music_w2

# w3: only the last note. Each step looks normal, but no memory of the
# bar, so bar lengths are uneven and the melody wanders.
w3 <- c(0, 0, 0, 1)
music_w3 <- generate_music(start = start_note, M = M, w = w3, notes = notes)
music_w3