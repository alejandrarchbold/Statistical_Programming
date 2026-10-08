# Practical 1: Music simulator

# Group 3: Alejandra Campo Archbold (s3023700), Weihang Gu(s2932201), Weichen Zhao(s2929481)

# Contributions:

# - Alejandra Campo Archbold: preprocess, tokens, final comments and tidy-up (33.3%)
# - Weihang Gu: make_matrix and next_note (33.3%)
# - Weichen Zhao: generate_music, weight comparison (w1, w2, w3), testing (33.3%)
# All three reviewed and tested the whole script together.

# GitHub URL project: https://github.com/alejandrarchbold/Statistical_Programming.git

#-------------------------------------------------------------------------------
# Overview:

# This script simulates new music using Markov model in the style of Giovanni Palestrina.
# We clean the text (preprocess), give each note a number (tokens), store all short note runs
# in a matrix (make_matrix), and pick each new note from what followed the last notes (next_note).
# In the process, we write the function generate_music that call next_note repeatedly until 
# the piece ends, eight bars are made or 200 notes are reached.
# Finally we compare three choices of mixture weights.
#-------------------------------------------------------------------------------

# setwd("put/your/local/repo/location/here")
music_raw <- readLines("palestrina.txt")   # training data
head(music_raw, 15)                        # look at the first 15 lines

# - preprocess:

preprocess <- function(music) {
  # Description: cleans the file and marks bar ends with "|" and piece ends with "||".
  # Input: music, the lines of the file.
  # Output: character vector of notes, rests, "|" (end of bar) and "||" (end of piece).
  # How it works: a bar is the last of its piece when the next line is blank.
  
  music <- music[-(1:3)]                         # drop the 3 header lines: remove composer, source and blank line
  music <- music[substr(music, 1, 1) != "#"]     # remove info lines, notes with #
  is_bar <- music != ""                          # lines that have notes (TRUE)
  last_bar <- c(music[-1] == "", TRUE)[is_bar]   # next line is blank, means last bar of a piece
  music[is_bar] <- paste(music[is_bar], ifelse(last_bar, "||", "|"))  # add || at the end of a piece, | otherwise
  unlist(strsplit(music[is_bar], " "))           # split bars into single notes
}

music_clean <- preprocess(music_raw)
length(music_clean) # 
tail(music_clean, 10) # 


# - tokens:
notes  <- unique(music_clean)          # list of uniques notes
tokens <- match(music_clean, notes)    # each note as its number in notes


# - make_matrix:

make_matrix <- function(tokens, end, mlag = 4) {
  # Description: stores every run of mlag+1 notes in a row, the data for the model.
  # Inputs: tokens (all notes as numbers), end (number of "||"), mlag (notes to look back).
  # Output: matrix M, with one row per run of notes.
  # How it works: builds M column by column and removes the rows that join two pieces.
  # The matrix M has n−mlag rows. If we use an i-loop from 1 to n−mlag, it will be very slow.
  # The logic of this matrix is that each column is the original tokens sequence shifted back.
  
  n <- length(tokens)                                      # number of notes
  M <- matrix(NA_integer_, n - mlag, mlag + 1)             # empty matrix
  for (j in 1:(mlag + 1)) M[, j] <- tokens[j:(n - mlag + j - 1)]   # each column starts one note later
  across <- rowSums(M[, 1:mlag, drop = FALSE] == end) > 0  # rows with "||" before the last column
  M[!across, , drop = FALSE]                               # remove them, they join or mix two pieces
}

M <- make_matrix(tokens = tokens, end = which(notes == "||"), mlag = 4)

# - next_note:

next_note <- function(key, M, w = rep(1, ncol(M) - 1)) {
  # Description: chooses the next note from the notes that followed key in the data.
  # Inputs: key (notes so far), M (from make_matrix), w (weight of each key length).
  # Output: one note, as a number.
  # How it works: finds key and its shorter endings in M, then samples one note using w.
  # picks the next token after key. w[1] is for the longest key,
  # w[mlag] is for the last note only.
  
  mlag <- ncol(M) - 1                           # notes to look back
  key <- tail(key, mlag)                        # keep only the latest notes
  k <- length(key)
  u <- p <- c()                                 # possible next notes and their probabilities
  for (i in 1:k) {                              # longest key first, then shorter
    sub_key <- key[i:k]
    mc <- mlag - length(sub_key) + 1            # columns to compare and weight to use
    if (w[mc] == 0) next                        # weight 0, skip
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))  # count differences in each row
    hit <- which(ii == 0)                       # rows that match exactly
    if (length(hit) > 0) {
      u <- c(u, M[hit, mlag + 1])               # notes that came next
      p <- c(p, rep(w[mc] / length(hit), length(hit)))  # split the weight among them
    }
  }
  if (length(u) == 0) return(sample(M[, mlag + 1], 1))  # no match: take any note
  if (length(u) == 1) return(u)                         # one match: no need to sample
  sample(u, 1, prob = p)                                # pick one by weight
}


# - generate_music:

generate_music <- function(start, M, w = rep(1, ncol(M) - 1), notes) {
  # Description: writes a new piece of music, one note at a time.
  # Inputs: start (first note), M (from make_matrix), w (weights), notes (note names).
  # Output: the piece as one text, notes separated by spaces.
  # How it works: calls next_note until "||", 8 bars or 200 notes.
  # makes a piece from start and returns it as text.

  end_tok <- which(notes == "||")     # number for ||
  bar_tok <- which(notes == "|")      # number for |
  piece <- start
  bars <- 0                           # bars made so far
  repeat {
    nxt <- next_note(piece, M, w)     # new note
    piece <- c(piece, nxt)            # add it to the piece
    if (nxt == bar_tok) bars <- bars + 1
    if (nxt == end_tok || bars == 8 || length(piece) >= 200) break  # stop rules
  }
  paste(notes[piece], collapse = " ") # tokens back to note names
}


# - compare the weights:

set.seed(2026) #useful for same random numbers every run.
start <- which(notes == "F#5w") # starting from F#5w.

# w1: uses all key lengths. Long keys keep the style and the bar length,
# short keys help when there are no long matches. Sounds most like the data.
sim_w1 <- generate_music(start, M, w = c(1, 1, 1, 1), notes = notes)
sim_w1

# w2: only keys of 4 notes. The first notes have no 4-note key, so they are
# random, and random keys almost never match. Big jumps and empty bars.
sim_w2 <- generate_music(start, M, w = c(1, 0, 0, 0), notes = notes)
sim_w2

# w3: only the last note. Each step looks fine, but the model does not
# remember the bar, so bar lengths change a lot and the melody wanders.
sim_w3 <- generate_music(start, M, w = c(0, 0, 0, 1), notes = notes)
sim_w3
