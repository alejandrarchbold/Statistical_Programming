# Practical 1: Music simulator
# Group: Alejandra Campo Archbold (s3023700), Weihang Gu(s2932201), Weichen Zhao(s2929481)
# Contributions:
# Alejandra Campo Archbold: preprocess, tokens, final comments and tidy-up (33.3%)
# Weihang Gu: make_matrix and next_note (33.3%)
# Weichen Zhao: generate_music, weight comparison (w1, w2, w3), testing (33.3%)
# All three reviewed and tested the whole script together.

#-------------------------------------------------------------------------------
# We learn from Palestrina's music and make new music with a Markov model.
# We clean the text, give each note a number, store all short note runs
# in a matrix, and pick each new note from what followed the last notes.
#-------------------------------------------------------------------------------

# setwd("put/your/local/repo/location/here")
music_raw <- readLines("palestrina.txt")   # one line of the file per element
head(music_raw, 20)                        # look at the first lines

# preprocess: raw lines in, vector of notes with "|" and "||" out.
preprocess <- function(music) {
  music <- music[-(1:3)]                         # remove composer, source and blank line
  music <- music[substr(music, 1, 1) != "#"]     # remove info lines (only 1st char, because F#4 has # too)
  is_bar <- music != ""                          # TRUE for lines with notes
  last_bar <- c(music[-1] == "", TRUE)[is_bar]   # next line is blank, so the piece ends here
  music[is_bar] <- paste(music[is_bar], ifelse(last_bar, "||", "|"))  # add "||" at the end of a piece, "|" otherwise
  unlist(strsplit(music[is_bar], " "))           # drop blank lines and split bars into notes
}

music_clean <- preprocess(music_raw)
length(music_clean)      # over 1 million notes
tail(music_clean, 10)    # should end with "G3b" "||"


notes  <- unique(music_clean)          # list of different notes
tokens <- match(music_clean, notes)    # each note becomes its position in notes


# make_matrix: each row has mlag notes and the note that came after them.
make_matrix <- function(tokens, end, mlag = 4) {
  n <- length(tokens)                                      # number of notes
  M <- matrix(NA_integer_, n - mlag, mlag + 1)             # empty matrix
  for (j in 1:(mlag + 1)) M[, j] <- tokens[j:(n - mlag + j - 1)]   # each column moves one note forward
  across <- rowSums(M[, 1:mlag, drop = FALSE] == end) > 0  # rows with "||" before the last column
  M[!across, , drop = FALSE]                               # remove them, they join two pieces
}

M <- make_matrix(tokens = tokens, end = which(notes == "||"), mlag = 4)


# next_note: picks the next token after key. w[1] is for the longest key,
# w[mlag] is for the last note only.
next_note <- function(key, M, w = rep(1, ncol(M) - 1)) {
  mlag <- ncol(M) - 1                           # how many notes we look back
  key <- tail(key, mlag)                        # use only the last mlag notes
  k <- length(key)
  u <- p <- c()                                 # possible next notes and their probabilities
  for (i in 1:k) {                              # longest key first, then shorter
    sub_key <- key[i:k]
    mc <- mlag - length(sub_key) + 1            # first column to compare (also the weight to use)
    if (w[mc] == 0) next                        # weight 0, so no need to search
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))  # count differences in each row
    hit <- which(ii == 0)                       # rows with no differences match
    if (length(hit) > 0) {
      u <- c(u, M[hit, mlag + 1])               # notes that came next
      p <- c(p, rep(w[mc] / length(hit), length(hit)))  # share the weight between them
    }
  }
  if (length(u) == 0) return(sample(M[, mlag + 1], 1))  # no match: take any note
  if (length(u) == 1) return(u)                         # one match: no need to sample
  sample(u, 1, prob = p)                                # pick one using the weights
}


# generate_music: makes a piece from start and returns it as text.
generate_music <- function(start, M, w = rep(1, ncol(M) - 1), notes) {
  end_tok <- which(notes == "||")     # token for end of piece
  bar_tok <- which(notes == "|")      # token for end of bar
  piece <- start
  bars <- 0                           # bars made so far
  repeat {
    nxt <- next_note(piece, M, w)     # new note from all notes so far
    piece <- c(piece, nxt)            # add it to the piece
    if (nxt == bar_tok) bars <- bars + 1
    if (nxt == end_tok || bars == 8 || length(piece) >= 200) break  # stop rules
  }
  paste(notes[piece], collapse = " ") # tokens back to note names
}


# Compare the weights, all starting from F#5w.
set.seed(2026)
start <- which(notes == "F#5w")

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
