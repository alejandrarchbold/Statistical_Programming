# Group members:3
# Names: Alejandra, Weichen, Hang


# 1. Read the data

music_raw <- readLines("palestrina.txt")
head(music_raw, 20) # Check


# 2. Clean the music data

# This function removes the extra information in the text file and turns the music into a single vector of notes and bar-line symbols.
#
# music: raw character vector read from palestrina.txt
#
# The blank lines are useful here because they show where one piece finishes and the next one begins.

preprocess <- function(music) {
  
  # Remove the composer, source and first blank line.
  music <- music[-c(1:3)]
  
  # Remove information lines such as the title, time signature and key.
  # We only remove lines which start with "#", so sharp notes such as C#5h are kept.
  music <- music[substr(music, 1, 1) != "#"]
  
  # Add a normal bar-line after each measure. If the next line is blank,
  # the current measure is the end of a piece, so use a double bar-line.
  for (i in 1:(length(music) - 1)) {
    
    if (music[i] != "") {
      
      if (music[i + 1] == "") {
        music[i] <- paste0(music[i], " ||")
      } else {
        music[i] <- paste0(music[i], " |")
      }
    }
  }
  
  # Blank lines are no longer needed once the piece endings are marked.
  music <- music[music != ""]
  
  # Split each line into separate notes, rests and bar-line symbols.
  music <- unlist(strsplit(music, " "))
  
  return(music)
}


music_clean <- preprocess(music_raw)


# -------------------------------------------------------------------
# 3. Convert notes to tokens
# -------------------------------------------------------------------

# notes contains each unique note or symbol in the cleaned data.
# tokens gives the position of each note in this vector.

notes <- unique(music_clean)

tokens <- match(music_clean, notes)


# -------------------------------------------------------------------
# 4. Build the sequence matrix
# -------------------------------------------------------------------

# This function creates the matrix used by the Markov model.
#
# With mlag = 4, each row contains five adjacent tokens:
# four previous tokens followed by the next token.
#
# tokens: integer vector representing the music
# end: token number for "||"
# mlag: maximum number of previous notes used by the model

make_matrix <- function(tokens, end, mlag = 4) {
  
  n <- length(tokens)
  
  # Create an empty matrix with one extra column for the next token.
  M <- matrix(NA,
              nrow = n - mlag,
              ncol = mlag + 1)
  
  # Shift the token vector by one position for each new column.
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # Remove rows where "||" occurs inside the first mlag columns.
  # These rows would otherwise contain a sequence crossing two pieces.
  bad_rows <- rowSums(
    M[, 1:mlag, drop = FALSE] == end
  ) > 0
  
  M <- M[!bad_rows, , drop = FALSE]
  
  return(M)
}


M <- make_matrix(
  tokens = tokens,
  end = which(notes == "||"),
  mlag = 4
)


# -------------------------------------------------------------------
# 5. Predict the next token
# -------------------------------------------------------------------

# next_note() chooses the next token using matches found in M.
#
# The model uses different amounts of recent history. For mlag = 4,
# it considers the last 4, 3, 2 and 1 tokens. The vector w controls how much weight is given to each of these sequence lengths.

next_note <- function(key, M, w = rep(1, ncol(M) - 1)) {
  
  mlag <- ncol(M) - 1
  
  # Only the most recent tokens are needed if key is longer than mlag.
  if (length(key) > mlag) {
    key <- tail(key, mlag)
  }
  
  next_tokens <- c()
  probs <- c()
  
  key_length <- length(key)
  
  # Start with the longest available key and then shorten it one step at a time.
  for (i in 1:key_length) {
    
    current_key <- key[i:key_length]
    
    # Find which columns of M should be compared with the current key.
    mc <- mlag - length(current_key) + 1
    
    # A value of zero means that the row matches the key exactly.
    ii <- colSums(
      !(t(M[, mc:mlag, drop = FALSE]) == current_key)
    )
    
    matches <- which(ii == 0)
    
    if (length(matches) > 0) {
      
      # The last column contains the tokens which followed each match.
      found_tokens <- M[matches, mlag + 1]
      
      # Pick the corresponding weight for this sequence length.
      weight_pos <- mlag - length(current_key) + 1
      
      # Split this weight equally across all matches.
      found_probs <- rep(
        w[weight_pos] / length(matches),
        length(matches)
      )
      
      next_tokens <- c(next_tokens, found_tokens)
      probs <- c(probs, found_probs)
    }
  }
  
  # If there are no usable matches, sample from all tokens that appear in the final column of M.
  if (length(next_tokens) == 0 || sum(probs) == 0) {
    
    result <- sample(M[, mlag + 1], 1)
    
  } else if (length(next_tokens) == 1) {
    
    # No sampling is needed if there is only one possible token.
    result <- next_tokens[1]
    
  } else {
    
    # Otherwise choose one token using the calculated weights.
    result <- sample(next_tokens, 1, prob = probs)
  }
  
  return(result)
}


# -------------------------------------------------------------------
# 6. Generate a piece of music
# -------------------------------------------------------------------

# This function repeatedly calls next_note() and adds each new token to the generated sequence.
# Generation stops if:
# - "||" is produced,
# - eight bar-lines have been produced, or
# - the sequence reaches 200 tokens.

generate_music <- function(start, M,
                           w = rep(1, ncol(M) - 1),
                           notes) {
  
  generated <- start
  
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  
  while (TRUE) {
    
    new_note <- next_note(
      key = generated,
      M = M,
      w = w
    )
    
    generated <- c(generated, new_note)
    
    if (new_note == end_token) {
      break
    }
    
    if (sum(generated == bar_token) >= 8) {
      break
    }
    
    if (length(generated) >= 200) {
      break
    }
  }
  
  # Convert the token numbers back into note names and combine them into one character string.
  result <- paste(
    notes[generated],
    collapse = " "
  )
  
  return(result)
}


# -------------------------------------------------------------------
# 7. Compare different lag weights
# -------------------------------------------------------------------

# w1 uses all four sequence lengths. This combines the extra context from longer sequences with the larger number of matches usually available from shorter sequences.

w1 <- c(1, 1, 1, 1)


# w2 only uses the longest four-token sequence. The prediction therefore depends on more recent context, but exact matches may be less common.

w2 <- c(1, 0, 0, 0)


# w3 only uses the most recent token. This gives more possible matches, but the model has less information about the notes before it.

w3 <- c(0, 0, 0, 1)


music_w1 <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = w1
)

music_w2 <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = w2
)

music_w3 <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = w3
)


# Print the three examples for comparison. The results will change between runs because the model includes random sampling.

music_w1
music_w2
music_w3

# Simple checks for the generated music

check_music <- function(x) {
  
  x <- unlist(strsplit(x, " "))
  
  cat("Number of tokens:", length(x), "\n")
  cat("Number of bars:", sum(x == "|"), "\n")
  cat("Last symbol:", tail(x, 1), "\n")
  cat("Within 200 tokens:", length(x) <= 200, "\n")
  cat("Within 8 bars:", sum(x == "|") <= 8, "\n")
}

check_music(music_w1)
check_music(music_w2)
check_music(music_w3)