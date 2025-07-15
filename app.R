library(shiny)
library(here) # for tidy and stable file paths
library(tidyverse) # for tidy data
library(quanteda.textstats) # for syllable counting (more accurate than sylcount)
library(hunspell) # for spell check (is this word in the dictionary)


# loading dhakal sentences
dhakal_typs <- read_tsv(here("data", "dhakal_typabilities_post_both_filters.txt")) %>%
 select(sentence, typability_z) %>%
  rename(text = sentence,
         typability = typability_z)

# function to standardise symbols
standardiseSymbols <- function(sentencesDF) {

  sentencesDF <- as.data.frame(sentencesDF)

  sentencesDF <- sentencesDF %>%
    mutate(across(1, ~ .x %>%
                    str_replace_all("\u00A0|\u202F", " ") %>%           # non-breaking & narrow no-break spaces -> space
                    str_replace_all("[\u2013\u2014\u2212]", "-") %>%    # en/em dash, minus sign -> hyphen
                    str_replace_all("[\u2018\u2019\u02BC]", "'") %>%    # curly apostrophes -> straight
                    str_replace_all("[\u201C\u201D\u00AB\u00BB]", "\"") %>% # curly and angled quotes -> straight
                    str_replace_all("\u2026", "...") %>%                # ellipsis -> ...
                    str_replace_all("\u00D7", "x")                      # multiplication sign -> 'x'
    ))
}

# function to calculate the predictor variables
calculatePredictorVariables <- function(sentencesDF) {

  # wants only a 1d df with the first column as sentences

  options(dplyr.summarise.inform = FALSE)

  sentencesDF <- as.data.frame(sentencesDF) %>%
    rename(sentence = 1) %>%
    mutate(sentNum = row_number()) %>%
    select(sentNum, sentence)

  # load in the helper data, e.g. bigram frequencies, hand categorisations

  # load bigram frequencies
  load(here("data", "predictor_calculation_aids", "GutenbergFrequencyTables.RData"))
  rm(AllTrigramsMLE, AllLetters)

  ## calc percentile ranks of bigrams
  AllBigramsMLE <- AllBigramsMLE %>%
    mutate(rank = rank(-Frequency, ties = "min"),
           # inflection point after first 15 top frequency bigrams in terms of their frequency (and frequency percent) - these bigrams account for 28.09% of all bigrams in English (in the corpus that Crump used)
           freqPerc = (Frequency / sum(Frequency)) * 100) %>%
    arrange(rank) %>%
    mutate(freqPercCumu = cumsum(freqPerc),
           top15bigram = if_else(rank <= 15, 1, 0))


  # read in the keyboard info
  keyboardInfo <- read_tsv(here("data", "predictor_calculation_aids", "keyboardInfo.txt"),
                           trim_ws = FALSE, quote = "", show_col_types = FALSE) %>%
    # for easier coding later, change the shifts to one-character symbols which are not keyboard characters in themselves
    mutate(character = case_when(character == "RSHIFT" ~ "ℛ",
                                 character == "LSHIFT" ~ "ℒ",
                                 TRUE ~ character))




  print("Calculating the proportion of each character type (e.g. uppercase, numbers) for each sentence.")

  # perform the simple calculations
  data_1_simpleCalcs <- sentencesDF %>%
    # for each unique sentence, calculate these things
    mutate(numActualWords = str_count(sentence, "\\S+"),
           numChars = nchar(sentence),
           meanWordLength = numChars / numActualWords,
           # what proportion the average word takes up of the sentence
           meanWordLengthPropSent = 1 / numActualWords,
           uppercase = str_count(sentence, "[A-Z]"),
           lowercase = str_count(sentence, "[a-z]"),
           letters = uppercase + lowercase,
           numbers = str_count(sentence, "[0-9]"),
           spaces = str_count(sentence, " "),
           symbols = numChars - (letters + numbers + spaces)) %>%
    # convert these counts into proportions
    mutate(across(c(uppercase, lowercase, letters, numbers,
                    spaces, symbols),
                  ~ .x / numChars,
                  .names = "{.col}Prop"),
           lowercasePropNonSpace = lowercase / (numChars - spaces),
           uppercasePropNonSpace = uppercase / (numChars - spaces),
           lettersPropNonSpace = letters / (numChars - spaces),
           symbolsPropNonSpace = symbols / (numChars - spaces),
           numbersPropNonSpace = numbers / (numChars - spaces),
           # create new columns converting the sentence to lowercase then removing the spaces (for some later calculations)
           sent_lower = tolower(sentence),
           sent_lower_no_space = str_replace_all(sent_lower, fixed(" "), ""))

  print("---Done.")


  print("Calculating the frequency that each bigram appears in the English language and averaging for each sentence.")

  data_2b_bigramFreqs <- data_1_simpleCalcs %>%
    # keep only these columns
    select(sentNum, sentence, sent_lower_no_space) %>%
    # separate out the  bigrams in each word
    mutate(bigram = lapply(sent_lower_no_space, str_split_bigrams)) %>%
    unnest(bigram) %>%
    # give each bigram a bigram number
    group_by(sentNum) %>%
    mutate(bigramNum = row_number()) %>%
    # add on the frequencies
    left_join(AllBigramsMLE,  by = c("bigram" = "Bigrams")) %>%
    mutate(n_bigrams = n(),
           n_unmatched_bigrams = sum(is.na(Frequency)),
           perc_unmatched_bigrams = n_unmatched_bigrams / n_bigrams,
    ) %>%
    ungroup()

  # average the bigram frequencies
  data_2_aveBigramFreqs <- data_2b_bigramFreqs %>%
    group_by(sentNum, sentence) %>%
    # summarise to collapse all individual bigram rows
    summarise(biFreqMean = mean(Frequency, na.rm = TRUE),
              propBiTop15 = mean(top15bigram, na.rm = TRUE)) %>%
    ungroup() %>%
    # join with an earlier dataset to keep the wpm info
    left_join(data_1_simpleCalcs, by = c("sentNum", "sentence"))

  print("---Done.")


  print("Calculating the proportion of bigrams within each hand category, e.g. hand alternation for each sentence.")

  # add hands
  data_3b_bigramFreqsAndHands <- data_2b_bigramFreqs %>%
    # the bigram frequency table used only has letters, so numbers and symbols will be NA in these 2 columns - this makes sure they are filled in for all character types (Pred = predecessor, first bigram character; Succ = successor, second bigram character)
    mutate(Pred = substr(bigram, 1, 1),
           Succ = substr(bigram, 2, 2)) %>%
    # separate the bigram characters onto different lines
    pivot_longer(cols = c("Pred", "Succ"),
                 names_to = "characterNum",
                 values_to = "character") %>%
    mutate(characterNum = recode(characterNum, "Pred" = 1, "Succ" = 2)) %>%
    # join with the hand categorisations
    left_join(select(keyboardInfo, character,
                     standardHand, standardFinger),
              by = "character") %>%
    # determine same/diff hand
    group_by(sentNum, bigramNum) %>%
    mutate(bigramHandCateg = case_when(
      # same character
      n_distinct(character) == 1 ~ "charRepetition",
      # same finger, different character
      (n_distinct(standardFinger) == 1 &&
         n_distinct(character) == 2) ~ "fingerRepetition",
      # same hand, different finger (and character)
      (n_distinct(standardHand) == 1 &&
         n_distinct(standardFinger) == 2) ~ "handRepetition",
      # different hand (and finger and character)
      n_distinct(standardHand) == 2 ~ "handAlternation"),
      # in case I want to add bigram hand as a predictor
      bigramHand = case_when(
        # same hand / finger / key
        n_distinct(standardHand) == 1 ~
          str_c(standardHand, "HandBigrams"),
        # both hands
        n_distinct(standardHand) == 2 ~ "bothHandBigrams")) %>%
    ungroup() %>%
    group_by(sentence) %>%
    ungroup() %>%
    # get rid of some columns for now
    select(-c(characterNum, character, standardHand, standardFinger)) %>%
    # go back to one row per bigram
    unique()

  # calc num of character reps, finger reps, hand reps and hand alts bigrams
  data_3c_bigramHandCategTally <- data_3b_bigramFreqsAndHands %>%
    group_by(sentNum, sentence, n_bigrams, bigramHandCateg) %>%
    tally() %>%
    ungroup() %>%
    pivot_wider(names_from = bigramHandCateg,
                values_from = n,
                values_fill = 0) %>%
    # ensure all 4 columns are created even if none of the input text contains
    # them
    add_cols(c("handAlternation", "handRepetition", "fingerRepetition", "charRepetition"))


  # put together the bigram frequency results and the hand alternation results
  data_3_sentFreqsAlts <- left_join(data_2_aveBigramFreqs,
                                    data_3c_bigramHandCategTally,
                                    by = c("sentNum", "sentence")) %>%
    # calculate the bigram hand categories as proportion of the number of bigrams
    mutate(across(matches(c("charRepetition", "fingerRepetition",
                            "handRepetition", "handAlternation")),
                  ~ .x / n_bigrams,
                  .names = "{.col}Prop")) %>%
    select(-c("handAlternation", "handRepetition", "fingerRepetition", "charRepetition"))

  print("---Done.")



  print("Calculating the proportion of characters that are on the right side of the keyboard for each sentence.")

  # separate out all the characters into rows
  data_4b_sentPropRightHand <- data_3_sentFreqsAlts %>%
    select(sentNum, sentence) %>%
    # split the characters, and make them lowercase
    mutate(character = str_split(sentence, ""),
           character_lower = str_split(tolower(sentence), "")) %>%
    unnest(c(character, character_lower)) %>%
    # join on the standard hand (side of the keyboard) and distance from home row
    left_join(select(keyboardInfo, character, standardHand,
                     homeRowDist), by = c("character_lower" = "character")) %>%
    # save the intermediate result in a dataframe called byCharacter
    {. ->> data_4c_byCharacter } %>%
    # filter out the spaces
    filter(character != " ") %>%
    # count how many there are for each hand
    group_by(sentNum, as.factor(standardHand), .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    rename("standardHand" = "as.factor(standardHand)") %>%
    # get the proportion of each hand
    group_by(sentNum) %>%
    mutate(totNonSpace = sum(n),
           prop = n / totNonSpace) %>%
    ungroup() %>%
    # filter to right hand
    filter(standardHand == "right") %>%
    # keep only these columns
    select(sentNum, prop) %>%
    # rename the column
    rename(propRightHand = prop)

  print("---Done.")


  print("Calculating the distance each key is from the 'home row' and averaging for each sentence.")

  # average distance from home row
  data_4d_aveHomeRowDist <- data_4c_byCharacter %>%
    # filter out the spaces
    filter(character != " ") %>%
    group_by(sentNum) %>%
    summarise(aveDistFromHomeRow = mean(homeRowDist))

  print("---Done.")

  print("Calculating the average number of keystrokes required per word for each sentence.")

  data_4e_minStrokesCalcs <- data_4c_byCharacter %>%
    group_by(sentNum) %>%
    mutate(isUppercase = if_else(str_count(character, "[A-Z]") == 1, TRUE, FALSE),
           isShiftedPunct = if_else(character %in%
                                      c("!", "\"", "£", "$", "%", "^",
                                        "&", "*", "(", ")", "_", "+",
                                        "{", "}", ":", "@", "~", "<",
                                        ">", "?"), TRUE, FALSE),
           # within each train of capital letters, number each capital consecutively
           upperTrainSeq = if_else(isUppercase,
                                   sequence(rle(isUppercase)$lengths),
                                   as.integer(0)),
           # give each uppercase train an ID
           upperTrainID = if_else(isUppercase,
                                  cumsum(upperTrainSeq == 1),
                                  as.integer(0))) %>%
    # calculate the total number of uppercase per uppercase train
    group_by(sentNum, upperTrainID) %>%
    mutate(upperTrainLen = max(upperTrainSeq)) %>%
    ungroup() %>%
    group_by(sentNum) %>%
    # assuming here that people will use shift if it's 1-2 caps in a row and
    # caps lock of 3+ caps in a row
    mutate(numStrokesForThisChar = case_when(
      # if it's a shifted punctuation, it's 2 keystrokes (key plus shift)
      # assumption: only one in row, not !!!!, ????, :::::
      isShiftedPunct ~ 2,
      # is first of an uppercase train of any length, it's 2 keystrokes
      # (letter plus shift/caps)
      isUppercase & upperTrainSeq == 1 ~ 2,
      # is the second in an uppercase train, it's 1 keystroke
      # (shift assumed already held down)
      isUppercase & upperTrainSeq == 2 ~ 1,
      # is the non-last in an uppercase train of 3+, it's 1 keystroke
      # (shift assumed already held down)
      isUppercase & upperTrainLen >= 3 &
        upperTrainSeq >= 3 & upperTrainSeq < upperTrainLen ~ 1,
      # if it's the last in an uppercase train of 3+, it's 2 keystrokes
      # (assumed to turn off caps lock)
      isUppercase & upperTrainLen >= 3 & (upperTrainSeq == upperTrainLen) ~ 2,
      # otherwise, 1
      TRUE ~ 1
    )) %>%
    ungroup()

  data_4f_minStrokes <- data_4e_minStrokesCalcs %>%
    group_by(sentNum, sentence) %>%
    summarise(minStrokes = sum(numStrokesForThisChar))


  # join onto the other variables
  data_4_sentIntermediate <- data_3_sentFreqsAlts %>%
    left_join(data_4b_sentPropRightHand, by = "sentNum") %>%
    left_join(data_4d_aveHomeRowDist, by = "sentNum") %>%
    left_join(data_4f_minStrokes, by = c("sentNum", "sentence")) %>%
    mutate(meanStrokesPerWord = minStrokes / numActualWords,
           meanStrokesPerChar = minStrokes / numChars)


  print("---Done.")

  print("Calculating the mean number of syllables per word for each sentence.")

  data_5_withSylls <- data_4_sentIntermediate %>%
    mutate(meanSyllsPerWord = textstat_readability(sentence,
                                                   "meanWordSyllables")$meanWordSyllables)

  print("---Done.")


  print("Calculating the proportion of words that are in the top 1000 most frequent.")


  top1000WordsLemma <- read_csv(here("data", "predictor_calculation_aids", "wordFrequencyLemma.csv"),
                                show_col_types = FALSE) %>%
    filter(lemRank <= 1000)

  # split the sentences into words
  data_6b_byWord <- select(data_5_withSylls, sentNum, sentence) %>%
    # words are things between spaces
    mutate(word = strsplit(sentence, " ")) %>%
    unnest(word) %>%
    # do it again but with commas as sometimes there isn't a space after a comma
    mutate(word = strsplit(word, ",")) %>%
    unnest(word) %>%
    # if the "word" doesn't end in a letter or number, remove that character(s)
    mutate(word_trimmedPunct = str_replace(word,
                                           "[^A-Za-z0-9]+$",
                                           ""),
           # if the "word" doesn't start with a number of letter, remove that character(s)
           word_trimmedPunct = str_replace(word_trimmedPunct,
                                           "^[^A-Za-z0-9]+",
                                           ""),
           # lowercase word (will search for all cases - e.g. if the word is at the
           # start of a sentence so has a capital, don't only want to get the word
           # frequency of when the word has a capital- want general word frequency
           word_trimmedPunct_lower = tolower(word_trimmedPunct)) %>%
    # get rid of rows that are empty (e.g. was a double space or a comma in
    # the wrong place
    filter(word_trimmedPunct_lower != "") %>%
    group_by(sentNum) %>%
    mutate(wordNum = 1:n()) %>%
    ungroup()

  data_6c_wordFreqs <- data_6b_byWord %>%
    mutate(word_trimmedPunct_lower_apostOff = case_when(
      # if the last two chars are in this list, remove them
      str_sub(word_trimmedPunct_lower, start = -2) %in% c("'s", "'m", "'d") ~
        str_sub(word_trimmedPunct_lower, end = -3),
      # if the last three chars are in this list, remove them
      str_sub(word_trimmedPunct_lower, start = -3) %in% c("'re", "'ve", "'ll", "n't") ~
        str_sub(word_trimmedPunct_lower, end = -4),
      # otherwise, leave it alone
      TRUE ~ word_trimmedPunct_lower),
      # is it in the top 1000 words (lemmas allowed)
      top1000lemm = if_else(word_trimmedPunct_lower_apostOff %in% top1000WordsLemma$word, 1, 0),
      # how many chars for the words in the top 1000?
      freqWordChars = if_else(top1000lemm == 1, nchar(word_trimmedPunct_lower), NA_real_))

  data_6d_sentWordFreqs <- data_6c_wordFreqs %>%
    group_by(sentNum, sentence) %>%
    summarise(highFreqWordProp = mean(top1000lemm),
              highFreqWordCharProp = sum(freqWordChars, na.rm = TRUE) / nchar(first(sentence))) %>%
    ungroup()

  data_6_withWordFreqs <- left_join(data_5_withSylls, data_6d_sentWordFreqs)

  print("---Done.")


  ## new ##

  print("Calculating the mean word frequency for each sentence.")

  # read in the word freq info (from reviewer)
  subtlex <- read_tsv(here("data", "predictor_calculation_aids",
                           "SUBTLEXus74286wordstextversion.txt"),
                      show_col_types = FALSE) %>%
    # convert to lowercase for case insensitivity (sentence words also converted to lower case)
    mutate(Word = tolower(Word))

  # split the sentences into words
  data_6e_byWord <- select(data_5_withSylls, sentNum, sentence) %>%
    # words are things between spaces
    mutate(word = strsplit(sentence, " ")) %>%
    unnest(word) %>%
    # do it again but with commas as sometimes there isn't a space after a comma
    mutate(word = strsplit(word, ",")) %>%
    unnest(word) %>%
    # clean the apostrophes
    mutate(word = str_replace_all(word,
                                  pattern = "[’‘`ʹʼ＇]",  replacement = "'")) %>%
    # do it again but now split along the apostrophe, to match with the word freq database
    mutate(word = strsplit(word, "'")) %>%
    unnest(word) %>%
    # if the "word" doesn't end in a letter or number, remove that character(s)
    mutate(word_trimmedPunct = str_replace(word,
                                           "[^A-Za-z0-9]+$",
                                           ""),
           # if the "word" doesn't start with a number of letter, remove that character(s)
           word_trimmedPunct = str_replace(word_trimmedPunct,
                                           "^[^A-Za-z0-9]+",
                                           ""),
           # lowercase word (will search for all cases - e.g. if the word is at the
           # start of a sentence so has a capital, don't only want to get the word
           # frequency of when the word has a capital- want general word frequency
           word_trimmedPunct_lower = tolower(word_trimmedPunct)) %>%
    # get rid of rows that are empty (e.g. was a double space or a comma in
    # the wrong place
    filter(word_trimmedPunct_lower != "") %>%
    group_by(sentNum) %>%
    mutate(wordNum = 1:n()) %>%
    ungroup()

  # add on the word frequencies
  data_6f_wordFreqs <- data_6e_byWord %>%
    left_join(select(subtlex, Word, FREQcount),
              by = c("word_trimmedPunct_lower" = "Word")) %>%
    mutate(FREQcount = if_else(is.na(FREQcount), 0, FREQcount))

  data_6g_sentWordFreqs <- data_6f_wordFreqs %>%
    group_by(sentNum, sentence) %>%
    summarise(meanWordFreq = mean(FREQcount)) %>%
    ungroup()

  data_6_withWordFreqs <- left_join(data_6_withWordFreqs, data_6g_sentWordFreqs)

  print("---Done.")


  print("Calculating the proportion of words that are not recognised as real words.")

  data_7b_nonDictWords <- data_6b_byWord %>%
    mutate(inDictUS = hunspell_check(word_trimmedPunct, dict = dictionary("en_US")),
           inDictGB = hunspell_check(word_trimmedPunct, dict = dictionary("en_GB")),
           inDictCA = hunspell_check(word_trimmedPunct, dict = dictionary("en_CA")),
           inDictAU = hunspell_check(word_trimmedPunct, dict = dictionary("en_AU"))) %>%
    rowwise() %>%
    mutate(nonDictWord = if_else(any(inDictUS, inDictGB, inDictCA, inDictAU), 0, 1)) %>%
    ungroup() %>%
    mutate(nonDictWordChars = if_else(nonDictWord == 1, nchar(word_trimmedPunct), NA_real_))

  data_7c_sentNonDictWords <- data_7b_nonDictWords %>%
    group_by(sentence) %>%
    summarise(propWordsNonDictWords = sum(nonDictWord) / n(),
              propCharsNonDictWords = sum(nonDictWordChars, na.rm = TRUE) / nchar(first(sentence)))

  data_7_withNonDictWords <- left_join(data_6_withWordFreqs, data_7c_sentNonDictWords)


  print("---Done.")


  options(dplyr.summarise.inform = TRUE)

  data_8_predictorVariables <- select(data_7_withNonDictWords,
                                      sentNum,
                                      sentence,
                                      numChars,
                                      minStrokes,
                                      numActualWords,
                                      meanStrokesPerWord,
                                      meanWordLength,
                                      meanWordLengthPropSent,
                                      highFreqWordProp,
                                      highFreqWordCharProp,
                                      meanWordFreq,
                                      propWordsNonDictWords,
                                      propCharsNonDictWords,
                                      meanSyllsPerWord,
                                      biFreqMean,
                                      propBiTop15,
                                      charRepetitionProp,
                                      fingerRepetitionProp,
                                      handRepetitionProp,
                                      handAlternationProp,
                                      lowercaseProp,
                                      uppercaseProp,
                                      numbersProp,
                                      symbolsProp,
                                      spacesProp,
                                      lowercasePropNonSpace,
                                      uppercasePropNonSpace,
                                      numbersPropNonSpace,
                                      symbolsPropNonSpace,
                                      meanStrokesPerChar,
                                      propRightHand,
                                      aveDistFromHomeRow)

  return(data_8_predictorVariables)

  beep()

}

# function to calculate typability
calculate_typability <- function(sentencesAndPredictors) {
  load(here("data", "typability-index-model.RData"))

  data_with_typability <- sentencesAndPredictors %>%
    mutate(typability = predict(model2, sentencesAndPredictors))

  data_with_typability
}

# function to split word/sentence into bigrams
str_split_bigrams = function(x) {
  substring(x, first = 1:(nchar(x) - 1), last = 2:nchar(x))
}

# function to add the stated columns if they don't exist
add_cols <- function(df, cols) {
  add <- cols[!cols %in% names(df)]
  if (length(add) != 0 ) df[add] <- 0
  return(df)
}


# function to suggest groups based on requested group type
suggest_groups <- function(text_with_scores,
                           grp_type,
                           num_grps,
                           grp_size = floor(nrow(text_with_scores) / num_grps)) {

  total_itms <- nrow(text_with_scores)

  # validate the input combination
  if (num_grps * grp_size > total_itms) {
    showModal(modalDialog(
      title = "Error",
      paste("The specified combination of ", num_grps, " groups with a group size of ",
            grp_size, " exceeds the total number of items (", total_itms, ").",
            "Please refresh the application and try again."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    return()  # stop the execution here
  }


  if (grp_type == "matched") {

    # assign a group to every item (temporarily), based on sorted typability
    suggested_groups <- text_with_scores %>%
      arrange(typability) %>%
      mutate(temp_group = (row_number() - 1) %% num_grps + 1)

    # calculate start and end positions of the final items that will be selected
    positions <- calculate_start_end_matched(total_itms, num_grps, grp_size)

    # add 'selected' column based on start and end
    suggested_groups <- suggested_groups %>%
      mutate(selected = row_number() >= positions$start &
               row_number() <= positions$end)

    # adjust group assignment to ensure the first group starts with 1
    suggested_groups <- suggested_groups %>%
      mutate(group = if_else(selected, temp_group, NA_integer_),
             # calculate adjustment for the first group
             adjustment = if_else(selected, (group[which(selected)[1]] - 1) %% num_grps, NA_integer_),
             # apply adjustment to group
             group = if_else(selected, ((group - adjustment - 1) %% num_grps + 1), NA_integer_)) %>%
      select(-adjustment) %>%
      mutate(group = factor(group))


  } else if (grp_type == "divergent") {

    # calculate start and end positions for divergent groups
    positions <- calculate_start_end_divergent(total_itms, num_grps, grp_size)

    # initialise the selected column in the suggested_groups dataframe
    suggested_groups <- text_with_scores %>%
      arrange(typability) %>%
      mutate(temp_group = ntile(typability, num_grps),
             selected = FALSE)

    # loop through each group and set selected based on positions
    for (i in 1:num_grps) {
      group_start <- positions$start[i]
      group_end <- positions$end[i]

      # if gap size is zero, handle start and end adjustments accordingly
      if (i > 1 && positions$start[i] == positions$start[i - 1]) {
        group_start <- positions$end[i - 1] + 1
        group_end <- group_start + grp_size - 1
      }

      # update the selected column for the current group
      suggested_groups$selected[seq_len(nrow(suggested_groups)) >= group_start &
                                  seq_len(nrow(suggested_groups)) <= group_end] <- TRUE
    }

    # assign groups
    suggested_groups <- suggested_groups %>%
      mutate(group = factor(if_else(selected, temp_group, NA_integer_)))
  }

  return(suggested_groups)

}




calculate_start_end_matched <- function(total_items, num_groups, group_size) {

  # calculate the position of the median-typability item
  median_pos <- floor(total_items / 2)

  # total number of selected items is the minimum of total items or num_groups * group_size
  total_selected <- min(total_items, num_groups * group_size)

  # calculate start based on total_selected being even or odd
  if (total_selected %% 2 == 0) {
    start <- median_pos - (total_selected / 2) + 1 # even
  } else {
    start <- median_pos - floor(total_selected / 2) # odd
  }

  # calculate end
  end <- start + total_selected - 1

  # return start and end
  return(list(start = start, end = end))
}




calculate_start_end_divergent <- function(total_items, num_groups, group_size) {

  # make an empty dataframe to store the positions
  positions <- data.frame(group = integer(num_groups),
                          start = integer(num_groups),
                          end = integer(num_groups))

  # calculate the start and end for the first group
  positions$group[1] <- 1
  positions$start[1] <- 1
  positions$end[1] <- group_size

  # calculate the start and end for the last group
  positions$group[num_groups] <- num_groups
  positions$end[num_groups] <- total_items
  positions$start[num_groups] <- total_items - group_size + 1

  # for the middle group(s)
  if(num_groups > 2) {

    # figure out the gap size between groups (set middle groups based on gap size between them)
    total_to_select <- num_groups * group_size
    total_unselected <- total_items - total_to_select
    gap_size <- total_unselected / (num_groups - 1)

    # calculate the start and end positions for the middle groups
    for (i in 2:(num_groups - 1)) {

      # set the group number
      positions$group[i] <- i

      # calculate start for current group based on previous group end and gap size
      positions$start[i] <- floor(positions$end[i - 1] + gap_size + 1)
      # previously ceiling(positions$end[i - 1] + gap_size)

      # calculate end position for current group
      positions$end[i] <- positions$start[i] + group_size - 1

    }
  }

  return(positions)
}




# Shiny UI
ui <- fluidPage(

  titlePanel("The Typability Index Web App"),

  sidebarLayout(
    sidebarPanel(
      # radio buttons to select input source - upload / dhakal
      radioButtons("input_choice", "Choose Input Source",
                   choices = c("Upload my own text" = "upload",
                               "Use Dhakal sentences" = "dhakal"),
                   selected = "upload"),

      # conditional UI to reset file input when "Upload my own text"  selected again
      uiOutput("file_input_ui"),

      # calculate typability button (visible for both upload and dhakal)
      actionButton("calculate_button", "Calculate Typability"),

      # UI to input number of groups and grouping type, initially hidden
      conditionalPanel(

        condition = "output.showGroupingOptions",

        h4("Suggest Groupings"),

        # radio buttons to select grouping type
        radioButtons("grouping_type", "Grouping Type",
                     choices = c("Matched" = "matched",
                                 "Divergent" = "divergent"),
                     selected = "matched"),

        numericInput("num_groups", "Number of Groups", value = 2, min = 2, step = 1),

        # group size input
        conditionalPanel(

          # show group size only if assign_all is unchecked
          condition = "input.assign_all == false",

          numericInput("group_size", "Group Size", value = 5, min = 1, step = 1)),

        # checkbox to assign all items to groups
        checkboxInput("assign_all",
                      label = div(
                        "Assign all* items to groups",
                        title = "Groups will be of equal size, with a remainder of 1 to (group_size - 1) items potentially unassigned."
                      ),
                      value = FALSE),

        # button to assign text to groups
        actionButton("suggest_groups_button", "Suggest Groupings"),

        # checkbox to hide items not assigned to groups. Only shown after suggesting groupings
        conditionalPanel(

          condition = "output.groupsSuggested",

          checkboxInput("hide_unassigned", "Hide items not assigned to groups", value = FALSE)),

        # button to download results (initially hidden)
        conditionalPanel(

          condition = "output.showGroupingOptions",
          downloadButton("download_button", "Download Results"))),

      # horizontal line
      hr()

    ),

    mainPanel(

      # conditional description
      conditionalPanel(
        condition = "!output.showGroupingOptions",

        p("Welcome to the Typability Index Web App!"),

        tags$div(
          "The accompanying manuscript is currently under review but can read the pre-print ",
          tags$a(href="https://osf.io/preprints/psyarxiv/qxuv5",
                 target="_blank",
                 "here.")),
        p(),
        p("This app calculates typability scores for text data and allows you to group items based on different criteria. Start by choosing an input source, then calculate typability."),
        p(),
        p("This version of the tool is for non-commercial research purposes only."),
        p(),
        p("Note: If your input text contains common symbols not insertable by a US QWERTY keyboard, these will be converted to their keyboard variants. Minus signs, en-dashes, and em-dashes will be converted to hyphens; curly apostrophes will be converted to straight apostrophes, curly and angled quotes will be converted to straight quotes, ellipsis will be converted to \"...\", and multiplication signs (×) will be converted to \"x\".")
      ),

      # density plot of groupings, only when groups have been suggested
      conditionalPanel(

        condition = "output.groupsSuggested",

        plotOutput("density_plot")),

      # table for typability scores and group assignments (only after Calculate Typability)
      tableOutput("typability_group_table")
    )
  )
)


# Server logic
server <- function(input, output, session) {

  # reactive to store the dataframe with text, typability, group columns
  scored_text <- reactiveVal(data.frame(text = character(), typability = numeric()))

  # create a fileInput ui that will reset when switching input sources
  output$file_input_ui <- renderUI({
    if (input$input_choice == "upload") {
      fileInput("input_file", "Choose TXT File", accept = ".txt")
    }
  })

  # reactive to control visibility of grouping options
  showGroupingOptions <- reactiveVal(FALSE)

  # reactive to track whether groupings suggested
  groupsSuggested <- reactiveVal(FALSE)

  # function to reset all inputs
  reset_inputs <- function() {
    updateNumericInput(session, "num_groups", value = 2)
    updateNumericInput(session, "group_size", value = 5)
    updateCheckboxInput(session, "assign_all", value = FALSE)
    updateRadioButtons(session, "grouping_type", selected = "matched")
    groupsSuggested(FALSE)
  }

  # when input source changes --------------------------------------

  observeEvent(input$input_choice, {

    # reset the reactive dataframe and clear displayed tables when switching sources
    scored_text(data.frame(text = character(), typability = numeric()))

    # hide the table
    output$typability_group_table <- renderTable(NULL)

    # reset show grouping options
    showGroupingOptions(FALSE)

    # reset all inputs to their default states
    reset_inputs()
  })

  # when calculate button is clicked -------------------------------------------

  observeEvent(input$calculate_button, {

    if (input$input_choice == "dhakal") {

      # load Dhakal data when Calculate Typability is clicked
      scored_text(data.frame(text = dhakal_typs$text, typability = dhakal_typs$typability))

    } else {

      # for uploaded text
      req(input$input_file)
      text_data <- readLines(input$input_file$datapath)
      text_data <- standardiseSymbols(text_data)
      predictors <- calculatePredictorVariables(text_data)
      typability_scores <- calculate_typability(predictors)$typability
      scored_text(data.frame(text = text_data, typability = typability_scores))
    }

    # display typability and groups in a single table without group column
    output$typability_group_table <- renderTable({
      scored_text()
    })

    # show grouping options after calculating typability
    showGroupingOptions(TRUE)
  })

  # bind the reactive value to the output
  output$showGroupingOptions <- reactive({
    showGroupingOptions()
  })

  # update  output binding
  outputOptions(output, "showGroupingOptions", suspendWhenHidden = FALSE)

  # suggest groupings when Suggest Groupings button is clicked --------------------------------------
  observeEvent(input$suggest_groups_button, {

    req(scored_text())
    current_data <- scored_text()

    # determine group size based on checkbox
    grp_size <- if (input$assign_all) {
      floor(nrow(current_data) / input$num_groups)
    } else {
      input$group_size
    }

    suggested_groups <- suggest_groups(current_data, input$grouping_type, input$num_groups, grp_size) %>%
      select(-temp_group, -selected)
    scored_text(suggested_groups)

    groupsSuggested(TRUE)

    # display updated table with typability and assigned groups
    output$typability_group_table <- renderTable({

      # filter out unassigned groups if the checkbox is checked
      if (input$hide_unassigned) {
        suggested_groups <- suggested_groups %>%
          filter(!is.na(group))
      }
      suggested_groups
    })

    # density plot ---------------------------------------
    output$density_plot <- renderPlot({

      req(scored_text()$group)
      plot_data <- scored_text()

      # filter out rows with NA values if hide_unassigned  checked
      if (input$hide_unassigned) {
        plot_data <- plot_data %>% filter(!is.na(group))
      }

      ggplot(plot_data, aes(x = typability, color = group, fill = group)) +
        geom_density(alpha = 0.3) +
        labs(
          title = "Density Plot of Typability by Group",
          x = "Typability Score",
          y = "Density") +
        theme_minimal() +
        scale_color_brewer(palette = "Set1")
    })
  })

  # download results as CSV -----------------------------------------
  output$download_button <- downloadHandler(
    filename = function() { "typability_scores.csv" },
    content = function(file) {
      write.csv(scored_text(), file, row.names = FALSE)
    }
  )

  # bind the reactive value to track if groups have been suggested
  output$groupsSuggested <- reactive({
    groupsSuggested()
  })

  outputOptions(output, "groupsSuggested", suspendWhenHidden = FALSE)
}

# run
shinyApp(ui = ui, server = server)
