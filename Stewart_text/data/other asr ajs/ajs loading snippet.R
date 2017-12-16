# the word count filenames are of the form: wordcount_XX.XXXX_XXXXXX.CSV
readCSVWithId <- function(id) {
  str.parts <- str_split(id, "/")[[1]]
  file.name <- paste0("wordcounts/wordcounts_", str.parts[1], "_", str.parts[2], ".CSV")
  word.counts <- read.csv(file.name)
  word.counts$id <- id
  return(word.counts)
}

# read in all the CSV files to a list of data frames
full.data <-  lapply(cites$id, readCSVWithId)
# unlist the data to form one data.frame
full.data.unlisted <- do.call(rbind.data.frame, full.data) 
# and then spread the terms out over the columns of a matrix using tidyr's spread()
doc.term <- spread(full.data.unlisted, key=WORDCOUNTS, value=WEIGHT, fill=0)
doc.ids <- names(doc.term[,1])
doc.term.mat <- as.matrix(doc.term[, -1])