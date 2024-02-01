# Download any data file from a private GitHub repo
# regardless of how large it is.
# Returns a string that will then need to be parsed
# by read_csv or the like to turn it into a data frame.

# Dependencies
require(tidyverse)
require(httr) 
require(rlist)
require(jsonlite)
library(rlist)

Sys.setenv(GITHUB_PAT = "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")


# Your GitHub username or team account name
gh_account <- 'bioth'

# This function accepts two arguments — the name of the repo
# and the path to the file of interest from the main
# directory (including the filename)
fetchGHdata <- function(repo, path) {

  # First you have to authenticate.
  # Store a personal access token in .Renviron
  # See https://blog.exploratory.io/extract-data-from-private-github-repository-with-rest-api-db804fa43d84
  auth <- authenticate(Sys.getenv("GITHUB_PAT"), "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")
  
  # Seperate the filename from the directory
  match <- regexpr("^(.*[\\/])", path)
  if (match[1] > 0) {
    dir <- path %>% substring(match[1], attr(match, "match.length"))
    file <- path %>% substring(attr(match, "match.length") + 1, nchar(path))
  } else {
    dir <- ""
    file <- path
  }
  
  # To handle files larger than 1MB, use this trick:
  # https://medium.com/@caludio/how-to-download-large-files-from-github-4863a2dbba3b
  req_meta <- 
    content(
      GET(
        paste("https://api.github.com/repos", gh_account, repo, "contents", dir, sep="/"), 
        auth
      )
    )
  
  entry <- req_meta %>% list.filter(name == file)
  sha <- entry[1][[1]]$sha
  
  # Grab contents, using sha as a reference
  req_blob <- GET(
    paste("https://api.github.com/repos", gh_account, repo, "git/blobs", sha, sep="/"), 
    auth
  )
  
  # Need to decode the contents, which are returned in base64
  d <- content(req_blob)$content %>%
    base64_dec() %>%
    rawToChar()
  
  return(d)
}


#Select file to load into R environment
github_repo <- "gut-microbiota-iron"
file <- "adult_dss_weight_measurement.csv"

#Transforming the data into a data.frame usable by R
gh_data <- read.table(text = fetchGHdata(github_repo,file), header = TRUE, sep = "\t")
