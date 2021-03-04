# I am using https://uf-repro.github.io/writing-R-packages/notes.html
# for reference
#' This is a random silly function
#' @param df A data.frame
#' @return a character vector
#' @export
f <- function(df)
{
  names(df)
  print('Hello World')
}
