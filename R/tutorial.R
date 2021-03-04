# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#
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
