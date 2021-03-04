# For FCS Files
#' This QCs all of the files in your folder
#' @param folder The location of your FCS files
#' @return A flowset
#' @export
generate_flowset <- function(folder) {
  fcs.dir <- dir(folder, full.names=TRUE)
  frames <- lapply(fcs.dir, read.FCS)
  as(frames, "flowSet")
}

#' This cleans the data and exports them into a new folder
#' @param fs A flowset
#' @return Cleaned data
#' @examples
#' clean_data(fs)
#' # When you need to re-initialize your clean data:
#' resQC <- read.flowSet(path = "resQC/", pattern = ".*fcs")
#' @export
clean_data <- function(fs) {
  require(flowAI)
  resQC <- flow_auto_qc(fs)
  write.flowSet(resQC, "resQC")
}

