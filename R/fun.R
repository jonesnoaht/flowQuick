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

#' Compensate the Data
#' @description It will attempt to find the compensation
#' controls using the regex '.*Setup.*fcs'
#' @param comp_controls the folder containing single-stained
#' compensation controls
#' @param channels the channels to be compensated in the
#' pattern 'BV421-A|BV570-A'
#' @param us_num the number of the unstained control tube in
#' the set
#' @return a list with a summary of the comp controls and a
#' spillover matrix
#' @examples
#' gs.uncomped <- GatingSet(resQC)
#' spillover_matrix <- compensate_data(comp_controls, 'BV421-A|BV570-A', 3)[[2]]
#' gs.comped <- compensate(gs.uncomped, spillover_matrix)
#' @export
compensate_data <- function(comp_controls, channels, us_num) {
  controls <- read.flowSet(path = comp_controls,
                           pattern = '.*Setup.*fcs')
  controls_summary <- summary(controls)

  spillover_matrix <- spillover(controls,
                                unstained = us_num,
                                fsc = 'FSC-A',
                                ssc = 'SSC-A',
                                patt = channels,
                                method = 'median',
                                stain_match = 'intensity',
                                useNormFilt = FALSE
                                )
  list(controls_summary,
       spillover_matrix
       )
}

#' Quickly gate your populations
#' @param gs a gatingset
#' @return the paths generated and a gatingset with gates
#' @export
quick_gate <- function(gs) {
  if (length(gs_get_pop_paths(gs)) > 1){
    gs_pop_remove(gs,paths[1])
    gs_pop_remove(gs,paths[2])
  }
  gs_get_pop_paths(gs)

  pgon <- matrix(c(1.2*10^4,1,
                   2*10^4,4*10^4,
                   4*10^4,1.2*10^5,
                   2.2*10^5,1.2*10^5,
                   2.2*10^5,1),
                 ncol=2, byrow = TRUE)
  colnames(pgon) <- c("FSC-A","SSC-A")
  cells <- polygonGate(.gate = pgon, filterId="Cells")
  gs_pop_add(gs, cells)
  FSC_singlets <- rectangleGate("FSC-A"=c(50,15000),
                                "SSC-A"=c(50000,100000),
                                filterId="FSC-singlet")
  gs_pop_add(gs, FSC_singlets, parent = "Cells")
  SSC_singlets <- rectangleGate("FSC-A"=c(50,15000),
                                "SSC-A"=c(50000,100000),
                                filterId="SSC-singlet")
  gs_pop_add(gs, SSC_singlets, parent = "FSC-singlet")
  recompute(gs)
  paths <- gs_get_pop_paths(gs[[1]])[-1] # what was this?
  list(paths,
       gs
  )
}
