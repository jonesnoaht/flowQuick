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
#' resQC <- flowCore::read.flowSet(path = "resQC/", pattern = ".*fcs")
#' @export
clean_data <- function(fs) {
  .Deprecated("get_fcs_resultsQC")
  resQC <- flowAI::flow_auto_qc(fs)
  flowCore::write.flowSet(resQC, "resQC")
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
#' @import flowCore
#' @examples
#' gs.uncomped <- flowWorkspace::GatingSet(resQC)
#' spillover_matrix <- compensate_data(comp_controls, 'BV421-A|BV570-A', 3)[[2]]
#' gs.comped <- compensate(gs.uncomped, spillover_matrix)
#' @export
compensate_data <- function(comp_controls, channels, us_num) {
  controls <- flowCore::read.flowSet(path = comp_controls,
                                     pattern = '.*Setup.*fcs')
  controls_summary <- summary(controls)

  spillover_matrix <- flowCore::spillover(controls,
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
  requireNamespace("flowCore")
  requireNamespace("flowWorkspace")
  if (length(gs_get_pop_paths(gs)) > 1){
    flowWorkspace::gs_pop_remove(gs,paths[1])
    flowWorkspace::gs_pop_remove(gs,paths[2])
  }
  flowWorkspace::gs_get_pop_paths(gs)

  pgon <- matrix(c(10^4, 1,
                   10^4, 10^4,
                   10^4, 10^5,
                   10^5, 10^5,
                   10^5, 1),
                 ncol=2, byrow = TRUE)
  colnames(pgon) <- c("FSC-A","SSC-A")

  cells <- flowCore::polygonGate(.gate = pgon, filterId="Cells")
  flowWorkspace::gs_pop_add(gs, cells)
  FSC_singlets <- flowCore::rectangleGate("FSC-H"=c(50, 15000),
                                "FSC-W"=c(50000, 100000),
                                filterId="FSC-singlet")
  flowWorkspace::gs_pop_add(gs, FSC_singlets, parent = "Cells")
  SSC_singlets <- flowCore::rectangleGate("SSC-H"=c(50, 15000),
                                "SSC-W"=c(50000, 100000),
                                filterId="SSC-singlet")
  flowWorkspace::gs_pop_add(gs, SSC_singlets, parent = "FSC-singlet")

  flowWorkspace::recompute(gs)

  paths <- gs_get_pop_paths(gs[[1]])[-1]
  list(paths,
       gs
       )
}
#' Set a polygon gate
#'
#' @param gating_set a gating set
#' @param gate_coordinates a vector of points in the format "c(x1, y1, x2, y2 . . . " and so on
#' @param gate_name the name by which you will call the gate
#' @param dimensions the X and Y parameters of the gate (in that order) as strings in a list
#' @param parent a string; the name of the parent gate
#' @importFrom magrittr "%>%"
#' @return modified gating set
#' @export
#'
#' @examples
#' \dontrun{
#' set_scatter_gate(gating_set = gs,
#'                  gate_coordinates = c(2.0e+04, 2.0e+04,
#'                                       1.4e+04, 1.5e+04,
#'                                       1.4e+04, 0.2e+04,
#'                                       3.0e+04, 0.2e+04,
#'                                       3.0e+04, 2.0e+04),
#'                  gate_name = "scatter")
#' }
set_scatter_gate <- function(gating_set = gs,
                             gate_coordinates = c(0.01e+05, 1.5e+05,
                                                  0.01e+05, 1.0e+05,
                                                  0.12e+05, 1.0e+05,
                                                  0.12e+05, 1.5e+05),
                             gate_name = "scatter",
                             dimensions = c("FSC-A","SSC-A"),
                             parent = "root") {
  scatter <- matrix(byrow = TRUE,
                    ncol = 2,
                    data = gate_coordinates)
  colnames(scatter) <- dimensions
  scatter <- flowCore::polygonGate(scatter, filterId = gate_name)
  flowWorkspace::gs_pop_add(gating_set, scatter, name = gate_name, parent = parent)
  flowWorkspace::recompute(gs)
  gating_set
}

#' Try your gate before setting your gate
#'
#' @param gating_set a gating set
#' @param gate_coordinates a vector of points in the format "c(x1, y1, x2, y2 . . . " and so on
#' @param parent a string; the name of the parent gate
#' @param dimensions the X and Y parameters of the gate (in that order) as strings in a list
#'
#' @return a plot
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
#' \dontrun{
#' set_scatter_gate(gating_set = gs,
#' gate_coordinates = c(2.0e+04, 2.0e+04,
#'                      1.4e+04, 1.5e+04,
#'                      1.4e+04, 0.2e+04,
#'                      3.0e+04, 0.2e+04,
#'                      3.0e+04, 2.0e+04),
#' gate_name = "scatter")
#' }
test_plot <- function(gating_set = gs,
                      gate_coordinates = c(0.01e+05, 1.5e+05,
                                           0.01e+05, 1.0e+05,
                                           0.12e+05, 1.0e+05,
                                           0.12e+05, 1.5e+05),
                      parent = "root",
                      dimensions = list("FSC-A",
                                        "SSC-A")) {
  flowWorkspace::recompute(gating_set)
  matrix(byrow = TRUE,
         ncol = 2,
         data = gate_coordinates) %>%
    data.frame() -> gate_coordinates_df
  colnames(gate_coordinates_df) <- c("x","y")
  flowWorkspace::gs_pop_get_data(gating_set, parent) %>%
    ggcyto::autoplot(x = dimensions[[1]],
             y = dimensions[[2]],
             bins = 256) +
    ggplot2::geom_shape(data = gate_coordinates_df,
               ggplot2::aes(x = x,
                   y = y,
                   alpha = 0.1))
}
#' Load your QC results or, if none, run QC on your data
#'
#' @param QC_folder the path to the folder with the QCed data
#' @param raw_folder the path to the folder with the raw files
#' @importFrom magrittr "%>%"
#' @return a cytoset with QCed data
#' @export
#'
#' @examples
#' \dontrun{
#' get_fcs_resultsQC(QC_folder = "./RG Microbiome -APCs panel/resultsQC/",
#'                   raw_folder = "../20200810NJ-ACEFlow/good/")
#' }
get_fcs_resultsQC <- function(QC_folder = "./resultsQC/",
                              raw_folder = "./rawData/") {
  cs <- c()
  try(cs <- flowCore::read.flowSet(path = QC_folder,
                         pattern = ".fcs"))
  if (length(cs) == 0) {
    ggplot2::theme_set(ggthemes::theme_clean())
    flowCore::read.flowSet(path = raw_folder,
                 pattern = ".fcs") %>%
      flowAI::flow_auto_qc() %>%
      flowWorkspace::flowSet_to_cytoset() -> cs
  }
  new_sample_names <- c()
  for (i in 1:length(cs)) {
    new_sample_names <- c(new_sample_names,
                          flowWorkspace::keyword(cs[[i]])$`TUBE NAME`)
  }
  flowWorkspace::sampleNames(cs) <- new_sample_names
  cs
}
#' Make Gating Set
#'
#' @param cyto_set a cytoset
#'
#' @return a gating set
#' @export
make_gs <- function(cyto_set) {
  gating_set <- flowWorkspace::GatingSet(cyto_set)
  Biobase::pData(gating_set) <- cbind(Biobase::pData(gating_set),row.names(Biobase::pData(gating_set)))
  colnames(Biobase::pData(gating_set))[2] <- "tube names"
  for (i in 1:length(gating_set)) {
    Biobase::pData(gating_set[[i]])$`tube names` <- flowWorkspace::keyword(gating_set[[i]])$`TUBE NAME`
  }
  gating_set
}
#' Set gates based on fluorescence-minus-one (FMO) controls
#'
#' @param gating_set a gating set
#' @param parent the name of the parent population
#' @param trim how many to remove from the end of the markernames list
#' @importFrom magrittr "%>%"
#' @return a gating set
#' @export
#'
#' @examples
#' \dontrun{
#' make_gates_from_fmos <- function(gating_set = gs, parent = "SSC-singlet", trim = 1)
#' }
make_gates_from_fmos <- function(gating_set = gs,
                                 parent = "SSC-singlet",
                                 trim = 0){
  l <- list()
  mn <- flowWorkspace::markernames(gating_set[[1]])
  for (i in flowWorkspace::sampleNames(gating_set)) {
    if (grepl("FMO", i)) {
      for (a in 1:(length(mn) - trim)) {
        if (grepl(mn[[a]], i)) {
          flowWorkspace::gs_pop_get_data(gating_set[[i]], parent) %>%
            ggplot2::fortify() %>%
            dplyr::pull(names(mn[a])) %>%
            quantile(.99) -> max
          l[[ names(mn[a]) ]] <- max[[1]]
        }
      }
    }
  }
  have_ld <- FALSE
  for (i in 1:length(l)) {
    if (grepl("dead", mn[[ names(l[i]) ]]) | grepl("Dead", mn[[ names(l[i]) ]])) {
      gate_name <- paste(mn[[ names(l[i]) ]], "neg")
      min_max <- list()
      min_max[[ names(l[i]) ]] <-  c(l[[i]], -Inf)
      gate <- flowCore::rectangleGate(filterId = gate_name, min_max)
      flowWorkspace::gs_pop_add(gating_set, gate, name = gate_name, parent = parent)
      have_ld <- TRUE
      ld_gate_name <- gate_name
    }
  }
  for (i in 1:length(l)) {
    if (have_ld & ! (grepl("dead", mn[[ names(l[i]) ]]) | grepl("Dead", mn[[ names(l[i]) ]]))) {
      gate_name <- paste(mn[[ names(l[i]) ]], "pos")
      min_max <- list()
      min_max[[ names(l[i]) ]] <-  c(l[[i]], Inf)
      gate <- flowCore::rectangleGate(filterId = gate_name, min_max)
      flowWorkspace::gs_pop_add(gating_set, gate, name = gate_name, parent = ld_gate_name)
    } else {
      gate_name <- paste(mn[[ names(l[i]) ]], "pos")
      min_max <- list()
      min_max[[ names(l[i]) ]] <-  c(l[[i]], Inf)
      gate <- flowCore::rectangleGate(filterId = gate_name, min_max)
      flowWorkspace::gs_pop_add(gating_set, gate, name = gate_name, parent = parent)
    }
  }
  flowWorkspace::recompute(gating_set)
  gating_set
}
