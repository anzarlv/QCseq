# Purpose: Script to launch Shiny App for QCseq
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

#' Launch Shiny App for QCseq
#'
#' This functions launches the Shiny app for QCseq The app is intended for
#' those with less coding experience and performs the entire QCseq workflow.
#'
#' @return launches the Shiny app
#'
#' @examples
#' \dontrun{
#' QCseq::runQCseq()
#' }
#'
#' @references
#' Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,
#' Allen J, McPherson J, Dipert A, Borges B (2025). shiny: Web Application
#' Framework for R. R package version 1.12.0, https://shiny.posit.co/.
#'
#' R Core Team (2025). _R: A Language and Environment for Statistical
#' Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.
#'
#' Shiny - File Upload. (2014, July 29). Shiny.
#' https://shiny.posit.co/r/gallery/widgets/file-upload/
#'
#' Shiny - Including HTML, text, and Markdown files. (2014, July 29). Shiny.
#' https://shiny.posit.co/r/gallery/application-layout/including-html-text-and-markdown-files/
#'
#' Shiny - Tabsets. (2014, July 30). Shiny.
#' https://shiny.posit.co/r/gallery/application-layout/tabsets/
#'
#' @export
#' @importFrom shiny runApp

runQCseq <- function(){
  app_dir <- system.file("shiny-scripts", package = "QCseq")
  action_shiny <- shiny::runApp(app_dir, display.mode = "normal")

  return(action_shiny)
}

#[END]
