# setup.R: packages, alliance definition, and the analysis panel `d`. Sourced at the top of every analysis script.
# The panel is built once per R session, so running all scripts through run_all.R builds it only once.
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr); library(sandwich); library(lmtest); library(fixest)
  library(splines); library(zoo); library(ggplot2); library(patchwork); library(sf); library(cshapes); library(rnaturalearth)
})
# The article uses the Correlates of War defence pact (Type I) as its alliance definition.
# Set ALLIANCE_DEF=all before starting R to count all four commitment types instead.
if (Sys.getenv("ALLIANCE_DEF") == "") Sys.setenv(ALLIANCE_DEF = "defense")
if (!identical(get0(".panel_def", envir = globalenv()), Sys.getenv("ALLIANCE_DEF"))) {
  eval(parse(text = readLines("scripts/panel.R")), envir = globalenv())
  assign(".panel_def", Sys.getenv("ALLIANCE_DEF"), envir = globalenv())
}
