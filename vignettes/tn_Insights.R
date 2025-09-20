## ----setup, include = FALSE---------------------------------------------------
#file.edit(normalizePath("~/.Renviron"))
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=TRUE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 3, 
  fig.width = 5
)

