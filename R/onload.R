.onAttach <- function(...) {
  if (!interactive()) return()
  
  tips = c(
    "Need help? E-mail the SMAC-Group: contact@smac-group.com",
    "Find out what's changed in the gmwm package at http://github.com/SMAC-Group/gmwm",
    "Use suppressPackageStartupMessages() to eliminate package startup messages."
  )
  
  tip = sample(tips, 1)
  
  tip = c(tip,"Please note, the optimization algorithm is a work in progress!")
  packageStartupMessage(paste(strwrap(tip), collapse = "\n"))
}