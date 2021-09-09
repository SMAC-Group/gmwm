# .onAttach <- function(...){
#   if (!interactive()) return()
#   
#   local_version = utils::packageDescription('gmwm')
#   cran_version = packageVersionCRAN("gmwm")
#   
#   if(!is.null(cran_version) && length(cran_version) != 0L){
#     latest_version = utils::compareVersion(cran_version[1],local_version$Version)
#     
#     d = if(latest_version == 0){
#      'CURRENT'
#     }else if(latest_version == 1){
#       'OUT OF DATE'
#     }else{
#       'DEVELOPMENT'
#     }
#   
#   }else{ # Gracefully fail.
#    d = "ERROR IN OBTAINING REMOTE VERSION INFO"
#    latest_version = 0
#   }
#   
#   packageStartupMessage('Version: ', local_version$Version, ' (', d,') built on ', local_version$Date)
#     if(latest_version == 1){
#       packageStartupMessage('\n!!! NEW VERSION ', cran_version[1] ,' RELEASED ON ', cran_version[2],' !!!')
#       packageStartupMessage('Download the latest version: ',cran_version[1],' from CRAN via `install.packages("gmwm")`\n')
#     }
#   packageStartupMessage('Check for updates and report bugs at https://github.com/SMAC-Group/gmwm')
#   packageStartupMessage('To see the user guides use browseVignettes("gmwm").')
# }