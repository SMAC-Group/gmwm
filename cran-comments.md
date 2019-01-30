This is the second submission of `gmwm` 

## Test environments
* local OS X install, R 3.2.3 
* win-builder (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs in R 3.2.2 with --as-cran or under win-devel.

We did receive a note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Stephane Guerrier <stef.guerrier@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  ARIMA (19:16)
  GMWM (14:53, 18:23)
  decompositions (20:58)

Suggests or Enhances not in mainstream repositories:
  imudata

Availability using Additional_repositories specification:
  imudata   yes   http://smac-group.com/datarepo
  
  This is because the package we reference a data package that is 40 mb in size and does not adhere to CRAN's policy.
  
## Downstream dependencies

There are no downstream dependencies.