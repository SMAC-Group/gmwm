# Contributing to a `SMAC-Group` repository

`SMAC-Group` repositories are able to receive community feedback in three ways:

1. Submit a bug report via an issue.
1. Request a feature request in an issue.
1. Providing a new feature or change via a pull request.

## Bug Reports

When creating a bug report, you must provide a minimal working example (MWE). Doing this enables us to verify the bug exists and provides us with an ability to test our solution against. We define a minimal working example to include:

  * (Data) minimal dataset, necessary to reproduce the error
    * Consider using `dput(dataset)` to obtain an R data structure export
    * Or use a very simple data.frame
  * (Code) The minimal **runnable** code necessary to reproduce the error, which can be run on the given dataset.
    * Within the code, please make sure it is **readable** by having approriate spacing, short variable names, comments
    * Emit anything that distracts from the key issue / bug.
  * (Platform) the necessary information on the used packages, R version and system it is run on.
    * Do *not* submit `sessionInfo()` without being asked to.
  * (RNG) In case of random processes, a seed (set by `set.seed()`) for reproducibility.

## Feature Requests

Currently, feature requests are taken under advisement, however, we do have our own internal plan for how the development of the `SMAC-Group` package will go. If enough community members need a given feature, we may elevate the priority of it.

## Pull requests

To submit a new feature or change done yourself, do the following:

1. Create a fork of the repository in git.
1. Make your changes in this fork of the repository.
1. Push fork to github.
1. Submit a pull request (e.g. PR).
1. Discuss the pull request.
1. Rinse'n'Repeat 4-5 until the change is accepted or clear it will not be merged in.
