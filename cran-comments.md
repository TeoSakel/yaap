## Resubmission

This is a resubmission. In this version I have:

* Updated the Title and Description fields so software and package names are
  formatted according to CRAN policy.
* Added method references to the Description field using CRAN's auto-linking format.
* Replaced commented-out examples in `archetypes.Rd` and
  `fitted.archetypes.Rd` with executable toy examples.
* Removed unnecessary `\dontrun{}` wrappers from examples, using small
  executable examples guarded for suggested packages where needed.
* Regenerated the Rd files.

## R CMD check results

Local Linux, R 4.6.0, `R CMD check --as-cran yaap_1.0.0.tar.gz`:

0 errors | 0 warnings | 1 note

* This is a new submission.

Windows R-devel, win-builder, <https://win-builder.r-project.org/66AID23ZP4cn/>:

0 errors | 0 warnings | 1 note

* This is a new submission.

macOS release, macOS Builder, <https://mac.R-project.org/macbuilder/results/1779894638-816890347dafa8d1/>:

0 errors | 0 warnings | 0 notes

## revdepcheck results

There are currently no downstream dependencies for this package because it is not yet on CRAN.
