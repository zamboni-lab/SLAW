# Checks requested by devtools::release()
library(devtools)

spell_check()

# R CMD check
check(cran=TRUE, remote=TRUE, manual=TRUE, run_dont_test=TRUE)

# devtools checks
release_checks()

check_rhub()

check_win_devel()

# Update NEWS
# Update DESCRIPTION
# Update cran-comments.md
# Update website

devtools:::git_checks()

# If pass
submit_cran()
