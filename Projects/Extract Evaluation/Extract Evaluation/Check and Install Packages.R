# List of packages you want to check
packages <- c("dplyr", "ggplot2", "tidyr")

# Function to check and install packages
check_and_install <- function(package) {
  if (!require(package, quietly = TRUE)) {
    cat(paste("Package", package, "is not installed. Installing now...\n"))
    install.packages(package)
    library(package, character.only = TRUE)
  } else {
    cat(paste("Package", package, "is already installed.\n"))
  }
}

# Loop through the list of packages and check/install
for (pkg in packages) {
  check_and_install(pkg)
}
