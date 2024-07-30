library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-p", "--param"), type = "character", default = "default",
              help = "parameter to specify", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Main function
main <- function(param) {
  message <- paste("The specified parameter is:", param)
  print(message)
}

# Run the main function with the specified parameter
main(opt$param)