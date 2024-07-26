library(pak)

# Get the list of installed packages
installed_packages <- pak::pkg_list()

# Extract the package names
package_names <- installed_packages$package

# Write the package names to packages.txt
writeLines(package_names, "/home/myuser/work/packages.txt")
