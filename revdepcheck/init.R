if (!require(xfun)) install.package("xfun")
pkg <- "TMB"; which <- "all"
db = available.packages(type = "source")
res = xfun:::check_deps(pkg, db, which)
message("Installing dependencies of reverse dependencies")
print(system.time(install.packages(res$install)))
pkgs <- res$check
tars = xfun:::download_tarball(pkgs, db, dir = ".")
oldname <- dir(pattern="tar.gz")
newname <- sub("_.*",".tar.gz",oldname)
file.rename(oldname, newname)
