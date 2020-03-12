
Rcpp::loadModule("temulator_module", TRUE)


# .onLoad <- function(libname, pkgname) {
#   setMethod("show", TEMULATOR:::TEMULATOR_object, function(object) {
#     object$print()
#   })
# }