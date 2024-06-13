.onAttach <- function(...) {
  msg <- create_startup_message()
  packageStartupMessage(msg)
}

# Create package startup message
create_startup_message <- function() {
  v_pkg <- create_desc("sfgp")
  msg <- paste0(
    "Attached sfgp", v_pkg,
    ". Type ?sfgp to get started."
  )
  return(msg)
}

# Create package description
create_desc <- function(pkg_name) {
  Lib <- dirname(system.file(package = pkg_name))
  pkgdesc <- suppressWarnings(
    utils::packageDescription(pkg_name, lib.loc = Lib)
  )
  if (length(pkgdesc) > 1) {
    out <- paste0(" ", pkgdesc$Version)
  } else {
    out <- ""
  }
  return(out)
}
