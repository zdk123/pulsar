#' find config file
#'
#' Find a default config file. First calls \code{batchtools::findConfFile} and then find a pulsar default.
#'
#' @param name name of default config or path to config file.
#' @examples
#' ## Default config file provided by pulsar runs code in interactive mode
#' ## This is for testing purposes and executes serially.
#' findConfFile()
#' ## Use the parallel package
#' ## slower than providing the 'ncores' argument to pulsar function, due to
#' ## the overhead of creating the batchtools registry.
#' findConfFile('parallel')
#'
#' ## Use the snow package to register/execute batch jobs on socket clusters.
#' findConfFile('snow')
#' ## Use a TORQUE / PBS queing system. Requires brew template file.
#' findConfFile('torque')
#' findTemplateFile('simpletorque')
#'
#' @details
#' See the batchtools functions \code{batchtools::findConfFile} and \code{batchtools::makeRegistry}. When calling \code{batch.pulsar}, we attempt to use batchtool's default lookup for a config file before calling \code{pulsar::findConfFile}.
#'
#' For clusters with a queuing submission system, a template file, for
#' defining worker node resources and executing the batch R code, will need to
#' be defined somewhere on the system. See \code{\link{findTemplateFile}}.
#' @seealso \code{\link{findTemplateFile}}
#' @export
findConfFile <- function(name='') {
 ## if x is not a file
 ## look for config file using batchtools rules,
 ## otherwise, look in the pulsar system package

  conffile <- batchtools::findConfFile()
  if (!is.na(conffile)) return(conffile)

  if (checkmate::testFileExists(name, access = "r"))
    return(fs::path_real(name))

  ## append type to file extension for default config files
  if (nchar(name)==0) name <- '.R'
  else name <- paste0('.', tools::file_path_sans_ext(name), '.R')

  conffile <- fs::path_real(system.file('config',
                  sprintf('batchtools.conf%s', name), package='pulsar'))
  # }
  if (checkmate::testFileExists(conffile, access = "r")) return(conffile)
  else return(character(0))
}

#' find template file
#'
#' Find a config file from batchtools or default file from pulsar
#'
#' @param name name of default template or path to template file.
#' @examples
#'  \dontrun{
#'  cluster.functions = batchtools::makeClusterFunctionsTORQUE(
#'                      template=pulsar::findTemplateFile('simpletorque'))
#'  }
#' @seealso findConfFile
#' @details
#' See the batchtools functions \code{batchtools::findTemplateFile}, \code{batchtools::makeClusterFunctionsTORQUE}, \code{batchtools::makeClusterFunctionsSGE}, etc, to employ batchtools' default lookup scheme for template files. Supply the output of this function to the \code{template} argument to override batchtools' default.
#'
#' In this case we look for "[name].tmpl" in the pulsar installation directory in the subfolder "templates".
#' @export
findTemplateFile <- function(name) {
  ## get non exported function
#  x   <- tryCatch(.batchtools_findTemplateFile(name), error=function(e) '')
#  if (checkmate::testFileExists(x, access = "r")) return(fs::path_real(x))
#  else {
    x <- system.file("templates", sprintf("%s.tmpl", name), package = "pulsar")
    if (checkmate::testFileExists(x, access = "r")) return(fs::path_real(x))
    else {
      stop(sprintf('Argument \'template\' (=\\"%s\\") must point to a template file or contain the template itself as string (containing at least one newline', name))
    }
#  }
}
