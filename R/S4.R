setOldClass("python.builtin.object")
setOldClass("python.builtin.method")

#' mirror OpenCRAVAT module metadata to R
#' @importClassesFrom S4Vectors character_OR_NULL
#' @export
setClass("OpenCRAVATModule",
 representation(code_size = "integer", data_size = "numeric", data_sources = "list", 
   data_versions = "list", datasource = "ANY", description = "character", 
   developer = "list", has_version = "python.builtin.method", 
   hidden = "logical", latest_version = "character", name = "character", 
   publish_time = "character", size = "integer", tags = "character_OR_NULL", 
     title = "character", type = "character", versions = "character"))

#' display OpenCRAVAT module metadata in R
#' @importFrom stringr str_wrap
#' @export
setMethod("show", "OpenCRAVATModule", function(object) {
 cat("OpenCRAVAT module ", object@name, "\n")
 cat("---\n")
 cat(stringr::str_wrap(object@description, 50), "\n")
})

#' define set of R representations of OpenCRAVAT module metadata
#' @export
setClass("OpenCRAVATModuleSet", representation(modset="list", created="ANY"))

#' display set of R representations of OpenCRAVAT module metadata
#' @export
setMethod("show", "OpenCRAVATModuleSet", function(object) {
cat("OpenCRAVATModuleSet, created ", object@created, "\n")
cat("There are", length(object@modset), "modules.\n")
cat("Module types/counts:")
print(table(sapply(object@modset, function(x) x@type)))
})

#' use reticulate to query OpenCRAVAT module store and populate an R object with metadata about available resources
#' @import reticulate
#' @import methods
#' @export
populate_module_set = function() {
 cra = reticulate::import("cravat", convert=TRUE)
 modnames = cra$admin_util$search_remote(".*")
 allinfo = lapply(modnames, function(x)
   cra$admin_util$get_remote_module_info(x))
 atts = names(allinfo[[1]])
 vals = lapply(allinfo, function(x) {
    tmp = lapply(atts, function(y) x[[y]])
    names(tmp) = atts
    do.call("new", c("OpenCRAVATModule", tmp))
    })
 names(vals) = modnames
 new("OpenCRAVATModuleSet", modset=vals, created=date())
}

#' convert some of the fields of module metadata to a data.frame
#' @importMethodsFrom BiocGenerics as.data.frame
#' @export
setMethod("as.data.frame", "OpenCRAVATModuleSet", function(x) {
   nna = function(x) {
          if (is.null(x)) return(NA)
          return(x)
          }
   tmp = lapply(x@modset, function(x) {
      data.frame(title=x@title, type=x@type, developer=nna(x@developer$email), latest_version=x@latest_version,
         name=x@name, description=x@description, published=x@publish_time)
      })
   do.call(rbind, tmp)
})
  
 

