 

#' use reticulate to determine where locally installed OpenCRAVAT modules reside
#' @export
local_oc_path = function() {
 proc = basilisk::basiliskStart(ocbenv) # avoid package-specific import
 on.exit(basilisk::basiliskStop(proc))
 basilisk::basiliskRun(proc, function() {
 sit = reticulate::import("site", convert=TRUE)
 sit$USER_SITE
 })
}

#' use reticulate/basilisk to list locally installed OpenCRAVAT annotators
#' @export
list_local_annotators = function() 
   dir(paste0(local_oc_path(), "/cravat/modules/annotators"))

#' connect to local annotator SQLite database
#' @import RSQLite
#' @param name character(1) OpenCRAVAT annotator name
#' @export
connect_local_annotator = function(name="clinvar") {
  avail = list_local_annotators()
  if (!(name %in% avail)) stop(paste(name, "not in list of local annotators"))
  pa = local_oc_path()
  dpath = paste0(pa, "/cravat/modules/annotators/", name, "/data/", name, ".sqlite")
  RSQLite::dbConnect(RSQLite::SQLite(), dpath)
}
