# load_dependencies_old <- function(pkg_list) {
#   # Load or install packages from list
#   for (pkg in pkg_list) {
#     if (substr(pkg, 1, 11) == "BiocManager") {
#       pkg = substr(pkg, 14, nchar(pkg))
#       if (!requireNamespace("BiocManager", quietly = TRUE)) {
#         install.packages("BiocManager")
#       }
#       if (!require(pkg, character.only = TRUE)) {
#         BiocManager::install(substr(pkg, 14, nchar(pkg)))
#       }
#     }
#     else{
#       if (!require(pkg, character.only = TRUE)) {
#           install.packages(pkg)
#         }
#       }
#       library(pkg, character.only = TRUE)
#     }
#   }

# load_dependencies <- function(pkg_list) {
#   # Load or install packages from list
#   for (pkg in pkg_list) {
#     if (!require(pkg, character.only = TRUE)) {
#       if (substr(pkg, 1, 11) == "BiocManager") {
#         # Install Bioconductor packages
#         if (!requireNamespace("BiocManager", quietly = TRUE)) {
#           install.packages("BiocManager")
#         }
#         pkg = substr(pkg, 14, nchar(pkg))
#         BiocManager::install(pkg)
#       } else {
#         # Install other packages
#         if (!require(pkg, character.only = TRUE)) {
#           install.packages(pkg)
#         }
#       }
#       library(pkg, character.only = TRUE)
#     }
#   }
# }
# 
# export_excel_old <-function(dflist, name = deparse(substitute(dflist))){
#   file_path <- file.path(getwd(), "stats", dataset, paste0(name, ".xlsx"))
#   # Check if the file already exists, and delete it if it does
#   if (file.exists(file_path)) {
#     file.remove(file_path)
#     message("File deleted.")
#   }
#   lapply(1:length(dflist), function(setid){
#     sheet_name = paste("Dataset_", setid, sep="")
#     print(c(setid, file_path, sheet_name))
#     if (!file.exists(file_path)){
#       print("File does not exist")
#       write.xlsx(x=dflist[[setid]], file=file_path, sheetName = paste("Dataset_", setid, sep=""), append=FALSE)
#     }
#     else {
#       print("File Exists")
#       write.xlsx(x=dflist[[setid]], file=file_path, sheetName = paste("Dataset_", setid, sep=""), append=TRUE)
#     }
#   })
# }

