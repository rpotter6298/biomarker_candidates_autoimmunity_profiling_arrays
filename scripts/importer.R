import_libs <- function() {
  function_files <- list.files(file.path("scripts", "backend"), full.names = TRUE)
  
  for (file in function_files) {
    import_env <- new.env()
    source(file, local = import_env)
    function_names <- ls(import_env)
    function_list <- mget(function_names, import_env)
    assign(basename(file), function_list, envir = .GlobalEnv)
  }
}



worklib = import_libs()