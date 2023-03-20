#Transformers functions
##Basic-ish Maths
geometric_mean <- function(numbers){
  gm = prod(numbers)^(1/length(numbers))
  return(gm)
}
#BoxCox
transformer_boxcox <- function(df, weighted = TRUE){
  t_df = df
  for (colname in colnames(df)[-1:-2]){
    #print(colname)
    lambda = determine_lambda(colname, df)
    if (weighted == TRUE){
      t_df[[colname]]=bc_weighted_transform(df[[colname]],lambda)        
    }
    else{
      t_df[[colname]]=bc_transform(df[[colname]],lambda)
    }
    
  }
  rm(dataf,column,envir=globalenv())
  return(t_df)
}
##Boxcox Modules
  determine_lambda <- function(colname, df) {
    dataf <<- as.data.frame(df)
    column <<- df[,colname]
    # column <<- column
    model <- lm(column ~ group, data = dataf)
    bc <- boxcox(model, lambda = seq(-5, 5))
    lambda <- bc$x[which(bc$y==max(bc$y))]
    return(lambda)
  }
  bc_transform <- function(y, lambda=0) {
    if (lambda == 0L) { log(y) }
    else { (y^lambda - 1) / lambda }
  }
  bc_weighted_transform <- function(y, lambda=0) {
    geom = geometric_mean(y)
    if (lambda == 0L) { log(y) }
    else { (y^lambda - 1) / lambda*geom^(lambda-1)}
  }
#RSN
transformer_rsn <- function(df){
  subset = as.matrix(df[-1:-2])
  sink(nullfile <- tempfile())
  rsn_transform = lumiN(subset, method= "rsn")
  sink(NULL)
  na_cols <- which(colSums(is.na(rsn_transform)) > 0)
  if (length(na_cols) > 0) {
    cat("Columns with NAs:", colnames(rsn_transform)[na_cols], "\n")
    rsn_transform <- rsn_transform[, -na_cols]
  }
  return(cbind(df[1:2],rsn_transform))
}
