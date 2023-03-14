# These functions are used to calculate the optimal number of columns and rows for a display of n items. 
# The display_division function takes an argument n 
# which represents the total number of items to be displayed.
# It then calls two other functions h_display_division and v_display_division 
# to calculate the optimal number of columns (h) and rows (v) for displaying the items.

h_display_division <- function(n, max_div = 5){
  divisors = 4:max_div
  if (n>max_div){
    valid_divisors = divisors[n %% divisors ==0]
    if (length(valid_divisors) > 0){
      return(max(valid_divisors))
    }
    else{
      remains = list()
      for(i in seq_along(divisors)){
        remains[i] = n %% divisors[i]}
      return(divisors[length(divisors)-which.max(rev(remains))+1])
    }
  }
  else{
    return(n)
  }
}

v_display_division <- function(n, h, max_div = 4){
  divisors = 3:max_div
  if (n/h< max_div){
    return(ceiling(n/h))
  }
  else {
    valid_divisors = divisors[n %% divisors*h ==0]
    if (length(valid_divisors) > 0){
      return(max(valid_divisors))
    }
    else{
      remains = list()
      for(i in seq_along(divisors)){
        remains[i] = n %% divisors[i]*h}
      return(divisors[length(divisors)-which.max(rev(remains))+1])
    }
  }
}



display_division <- function(n){
  h = h_display_division(n)
  v = v_display_division(n,h)
  d = ceiling(n/(h*v))
  return(c(h,v,d))
}


