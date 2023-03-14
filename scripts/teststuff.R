X02 = set_colname_adapter(X01_bg_purge)

cutoff = 45


library(ggplot2)

# Create example data
set.seed(123)
data <- data.frame(
  group = sample(c(0,1), 100, replace = TRUE),
  var1 = rnorm(100),
  var2 = rnorm(100),
  var3 = rnorm(100),
  var4 = rnorm(100),
  var5 = rnorm(100)
)
type(X02)
type(X02[2])
str(X02[[2]])
data = X02[[1]]
type(data)
# Reshape data to long format
data_long <- tidyr::pivot_longer(data, cols = -(1:2), names_to = "variable")

# Create violin plot with ggplot2
ggplot(data_long, aes(x = variable, y = value, fill = factor(group))) +
  geom_violin(scale = "count", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "Variable", y = "Value", fill = "Group") +
  theme_bw()


# Create box plots for each variable by group
par(mfrow = c(1, 5)) # Set up a 1x3 grid of plots
for (i in 2:3) {
  boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group")
}

# Generate example data
data <- data.frame(
  group = rep(LETTERS[1:2], each = 10),
  var1 = rnorm(20),
  var2 = rnorm(20),
  var3 = rnorm(20)
)


data = P01_Preprocessed[[2]]

pp_box_plot_multi <- function(set, name){
  # Set up plot area and device
  dim = display_division(ncol(data)-2)
  if ((ncol(data)-2)<=(dim[1]*dim[2])){
    filename = name
    col_range = 3:ncol(data)
    png(filename = paste0(filename,".png"), width = 1200+100*dim[1], height = 1200+150*dim[2], res=250)
    par(mfrow = c(dim[2],dim[1]))
    for (i in col_range) {
      plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group")
    }
    dev.off()
  }
  else{
    modifier = dim[1]*dim[2]
    for (d in seq(dim[3])){
      filename = paste(name,d,sep="_")
      png(filename = paste0(filename,".png"), width = 1200+300*dim[1], height = 900+300*dim[2], res=250)
      col_range = 3:(modifier+2)+(modifier*(d-1))
      if (d == dim[3]){
        #dim = display_division((ncol(data)-2)%%(dim[3]-1))
        col_range = col_range[col_range<=ncol(data)]
      }
      par(mfrow = c(dim[2],dim[1]))
      print(col_range)
      for (i in col_range) {
        plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group")
      }
      dev.off()
    }
  }
}
graphics.off()
pp_box_plot_multi(data, "pp_test/file")
# Create boxplots

PP_box_plot_multi <- function(dflist){
  name = deparse(substitute(dflist))
}


data = P01_Preprocessed[[1]]
PP_box_plot_multi(P01_Preprocessed[1], "Test_1")

pp_box_plot <- function(set){
  melted = data_long
  #melted = melt(set)
  melted$variable = as.factor(melted$variable)
  pplot = ggplot(melted, aes(x=variable, y=value)) + 
    geom_boxplot() + 
    theme(axis.text.x=element_blank())
  #  ggtitle("Box plot of patient antibody expression levels by analyte") +
 #   scale_y_continuous(limits=c(1,2500)) +
 #   geom_hline(yintercept = 250, colour = "red") +
 #   xlab("Analyte") + ylab("Median Flourescent Intensity (MFI)")
  pplot
}
  
pp_box_plot(X02[1])



isprime <- function(n) {
  if (n < 2) {
    return(FALSE)
  }
  for (i in 2:(floor(sqrt(n)))) {
    if (n %% i == 0) {
      return(FALSE)
    }
  }
  return(TRUE)
}

get_prime_factors <- function(n) {
  primes <- c() # Vector to store the prime factors
  d <- 2 # Start with the smallest prime factor
  while (n > 1) {
    if (n %% d == 0) { # Check if d is a factor of n
      primes <- c(primes, d) # Add d to the list of prime factors
      n <- n / d # Divide n by d
    } else {
      d <- d + 1 # Increment d to the next prime number
      while (!isprime(d)) { # Check if d is prime
        d <- d + 1
      }
    }
  }
  return(primes)
}

# Example usage
get_prime_factors(57)

isprime(57)

get_display_division <- function(n, max_div) {
  if (n >= max_div) {
    div <- floor(n / (n %/% max_div))
    if (div > max_div) {
      factors <- unique(c(1L, sort(as.integer(sapply(3L:max_div, function(x) {
        if (max_div %% x == 0L) x else NULL
      }))), decreasing = TRUE))
      div <- max(factors[factors <= n])
      if (is.na(div)) {
        div <- max_div
      }
    }
  } else {
    factors <- unique(c(1L, sort(as.integer(sapply(1L:(n %/% 2 + 1L), function(x) {
      if (n %% x == 0L) c(x, n %/% x) else NULL
    })))))
    div <- min(factors[factors >= max_div])
    if (is.na(div)) {
      div <- max_div
    }
  }
  return(div)
}

# Example usage
get_display_division(12, 5) # Returns 4
get_display_division(25, 5) # Returns 5

n = 12
max_div = 5

h_display_division <- function(n, max_div = 5){
  divisors = 3:max_div
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
  return(c(h,v))
}

display_division(20)
display_division(40)
display_division(31)

h_display_division(n)
v_display_division(31)
max_div
n=31

my_list <- list("apple", "banana", "orange")
last_index <- length(my_list) - 1
last_element <- my_list[[last_index + 1]]
