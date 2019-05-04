
#private function to check vector of probabilities
check_prob <- function(prob) {
  if (is.na(prob)) {
    stop("invalid prob value")
  }
  if (prob < 0 | prob > 1) {
    stop("prob has to be a number between 0 and 1")
  }
  TRUE
}

#private function to check if trials input is a value value for number of trials
check_trials <- function(trials) {
  if (is.na(trials)) {
    stop("invalid trials value")
  }
  if (!is.numeric(trials)) {
    stop("trials must be numeric")
  }
  if (trials < 0) {
    stop("trials must be non-negative")
  }
  TRUE
}

#private function to check if success input is a valid value for number of successes
check_success <- function(success, trials) {
  if (success > trials | success < 0) {
    stop("invalid success value")
  }
  TRUE
}

#private auxiliary functions that return corresponding value from computed summary measure
aux_mean <- function(trials, prob) {
  trials*prob
}
aux_variance <- function(trials, prob) {
  trials*prob*(1-prob)
}
aux_mode <- function(trials, prob) {
  integer((trials*prob) + prob)
}
aux_skewness <- function(trials, prob) {
  (1-(2*prob))/((trials*prob*(1-prob))^0.5)
}
aux_kurtosis <- function(trials, prob) {
  (1-(6*prob*(1-prob)))/(trials*prob*(1-prob))
}

#' @title Binomial Choosing
#' @description calculates the number of combinations in which k successes can occur in n trials
#' @param n vector of number of trials (numeric)
#' @param k vector of number of successes (numeric)
#' @return number of calculated combinations
#' @export
#' @examples
#' 
#' bin_choose(n = 10, k = 5)
#' 
bin_choose <- function(n, k) {
  check_success(k,n)
  
  combo <- (factorial(n))/(factorial(k)*factorial(n-k))
  combo
}

#' @title Binomial Probability
#' @description calculates the probability of getting k successes in n trials
#' @param trials vector of number of trials (numeric)
#' @param success vector of number of successes (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated probability
#' @export
#' @examples
#' 
#' bin_probability(success = 2, trials = 5, prob = 0.5)
#' 
bin_probability <- function(success, trials, prob) {
  check_trials(trials)
  check_prob(prob)
  check_success(success, trials)
  
  binprob <- bin_choose(trials, success)*(prob^success)*((1-prob)^(trials-success))
  binprob
}

#' @title Binomial Distribution
#' @description calculates the binomial distribution with given params
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return data frame containing successes in one column and probabilities in another
#' @export
#' @examples
#' 
#' bin_distribution(trials = 5, prob = 0.5)
#' 
#' dis11 <- bin_distribution(trials = 5, prob = 0.5)
#' plot(dis1)
#' 
bin_distribution <- function(trials, prob) {
  kth <- rep(0, trials+1)
  time <- c(0:trials)
  for (i in time) {
    kth[i+1] <- bin_probability(i, trials, prob)
  }
  
  dist <- data.frame(success = c(0:trials), probability = kth)
  class(dist) <- c("bindis", "data.frame")
  dist
}

#' @export
plot.bindis <- function(x) {
  plot1 <- barplot(height = x$probability, names.arg = x$success)
  plot1
}

#' @title Binomial Cumulative
#' @description calculates the cumulative binomial distribution with given params
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return data frame containing successes, probability, and cumulative prob
#' @export
#' @examples
#' 
#' bin_cumulative(trials = 5, prob = 0.5)
#' 
#' dis2 <- bin_cumulative(trials = 5, prob = 0.5)
#' plot(dis2)
#' 
bin_cumulative <- function(trials, prob) {
  cumulative <- rep(0, trials+1)
  time <- c(0:trials)
  distdat <- bin_distribution(trials, prob)
  for (i in time) {
    cumulative[i+1] <- sum(cumulative[i]) + distdat$probability[i+1]
  }
  cumdist <- cbind(distdat, cumulative)
  class(cumdist) <- c("bincum", "data.frame")
  cumdist
}

#' @export
plot.bincum <- function(x) {
  cumplot <- plot(x$success, x$cumulative, type = "b")
  cumplot
}

#' @title Binomial Variable
#' @description creates a list of binomial random variables
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return named list
#' @export
#' @examples
#' 
#' bin_variable(trials = 10, prob = 0.3)
#' 
#' bin1 <- bin_variable(trials = 10, prob = 0.3)
#' print(bin1)
#' summary(bin1)
#' 
#' binsum1 <- summary(bin1)
#' print(binsum1)
#' 
bin_variable <- function(trials, prob) {
  check_trials(trials)
  check_prob(prob)
  
  trialsv <- trials
  probv <- prob
  names(trialsv) <- "trials"
  names(probv) <- "prob"
  listvar <- list(trialsv, probv)
  class(listvar) <- "binvar"
  
  listvar
}

#' @export
print.binvar <- function(x) {
  cat('\n\n"Binomial variable"\n\n')
  
  cat('Parameters\n')
  cat(paste0("- number of trials:", " ", x[[1]]), "\n")
  cat(paste0("- prob of success:", " ", x[[2]]), "\n")
}

#' @export
summary.binvar <- function(x) {
  n <- x[[1]]
  p <- x[[2]]
  sumvar <- list(n,
                 p,
                 aux_mean(n,p),
                 aux_variance(n,p),
                 aux_mode(n,p),
                 aux_skewness(n,p),
                 aux_kurtosis(n,p))
  names(sumvar) <- c("trials",
                     "prob",
                     "mean",
                     "variance",
                     "mode",
                     "skewness",
                     "kurtosis")
  class(sumvar) <- "summary.binvar"
  
  sumvar
}

#' @export
print.summary.binvar <- function(x) {
  cat('\n\n"Binomial variable"\n\n')
  
  cat('Parameters\n')
  cat(paste0("- number of trials:", " ", x[[1]]), "\n")
  cat(paste0("- prob of success:", " ", x[[2]]), "\n\n")
  
  cat('Measures\n')
  cat(paste0("- mean:", " ", x[[3]]), "\n")
  cat(paste0("- variance:", " ", x[[4]]), "\n")
  cat(paste0("- mode:", " ", x[[5]]), "\n")
  cat(paste0("- skewness:", " ", x[[6]]), "\n")
  cat(paste0("- kurtosis:", " ", x[[7]]), "\n")
}

#' @title Binomial Mean
#' @description lists binomial mean
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated mean
#' @export
#' 
bin_mean <- function(trials, prob) {
  aux_mean(trials, prob)
}

#' @title Binomial Var
#' @description lists binomial variance
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated variance
#' @export
#' 
bin_variance <- function(trials, prob) {
  aux_variance(trials, prob)
}

#' @title Binomial Mode
#' @description lists binomial mode
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated mode
#' @export
#' 
bin_mode <- function(trials, prob) {
  aux_mode(trials, prob)
}

#' @title Binomial Skew
#' @description lists binomial skewness
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated skewness
#' @export
#' 
bin_skewness <- function(trials, prob) {
  aux_skewness(trials, prob)
}

#' @title Binomial Kurtosis
#' @description lists binomial kurtosis
#' @param trials vector of number of trials (numeric)
#' @param prob vector of probabilities (numeric)
#' @return calculated kurtosis
#' @export
#' 
bin_kurtosis <- function(trials, prob) {
  aux_kurtosis(trials, prob)
}


