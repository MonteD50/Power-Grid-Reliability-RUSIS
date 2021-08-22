library(polynom)
library(gtools)
library(stringr)
start_time <- Sys.time()

df <- data.frame(CWP=integer(), 
                 prob=numeric())

PDF.real = data.frame(CWP=integer(), 
                      prob=numeric())

weibullFailureRate <- function(componentAge, initialF = 0.07, beta = 0.15) {
  return(1 - (initialF*exp(beta*log(componentAge))))
}

t <- function(t, p) {
  return(1 - ((p*t / (1 + t))))
}

CWP.individual <- function(customers, componentAge, branchProb) {
  # get the probabilities with respect to time
  probs = weibullFailureRate(componentAge)
  branchProb = weibullFailureRate(branchProb)
  
  x = list()
  for (i in 1:length(customers)) {
    x[[i]] <- c(customers[i], i)
  }
  allProbs = numeric()
  # calculate the zero prob
  zeroProb = (1 - branchProb) + (branchProb*prod(1 - probs))
  df[nrow(df) + 1,] = list(CWP=0, prob=zeroProb)
  allProbs = append(allProbs, zeroProb)
  allCombos = integer()
  allCombos = append(allCombos, 0)
  
  for (i in 1:length(x)) {
    # create all i combos of the fork
    combo = combn(x, i)
    for (j in 1:ncol(combo)) {
      comb = combo[,j]
      # get the indexes of the combinations
      indexes = which(x %in% comb)
      s = 0
      for (c in comb) {
        s = s + c[1]
      }
      # change the indices of prob of failure to prob of 
      # success
      failure = 1 - probs
      for (ind in indexes) {
        failure[ind] = 1 - failure[ind]
      }
      # multiply the probs
      prob = branchProb*prod(failure)
      if (s %in% allCombos) {
        prior = df$prob[df$CWP == s]
        df$prob[df$CWP == s] <- prior + prob
      } else {
        allCombos = append(allCombos, s)
        df[nrow(df) + 1,] = list(CWP=s, prob=prob)
        allProbs = append(allProbs, prob)
      }
    }
  }
  # check for duplicates. If found add the probs up
  CWPSeen = integer()
  for (i in 1:nrow(df)) {
    cwp = df[i,1]
    prob = df[i, 2]
    if (cwp %in% CWPSeen) {
      prior = PDF.real$prob[PDF.real$CWP == cwp]
      PDF.real$prob[PDF.real$CWP == cwp] <- prior + prob
    } else {
      CWPSeen = append(CWPSeen, cwp)
      PDF.real[nrow(PDF.real) + 1,] = list(CWP=cwp, prob=prob)
    }
  }
  return(PDF.real)
}

CWP.all <- function(...) {
  powerGrid = list(...)
  lPowerGrid = length(powerGrid)
  firstLoopCount = lPowerGrid / 3
  print(paste("Running for", firstLoopCount, "forks"))
  # start with polynomial 1
  all.polynomial = polynomial(1)
  count = 1
  for (i in 1:firstLoopCount) {
    if (count %% 1000 == 0) {
      print(paste("On Loop:", i))
    }
    # get the PMF of the fork_i
    if (i == 1) {
      individual.pdf = CWP.individual(powerGrid[i][[1]], 
                                      powerGrid[i + 1][[1]],
                                      powerGrid[i + 2][[1]])
    } else {
      individual.pdf = CWP.individual(powerGrid[(i*3)-2][[1]], 
                                      powerGrid[(i*3)-1][[1]],
                                      powerGrid[(i*3)][[1]])
    }
    print("done")
    # create a vector of zeros
    x.polynomial = integer(max(individual.pdf$CWP)+1)
    # assign CWP indexes + 1 to the probabilities
    x.polynomial[individual.pdf$CWP+1] <- individual.pdf$prob
    # convert the vector of probs to a polynomial
    x.polynomial <- as.polynomial(x.polynomial)
    # multiply the previous polynomial to the current one
    all.polynomial = all.polynomial*x.polynomial
    count = count + 1
  }
  # extract the exponent polynomial (is the CWP)
  expression = as.character(all.polynomial)
  allValues = str_split(expression, "\\+")[[1]]
  zeroMet = FALSE
  allCombos = integer()
  for (exp in allValues) {
    val = str_trim(str_split(exp, "\\^")[[1]][2])
    if (is.na(val)) {
      if (zeroMet == FALSE) {
        zeroMet = TRUE
        allCombos = append(allCombos, 0)
      } else {
        allCombos = append(allCombos, 1)
      }
    } else {
      allCombos = append(allCombos, strtoi(val))
    }
  }
  # get the polynomial coefficients (which is the probabilities)
  allProbs = coef(all.polynomial)
  # removes all zero coefficients from coefficients
  allProbs = unlist(lapply(allProbs, function(x) {x[x!=0]}))
  finalPDF = data.frame(allCombos, allProbs)
  
  print(format(tail(finalPDF), scientific=F))
  print(paste("Sum of Probs is:", sum(allProbs)))
  
  plot(allCombos, allProbs, pch = ".", xlab = "CWP", ylab = "Probability",
       main = "PMF of Power Grid", col = "red")
  print(skewness(allProbs))
}

# Main Function Call: 
# One Fork is: c(customers), c(component ages), fork Age
# This is an example of 6 forks
#CWP.all(c(5, 15, 16, 1, 2, 1), c(1, 2, 3, 3, 2, 1), 10,
#        c(3, 10, 11, 3, 5, 3), c(1, 2, 3, 3, 2, 1), 10,
#        c(5, 7, 9, 1, 2), c(1, 2, 1, 1, 1), 1,
#        c(3, 10, 8, 9, 3, 4), c(1, 2, 3, 3, 2, 1), 10,
#        c(4, 2, 14, 3, 3, 3), c(1, 2, 3, 3, 2, 1), 10,
#        c(3, 15, 1), c(1, 2, 3, 3, 2, 1), 10)
CWP.all(rep(4, 20), rep(1, 20), 1)

end_time <- Sys.time()
print(paste("This took:", end_time - start_time))