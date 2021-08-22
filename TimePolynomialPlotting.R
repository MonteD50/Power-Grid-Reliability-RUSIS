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

CWP.individual <- function(customers, yearNum) {
  probs = weibullFailureRate(rep(yearNum, length(customers)))
  branchProb = weibullFailureRate(yearNum)
  x = list()
  for (i in 1:length(customers)) {
    x[[i]] <- c(customers[i], i)
  }
  allProbs = numeric()
  zeroProb = (1 - branchProb) + (branchProb*prod(1 - probs))
  df[nrow(df) + 1,] = list(CWP=0, prob=zeroProb)
  allProbs = append(allProbs, zeroProb)
  allCombos = integer()
  allCombos = append(allCombos, 0)
  
  for (i in 1:length(x)) {
    combo = combn(x, i)
    for (j in 1:ncol(combo)) {
      comb = combo[,j]
      indexes = which(x %in% comb)
      s = 0
      for (c in comb) {
        s = s + c[1]
      }
      failure = 1 - probs
      for (ind in indexes) {
        failure[ind] = 1 - failure[ind]
      }
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

CWP.all <- function(yearNum, ...) {
  powerGrid = list(...)
  lPowerGrid = length(powerGrid)
  print(paste("Running for", lPowerGrid, "forks"))
  all.polynomial = polynomial(1)
  count = 1
  for (i in 1:lPowerGrid) {
    if (count %% 1000 == 0) {
      print(paste("On Loop:", i))
    }
    individual.pdf = CWP.individual(powerGrid[i][[1]], 
                                      yearNum)
    x.polynomial = integer(max(individual.pdf$CWP)+1)
    x.polynomial[individual.pdf$CWP+1] <- individual.pdf$prob
    x.polynomial <- as.polynomial(x.polynomial)
    all.polynomial = all.polynomial*x.polynomial
    count = count + 1
  }
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
  allProbs = coef(all.polynomial)
  allProbs = unlist(lapply(allProbs, function(x) {x[x!=0]}))
  finalPDF = data.frame(allCombos, allProbs)
  
  print(format(tail(finalPDF), scientific=F))
  print(paste("Sum of Probs is:", sum(allProbs)))
  
  return(finalPDF)
}

first = CWP.all(0.01, c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1))
second = CWP.all(1, c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1))
third = CWP.all(5, c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1),
                c(3, 10, 11, 3, 3, 3),
                c(5, 19, 3, 5, 2),
                c(1, 6, 8),
                c(19, 22, 1, 5, 4),
                c(7, 4, 3, 1))
fourth = CWP.all(20, c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1),
                 c(3, 10, 11, 3, 3, 3),
                 c(5, 19, 3, 5, 2),
                 c(1, 6, 8),
                 c(19, 22, 1, 5, 4),
                 c(7, 4, 3, 1))

plot(first$allCombos, first$allProbs, pch = ".", xlab = "CWP", ylab = "Probability",
     main = "PMF of Power Grid", col = "red")
points(second$allCombos, second$allProbs, pch = ".", 
       col = "green")
points(third$allCombos, third$allProbs, pch = ".", 
       col = "blue")
points(fourth$allCombos, fourth$allProbs, pch = ".",
       col = "pink")
legend(x = "topleft", y = 0.92, 
       legend=c("0.01 Years", "1 Year", "5 Years", "20 Years"),
       col = c("red", "green", "blue", "pink"),
       cex=0.5, lty = 1:2)
end_time <- Sys.time()
print(paste("This took:", end_time - start_time))