library(gtools)

start_time <- Sys.time()

df <- data.frame(CWP=integer(), 
                 prob=numeric())

PDF.real = data.frame(CWP=integer(), 
                      prob=numeric())

CWP.individual <- function(customers, probs, branchProb) {
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
  return(df)
}

combinedMGF <- function(CWPCombos, ProbCombos) {
  combos = expand.grid(CWPCombos)
  combosProb = expand.grid(ProbCombos)
  columns = ncol(combos)
  if (columns == 1) {
    CWP = with(combos, Var1)
    prob = with(combosProb, Var1)
  }else if (columns == 2) {
    CWP = with(combos, Var1+Var2)
    prob = with(combosProb, Var1*Var2)
  } else if (columns == 3) {
    CWP = with(combos, Var1+Var2+Var3)
    prob = with(combosProb, Var1*Var2*Var3)
  } else if (columns == 4) {
    CWP = with(combos, Var1+Var2+Var3+Var4)
    prob = with(combosProb, Var1*Var2*Var3*Var4)
  } else if (columns == 5) {
    CWP = with(combos, Var1+Var2+Var3+Var4+Var5)
    prob = with(combosProb, Var1*Var2*Var3*Var4*Var5)
  } else if (columns == 6) {
    CWP = with(combos, Var1+Var2+Var3+Var4+Var5+Var6)
    prob = with(combosProb, Var1*Var2*Var3*Var4*Var5*Var6)
  } else if (columns == 7) {
    CWP = with(combos, Var1+Var2+Var3+Var4+Var5+Var6+Var7)
    prob = with(combosProb, Var1*Var2*Var3*Var4*Var5*Var6*Var7)
  } else if (columns == 8) {
    CWP = with(combos, Var1+Var2+Var3+Var4+Var5+Var6+Var7+Var8)
    prob = with(combosProb, Var1*Var2*Var3*Var4*Var5*Var6*Var7*Var8)
  } else if (columns == 9) {
    CWP = with(combos, Var1+Var2+Var3+Var4+Var5+Var6+Var7+Var8+Var9)
    prob = with(combosProb, Var1*Var2*Var3*Var4*Var5*Var6*Var7*Var8*Var9)
  }
  
  PDF = data.frame(cbind(CWP, prob))
  
  CWPSeen = integer()
  for (i in 1:nrow(PDF)) {
    cwp = PDF[i,1]
    prob = PDF[i, 2]
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

CWP.all <- function(numSplit, ...) {
  powerGrid = list(...)
  lPowerGrid = length(powerGrid)
  firstLoopCount = lPowerGrid / 3
  print(paste("Running for", firstLoopCount, "forks"))
  secondLoopCount = round(firstLoopCount / numSplit, 0)
  stored.combos = list()
  stored.probs = list()
  for (i in 1:firstLoopCount) {
    if (i == 1) {
      individual.pdf = CWP.individual(powerGrid[i][[1]], 
                                      powerGrid[i + 1][[1]],
                                      powerGrid[i + 2][[1]])
    } else {
      individual.pdf = CWP.individual(powerGrid[(i*3)-2][[1]], 
                                      powerGrid[(i*3)-1][[1]],
                                      powerGrid[(i*3)][[1]])
    }
    stored.combos = append(stored.combos, list(individual.pdf$CWP))
    stored.probs = append(stored.probs, list(individual.pdf$prob))
  }
  final.combos = list()
  final.probs = list()
  for (i in 1:secondLoopCount) {
    print(paste("On Loop:", i))
    if (i == 1) {
      allSubsetCombos = stored.combos[i:numSplit]
      allSubsetProbs = stored.probs[i:numSplit]
    } else {
      if (i*numSplit <= firstLoopCount) {
        endIndex = i*numSplit
      } else {
        endIndex = firstLoopCount
      }
      startIndex = ((i - 1)*numSplit) + 1
      allSubsetCombos = stored.combos[startIndex:endIndex]
      allSubsetProbs = stored.probs[startIndex:endIndex]
    }
    subsetMGF = combinedMGF(allSubsetCombos, allSubsetProbs)
    final.combos = append(final.combos, list(subsetMGF$CWP))
    final.probs = append(final.probs, list(subsetMGF$prob))
  }
  PDF.real = combinedMGF(final.combos, final.probs)
  print(format(PDF.real, scientific=F))
  print(sum(PDF.real$prob))
}

CWP.all(4, 
        c(3, 10, 11, 3, 3), c(0.93, 0.93, 0.93, 0.93, 0.93), 0.93,
        c(2, 5, 3, 3), c(0.93, 0.93, 0.93, 0.93), 0.93,
        c(11, 21, 3, 3), c(0.93, 0.93, 0.93, 0.93), 0.93,
        c(4, 6, 13, 3), c(0.93, 0.93, 0.93, 0.93), 0.93,
        c(8, 12, 11, 3), c(0.93, 0.93, 0.93, 0.93), 0.93,
        c(3, 10, 11, 4, 3), c(0.93, 0.93, 0.93, 0.93, 0.93), 0.93,
        c(2, 5, 3, 3), c(0.93, 0.93, 0.93, 0.93), 0.93,
        c(11, 21, 3), c(0.93, 0.93, 0.93), 0.93,
        c(3, 10, 11, 3, 3), c(0.93, 0.93, 0.93, 0.93, 0.93), 0.93,
        c(2, 5, 3), c(0.93, 0.93, 0.93), 0.93,
        c(2, 5, 2), c(0.93, 0.93, 0.93), 0.93,
        c(11, 21, 5), c(0.93, 0.93, 0.93), 0.93)

end_time <- Sys.time()
print(paste("This took:", end_time - start_time))