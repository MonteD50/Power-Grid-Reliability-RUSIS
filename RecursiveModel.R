library(gtools)

start_time <- Sys.time()

df <- data.frame(CWP=integer(), 
                 prob=numeric(),
                 forkNum=integer())

CWP.individual <- function(customers, probs, branchProb, forkNumber) {
  x = list()
  for (i in 1:length(customers)) {
    x[[i]] <- c(customers[i], i)
  }
  allProbs = numeric()
  zeroProb = (1 - branchProb) + (branchProb*prod(1 - probs))
  print(paste("CWP=", 0, "with Prob=", zeroProb))
  df[nrow(df) + 1,] = list(CWP=0, prob=zeroProb, forkNum=forkNumber)
  allProbs = append(allProbs, zeroProb)
  allCombos = integer()
  
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
      allCombos = append(allCombos, s)
      print(paste("CWP=", s, "with Prob=", prob))
      df[nrow(df) + 1,] = list(CWP=s, prob=prob, forkNum=forkNumber)
      allProbs = append(allProbs, prob)
    }
  }
  print(paste("Sum of probs=", sum(allProbs)))
  allCombos = append(allCombos, 0)
  return(list(df=df, combos=unique(allCombos)))
}

x = c(3, 10, 11, 3, 3, 3)
x.p = c(0.93, 0.93, 0.93, 0.93, 0.93, 0.93)
forks.x = 0.93

y = c(2, 5, 3, 3)
y.p = c(0.93, 0.93, 0.93, 0.93)
forks.y = 0.93

z = c(11, 21, 3, 3)
z.p = c(0.93, 0.93, 0.93, 0.93)
forks.z = 0.93

w = c(32, 44, 77, 3)
w.p = c(0.93, 0.93, 0.93, 0.93)
forks.w = 0.93

result.x = CWP.individual(x, x.p, forks.x, 1)
df = result.x$df
result.y = CWP.individual(y, y.p, forks.y, 2)
df = result.y$df
result.z = CWP.individual(z, z.p, forks.z, 3)
df = result.z$df
result.w = CWP.individual(w, w.p, forks.w, 4)
df = result.w$df
result.s = CWP.individual(w, w.p, forks.w, 5)
df = result.s$df
result.p = CWP.individual(w, w.p, forks.w, 6)
df = result.p$df

allTheCombos = list(result.x$combos,
                    result.y$combos,
                    result.z$combos,
                    result.w$combos,
                    result.s$combos,
                    result.p$combos)

combos = expand.grid(allTheCombos)
combos$value = apply(combos, 1, sum)
combinedCombos = unique(combos$value)

nestedProb <- function(all.combos, specificAllCombo, forkNumber) {
  if (forkNumber == 2) {
    prob = 0
    for (j in all.combos[[2]]) {
      first.prob = df[df$CWP == (specificAllCombo - j) & df$forkNum == 1,]
      if (nrow(first.prob) != 0) {
        firstProb = 0
        for (w in 1:nrow(first.prob)) {
          firstProb = firstProb + first.prob[w, "prob"]
        }
        second.prob = df[df$CWP == j & df$forkNum == 2,]
        secondProb = 0
        for (w in 1:nrow(second.prob)) {
          secondProb = secondProb + second.prob[w, "prob"]
        }
        prob = prob + firstProb*secondProb
      }
    }
    return(prob)
  }
  prob = 0
  for (j in all.combos[[forkNumber]]) {
    second.prob = df[df$CWP == j & df$forkNum == forkNumber,]
    secondProb = 0
    for (w in 1:nrow(second.prob)) {
      secondProb = secondProb + second.prob[w, "prob"]
    }
    prob = prob + (nestedProb(allTheCombos, specificAllCombo - j, 
                              forkNumber - 1)*secondProb)
  }
  return(prob)
}
allProbs = numeric()
forkN = 6
combinedCombos = sort(combinedCombos)
for (i in combinedCombos) {
  prob = nestedProb(allTheCombos, i, forkN)
  allProbs = append(allProbs, prob)
  print(paste("CWP=", i, "with Prob=", prob))
  
}
print(paste("Sum of all Probs=", sum(allProbs)))

plot(combinedCombos, allProbs, main="PDF of Power Grid",
     xlab="CWP", ylab="Probability")

end_time <- Sys.time()
print(paste("This took:", end_time - start_time))