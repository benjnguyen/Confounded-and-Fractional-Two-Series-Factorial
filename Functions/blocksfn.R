library(cfcdae)
library(conf.design)

blocks <- function(k, q, x, binaryCol = FALSE, table = FALSE, printset = FALSE)
{
  # ToDo: Allow x to be a string vector, e.g. x = c(BCD, ACD)
  #       This is done by finding the treatment terms in the same
  #       row as the correct binary string.
  
  
  # ToDo: Find the generalized interaction confounded variable
  #       This is done by taking confounders, say x = c("01110", "10110")
  #       adding them component wise
  #       and doing entrywise mod 2
  #       then matching the result to a sequence that appears
  #       in the list of standard orders (str.input)
  #       and finding the treatment letter in the (term) column.
  
  binaryseq = as.matrix(strsplit(x, ""))
  parsebin = do.call(rbind, binaryseq)
  matrix(as.numeric(parsebin), nrow = nrow(parsebin))
  testgen = matrix(as.numeric(parsebin), nrow = nrow(parsebin))
  if (q == 1) {
    set = conf.set(testgen, 2)
    str.set = matrix(0, nrow = 1, ncol = 1)
    str.set[1] = paste(set, collapse = "")
  }
  else if (q >= 2) {
    set = conf.set(testgen, 2)
    str.set = matrix(rep(0, nrow(set)), nrow = nrow(set), ncol = 1)
    for (i in 1:nrow(set)) {
      str.set[i] = paste(set[i,], collapse = "")
    }
  }
  # ToDo: Add some stopping conditions that make sure we have valid inputs
  #       For example, for q defining contrasts, we need x of length q
  if (!q == length(x))
    stop("The number of defining contrasts specified in 'x' is not equal to the number of required contrasts q; check length(x).")
  for(i in 1:length(x)) {
    if(!k == nchar(x[i]))
      stop("At least one binary sequence in 'x' is not equal to the required length k; check nchar(x) == k.")
  }
  # Enumerating standard order for k factors
  ntrt = 2^k
  nblock = 2^q
  input = expand.grid(replicate(k, 0:1, simplify = FALSE))
  str.input = matrix(rep(0, nrow(input)), nrow = nrow(input), ncol = 1)
  for (i in 1:nrow(input)) {
    str.input[i] =  paste(input[i,], collapse = "")
  }
  # Enumerating inputted confounders
  confounder <- rep(list(NULL), q)
  names(confounder) <- paste("confounder", 1:q, sep = "")
  for (i in 1:q) {
    confounder[[i]] <- x[i]
  }
  # ToDo: Find the generalized interaction confounded variable
  len = nchar(confounder[[1]])
  # Taking subset of standard ordered treatments
  # Subset is defined by confounding variable
  # E.g. if we're in 2^4 design {A,B,C,D} main effects
  #      and we want to confound BCD
  #      We should only inspect columns B,C,D when doing addition mod 2
  # Three is an implausible value in arithmetic modulo 
  matching <- rep(list(NULL), q)
  names(matching) = paste("match", 1:q, sep = "")
  for (i in 1:q) {
    matching[[i]] =  matrix(rep(rep(3, len), nrow(input)), nrow = nrow(input))
    for (j in 1:(ntrt)) {
      for (b in 1:len) {
        matching[[i]][j, b] = (substr(str.input[j],b,b) == substr(confounder[[i]],b,b))
      }
    }
  }
  for (i in 1:q) {
    matching[[i]] <- matching[[i]][, which(strsplit(x[i], "")[[1]]=="1")]
    matching[[i]] <- as.matrix(apply(matching[[i]], 1, sum))
    matching[[i]] <- matching[[i]] %% 2
  }
  # Re-naming 0 to even, and 1 to odd
  for (i in 1:q) {
    for (m in 1:ntrt) {
      if (matching[[i]][m] == 1) {
        matching[[i]][m] <- "o"
      } 
      else if (matching[[i]][m] == 0) {
        matching[[i]][m] <- "e"
      }
    }
  }
  # Giving treatment names to standard ordered treatments
  term <- matrix((rep("", ntrt)), nrow = ntrt)
  term[1] <- "(1)"
  for (i in 1:ntrt) {
    for(j in 1:k) {
      if (substr(str.input[i], j, j) == 1) {
        term[i] <- paste(term[i], letters[j], sep = "")
      }
    }
  }
  # Declare and initialize 2^q blocks
  block <- rep(list(NULL), nblock)
  names(block) <- paste("block", 1:nblock, sep = "")
  parity = expand.grid(replicate(q, c("e", "o"), simplify = FALSE))
  str.parity = matrix(rep("", nrow(parity)), nrow = nrow(parity), ncol = 1)
  for (i in 1:nrow(parity)) {
    # Generalize it to 1,...,q
    str.parity[i] = paste(as.matrix(parity[i,1:q]), collapse = "")
  }
  constructblock <- cbind(term, str.input)
  for (i in 1:q) {
    constructblock <- cbind(constructblock, matching[[i]])
  }
  if (q == 1) { 
    l = constructblock[,3]
  } else {
    l = apply(subset(constructblock[,3:(2+q)]), 1, paste, collapse = "")
  }
  
  for (i in 1:nblock){
    #block[[i]] = which(paste(constructblock[1:32,3:(2+q)], collapse = "") == str.parity[i])
    block[[i]] = which(l == str.parity[[i]])
  }
  blocknumber <- matrix((rep(0, dim(constructblock)[1])), nrow = dim(constructblock)[1])
  for(i in 1:nblock) {
    blocknumber[block[[i]]] <- i
  }
  constructblock <- cbind(constructblock, blocknumber)
  # Add the generalized interaction (confounded interaction terms)
  # This would correspond to taking the confounders, summing entrywise, taking mod 2
  # then associating the result with an input in the list of standard ordered treatments.
  # names(constructblock) <- c("Treatment", "Std. Order", ""....,, blocknumber)
  # this needs to spell out the confounder over the .... columns
  results <- rep(list(NULL), nblock)
  names(results) <- paste(str.parity, "block", sep = "")
  for (i in 1:length(results)) {
    results[[i]] <- constructblock[block[[i]], 1]
  }
  # Print confounded set (defn contrasts + generalized interaction)
  ind <- matrix(rep(0, length(str.set)), nrow = length(str.set))
  for (i in 1:length(str.set)) {
    ind[i] = which(str.set[i] == constructblock[,2])
  }
  confset = constructblock[ind, 1]
  if (printset == TRUE) {
    cat("The confounded set is \n")
    cat("{", toupper(as.character(confset)),"}\n")
  }
  constructblock <- as.data.frame(constructblock)
  if (binaryCol == TRUE) {
    colnames(constructblock) <- c("treatment", "binary trt", 
                                  x, "blocknumber")
  }
  colnames(constructblock) <- c("treatment", "binary trt", 
                                toupper(as.character(constructblock[match(x, constructblock[,2]),1])),
                                "blocknumber")
  if (table == TRUE) {
    print(constructblock)
  }
  results <- as.data.frame(results)
  results
}
