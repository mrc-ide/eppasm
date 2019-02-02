FOIs <- function(pop, hivpop, artpop, i, r_ts, mx, fp) {
# NOTE: definitely slightly overestimate FOI in case of counting both?
# browser()
nRg    <- length(mx$gamma)
cArray <- array(0, c(fp$ss$pAG, fp$ss$pAG, nRg^2, 2, 2))
                   # age a, age a', group combinations, sex, and
                   # 2 matrix (unweight, weight) replace unweight with weighted

ages <- fp$ss$AGE_START:fp$ss$AGE_END 
K <- expand.grid(1:nRg, 1:nRg) # groups combination
# index to lookup which the opposite sex group position in group combinations
iid <- rep(1:2, each=nRg) + nRg * rep(0:(nRg-1), 2) 
# Calculate CsagSAG
# cArray <- mixG(cArray, pop, i, mx, fp)
# -----------------------------------------------------------------------------
for (d in 1:2) {
  for (j in 1:(nRg^2)) {
    cArray[,,j,d,1] <- outer(ages, ages, 
      function(x, y) Cmix(d, x, y, K[j,2], K[j,1], pop, i, mx, fp))
  }   
}
# Weights.
# cArray <- mixG.w(cArray, i, pop, fp, mx)
# -----------------------------------------------------------------------------
for (d in 1:2) { # loop through sex d <> b
	for (j in 1:(nRg^2)) {
		b <- ifelse(d==1, 2, 1) # with the opposite sex
		popA <- rowSums(pop[, sid(d,fp)[K[j,2]], , i]) # myself
		popB <- rowSums(pop[, sid(b,fp)[K[j,1]], , i]) # other
		cArray[,,j,d,2] <- sweepX(cArray[,,j,d,1], 1, popA) / 
					     t(sweepX(cArray[,,iid[j],b,1], 1, popB))
	}
}
cArray[,,,,2] <- na2num(cArray[,,,,2], 1)
for (b in 1:(nRg^2)) cArray[,,b,1,1] %<>% '*'(sqrt(cArray[,,b,1,2])) # Male
for (b in 1:(nRg^2)) cArray[,,b,2,1] %<>% '/'(sqrt(cArray[,,b,2,2])) # Female

# New prevalence with different transmision rate
# -----------------------------------------------------------------------------
pF <- array(0, c(fp$ss$pAG, fp$ss$NG))

for (j in 1:(nRg*2)) { # in each group and sex, count H+ and on art * reduction
	if (any(pF > 1)) message("Risk > 1, checkpF")
	Hstar <- hivpop[,,j,i] + colSums(artpop[,,,j,i]) * mx$phi
	Hsag  <- colSums(sweepX(Hstar, 1, mx$wD[,3])) # weighted by CD4, ∑ over
	nAlpha <- rowSums(pop[,j,,i]) / rep(sumByAG(rowSums(pop[,j,,i])), 
										times = fp$ss$h.ag.span) # expand age
	pF[,j] <- nAlpha * rep(Hsag, times = fp$ss$h.ag.span) # split to each age
	pF[,j] <- pF[,j] / rowSums(pop[,j,,i])
}

pF <- pF * r_ts/(max(fp_par$rvec))

# FOI.
# cArray <- mixRisk(cArray, pF, fp, mx)
# -----------------------------------------------------------------------------
for (d in 1:2) { # loop through sex 
	for (j in 1:(nRg^2)) { # through all sex * group
		b <- ifelse(d==1, 2, 1)
		cArray[,,j,d,1] <- sweepX(cArray[,,j,d,1], 2, 
			pF[,sid(b, fp)[K[j,1]]]) * (1 - mx$cEf * mx$condom[,,d])
	}
}

# Summing over other ages, groups
FOI <- array(0, c(fp$ss$pAG, fp$ss$NG))
if (length(mx$gamma) == 1) { # TODO: find way for auto indexing
	FOI[, 1] <- rowSums(cArray[,,1,fp$ss$m.idx,1]) # f_id(1, g) 1
	FOI[, 2] <- rowSums(cArray[,,1,fp$ss$f.idx,1]) # f_id(2, g) 2
} else if (length(mx$gamma)==2) {
	FOI[, 1] <- rowSums(cArray[,,1:2,fp$ss$m.idx,1]) # f_id(1, g) 1
	FOI[, 2] <- rowSums(cArray[,,3:4,fp$ss$m.idx,1]) # f_id(2, g) 1
	FOI[, 3] <- rowSums(cArray[,,1:2,fp$ss$f.idx,1]) # f_id(1, g) 2
	FOI[, 4] <- rowSums(cArray[,,3:4,fp$ss$f.idx,1]) # f_id(2, g) 2	
} else {
	stop("> 2 groups is not supported")
}
return(FOI)
}


# Not use.
# -----------------------------------------------------------------------------
# Calculate the unweighted group mixing matix for each risk group in each sex.
mixG <- function(cArray, pop, i, mx, fp) {
  ages <- fp$ss$AGE_START:fp$ss$AGE_END 
  nRg  <- length(mx$gamma)
	kid <- expand.grid(1:nRg, 1:nRg)
  for (k in 1:2) {
    for (j in 1:(nRg^2)) {
      cArray[,,j,k,1] <- outer(ages, ages, 
        function(x, y) Cmix(k, x, y, kid[j,1], kid[j,2], pop, i, mx, fp))
    }   
  }
  return(cArray)
}

# Weighting the mix G matrices, seemed expandable to more than 2 groups
# See BackupFOI for step by step calculation process for two groups
mixG.w <- function(cArray, i, pop, fp, mx) {
	nRg <- length(mx$gamma)
	# special index to lookup which the opposite sex group position
	iid <- rep(1:2, each=nRg) + nRg * rep(0:(nRg-1), 2) 
	kid <- expand.grid(1:nRg, 1:nRg)
	for (q in 1:2) { # loop through sex
		for (j in 1:(nRg^2)) {
			p <- ifelse(q==1, 2, 1)
			popA <- rowSums(pop[, sid(q,fp)[kid[j,2]], , i])
			popB <- rowSums(pop[, sid(p,fp)[kid[j,1]], , i])
			cArray[,,j,q,2] <- sweepX(cArray[,,j,q,1], 1, popA) / 
						     t(sweepX(cArray[,,iid[j],p,1], 1, popB))
		}
	}
	return(cArray)
}

# Repeat some calculation, should be simplified
mixRisk <- function(cArray, pF, fp, mx) {
	nRg <- length(mx$gamma)
	# special index to lookup which the opposite sex group position
	iid <- rep(1:2, each=nRg) + nRg * rep(0:(nRg-1), 2) 
	kid <- expand.grid(1:nRg, 1:nRg)
	for (q in 1:2) { # loop through sex 
		for (j in 1:(nRg^2)) { # through all sex * group
			p <- ifelse(q==1, 2, 1)
			cArray[,,j,q,1] <- sweepX(cArray[,,j,q,1], 2, pF[,sid(p, fp)[kid[j,1]]]) * (1 - mx$condom[,,q])
		}
	}
	return(cArray)
}