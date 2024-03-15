nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  Y <- -0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)] + rnorm(n, 0, 2)
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h")
  
  pval <- c(pval, gomp_lm(formula('expr ~ a'), NULL, dat_sub, character(0), c("f", "h"), character(0))[,3])

  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}

# These look right:
# 0, 1 with 1 constraint under null
# 0, 2 with 2 constraints under null
# 0, 3 with 3 constraints under null

nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  Y <- rbinom(n, 1, prob = plogis(-0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)]))
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h")
  w <- runif(n, 0.7, 1)
  
  pval <- c(pval, gomp_glm(formula('expr ~ a'), NULL, dat_sub, character(0), c("f", "h"), character(0))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}

# https://rpubs.com/bbolker/glmerconstr
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  
  ngroups <- 10
  X$group <- sample(letters[1:ngroups], n, replace = T)
  group_effects <- rnorm(ngroups, mean = 0, sd = 0.2)
  
  Y <- rbinom(n, 1, prob = plogis(-0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)]))
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h", "group")
  
  pval <- c(pval, gomp_glmer(formula('expr ~ a + (1|group)'), '(1|group)', dat_sub, character(0), c("f", "h"), character(0))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}

#LMER
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  
  ngroups <- 10
  X$group <- sample(letters[1:ngroups], n, replace = T)
  group_effects <- rnorm(ngroups, mean = 0, sd = 0.2)
  
  Y <- -0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)] + rnorm(n, 0, 2)
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h", "group")
  
  pval <- c(pval, gomp_lmer(formula('expr ~ a + (1|group)'), '(1|group)', dat_sub, character(0), c("f", "h"), character(0))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}

########
# OMPS #
########

# LM
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  Y <- -0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)] + rnorm(n, 0, 2)
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h")
  
  pval <- c(pval, omp_lm(formula('expr ~ a'), NULL, dat_sub, character(0), character(0), c("f", "h"))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval < 0.05), 4)))
    hist(pval, breaks = seq(0, 1, 0.02))
  }
}

# GLM
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 1000
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  Y <- rbinom(n, 1, prob = plogis(-0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)]))
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h")
  
  pval <- c(pval, omp_glm_augment(formula('expr ~ a'), NULL, dat_sub, character(0), character(0), c("f", "h"))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval < 0.05), 4)))
    hist(pval, breaks = seq(0, 1, 0.02))
  }
}

# LMER
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  
  ngroups <- 10
  X$group <- sample(letters[1:ngroups], n, replace = T)
  group_effects <- rnorm(ngroups, mean = 0, sd = 0.2)
  
  Y <- -0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)] + rnorm(n, 0, 2)
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h", "group")
  
  pval <- c(pval, omp_lmer(formula('expr ~ a + (1|group)'), '(1|group)', dat_sub, character(0), character(0), c("f", "h"))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}

# GLMER
nsims <- 2000
pval <- c()
for (i in 1:nsims) {
  n <- 40
  X <- data.frame(matrix(runif(1 * n), ncol=1))
  colnames(X) <- 'a'
  X$f <- factor(sample(letters[1:3], n, replace = T), levels = c('c', 'b', 'a'))
  X$h <- factor(sample(letters[1:4], n, replace = T))
  group_effects_f <- c(0,0,0)
  names(group_effects_f) <- letters[1:3]
  group_effects_h <- c(0,0,0,0)
  names(group_effects_h) <- letters[1:4]
  
  ngroups <- 10
  X$group <- sample(letters[1:ngroups], n, replace = T)
  group_effects <- rnorm(ngroups, mean = 0, sd = 0.2)
  
  Y <- rbinom(n, 1, prob = plogis(-0.5 + as.matrix(X[,-c(2:4), drop=F]) %*% c(1) + group_effects_h[as.character(X$h)] + group_effects_f[as.character(X$f)]))
  dat_sub <- cbind(Y, X)
  colnames(dat_sub) <- c("expr", "a", "f", "h", "group")
  
  pval <- c(pval, omp_glmer_augment(formula('expr ~ a + (1|group)'), '(1|group)', dat_sub, character(0), character(0), c("f", "h"))[,3])
  
  if (i %% 10 == 0) { 
    print(paste0(i, ": ", round(mean(pval[!is.na(pval) & pval > 0] < 0.05), 4)))
    hist(pval[pval > 0], breaks = seq(0, 1, 0.02))
  }
}





