############ funcoes #####################

envelope.poi <- function(fit.model){
par(mfrow=c(1,1))
X <- model.matrix(fit.model)
n <- nrow(X)
p <- ncol(X)
w <- fit.model$weights
W <- diag(w)
H <- solve(t(X)%*%W%*%X)
H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
h <- diag(H)
td <- resid(fit.model,type="deviance")/sqrt((1-h))
e <- matrix(0,n,100)
#
for(i in 1:100){
  nresp <- rpois(n, fitted(fit.model))
  fit <- glm(nresp ~ X, family=poisson)
  w <- fit$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  e[,i] <- sort(resid(fit,type="deviance")/sqrt(1-h))}
#
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:n){
  eo <- sort(e[i,])
  e1[i] <- (eo[2]+eo[3])/2
  e2[i] <- (eo[97]+eo[98])/2}
#
med <- apply(e,1,mean)
faixa <- range(td,e1,e2)
par(pty="s")
qqnorm(td,xlab="Percentil da N(0,1)",
       ylab="Componente do Desvio", ylim=faixa, pch=16)
par(new=T)
#
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2) }
envelope.bn <- function(fit.model){
  par(mfrow=c(1,1))
  X <- model.matrix(fit.model)
  n <- nrow(X)
  p <- ncol(X)
  fi <- fit.model$theta
  w <- fi*fitted(fit.model)/(fi + fitted(fit.model))
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  td <- resid(fit.model,type="deviance")/sqrt(1-h)
  fi <- fit.model$theta
  e <- matrix(0,n,100)
  #
  for(i in 1:100){
    resp <- rnegbin(n, fitted(fit.model),fi)
    fit <- glm.nb(resp ~ X, control = glm.control(maxit = 50))
    w <- fit$weights
    W <- diag(w)
    H <- solve(t(X)%*%W%*%X)
    H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h <- diag(H)
    e[,i] <- sort(resid(fit,type="deviance")/sqrt(1-h))}
  #
  e1 <- numeric(n)
  e2 <- numeric(n)
  #
  for(i in 1:n){
    eo <- sort(e[i,])
    e1[i] <- (eo[2]+eo[3])/2
    e2[i] <- (eo[97]+eo[98])/2}
  #
  med <- apply(e,1,mean)
  faixa <- range(td,e1,e2)
  par(pty="s")
  qqnorm(td,xlab="Percentil da N(0,1)",
         ylab="Componente do Desvio", ylim=faixa, pch=16)
  par(new=T)
  #
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2)
}

########## IMPORTACAO E PACOTES ##############

library(corrplot)
library(dplyr)
library(ggplot2)

##############################################

mat <- read.csv('dados/Maths.csv')
pot <- read.csv('dados/Portuguese.csv')
# selecao manual de variaveis #################
mat |> summary()
mat |> colnames()
mat[,mat |> sapply(is.character) | mat|>sapply(is.integer)] |> colnames()
mat[sapply(mat, is.character)] <- lapply(mat[sapply(mat, is.character)],
                                       as.factor)
y <- c('absences')
x <- c('school','sex','age','famsize','famsize','Medu','Fedu','reason','guardian','traveltime',
       'failures','schoolsup','famsup','activities','higher','internet','romantic','famrel','freetime','goout',
       'Dalc')
##  exploratoria ############################



mat |> ggplot() +
    aes(school, fill = school) + geom_bar(show.legend = F) + facet_wrap(vars(sex))

mat |> ggplot() +
  aes(age) + geom_bar() + facet_wrap(vars(failures))

mat |> ggplot() +
  aes(y, fill =guardian) + geom_bar(position = 'dodge')
mat |> ggplot() +
  aes(guardian) + geom_bar()
mat[,y] |> max()

mat <- mat[,c(x,y)]

### selecao de variaveis pelo step ###############################
glm(absences ~. ,mat, family='poisson') |> step(direction = 'backward')

fit1 <- glm(formula = absences ~ school + sex + age + famsize + Medu +
      reason + guardian + traveltime + schoolsup + higher + internet +
      romantic + famrel + freetime + Dalc, family = "poisson",
    data = mat)

envelope.poi(fit1)

######### modelo binomial negativo #############################
library(MASS)
fit2 <- glm.nb(formula = absences ~ school + sex + age + famsize + Medu +
      reason + guardian + traveltime + schoolsup + higher + internet +
      romantic + famrel + freetime + Dalc,
    data = mat)

fit2 |> summary()

glm.nb(absences ~ ., data = mat ) |> step(direction = 'backward')

fit3 <- glm.nb(absences ~ school + sex + age + Medu + reason +
         internet + famrel + Dalc, data = mat, init.theta = 0.7000362531,
       link = log)

fit3 |> summary()

fit3.inter <- glm.nb(absences ~ (school + sex + age + Medu + reason +
                 internet + famrel + Dalc)^2, data = mat, init.theta = 0.7000362531,
               link = log)
fit3 |> envelope.bn()
