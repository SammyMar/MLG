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
diagnostico.bn <- function(fit.model){
  X <- model.matrix(fit.model)
  n <- nrow(X)
  p <- ncol(X)
  fi <- fit.model$theta
  w <- fi*fitted(fit.model)/(fi + fitted(fit.model))
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  ts <- resid(fit.model,type="pearson")/sqrt(1-h)
  td <- resid(fit.model,type="deviance")/sqrt(1-h)
  di <- (h/(1-h))*(ts^2)
  par(mfrow=c(2,2))
  a <- max(td)
  b <- min(td)
  plot(fitted(fit.model),h,xlab="Valores Ajustados", ylab="Medida h",
       pch=16)
  identify(fitted(fit.model), h, n=5)
  title(sub="(a)")
  #
  plot(di,xlab="Indice", ylab="Distancia de Cook", pch=16)
  identify(di,n=3)
  title(sub="(b)")
  #
  plot(td,xlab="Indice", ylab="Residuo Componente do Desvio",
       ylim=c(b-1,a+1), pch=16)
  abline(2,0,lty=2)
  abline(-2,0,lty=2)
  identify(td,n=1)
  title(sub="(c)")
  #
  eta = predict(fit.model)
  z = eta + resid(fit.model, type="pearson")/sqrt(w)
  plot(predict(fit.model),z,xlab="Preditor Linear",
       ylab="Variavel z", pch=16)
  lines(smooth.spline(predict(fit.model), z, df=2))
  title(sub="(d)")
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
mat$reason |> unique()
mat[sapply(mat, is.character)] <- lapply(mat[sapply(mat, is.character)],
                                       as.factor)
y <- c('absences')
x <- c('school','sex','age','famsize','famsize','Medu','Fedu','reason','guardian','traveltime',
       'failures','schoolsup','famsup','activities','higher','internet','romantic','famrel','freetime','goout',
       'Dalc')
##  exploratoria ############################

str(mat)

mat |> ggplot() +
  aes(x = absences, y = internet) +
  geom_boxplot(fill = "darkblue") +
  theme_minimal()
mat |> ggplot() +
  aes(x = absences, y = internet) +
  geom_boxplot(fill = "darkblue") +
  theme_minimal() + facet_wrap(vars(sex))
ggplot(mat) +
  aes(x = internet, weight = absences) +
  geom_bar(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(sex))
mat |> ggplot() +
    aes(school, fill = school) + geom_bar(show.legend = F) + facet_wrap(vars(sex))

mat |> ggplot() +
  aes(age) + geom_bar() + facet_wrap(vars(failures))

mat |> ggplot() +
  aes(guardian, weight =absences) + geom_bar(position = 'dodge')

mat |> ggplot() +
  aes(guardian) + geom_bar()
mat |> ggplot() +
  aes(x = absences, y = guardian) +
  geom_boxplot(fill = "darkblue") +
  theme_minimal()
mat |> ggplot() +
  aes(x = reason, fill = absences) +
  geom_bar() +
  scale_fill_gradient() +
  theme_minimal()

 mat |> ggplot() +
  aes(absences, fill =school) + geom_bar(position = 'dodge')

table(mat$school) # NAO SEI SE VALE A PENA INCLUIR ESCOLA POR CAUDA DA DIFERENCA DE OBSERVACOES PRA CADA UMA
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
fit3.inter |> summary()
fit4 <- glm.nb(absences ~ school + sex + age + Medu + reason +
                 internet  + Dalc, data = mat, init.theta = 0.7000362531,
               link = log)
fit4 |> summary()
anova(fit4,fit3,test = 'LR')

fit5 <-  glm.nb(absences ~   sex + age + Medu + reason +
                  internet  + Dalc, data = mat, init.theta = 0.7000362531,
                link = log)
fit5 |> summary()
anova(fit5,fit4)
diagnostico.bn(fit4)

