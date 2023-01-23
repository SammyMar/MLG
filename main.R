########## IMPORTACAO E PACOTES ##############
library(corrplot)
library(dplyr)
################

mat <- read.csv('dados/Maths.csv')
pot <- read.csv('dados/Portuguese.csv')

mat['subject'] <- 'Math'
pot['subject'] <- 'Port'

df.complete <- rbind(mat,pot)

notas <- df.complete[,c('G1','G2','G3')]

corrplot(cor(notas), method = 'color', order = 'alphabet')


df.complete[sapply(df.complete, is.character)] <- lapply(df.complete[sapply(df.complete, is.character)],
                                       as.factor)
str(df.complete)
# ABORDAGEM 1 - 3 MODELOS, UM PARA CADA NOTA COM MATERIA COMO COVARIAVEL
# CONSIDERE UM MODELO BIN(N,P), ONDE N = NUMERO MAXIMO DE QUESTOES
# E UM P = NOTA / N
nrow(df.complete)
summary(df.complete)
#G1
g1 <- df.complete |>
  select(-c('G2',"G3")) |>
  mutate(
    n = 20,
    p = G1/n
  )
y <- cbind(g1$G1, 20 - g1$G1)

covariaveis <- g1 |> select(-c('n','p','G1'))

glm(y ~ . ,family = binomial,data = covariaveis) |> summary()




g2 <- df.complete |>
  select(-c('G1',"G3"))|>
  mutate(
    n = 20,
    p = G2/n
  )
g3 <- df.complete |>
  select(-c('G2',"G1"))|>
  mutate(
    n = 20,
    p = G3/n
  )
glm()
# ABORDAGEM 2 - 1 MODELO COM 3 LINHAS PARA CADA INDIVIDUO INDICANDO CADA UMA DAS NOTAS
# COM MATERIA COMO COVARIAVEL





