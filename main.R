########## IMPORTACAO E PACOTES ##############

mat <- read.csv('dados/Maths.csv')
pot <- read.csv('dados/Portuguese.csv')

mat['subject'] <- 'Math'
pot['subject'] <- 'Port'

df.complete <- rbind(mat,pot)
