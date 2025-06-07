library(tox)
library(dplyr)

datos  <- database_tox  # Se cargan datos de ejemplo del paquete tox

Lambda <- datos%>% filter(especie=="Boana pulchella",
               exposicion_horas==96,compuesto=="Lambdacialotrina")

Glufo <- datos%>% filter(especie=="Boana pulchella",
                exposicion_horas==96,compuesto=="Glufosinato")


modelo_lambda <- tox(Lambda)
modelo_glufo <- tox(Glufo)

plot_probit_modelos(modelo_glufo, modelo_lambda, model_names = c("Glufosinato", "Lambdacialotrina"))
