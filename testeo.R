library(tox)
library(dplyr)



base  <- database_tox   # Se cargan datos de ejemplo del paquete tox

datos <- base %>% filter(especie=="Boana pulchella",
               exposicion_horas==96,compuesto=="Glufosinato",
               fecha != "2022-05-31") %>%
               group_by(dosis) %>%
               summarise(n=sum(n), respuesta=sum(respuesta))

datos2 <- base %>% filter(especie=="Boana pulchella",
               exposicion_horas==96,compuesto=="Lambdacialotrina",
               fecha != "2022-05-31") %>%
              group_by(dosis) %>%
              summarise(n=sum(n), respuesta=sum(respuesta))


modelo1 <- probit(datos)
modelo2 <- probit(datos2)

plot_probit_modelos(modelo1, modelo2)

#plot(qnorm(resultados_df$Probabilidad), log10(resultados_df$DL), type = "l")
#points(qnorm(data$Presp), log10(data$dosis), pch = 19, col = "red")
