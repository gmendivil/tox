library(tox)
library(dplyr)

datos  <- database_tox  # Se cargan datos de ejemplo del paquete tox

datos <- datos%>% filter(especie=="Cnesterodon decemmaculatus",
               exposicion_horas==96,compuesto=="Glufosinato",
               fecha=="2022-10-24") %>%
               group_by(dosis) %>%
               summarise(n=sum(n), respuesta=sum(respuesta))



probit(datos)


#plot(qnorm(resultados_df$Probabilidad), log10(resultados_df$DL), type = "l")
#points(qnorm(data$Presp), log10(data$dosis), pch = 19, col = "red")
