agregar <- function(df){
  #' @title Genera un dataframe con los datos agregados.
  #'
  #' @description Suma los valores de respuesta de las réplicas y los n's.
  #'
  #' @param df Dataframe con los resultados de todas las rèplicas. El mismo debe contener las columnas dosis , respuesta y n.
  #'
  #' @return Dataframe con los datos agregados.
  #'
  #' @examples
  #' datos_agregados <- agregar(df)

  datos_agregados <- df %>%
    group_by(dosis) %>%
    summarise(respuesta = sum(respuesta),
              n = sum(n))
  return(datos_agregados)
}


corr_abbott <- function(df, control = 0){
  Pctrl <- df[df$dosis == control, ]
  Pctrl <- Pctrl$respuesta / Pctrl$n
  #Presp <- datos[datos$dosis != 0, ]
  df$Presp <- df$respuesta / df$n
  Pcorr <- (df$Presp-Pctrl)/(1- Pctrl )
  df$Pcorr <- Pcorr
  df$Pcorr <- ifelse(df$Pcorr < 0,  (0.5/(df$n+1)), df$Pcorr)
  return(df)
}

datos <- agregar(df_24_octubre)
datos <- corr_abbott(datos)



modelo <- glm(Pcorr ~ log(dosis),
              family = binomial(link = "probit"),
              weights = n,
              data = datos[datos$dosis > 0, ])

summary(modelo)

predichos <- predict(modelo, type = "response")
predichos
