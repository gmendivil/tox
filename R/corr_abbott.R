# corr_abbott <- function(df, control = 0){
#   Pctrl <- df[df$dosis == control, ]
#   Pctrl <- Pctrl$respuesta / Pctrl$n
#   #Presp <- datos[datos$dosis != 0, ]
#   df$Presp <- df$respuesta / df$n
#   Pcorr <- (df$Presp-Pctrl)/(1- Pctrl )
#   df$Pcorr <- Pcorr
#   df$Pcorr <- ifelse(df$Pcorr < 0,  (0.5/(df$n+1)), df$Pcorr)
#   return(df)
# }

corr_abbott <- function(df, control = 0){
  ctrl <- df[df$dosis == control, ]
  Pctrl <- ctrl$suavizados
  Presp <- df[df$dosis != 0, ]
  df$Presp <- df$respuesta / df$n
  Pabbott <- (df$Presp-Pctrl)/(1- Pctrl )
  df$Pabbott <- Pabbott
  df[1,ncol(df)] <- 0
  df$Pabbott <- ifelse(df$Pabbott < 0,  (0.5/(df$n+1)), df$Pabbott)
  return(df)
}
