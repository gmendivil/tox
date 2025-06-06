#' Correcion de Abbott para datos de mortalidad natural
#'
#' @param df data frame
#' @param control Variable que indica como se codifica la columna control (Por ej.: 0, "control", etc)
#'
#' @returns retorna el data frame corregido
#' @export
#'

corr_abbott <- function(df, control = 0){
  ctrl <- df[df$dosis == control, ]
  if (ctrl$respuesta > 0) { # La condición debe devolver TRUE o FALSE
    Pctrl <- ctrl$suavizados# Ejecuta un código
  } else {
    Pctrl <- 0# Ejecuta otro código
  }

  Presp <- df[df$dosis != 0, ]  #datos????
  df$Presp <- df$respuesta / df$n
  Pabbott <- (df$Presp-Pctrl)/(1- Pctrl )
  df$Pabbott <- Pabbott
  df[1,ncol(df)] <- 0
  df$Pabbott <- ifelse(df$Pabbott < 0,  (0.5/(df$n+1)), df$Pabbott)
  return(df)
}
