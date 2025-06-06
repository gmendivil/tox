agregar <- function(df){
  #' @title Genera un dataframe con los datos agregados.
  #'
  #' @description Suma los valores de respuesta de las réplicas y los n's.
  #'
  #' @param df Dataframe con los resultados de todas las rèplicas. El mismo debe contener las columnas dosis , respuesta y n.
  #'
  #' @return Dataframe con los datos agregados.
  #'
  # @examples
  # datos_agregados <- agregar(data) CORREGIR!

  datos_agregados <- df %>%
    dplyr::group_by(dosis) %>%
    dplyr::summarise(n=sum(n), respuesta=sum(respuesta))

  return(datos_agregados)
}

