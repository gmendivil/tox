suavizar <- function(x) {
  n <- nrow(x)
  z <- x$respuesta / x$n
  if (n <= 1) return(x)  # No hay nada que suavizar

  for (i in 2:n) {
    if (z[i] < z[i - 1]) {  # Verifica si hay un incremento no permitido
      # Calcula el promedio de todos los elementos desde 1 hasta i
      avg <- mean(z[1:i])
      # Asigna el promedio a todos los elementos en la ventana
      z[1:i] <- avg
    }
  }
  x$suavizados <- z
  return(x)
}

