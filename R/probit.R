probit <- function(data) {

  datos_suav<- suavizar(data)
  data_abbot <- corr_abbott(datos_suav)

  datos <- data_abbot %>% dplyr::filter(data_abbot$dosis!=0)

  model <- stats::glm(Pabbott ~ log10(datos$dosis),
               family = stats::quasibinomial(link = "probit"),
               weights = datos$n,
               data = datos)


  summary(model)
  # Calcular LD50 y su intervalo de confianza
  intercept <- stats::coef(model)[1]
  slope <- stats::coef(model)[2]
  log_ld50 <- (0 - intercept) / slope
  ld50 <- 10^log_ld50



  stats::confint(model)


  # Error estandar usando método delta
  vcov_matrix <- stats::vcov(model)
  var_log_ld50 <- (vcov_matrix[1,1] + vcov_matrix[2,2]*(log_ld50^2) +
                     2*log_ld50*vcov_matrix[1,2]) / (slope^2)
  se_log_ld50 <- sqrt(var_log_ld50)
  # Intervalo de confianza para LD50
  z <- stats::qnorm(1 - (1 - 0.95)/2)
  ci_lower <- 10^(log_ld50 - z * se_log_ld50)
  ci_upper <- 10^(log_ld50 + z * se_log_ld50)

  # Calcular pendiente en el punto medio (LD50)
  slope_midpoint <- slope * stats::dnorm(0)  # Derivada de la funcion Probit en 0
  se_slope <- sqrt(vcov_matrix[2,2]) * stats::dnorm(0)

  # Estadisticos de bondad de ajuste
  predicted <- stats::predict(model, type = "response") * datos$n


  # Estadisticos de bondad de ajuste

  null_deviance <- model$null.deviance
  df_null <- model$df.null
  deviance <- model$deviance
  df_deviance <- model$df.residual
  # https://www.statology.org/null-residual-deviance/
  chi_sq <- null_deviance - deviance
  df <- df_null - df_deviance
  p_value <- stats::pchisq(chi_sq, df, lower.tail = FALSE)

  chi_crit <-  stats::qchisq(p = 1 - 0.05, df = df_deviance)

  aic <- stats::AIC(model)



  # Coeficientes del modelo
  b0 <- stats::coef(model)[1]  # Intercepto
  b1 <- stats::coef(model)[2]  # Pendiente (coeficiente de log10(dosis))
  vcov_matrix <- stats::vcov(model)  # Matriz de varianza-covarianza

  # Secuencia de probabilidades
  probabilidades <- c(0.01, 0.05, 0.1, 0.15, 0.5, 0.85, 0.9, 0.95, 0.99)

  # Funcion para calcular DLp y su IC
  calcular_dlp_manual <- function(p) {
    # Valor de probit para la probabilidad p
    probit_p <- stats::qnorm(p)

    # Calcular log10(DLp)
    log_dlp <- (probit_p - b0) / b1

    # Aplicar metodo delta para el error estandar
    gradiente <- c(
      -1 / b1,                     # Derivada respecto al intercepto (b0)
      -(probit_p - b0) / (b1^2)    # Derivada respecto a la pendiente (b1)
    )

    # Calcular varianza y error estandar
    var_log_dlp <- t(gradiente) %*% vcov_matrix %*% gradiente
    se_log_dlp <- sqrt(var_log_dlp)

    # Convertir a escala original y crear dataframe
    data.frame(
      Punto = p,
      Concentracion = 10^log_dlp,
      LCI95 = 10^(log_dlp - 1.96 * se_log_dlp),
      LCS95 = 10^(log_dlp + 1.96 * se_log_dlp)
    )
  }

  # Calcular para todas las probabilidades
  resultados <- lapply(probabilidades, calcular_dlp_manual)
  resultados_df <- do.call(rbind, resultados)

  # Mostrar resultados
  rownames(resultados_df) <- NULL

  (resultados_df)



  #plot(qnorm(resultados_df$Probabilidad), log10(resultados_df$DL), type = "l")
  #points(qnorm(data$Presp), log10(data$dosis), pch = 19, col = "red")
}
