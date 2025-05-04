probit <- function(data, control_dosis = 0, conf.level = 0.95) {
#' Función para análisis Probit con correcciones
#'
#' @param data Dataframe con columnas "dosis" - "n" - "respuesta".
#' @param control_dosis Factor que representa el control si lo hay. Ej: 0, C, control
#' @param conf.level Nivel de confianza
#'
#' @return Una lista.
#' @export
#'
#' @examples
#' datos <- importar(data)
#' probit(datos)

        # Verificar y corregir datos
  if (!all(c("dosis", "n", "respuesta") %in% colnames(data))) {
    stop("El dataframe debe contener las columnas: dosis, n, respuesta")
  }


  # Corrección de Abbott (mortalidad natural)
  if (control_dosis %in% data$dosis) {
    control <- data[data$dosis == control_dosis, ]
    p_control <- control$respuesta / control$n
    data <- data[data$dosis != control_dosis, ]
    data$respuesta<- round((data$respuesta/data$n - p_control) / (1 - p_control) * data$n)
    data$respuesta <- pmax(0, pmin(data$respuesta, data$n))  # Asegurar valores entre 0 y n
  }

  data <- data %>% dplyr::filter(dosis != 0)

  # Corrección de Litchfield-Wilcoxon para 0 y 1
  data$r_adj <- data$respuesta+ 0.5
  data$n_adj <- data$n + 1
  data$p <- data$r_adj / data$n_adj

  # Ajustar modelo Probit
  model <- stats::glm(cbind(r_adj, n_adj - r_adj) ~ log10(dosis),
               family = stats::binomial(link = "probit"),
               data = data)

  # Calcular LD50 y su intervalo de confianza
  intercept <- stats::coef(model)[1]
  slope <- stats::coef(model)[2]
  log_ld50 <- (0 - intercept) / slope
  ld50 <- 10^log_ld50

  # Error estándar usando método delta
  vcov_matrix <- stats::vcov(model)
  var_log_ld50 <- (vcov_matrix[1,1] + vcov_matrix[2,2]*(log_ld50^2) +
                     2*log_ld50*vcov_matrix[1,2]) / (slope^2)
  se_log_ld50 <- sqrt(var_log_ld50)

  # Intervalo de confianza para LD50
  z <- stats::qnorm(1 - (1 - conf.level)/2)
  ci_lower <- 10^(log_ld50 - z * se_log_ld50)
  ci_upper <- 10^(log_ld50 + z * se_log_ld50)

  # Calcular pendiente en el punto medio (LD50)
  slope_midpoint <- slope * stats::dnorm(0)  # Derivada de la función Probit en 0
  se_slope <- sqrt(vcov_matrix[2,2]) * stats::dnorm(0)

  # Estadísticos de bondad de ajuste
  predicted <- stats::predict(model, type = "response") * data$n_adj


  # Estadísticos de bondad de ajuste

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

  # Resultados
  list(
    model = model,
    null_deviance = null_deviance,
    df_null = df_null,
    deviance = deviance,
    df_deviance = df_deviance,
    chi_crit = chi_crit,
    LD50 = ld50,
    LD50_CI = c(ci_lower, ci_upper),
    LD50_SE = se_log_ld50,
    Slope = slope_midpoint,
    Slope_SE = se_slope,
    ChiSq = c(Statistic = chi_sq,
              Df = df,
              Pvalue = p_value),
    AIC = aic,
    CorrectedData = data
  )
}
