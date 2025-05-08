corr_abbott <- function(df, control = 0){
  ctrl <- df[df$dosis == control, ]
  Pctrl <- ctrl$suavizados
  Presp <- datos[datos$dosis != 0, ]
  df$Presp <- df$respuesta / df$n
  Pabbott <- (df$Presp-Pctrl)/(1- Pctrl )
  df$Pabbott <- Pabbott
  df[1,ncol(df)] <- 0
  df$Pabbott <- ifelse(df$Pabbott < 0,  (0.5/(df$n+1)), df$Pabbott)
  return(df)
}

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

# Datos con log(conc) para tratamientos
datos <- data.frame(
  dosis = c(0, 4.0, 4.5, 5.0, 5.5, 6.0),
  respuesta = c(4,2, 12, 19, 20, 24),
  n = rep(30, 6)
)

datos_suav<- suavizar(datos)

data_abbot <- corr_abbott(datos_suav)
library(dplyr)

data <- data_abbot %>% filter(dosis!=0)

model <- glm(Pabbott ~ log10(dosis),
             family = binomial(link = "probit"),
             data = data)
summary(model)
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
z <- stats::qnorm(1 - (1 - 0.95)/2)
ci_lower <- 10^(log_ld50 - z * se_log_ld50)
ci_upper <- 10^(log_ld50 + z * se_log_ld50)

# Calcular pendiente en el punto medio (LD50)
slope_midpoint <- slope * stats::dnorm(0)  # Derivada de la función Probit en 0
se_slope <- sqrt(vcov_matrix[2,2]) * stats::dnorm(0)

# Estadísticos de bondad de ajuste
predicted <- stats::predict(model, type = "response") * data$n


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
