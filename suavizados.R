corr_abbott <- function(df, control = 0){
  ctrl <- df[df$dosis == control, ]
  Pctrl <- ctrl$suavizados
  Presp <- df[df$suavizados != 0, ]  #datos????
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
# Control     30          4      0.1333        0.0000         0.1006
# 4.0000      30          2      0.0667        -.0378         0.0971
# 4.5000      30         12      0.4000        0.3329         0.2565
# 5.0000      30         19      0.6333        0.5923         0.4688
# 5.5000      30         20      0.6667        0.6294         0.6711
# 6.0000      30         24      0.8000        0.7776         0.8209

datos_suav<- suavizar(datos)

data_abbot <- corr_abbott(datos_suav)


library(dplyr)

data <- data_abbot %>% filter(dosis!=0)

model <- suppressWarnings(glm(Pabbott ~ log10(dosis),
             family = binomial(link = "probit"),
             weights = n,
             data = data))


#use model to predict value of am
pred <- predict(model, type="response")
data
summary(lm(pred~data$Presp))
newdata = data.frame(respuesta=0)


#  Deviance Residuals:
quantile(residuals(model))


# Calculo del CHI2 con p-valor
vcdExtra::LRstats(model)
chi2 <- model$deviance
print(chi2)
pchisq(chi2, model$df.residual, lower.tail = FALSE)


summary(model)
# Calcular LD50 y su intervalo de confianza
intercept <- stats::coef(model)[1]
slope <- stats::coef(model)[2]
log_ld50 <- (0 - intercept) / slope
print(log_ld50)
ld50 <- 10^log_ld50
print(ld50)


confint(model)
# Asegúrate de cargar el paquete MASS
library(MASS)

# Generar secuencia de probabilidades del 1% al 99%
probabilidades <- seq(0.0, 0.99, by = 0.01)

# Función para calcular DLp y su intervalo de confianza
calcular_dlp <- function(model, p) {
  dl <- MASS::dose.p(model, p = p)        # Calcula en escala log10
  log_dl <- dl                   # Valor en log10(dosis)
  se <- attr(dl, "SE")             # Error estándar en log10

  # Convertir a escala original (dosis)
  dl_estimada <- 10^log_dl
  li <- 10^(log_dl - 1.96 * se)    # Límite inferior (95% CI)
  ls <- 10^(log_dl + 1.96 * se)    # Límite superior (95% CI)

  data.frame(
    Probabilidad = p,
    DL = dl_estimada,
    LI = li,
    LS = ls
  )
}

# Aplicar la función a todas las probabilidades
resultados <- lapply(probabilidades, function(p) calcular_dlp(model, p))
resultados_df <- do.call(rbind, resultados)

# Mostrar resultados (ejemplo)
(resultados_df)

# Error estándar usando método delta
vcov_matrix <- stats::vcov(model)
var_log_ld50 <- (vcov_matrix[1,1] + vcov_matrix[2,2]*(log_ld50^2) +
                   2*log_ld50*vcov_matrix[1,2]) / (slope^2)
se_log_ld50 <- sqrt(var_log_ld50)
print(se_log_ld50)
# Intervalo de confianza para LD50
z <- stats::qnorm(1 - (1 - 0.95)/2)
ci_lower <- 10^(log_ld50 - z * se_log_ld50)
print(ci_lower)
ci_upper <- 10^(log_ld50 + z * se_log_ld50)
print(ci_upper)

# Calcular pendiente en el punto medio (LD50)
slope_midpoint <- slope * stats::dnorm(0)  # Derivada de la función Probit en 0
print(slope_midpoint)
se_slope <- sqrt(vcov_matrix[2,2]) * stats::dnorm(0)
print(slope)

# Estadísticos de bondad de ajuste
predicted <- stats::predict(model, type = "response") #* data$n
predicted

# Estadísticos de bondad de ajuste

null_deviance <- model$null.deviance
df_null <- model$df.null
deviance <- model$deviance
print(deviance)
df_deviance <- model$df.residual
# https://www.statology.org/null-residual-deviance/
chi_sq <- null_deviance - deviance
print(chi_sq)
df <- df_null - df_deviance
p_value <- stats::pchisq(chi_sq, df, lower.tail = FALSE)
print(p_value)

chi_crit <-  stats::qchisq(p = 1 - 0.05, df = df_deviance)
print(chi_crit)
aic <- stats::AIC(model)



# Coeficientes del modelo
b0 <- coef(model)[1]  # Intercepto
b1 <- coef(model)[2]  # Pendiente (coeficiente de log10(dosis))
vcov_matrix <- vcov(model)  # Matriz de varianza-covarianza

# Secuencia de probabilidades (1% a 99%)
probabilidades <- c(0.01, 0.05, 0.1, 0.15, 0.5, 0.85, 0.9, 0.95, 0.99)

# Función para calcular DLp y su IC
calcular_dlp_manual <- function(p) {
  # Valor de probit para la probabilidad p
  probit_p <- stats::qnorm(p)

  # Calcular log10(DLp)
  log_dlp <- (probit_p - b0) / b1

  # Aplicar método delta para el error estándar
  gradiente <- c(
    -1 / b1,                     # Derivada respecto al intercepto (b0)
    -(probit_p - b0) / (b1^2)    # Derivada respecto a la pendiente (b1)
  )

  # Calcular varianza y error estándar
  var_log_dlp <- t(gradiente) %*% vcov_matrix %*% gradiente
  se_log_dlp <- sqrt(var_log_dlp)

  # Convertir a escala original y crear dataframe
  data.frame(
    Probabilidad = p,
    DL = 10^log_dlp,
    LI = 10^(log_dlp - 1.96 * se_log_dlp),
    LS = 10^(log_dlp + 1.96 * se_log_dlp)
  )
}

# Calcular para todas las probabilidades
resultados <- lapply(probabilidades, calcular_dlp_manual)
resultados_df <- do.call(rbind, resultados)

# Mostrar resultados
rownames(resultados_df) <- NULL
resultados_df$Probabilidad <- resultados_df$Probabilidad
(resultados_df)



plot(qnorm(resultados_df$Probabilidad), log10(resultados_df$DL), type = "l")
points(qnorm(data$Presp), log10(data$dosis), pch = 19, col = "red")
