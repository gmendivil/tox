# corr_abbott <- function(df, control = 0){
#   ctrl <- df[df$dosis == control, ]
#   Pctrl <- ctrl$suavizados
#   Presp <- datos[datos$dosis != 0, ]
#   df$Presp <- df$respuesta / df$n
#   Pabbott <- (df$Presp-Pctrl)/(1- Pctrl )
#   df$Pabbott <- Pabbott
#   df[1,ncol(df)] <- 0
#   df$Pabbott <- ifelse(df$Pabbott < 0,  (0.5/(df$n+1)), df$Pabbott)
#   return(df)
# }

# suavizar <- function(x) {
#   n <- nrow(x)
#   z <- x$respuesta / x$n
#   if (n <= 1) return(x)  # No hay nada que suavizar
#
#   for (i in 2:n) {
#     if (z[i] < z[i - 1]) {  # Verifica si hay un incremento no permitido
#       # Calcula el promedio de todos los elementos desde 1 hasta i
#       avg <- mean(z[1:i])
#       # Asigna el promedio a todos los elementos en la ventana
#       z[1:i] <- avg
#     }
#   }
#   x$suavizados <- z
#   return(x)
# }

dai <- data
datos <- data %>% dplyr::filter(especie=="Boana pulchella",
              exposicion_horas==96,compuesto=="Lambdacialotrina") %>%
              dplyr::group_by(dosis) %>%
              dplyr::summarise(n=sum(n), respuesta=sum(respuesta))


datos_suav<- suavizar(datos)

data_abbot <- corr_abbott(datos_suav)

probit(datos)


library(dplyr)

datos <- data_abbot %>% filter(dosis!=0)

model <- glm(Pabbott ~ log10(dosis),
             family = binomial(link = "probit"),
             weights = n,
             data = datos)


summary(model)
# Calcular LD50 y su intervalo de confianza
intercept <- stats::coef(model)[1]
slope <- stats::coef(model)[2]
log_ld50 <- (0 - intercept) / slope
ld50 <- 10^log_ld50



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



#plot(qnorm(resultados_df$Probabilidad), log10(resultados_df$DL), type = "l")
#points(qnorm(data$Presp), log10(data$dosis), pch = 19, col = "red")
