#' @export
plot_probit_modelos <- function(results1, results2, titulo ="Comparación de Modelos",
                                ld50_line = TRUE,
                                ci = TRUE,
                                model_names = c("Modelo 1", "Modelo 2"),
                                colors = c("#4E79A7", "#59A14F")) {

  # Función interna para procesar cada modelo
  prepare_model_data <- function(results, color, extended_dosis) {
    model <- results$model
    data <- results$CorrectedData
    ld50 <- results$LD50
    ci_ld50 <- results$LD50_CI

    # Dataframe para predicción con rango extendido
    new_data <- data.frame(dosis = extended_dosis)

    # Predecir probabilidades
    pred <- stats::predict(model,
                    newdata = new_data,
                    type = "response",
                    se.fit = TRUE)

    data.frame(
      dosis = new_data$dosis,
      prob = pred$fit,
      lower = pred$fit - 1.96 * pred$se.fit,
      upper = pred$fit + 1.96 * pred$se.fit,
      model = model_names[which(colors == color)],
      color = color,
      ld50 = ld50,
      ci_lower = ci_ld50[1],
      ci_upper = ci_ld50[2]
    )
  }

  # Calcular rango extendido para ambos modelos
  all_doses <- c(results1$CorrectedData$dosis, results2$CorrectedData$dosis)
  log_range <- log10(all_doses)
  extended_log <- c(min(log_range) - 1, max(log_range) + 1) # Extender 1 orden de magnitud
  extended_dosis <- 10^seq(extended_log[1], extended_log[2], length.out = 200)

  # Generar datos para ambos modelos
  plot_data1 <- prepare_model_data(results1, colors[1], extended_dosis)
  plot_data2 <- prepare_model_data(results2, colors[2], extended_dosis)
  plot_data <- rbind(plot_data1, plot_data2)

  # Combinar datos experimentales
  exp_data <- rbind(
    cbind(results1$CorrectedData, model = model_names[1]),
    cbind(results2$CorrectedData, model = model_names[2])
  )

  # Crear gráfico base
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dosis, y = prob, color = model)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(data = exp_data, ggplot2::aes(x = dosis, y = Presp), size = 3) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^extended_log[1], 10^extended_log[2])
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = expression(paste("Dose (", mu, "g/L)")),
         y = "Dead/exposed",
         title = titulo) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          legend.position = c(0.1,0.8),
          legend.title = ggplot2::element_blank())

  # Añadir intervalos de confianza
  if(ci) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = model),
                  alpha = 0.2,
                  color = NA) +
      ggplot2::scale_fill_manual(values = colors)
  }
  # Añadir líneas de DL50

  if(ld50_line) {
    p <- p +
      ggplot2::geom_vline(data = unique(plot_data[,c("model","ld50","ci_lower","ci_upper")]),
                ggplot2::aes(xintercept = ld50, color = model),
                 linetype = "dashed",
                 show.legend = FALSE) +
      ggplot2::geom_label(data = unique(plot_data[,c("model","ld50","ci_lower","ci_upper")]),
                 ggplot2::aes(x = ld50,
                     y = 0.5,
                     label = sprintf("LD50 = %.3f\n(%.3f - %.3f)",
                                     ld50, ci_lower, ci_upper),
                     color = model),
                 hjust = -0.1,
                 vjust = 0.5,
                 show.legend = FALSE) +
      ggplot2::annotation_logticks(sides = "b")
  }

  return(p)
}
