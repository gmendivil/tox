% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agregar.R
\name{agregar}
\alias{agregar}
\title{Genera un dataframe con los datos agregados.}
\usage{
agregar(df)
}
\arguments{
\item{df}{Dataframe con los resultados de todas las rèplicas. El mismo debe contener las columnas dosis , respuesta y n.}
}
\value{
Dataframe con los datos agregados.
}
\description{
Suma los valores de respuesta de las réplicas y los n's.
}
