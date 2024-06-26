% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_haz.R
\name{predict.haz}
\alias{predict.haz}
\title{wrapper for piece.haz und GG.haz for predction of hazard rates (Lam_o*exp(xb))}
\usage{
\method{predict}{haz}(object, newdata, time = 1, haz.strata)
}
\arguments{
\item{object}{ph.surv estimation object}

\item{newdata}{newdata: dataframe, first column is time, remaining colmns X}

\item{time}{gives position of time column or its name}

\item{haz.strata}{gives position of stratum variable of baseline hazard in the same way as time if baseline hazard is not stratified NA has to be passed!}
}
\value{
matrix of cumlative baseline hazards
}
\description{
wrapper for piece.haz und GG.haz for predction of hazard rates (Lam_o*exp(xb))
}
