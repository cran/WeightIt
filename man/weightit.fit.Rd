\name{weightit.fit}
\alias{weightit.fit}
\title{Generate Balancing Weights with Minimal Input Processing}
\description{
\code{weightit.fit()} dispatches one of the weight estimation methods determined by \code{method}. It is an internal function called by \fun{weightit} and should probably not be used except in special cases. Unlike \code{weightit()}, \code{weightit.fit()} does not accept a formula and data frame interface and instead requires the covariates and treatment to be supplied as a numeric matrix and atomic vector, respectively. In this way, \code{weightit.fit()} is to \code{weightit()} what \fun{lm.fit} is to \fun{lm} - a thinner, slightly faster interface that performs minimal argument checking.
}
\usage{
weightit.fit(covs,
             treat,
             method = "ps",
             s.weights = NULL,
             by.factor = NULL,
             estimand = "ATE",
             focal = NULL,
             stabilize = FALSE,
             ps = NULL,
             moments = NULL,
             int = FALSE,
             subclass = NULL,
             is.MSM.method = FALSE,
             missing = NULL,
             verbose = FALSE,
             include.obj = FALSE,
             ...)
}
\arguments{
  \item{covs}{
a numeric matrix of covariates.
}
  \item{treat}{
a vector of treatment statuses.
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See \fun{weightit} for allowable options. The default is \code{"ps"} for propensity score weighting.
}
  \item{s.weights}{
a numeric vector of sampling weights. See the individual pages for each method for information on whether sampling weights can be supplied.
}
  \item{by.factor}{
a factor variable for which weighting is to be done within levels. Corresponds to the \code{by} argument in \fun{weightit}.
}
  \item{estimand}{
the desired estimand. For binary and multi-category treatments, can be "ATE", "ATT", "ATC", and, for some methods, "ATO", "ATM", or "ATOS". The default for both is "ATE". This argument is ignored for continuous treatments. See the individual pages for each method for more information on which estimands are allowed with each method and what literature to read to interpret these estimands.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. For the methods that involve estimating propensity scores, this involves multiplying each unit's weight by the proportion of units in their treatment group. Default is \code{FALSE}.
}
  \item{focal}{
when multi-category treatments are used and ATT weights are requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. Must be non-\code{NULL} if \code{estimand = "ATT"} or \code{"ATC"}.
}
  \item{ps}{
a vector of propensity scores. If specified, \code{method} will be ignored and set to \code{"ps"}.
}
  \item{moments, int, subclass}{
arguments to customize the weight estimation. See \fun{weightit} for details.
}
  \item{is.MSM.method}{
see \fun{weightitMSM}. Typically can be ignored.
}
  \item{missing}{
\code{character}; how missing data should be handled. The options depend on the \code{method} used. If \code{NULL}, \code{covs} covs will be checked for \code{NA} values, and if present, \code{missing} will be set to \code{"ind"}. If \code{""}, \code{covs} covs will not be checked for \code{NA} values; this can be faster when it is known there are none.
}
  \item{verbose}{
whether to print additional information output by the fitting function.
}
  \item{include.obj}{
whether to include in the output any fit objects created in the process of estimating the weights. For example, with \code{method = "ps"}, the \code{glm} objects containing the propensity score model will be included. See the individual pages for each method for information on what object will be included if \code{TRUE}.
}
  \item{...}{
other arguments for functions called by \code{weightit.fit()} that control aspects of fitting that are not covered by the above arguments.
}

}
\value{
  A \code{weightit.fit} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when \code{method} is \code{"ps"}, \code{"gbm"}, \code{"cbps"}, \code{"super"}, or \code{"bart"}.}
  \item{fit.obj}{When \code{include.obj = TRUE}, the fit object.}
  \item{info}{Additional information about the fitting. See the individual methods pages for what is included.}

  The \code{weightit.fit} object does not have specialized \code{print()}, \code{summary()}, or \code{plot()} methods. It is simply a list containing the above components.
}
\details{
\code{weightit.fit()} is called by \code{weightit()} after the arguments to \code{weightit()} have been checked and processed. \code{weightit.fit()} dispatches the function used to actually estimate the weights, passing on the supplied arguments directly. \code{weightit.fit()} is not meant to be used by anyone other than experienced users who have a specific use case in mind. The returned object doesn't contain any information about the supplied arguments or details of the estimation method; all that is processed by \code{weightit}.

Less argument checking or processing occurs in \code{weightit.fit()} than does in \code{weightit()}, which means supplying incorrect arguments can result in errors, crashes, and invalid weights, and error and warning messages may not be helpful in diagnosing the problem. \code{weightit.fit()} does check to make sure weights were actually estimated, though.

\code{weightit.fit()} may be most useful in speeding up simulation simulation studies that use \code{weightit()} because the covariates can be supplied as a numeric matrix, which is often how they are generated in simulations, without having to go through the potentially slow process of extracting the covariates and treatment from a formula and data frame. If the user is certain the arguments are valid (e.g., by ensuring the estimated weights are consistent with those estimated from \code{weightit()} with the same arguments), less time needs to be spent on processing the arguments. Also, the returned object is much smaller than a \code{weightit} object because the covariates are not returned alongside the weights.
}
\seealso{
\fun{weightit}, which you should use for estimating weights unless you know better.
}

\author{
Noah Greifer
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
covs_mat <- as.matrix(lalonde[c("age", "educ", "married",
                                "nodegree", "re74", "re75")])
wf1 <- weightit.fit(covs_mat, treat = lalonde$treat,
                    method = "ps", estimand = "ATT")
str(wf1)
}
