#' Create nice labels for Factor Maps and other maps (must be ggplot2 graphs), for both the x and y axes at the same time.
#'
#' @param x_axis Dimension number for x axis (default to 1) ;
#' @param y_axis Dimension number for y axis (default to 2) ;
#' @param lambda Eigenvalues ;
#' @param tau Percentage of explained variances ;
#' @param digit4lambda How many digits should be displayed for eigenvalues ;
#' @param digit4tau How many digits should be displayed for percentages of explained variances ;
#' @param axisName Name for the axes (e.g. "Dimension", "Component").
#'
#' @return x and y axis labels
#' @export
#'
#' @examples
#' createxyLabels.gen.digit(lambda = 2:1, tau = 10*2:1)
createxyLabels.gen.digit <- function(
    x_axis = 1,
    y_axis = 2,
    lambda,
    tau,
    digit4lambda = 3,
    digit4tau = 0,
    axisName = "Dimension ") {
    xyLabels = ggplot2::labs(
      x = createLabel.gen.digit(
        zeAxis = x_axis,
        lambda = lambda[x_axis],
        tau = tau[x_axis],
        digit4lambda = digit4lambda,
        digit4tau = digit4tau,
        axisName
      ),
      y = createLabel.gen.digit(
        zeAxis = y_axis,
        lambda = lambda[y_axis],
        tau = tau[y_axis],
        digit4lambda = digit4lambda,
        digit4tau = digit4tau,
        axisName
      )
    )
    return(xyLabels)
  }

#' @export
#' @keywords internal

createLabel.gen.digit <- function (
  zeAxis,
  lambda,
  tau,
  digit4lambda = 3,
  digit4tau = 0,
  axisName = "Dimension ")
  {
    lambda <- round(lambda, digit4lambda)
    tau <- round(tau, digit4tau)
    genLabel <- bquote(.(axisName) * .(zeAxis) * .(". ") ~
                         ~ lambda == .(lambda) * . ~ ~ tau == .(tau) * .("%"))
    return(genLabel)
  }
