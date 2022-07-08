createLabel.gen.digit <- function (zeAxis, lambda, tau, digit4lambda = 3, digit4tau = 0, axisName = "Dimension ") 
{
  lambda <- round(lambda, digit4lambda)
  tau <- round(tau, digit4tau)
  genLabel <- bquote(.(axisName) * .(zeAxis) * .(". ") ~ 
                       ~lambda == .(lambda) * . ~ ~tau == .(tau) * .("%"))
  return(genLabel)
}