#' ODE solution
#'
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
model <- function(t, y, param){
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param["N"]
  beta <- param["beta"]
  mu <- param["mu"]
  gamma <- param["gamma"]
  lamda <- param["lamda"]
  epsilon <- param["epsilon"]
  dSt <- mu * (N - S) - beta * S * (I1 + I2)/N + epsilon * R
  dEt <- beta * S * (I1 + I2)/N - mu * E - 2 * lamda * E
  dI1t <-  -(mu + gamma) * I1 + lamda * E
  dI2t <- -(mu + gamma) * I2 + lamda * E
  dRt <- gamma * (I1 + I2) - R * (mu + epsilon)
  outcome <- c(dSt, dEt, dI1t, dI1t, dRt)
  list(outcome)
}
  times <- seq(0, 52, by =1/7)
  param <- c(mu = 0.00064, lamda = 0.1, beta = 5, gamma = 0.095, N = 1, epsilon = 0.1)
  init <- c(S = 0.8, E = 0.000032, I1 = 0.00003, I2 = 0.00002, R = 0.000082)
  result <- deSolve::ode(y = init, times = times, func = model, parms = param)
  result <- as.data.frame(result)
  tail(round(result, 3.6), 10)
  #' @export
  seirplot <- ggplot2::ggplot(data = result) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = S, col = "S"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = I1, col = "I1"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = I2, col = "I2"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = R, col = "R"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = E, col = "E"), lwd = 2) +
    ggplot2::labs(x = "Time", y = "Ratio") +
    ggplot2::scale_color_manual(name = "SEIR", values = c("S" = "orange", "E" = "purple", "I1" = "red", "I2" = "blue", "R" = "green"))
  seirplot
  ggplot2::ggsave(seirplot, file = "seir.pdf", width = 7, height = 6)
  ggplot2::ggsave(seirplot, file = "seir.svg", width = 7, height = 6)

