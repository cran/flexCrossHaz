bsplineSurv <-function(x, nterm, deg=3, deriv=0){
      degree<-deg
      rx <- range(x)
      nterm <- round(nterm)
      dx <- (rx[2] - rx[1])/nterm
      knots <- c(rx[1] + dx * ((-degree):(nterm - 1)), rx[2] + dx * (0:degree))
      newx <- spline.des(knots, x, degree + 1,derivs=rep(deriv,length(x)))$design
      newx
}

