make_trapezoid <- function(w_top, w_bot, H_bank, Z_fp, dx = 0.1){
  # CREATE A CHANNEL CROSS SECTION with a trapexzoidal geometry
  dx <- 0.1
  dw <- (w_top - w_bot) / 2
  dist <- c(0, dw, dw + w_bot, w_top)
  stad <- c(0, H_bank, H_bank, 0)
  z <- Z_fp - stad
  
  xs_data <- approx(x = dist,
                    y = z,
                    xout = seq(min(dist), max(dist), dx))
  xs_data <- as.data.frame(xs_data)
  colnames(xs_data) <- c("dist", "Z")
  return(xs_data)
}

estimate_hydraulics <- function(xs_data,
                                S,
                                d,
                                d84 = 2 * d,
                                ymax = 3.5,
                                dx = 0.05,
                                t_crit = t_crit_soul(d),
                                g = 9.81,
                                rho = 1000,
                                rho_s = 2650,
                                a1 = 6.5,
                                a2 = 2.5,
                                kv = 1e-6,
                                fun = "ec",
                                rp = FALSE){
  Gs <- (rho_s - rho) / rho
  #specify the range of target WSE, and calculate dx and dp
  y_range <- seq(0, ymax, dx)  #make calculations for a range of depths for each section
  dx <- (c(0,diff(xs_data$dist)) + c(diff(xs_data$dist),0) ) /2
  px <- sqrt(diff(xs_data$dist)^2 + diff(xs_data$Z)^2)
  dp <- (c(0, px) + c(px, 0)) / 2
  invert <- min(xs_data$Z, na.rm = T)
  thalweg <- which(xs_data$Z == min(xs_data$Z, na.rm = T))[1]
  xs_hyd <- data.frame(flow = array(NA, length(y_range)),
                       area = array(NA, length(y_range)),
                       width = array(NA, length(y_range)),
                       perimeter = array(NA, length(y_range)),
                       velocity = array(NA, length(y_range)),
                       transport = array(NA, length(y_range)),
                       wse = array(NA, length(y_range)),
                       t_star = array(NA, length(y_range)))
  
  for(j in seq_along(y_range)){
    wse <- invert + y_range[j]
    depths <- calc_depths(xs_data, wse)
    #plot(depths)
    A <- sum(depths * dx, na.rm = T)
    W <- sum(dx[!is.na(depths)], na.rm = T)
    P <- sum(dp[!is.na(depths)], na.rm = T)
    R <- A / P
    U <- u_ferg(R, S, d84, a1, a2, g)
    tau <- tau_1D(R, S, g, rho)
    t_star <- t_star(tau, d, g, rho, rho_s)
    if(fun == "eb"){
      q_star <- eb(R, S, d, g, rho, rho_s, kv)
    }
    if(fun == "ec"){
      q_star <- ec(R, S, d, d84, g, rho, rho_s, U, t_crit)
    }
    if(fun == "rk"){
      d84_transp <- d * 2  #recalculate based on d
      q_star <- rk(R, S, d, d84_transp, g, rho, rho_s, rp)
    }
    if(fun == "mpm"){
      q_star <- mpm(R, S, d, t_crit, g, rho, rho_s)
    }
    if(fun == "wp"){
      q_star <- wp(R, S, d, t_crit, g, rho, rho_s)
    }
    if(fun == "vr"){
      q_star <- vr(R, S, d, t_crit, g, rho, rho_s, kv)
    }
    if(fun == "yln"){
      q_star <- yln(R, S, d, t_crit, g, rho, rho_s, kv)
    }
    if(fun == "rk"){
      Qb <- qb_conversion(q_star, d84_transp, g, Gs) * W
    } else {
      Qb <- qb_conversion(q_star, d, g, Gs) * W
    }
    Q <- A * U
    xs_hyd[j,] <- c(Q, A, W, P, U, Qb, wse, t_star)
  }
  return(xs_hyd)
}



qb_conversion <- function(q_star, d, g = 9.81, Gs = 1.65){
  qb <- q_star * sqrt(Gs * g * d^3)
  return(qb)
}

tau_1D <- function(R, S, g = 9.81, rho = 1000){
  tau <- g * rho * R * S
  return(tau)
}

u_ferg <- function(R, S, d84, a1 = 6.5, a2 = 2.5, g = 9.81){
  Res = a1 * a2 * (R / d84) / (a1^2 + a2^2 * (R / d84)^(5/3))^(1/2)
  U = Res * sqrt(g * R * S)
  return(U)
}

t_star <- function(tau, d, g = 9.81, rho = 1000, rho_s = 2650){
  Gs <- (rho_s - rho) / rho
  t_star <- tau / (g * rho * Gs  * d)
  return(t_star)
}

t_crit_soul <- function(d, rho = 1000, rho_s = 2650, g = 9.81, kv = 1e-6){
  Gs <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)
  tc_star <- 0.3 / (1 + 1.2 * d_s) + 0.055*(1 - exp(-0.020 * d_s))
  return(tc_star)
}

t_crit_vr <- function(d, rho = 1000, rho_s = 2650, g = 9.81, kv = 1e-6){
  Gs <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)
  if(d_s <= 4){
    t_crit <- 0.24 * d_s^-1
  }
  if(d_s > 4 & d_s <= 10){
    t_crit <- 0.14 * d_s^-0.64
  }
  if(d_s > 10 & d_s <= 20){
    t_crit <- 0.04 * d_s^-0.1
  }
  if(d_s > 20 & d_s <= 150){
    t_crit <- 0.013 * d_s^0.29
  }
  if(d_s > 150){
    t_crit <- 0.056
  }
  return(t_crit)
}

d_star <- function(d, Gs = 1.65, g = 9.81, kv = 1e-6){
  d_star <- d * (Gs * g / kv^2)^(1/3)
  return(d_star)
}

mpm <- function(R, S, d, t_crit = 0.047, g = 9.81, rho = 1000, rho_s = 2650){
  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  
  if(t_s > t_crit) {
    q_star <- 8 * (t_s - t_crit)^(3/2)
  } else {
    q_star <- 0
  }
  
  return(q_star)
}

wp <- function(R, S, d, t_crit = 0.047, g = 9.81, rho = 1000, rho_s = 2650){
  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  
  if(t_s > t_crit) {
    q_star <- 4 * (t_s - t_crit)^(3/2)
  } else {
    q_star <- 0
  }
  
  return(q_star)
}

vr <- function(R, S, d, t_crit = t_crit_vr(d, rho, rho_s, g, kv), g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6){
  tau <- tau_1D(R,S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  Gs <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)
  q_star <- (0.053 / d_s^(0.3))*(t_s / t_crit - 1)^2.1
  return(q_star)
}

yln <- function(R, S, d, t_crit = t_crit_soul(d), g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6){
  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  Gs <- (rho_s - rho) / rho
  r <- t_s / t_crit - 1
  sigma <- 2.45 * sqrt(t_crit) /  (Gs + 1)^0.4
  
  if(t_s > t_crit) {
    q_star <- 0.635 * r * sqrt(t_s) * (1 - (1 /  (sigma * r)) * log(1 + sigma * r))
  } else {
    q_star <- 0
  }
  
  return(q_star)
}

eb <- function(R, S, d, g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6) {
  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  Gs <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)
  K <- sqrt(2/3 + 36 / d_s^3) - sqrt(36 / d_s^3)
  if(t_s >= 0.19) {
    q_star <- 40 * K * t_s^3
  } else {
    q_star <- (K * exp(-0.391 / t_s)) / 0.465
  }
  return(q_star)
}

ec <- function(R,
               S,
               d,
               d84 = 2 * d,
               g = 9.81,
               rho = 1000,
               rho_s = 2650,
               U = u_ferg(R, S, d84, a1, a2, g),
               a1 = 6.5,
               a2 = 2.5,
               kv = 1e-6,
               t_crit = t_crit_soul(d, rho, rho_s, g, kv)
){
  tau <- tau_1D(R,S, g, rho)
  Gs <- (rho_s - rho) / rho
  omega_s <- tau * U / (rho * (g * Gs * d)^(3/2))  #dimensionless stream power
  dcrit <- t_crit * Gs * d / S #depth at which transport is initiated
  Res <- u_ferg(dcrit, S, d84, a1, a2, g) / sqrt(g * dcrit * S) #reference resistance
  om_crit <- Res * t_crit^(3/2)  #reference dimensionless stream power
  if (omega_s > 0) {
    E_star <- (0.92 - 0.25 * sqrt(om_crit / omega_s))^9  #function fit by Eaton and Church 2011
  } else {
    E_star <- 0
  }
  q_star = E_star * omega_s  #translate efficiency to dimensionless bedload
  if(q_star < 0.00001){
    q_star <- 0  #return 0 for in the event that q_star is less than 0.00001
  }
  return(q_star)
}


load_section <- function(xs_nm, dx = NA, rev = FALSE, trim = FALSE, mkplt = FALSE){
  #load data and calculate distance along the profile
  xs_data <- utils::read.csv(file = xs_nm, skip = 1)
  colnames(xs_data) <- c("E", "N", "Z")
  xs_data$dist <- c(0, cumsum(sqrt(diff(xs_data$N)^2 + diff(xs_data$E)^2)))
  if(rev){
    xs_data$dist <- max(xs_data$dist, na.rm = T) - xs_data$dist
  }
  if(mkplt){
    graphics::plot(xs_data$dist, xs_data$Z,
                   type = "l",
                   col = "darkgreen",
                   asp = 2,
                   lty = 1,
                   lwd = 2,
                   xlab = "distance from left bank (m)",
                   ylab = "elevation above sea level (m)")
    graphics::legend("bottomright",
                     inset = c(0.01, 0.02),
                     legend = c("Original", "Trimmed", "Resampled"),
                     col = c("darkgreen", "darkorange", "darkred"),
                     lty = c(1, 1, NA),
                     pch = c(NA, NA, 19),
                     cex = 0.7)
  }
  if(trim){
    #find the thalweg of the channel, where long profile intersects the cross section
    center <- which(xs_data$Z == min(xs_data$Z, na.rm = T))[1]
    
    #find left bank & delete ground points left of bank
    filt <- xs_data$dist < xs_data$dist[center]
    leftbank <- which(xs_data$Z == max(xs_data$Z[filt], na.rm = T))[1]
    filt <- xs_data$dist < xs_data$dist[leftbank]
    xs_data$Z[filt] <- NA  #delete topography outside the main channel
    
    #find the right bank
    filt <- xs_data$dist > xs_data$dist[center]
    rightbank <- which(xs_data$Z == max(xs_data$Z[filt], na.rm = T))[1]
    filt <- xs_data$dist > xs_data$dist[rightbank]
    xs_data$Z[filt] <- NA #delete topography outside the main channel
    if(mkplt){
      graphics::lines(xs_data$dist, xs_data$Z, lty = 1, lwd = 2, col = "darkorange")
    }
  }
  #drop the data outside the banks
  xs_data <-xs_data[!is.na(xs_data$Z),]
  
  if(is.na(dx)){
    return(data.frame(E = xs_data$E,
                      N = xs_data$N,
                      Z = xs_data$Z,
                      dist = xs_data$dist))
  } else{
    sample_d <- seq(from = min(xs_data$dist),
                    to = max(xs_data$dist),
                    by = dx)
    sample_z <- stats::approx(x = xs_data$dist,
                              y = xs_data$Z,
                              xout = sample_d)[[2]]
    sample_x <- stats::approx(x = xs_data$dist,
                              y = xs_data$E,
                              xout = sample_d)[[2]]
    sample_y <- stats::approx(x = xs_data$dist,
                              y = xs_data$N,
                              xout = sample_d)[[2]]
    if(mkplt){
      graphics::points(sample_d, sample_z, pch = 19, cex = 0.6, col = "darkred")
    }
    return(data.frame(E = sample_x,
                      N = sample_y,
                      Z = sample_z,
                      dist = sample_d))
  }
}

calc_depths <- function(xs_data, wse, mkplt = F){
  thalweg <- which(xs_data$Z == min(xs_data$Z))[1]
  depth <- wse - xs_data$Z
  depth[depth < 0] <- NA
  for (j in seq_along(xs_data$dist)){
    if(j < thalweg){
      if(sum(is.na(depth[j:thalweg])) > 0){
        depth[j] <- NA  #set all depths separated from main channel by a dry patch to zero
      }
    }
    if(j > thalweg){
      if(sum(is.na(depth[thalweg:j])) > 0){
        depth[j] <- NA  #set all depths separated from main channel by a dry patch to zero
      }
    }
  }
  if(mkplt){
    graphics::plot(xs_data$dist,
                   xs_data$Z,
                   type = "l",
                   col = "darkgreen",
                   asp = 2,
                   main = "locations where water depths have been calculated",
                   xlab = "distance from left bank (m)",
                   ylab = "elevation above sea level (m)")
    graphics::abline(h = wse, lty = 2, col = "blue")
    # Filter out NA values for plotting points
    valid_points <- !is.na(depth)
    graphics::points(xs_data$dist[valid_points], xs_data$Z[valid_points],
                     pch = 19, col = "blue")
  }
  return(depth)
}


