# Helper function of BF

ifunc =  function(x, t_BF, tau, n1, n2){
  return (pt(t_BF * sqrt(1 + n2 / n1 * tau) + sqrt(n2 / n1 * tau)*x, n2-1) * dt(x, n1-1))
}

BF_pvalue_BF = function(t_BF, tau, n1, n2) {
  
  p.L = integrate(function(x) ifunc(x, t_BF, tau, n1, n2), -Inf, Inf)$value
  
  return (min(2*p.L, 2*(1-p.L)))
}

BF_pvalue_BF_target_function = function(t_BF, tau, n1, n2, target_p_value) {
  return(BF_pvalue_BF(t_BF, tau, n1, n2) - target_p_value)
}

BF_pvalue_BF_solver_positive = function(tau, n1, n2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_BF_target_function(t_BF, tau, n1, n2, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of VREPB

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

BF_pvalue_given_lambda_VREPB = function(t_BF, n1, n2, tau, lambda) {
  
  c = tau / (tau + n1/n2)
  gamma = lambda / (lambda + n1/n2)
  
  upper = t_BF * sqrt(n1 + n2 -2)
  lower = sqrt((n1-1) * c / gamma + (n2-1) * (1-c) / (1-gamma))
  
  test_stat = upper/lower
  
  p = pt(q = abs(test_stat), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  
  out = p
}

mass_given_tau_j = function(n1, n2, grid, mass, tau_j) {
  f_tau_j_given_grid = p_tau_j_given_lambda(n1, n2, tau_j, grid)
  f_tau_j_grid = f_tau_j_given_grid * mass
  post_mass = f_tau_j_grid / sum(f_tau_j_grid)
  return(post_mass)
}

BF_pvalue_VREPB = function(t_BF, n1, n2, grid, mass, tau) {
  P_value_joint = BF_pvalue_given_lambda_VREPB(t_BF, n1, n2, tau, grid)
  post_mass = mass_given_tau_j(n1, n2, grid, mass, tau)
  P_value = sum(P_value_joint * post_mass)
  return (P_value)
}

BF_pvalue_VREPB_target_function = function(t_BF, n1, n2, grid, mass, tau, target_p_value) {
  return(BF_pvalue_VREPB(t_BF, n1, n2, grid, mass, tau) - target_p_value)
}

BF_pvalue_VREPB_solver_positive = function(n1, n2, grid, mass, tau, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_VREPB_target_function(t_BF, n1, n2, grid, mass, tau, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of Welch

BF_pvalue_Welch = function(t_BF, tau, n1, n2) {
  
  df_estimate = (tau/n1 + 1/n2)^2 / (1/(n1-1) * (tau / n1)^2 + 1/(n2-1) * (1/n2)^2)
  p = pt(abs(t_BF), df = df_estimate, lower.tail = FALSE) * 2
  
  return (p)
}

BF_pvalue_Welch_target_function = function(t_BF, tau, n1, n2, target_p_value) {
  return(BF_pvalue_Welch(t_BF, tau, n1, n2) - target_p_value)
}

BF_pvalue_Welch_solver_positive = function(tau, n1, n2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_Welch_target_function(t_BF, tau, n1, n2, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of Pooled t-test 

BF_pvalue_EV = function(t_BF, tau, n1, n2) {
  
  test_statistics = t_BF * sqrt((tau / n1 + 1 / n2)/((1/n1 + 1/n2) * (((n1-1) * tau + (n2-1))/(n1+n2-2))))
  p = pt(abs(test_statistics), df = n1+n2-2, lower.tail = FALSE) * 2
  
  return (p)
}

BF_pvalue_EV_target_function = function(t_BF, tau, n1, n2, target_p_value) {
  return(BF_pvalue_EV(t_BF, tau, n1, n2) - target_p_value)
}

BF_pvalue_EV_solver_positive = function(tau, n1, n2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_EV_target_function(t_BF, tau, n1, n2, target_p_value), interval = c(0, 100))$root
  
}










