source('../EPB main/EPB.R')

# Power, FDR

Power = function(discovery, flag_list) {
  dis_count = sum(flag_list)
  true_dis_count = sum(flag_list[discovery])
  power = true_dis_count/dis_count
  return (power)
}

FDP = function(discovery, flag_list) {
  false_dis = length(discovery) - sum(flag_list[discovery])
  fdp = false_dis/max(length(discovery) ,1)
  return (fdp)
}

# Data Generator

data_generator = function(n1, n2, data_generation_parameter, equ_var) {
  
  n1=n1
  n2=n2
  k=data_generation_parameter$k
  d1=data_generation_parameter$d1
  d2=data_generation_parameter$d2
  m=data_generation_parameter$m
  mu1=data_generation_parameter$mu1
  mu2=data_generation_parameter$mu2
  mean_var2 = data_generation_parameter$mean_var2
  var_var2 = data_generation_parameter$var_var2
  pi0 = data_generation_parameter$pi0
  mu0 = data_generation_parameter$mu0
  
  null_count = as.integer(pi0 * m)
  dis_count = m - null_count
  flag_list = c(rep(0, null_count), rep(1, dis_count)) # 0 for null, 1 for truth
  
  lambda = k * rf(m, d1, d2)
  var2 = abs(rnorm(m, mean_var2, sqrt(var_var2)))
  var1 = lambda * var2
  if (equ_var) {
    var1 = var2
  }
  
  X1 = matrix(0, m, n1)
  X2 = matrix(0, m, n2)
  
  for (j in 1:m) {
    if (j <= null_count) {
      X1[j, ] = rnorm(n1, mu0, sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mu0, sqrt(var2[j]))
    }
    else {
      X1[j, ] = rnorm(n1, rnorm(1, mean = 0, sd = sqrt(mu1 * var1[j])), sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mu2, sqrt(var2[j]))
    }
  }
  
  output = list('X1' = X1, 'X2' = X2, 'flag_list' = flag_list)
  
  return (output)
}

dir_name = function(equ_var, n1, n2, data_generation_parameter, VR_parameter, DV_parameter, alpha) {
  
  equ = 'unequal'
  if (equ_var) {
    equ = 'equal'
  }
  
  n1=n1
  n2=n2
  k=data_generation_parameter$k
  d1=data_generation_parameter$d1
  d2=data_generation_parameter$d2
  m=data_generation_parameter$m
  mu1=data_generation_parameter$mu1
  mu2=data_generation_parameter$mu2
  mean_var2 = data_generation_parameter$mean_var2
  var_var2 = data_generation_parameter$var_var2
  pi0 = data_generation_parameter$pi0
  mu0 = data_generation_parameter$mu0
  
  B = VR_parameter[1]
  l1 = VR_parameter[2]
  u1 = VR_parameter[3]
  
  B1 = DV_parameter[1]
  B2 = DV_parameter[2]
  l2 = DV_parameter[3]
  u2 = DV_parameter[4]
  
  base_dir = paste('Simulation_result/', equ, sep = '')
  
  if (!dir.exists(base_dir)) {
    print('Create Base Directory')
    dir.create(base_dir)
  }
  
  head = paste(base_dir, '/(', paste(n1, n2, k, d1, d2,  m,  mu1,  mu2,  mean_var2,  var_var2,  pi0, mu0, B, l1, u1, B1, B2, l2, u2, alpha, sep = ','), ')', sep = '')
  
  return(head)
}

file_name = function(rounds, algorithm_list) {
  time = Sys.time()
  file = paste(time, ': R=', rounds, '', sep = '')
  for (code in algorithm_list) {
    file = paste(file, code, sep = ',')
  }
  return (paste('(', file, ')', sep = ''))
}

simulator = function(seed, data_generation_parameter, VR_parameter, DV_parameter, alpha, rounds, algorithm_list, equ_var) {
  
  if (!dir.exists('Simulation_result')) {
    print('Create Data Directory')
    dir.create('Simulation_result')
  }
  
  set.seed(seed)
  
  args = commandArgs(TRUE)
  
  n1 = as.integer(args[1])
  n2 = as.integer(args[2])
  
  dir = dir_name(equ_var, n1, n2, data_generation_parameter, VR_parameter, DV_parameter, alpha)
  file = file_name(rounds, algorithm_list)
  
  if (!dir.exists(dir)) {
    print('Create Directory')
    dir.create(dir)
  }
  
  algorithm_name = c('VREPB', 'DVEPB', 'Welch', 'B_F', 'EV_test')
  
  FDP_of_algorithms = matrix(0, 5, rounds)
  Power_of_algorithms = matrix(0, 5, rounds)
  
  print('Simulation Start')
  print(paste('Parameters:', file.path(dir, file)))
  
  r = 1
  
  while(r <= rounds) {
    
    print(paste('Start of round', r))
    
    Rerun = FALSE
    
    output_r = data_generator(n1, n2, data_generation_parameter, equ_var)
    X1 = output_r$X1
    X2 = output_r$X2
    flag_list = output_r$flag_list
    information = information_extractor(X1, X2)
    m = information$m
    
    for (code in algorithm_list) {
      
      if (code == 1) {
        print(paste('start of', algorithm_name[code]))
        
        P_list = P_value_VREPB(information, VR_parameter)
        
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 2) {
        print(paste('start of', algorithm_name[code]))
        
        DV_NPMLE_result = DV_NPMLE(information, B1 = DV_parameter[1], B2 = DV_parameter[2], lower_quantile = DV_parameter[3], upper_quantile = DV_parameter[4])
        P_list = P_value_DVEPB(information, DV_NPMLE_result$grid, DV_NPMLE_result$mass)
        
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        
        if (power == 0) {
          print(paste('Error detected and rerun round', r))
          Rerun = TRUE
          break
        }
        
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 3) {
        print(paste('start of', algorithm_name[code]))
        
        P_list = P_value_Welch(information)
        
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 4) {
        print(paste('start of', algorithm_name[code]))
        
        P_list = P_value_BF(m, X1, X2)
        
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 5) {
        print(paste('start of', algorithm_name[code]))
        
        P_list = P_value_EV_test(m, X1, X2)
        
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
    }
    
    if (Rerun) {
      next
    }
    
    print(paste('End of round', r))
    print('')
    
    print('Updating Round Data')
    
    power_df_r = data.frame('Round_Power' = c(1:r), 'VREPB' = Power_of_algorithms[1, 1:r], 'DVEPB' = Power_of_algorithms[2, 1:r], 'Welch' = Power_of_algorithms[3, 1:r], 'B_F' = Power_of_algorithms[4, 1:r], 'EV_test' = Power_of_algorithms[5, 1:r])
    write.csv(power_df_r, file = paste(file.path(dir, file), '_power.csv', sep = ''))
    fdp_df_r = data.frame('Round_FDP' = c(1:r), 'VREPB' = FDP_of_algorithms[1, 1:r], 'DVEPB' = FDP_of_algorithms[2, 1:r], 'Welch' = FDP_of_algorithms[3, 1:r], 'B_F' = FDP_of_algorithms[4, 1:r], 'EV_test' = FDP_of_algorithms[5, 1:r])
    write.csv(fdp_df_r, file = paste(file.path(dir, file), '_fdp.csv', sep = ''))
    
    r = r + 1
    
  }
  
  print('Simulation Over')
  print('Saving Simulation Result')
  
  Power_list = rowMeans(Power_of_algorithms)
  FDR_list = rowMeans(FDP_of_algorithms)
  simulation_result = data.frame('Algorithm' = algorithm_name[algorithm_list], 'Power' = Power_list[algorithm_list], 'FDR' = FDR_list[algorithm_list])
  write.csv(simulation_result, file = paste(file.path(dir, file), '_summary.csv', sep = ''))
  
  print(simulation_result)
  
  return (0)
}

#seed = Sys.time()
seed = 1
alpha = 0.1
rounds = 50
VR_parameter = c(1000, 0, 1.0)
DV_parameter = c(80, 80, 0.01, 1.0)
algorithm_list = c(1,2,3,4,5)
equ_var = FALSE
data_generation_parameter = data.frame('k' = 2, 'd1' = 8, 'd2' = 12, 'm' = 5000, 'mu1' = 12, 'mu2' = 0, 'mean_var2' = 6, 'var_var2' = 4, 'pi0' = 0.9, 'mu0' = 0)

simulator(seed, data_generation_parameter, VR_parameter, DV_parameter, alpha, rounds, algorithm_list, equ_var)
