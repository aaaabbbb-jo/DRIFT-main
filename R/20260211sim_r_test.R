# 1. set up ------------
library(mirt)
library(tidyverse)
library(CVXR)
library(ggplot2)

source("/n/home09/lyy/IRT/0908function.R")
source("/n/home09/lyy/IRT/1009function_add.R")
folder_name = "/n/home09/lyy/IRT/20260211_Result/"

###   Code Log: Main changes here: 

#1  opt_delta: opt.delta = select_delta(est.Z, simul.data$O, simul.data$Y) better  Complete

#2  gamma_Y distribution

seed_num = 42
nx = 5
nz = 3
ny = 30 # number of response items
nU = 1000 # number of evaluating vactor
is_exclude = T

##### Lambda List #####
set.seed(1)
simul.Lambda = list(matrix(round(rnorm(nx*nz, sd = 0.5),2), nrow = nx, ncol = nz), 
                    matrix(round(rnorm(nx*nz),2), nrow = nx, ncol = nz)) 

gamma_GEO = c(rep(1, nz), 0)
mu_vmf = gamma_GEO[1:nz]/norm(gamma_GEO[1:nz],"2")
kappa = 2

dis_alpha = 5
dis_beta = 1.5
scale_beta = 1
hist(rbeta(nU,dis_alpha,dis_beta))

set.seed(1)
alpha_angle <- rvmf_rt(ny, mu_vmf, kappa)
alpha_dis <- scale_beta * rbeta(ny,dis_alpha,dis_beta)
simul.alpha <- alpha_angle*alpha_dis + gamma_GEO[1:nz]

##### testing gamma_U for nU times #####
set.seed(1)
gamma_angle <- rvmf_rt(nU, mu_vmf, kappa)
gamma_dis <- scale_beta * rbeta(nU,dis_alpha,dis_beta)
gamma_mat <- gamma_angle*gamma_dis + gamma_GEO[1:nz]

#delta_adj = c(0,0.2,0.5,1,1.5)
quan_delta_select =1
Sigma <- diag(nz)     
method_list <- c("R-learner", "Observed Maximin", "Factorized GEO", "DRIFT","GEO_N", "DRIFT_N")
num_method<- length(method_list)

r_list = c(1,1.5)
# nsub_list = c(100,300,900)
#r_list = c(0.6,1)
nrep = 100
nsub = 300
#r = 0.6 
reset = T # whether re-generate alpha each time
plot_df = data.frame(c())
fig_n = paste0(folder_name,"change_r.jpg")
for (r in r_list){
  file_names=paste0(folder_name,"_nsub=",nsub,"_r=",r,"_kappa=",kappa,".RData")
  fig_names=paste0(folder_name,"_nsub=",nsub,"_r=",r,"_kappa=",kappa,".jpg")
  scale_beta  = r*norm(gamma_GEO[1:nz],"2") 
  simul.res = list()
  cat.y = rep(2, ny) # all binary, but note here it's y=1 and =2
  simul.tau = cbind(rep(0,ny),rep(Inf, ny)) # tau is the threshold for each category, the last category is always Inf 
  simul.res$acc_metric1_n200 = array(dim = c(num_method, nU, nrep))
  simul.res$acc_metric2_n200 = array(dim = c(num_method, nU, nrep))
  simul.res$acc_metric3_n200 = array(dim = c(num_method, nU, nrep))
  
  # LOAD PREVIOUS simul.res before running so that it can save all case result
  ix <- 1
  sd = seed_num
  set.seed(ix)
  alpha_dis <- scale_beta * rbeta(ny,dis_alpha,dis_beta)
  gamma_dis <- scale_beta * rbeta(nU,dis_alpha,dis_beta)
  
  while (ix <=nrep) {
    set.seed(sd)
    print(paste0("r = ", r,"; Seed = ",sd,"; Nrep = ",ix))
    sd = sd + 1
    sk=F
    
    ##### Generating alpha_Y for ny times #####
    alpha_angle <- rvmf_rt(ny, mu_vmf, kappa)
    gamma_angle <- rvmf_rt(nU, mu_vmf, kappa)
    
    simul.alpha <- alpha_angle*alpha_dis + gamma_GEO[1:nz]
    gamma_mat <- gamma_angle*gamma_dis + gamma_GEO[1:nz]
    
    simul.data = generate_FAS_0930(nsub, ny, cat.y, nz, nx, simul.Lambda, simul.tau, simul.alpha, gamma_GEO, centerZ = T)
    
    item.data = simul.data$Y
    colnames(item.data) <- paste0("item_", 1:ncol(item.data))
    design_matrix <- data.frame(X = simul.data$X, A = simul.data$A)
    x_vars <- colnames(design_matrix)[1:(nx-1)]
    main_terms <- paste(colnames(design_matrix[-nx]), collapse = "+")
    interaction_terms <- paste(paste0("A:", x_vars), collapse = "+")
    functext <- paste0("mod1b <- mixedmirt(item.data, design_matrix, model=",nz,", lr.fixed = ~ ",paste(main_terms, interaction_terms, sep = "+"), 
                       ", fixed = ~ 0, itemtype = 'graded', SE = F)")
    mod1b = eval(parse(text = functext))
    
    # factor parameter
    est.item = do.call(rbind,coef(mod1b, simplify=T)[1:ny])# if has CI [seq(1,60,3),]
    est.alpha = est.item[,1:nz] # our sign is different # maybe automatically let the last two triangle to be 0 to fix factor?? but factor loading not identifiable anyway
    est.tau = est.item[,(nz+1)] # intercept not zero, not identifiable
    # latent regression parameter
    lr.betas = coef(mod1b)$lr.betas[1,] # note the est.Lambda reading is specifically for nz=3, nx=4; for other nz, nx, need to change
    est.Lambda = getLambda(nx,nz,lr.betas)
    est.Z = t(sapply(1:nsub, function(ix) {f_w(simul.data$X[ix,], est.Lambda[[(simul.data$A + 1)[ix]]])}))  # mean to estimate
    test.data = tibble(response = simul.data$O, X = simul.data$X[,1:(nx-1)], A = simul.data$A, Z = simul.data$Z, post_Z = est.Z)
    
    tmp = generate_U(simul.data$Z, gamma_mat)
    simul.data$U = tmp$U
    simul.data$gamma_mat = tmp$gamma_mat
    
    ref.TE_true = simul.data$X %*% apply(simul.data$gamma_mat, 1, ref_TE_true, 
                                         Lambda1 = simul.Lambda[[2]], Lambda0 = simul.Lambda[[1]])
  
    fit.base = fit_baseline(test.data)  
    simul.res$acc_metric1_n200[1, , ix] = colMeans(fit.base$TE * ref.TE_true > 0)
    simul.res$acc_metric2_n200[1, , ix] = cor(fit.base$TE,ref.TE_true)
    simul.res$acc_metric3_n200[1, , ix] = colMeans((fit.base$TE - ref.TE_true)^2)
    
    fit.prep = fit_prep(simul.data$X, est.Lambda[[2]], est.Lambda[[1]])
    opt.delta = select_delta_new(est.Z,simul.data$O,simul.data$Y,quan=quan_delta_select, exclude_null = is_exclude,upper_b = 0.95)  
    
    ###### For NO GEO #######
    
    Cov_Z = as.matrix(t(est.Z)%*%est.Z/nrow(est.Z))
    L = chol(Cov_Z)
    get_GEO = geo_center(L, est.alpha) 
    logits = test.data$post_Z%*% get_GEO$gamma
    est_O=rbinom(n = nsub, size = 1, prob = 1 / (1 + exp(-logits)))
    opt.delta_no = select_delta_new_noGEO(est.Z,est_O,simul.data$Y,get_GEO$gamma, quan=quan_delta_select, exclude_null = is_exclude,upper_b = 0.95)  
  
    
    ### Group DRO
    DRO_fit = group_DRO(fit.prep, test.data, est.alpha, solver = c("ECOS", "SCS"))
    simul.res$acc_metric1_n200[2, , ix] = colMeans(DRO_fit$pred.te * ref.TE_true > 0)
    simul.res$acc_metric2_n200[2, , ix] = cor(DRO_fit$pred.te,ref.TE_true)
    simul.res$acc_metric3_n200[2, , ix] = colMeans((DRO_fit$pred.te-ref.TE_true )^2)
    
    id_list = c(0,1)
    for (id in 1:length(id_list)) {
      result <- tryCatch({
        fit.maximin = fit_maximin(id_list[id]*opt.delta, fit.prep, test.data)
        acc.maximin = colMeans(fit.maximin$pred.te * ref.TE_true > 0)
        simul.res$acc_metric1_n200[id + 2, , ix] = c(acc.maximin)
        simul.res$acc_metric2_n200[id + 2, , ix] = cor(fit.maximin$pred.te,ref.TE_true)
        simul.res$acc_metric3_n200[id + 2, , ix] = colMeans((fit.maximin$pred.te - ref.TE_true)^2)
        
        fit.maximin =  fit_maximin_new_noGEO(id_list[id]*opt.delta_no, fit.prep, get_GEO$gamma,est_O, test.data)
        acc.maximin = colMeans(fit.maximin$pred.te * ref.TE_true > 0)
        simul.res$acc_metric1_n200[id + 4, , ix] = c(acc.maximin)
        simul.res$acc_metric2_n200[id + 4, , ix] = cor(fit.maximin$pred.te,ref.TE_true)
        simul.res$acc_metric3_n200[id + 4, , ix] = colMeans((fit.maximin$pred.te - ref.TE_true)^2)

        TRUE
      }, error = function(e) {
        message(paste("Error in delta iteration", id, "of ix =", ix, ":", e$message))
        FALSE
      })
      
      if (!result){
        sk=T} 
      if (!result) next
    }
    if (!sk){
      save(simul.res, file = file_names)
      ix=ix+1
    }
  }
  
  # boxplot
  acc_data = data.frame(acc = rbind(apply(simul.res$acc_metric1_n200, c(1, 3), min, na.rm = T), 
                                    apply(simul.res$acc_metric1_n200, c(1, 3), quantile, probs = 0.01, na.rm = T),  # 1% quantile
                                    apply(simul.res$acc_metric1_n200, c(1, 3), quantile, probs = 0.02, na.rm = T)),   # 2% quantile
                        cor = rbind(apply(simul.res$acc_metric2_n200, c(1, 3), min, na.rm = T), 
                                     apply(simul.res$acc_metric2_n200, c(1, 3), quantile, probs = 0.01, na.rm = T),  # 1% quantile
                                     apply(simul.res$acc_metric2_n200, c(1, 3), quantile, probs = 0.02, na.rm = T)),
                        mse = rbind(apply(simul.res$acc_metric3_n200, c(1, 3), max, na.rm = T), 
                                    apply(simul.res$acc_metric3_n200, c(1, 3), quantile, probs = 0.99, na.rm = T),  # 1% quantile
                                    apply(simul.res$acc_metric3_n200, c(1, 3), quantile, probs = 0.98, na.rm = T)),
                        method = rep(method_list, 3), 
                        nsub = rep(c(rep(nsub, num_method)), 3),
                        alpha = rep(c(rep(dis_alpha, num_method)), 3),
                        beta = rep(c(rep(dis_beta, num_method)), 3),
                        kappa = rep(c(rep(kappa, num_method)), 3),
                        r = rep(c(rep(r, num_method)), 3), 
                        type = rep(c("worst", "lower 1% quantile", "lower 2% quantile"), each = num_method))
  
  acc_data$type <- factor(acc_data$type, levels = c("worst", "lower 1% quantile", "lower 2% quantile"))
  plot_df = rbind(plot_df, acc_data)
}

# Create the plot
acc_data_long = plot_df %>% pivot_longer(cols = starts_with("acc"), names_to = NULL, names_prefix = "acc.", values_to = "acc") %>%
  pivot_longer(cols = starts_with("cor"), names_to = NULL, names_prefix = "cor.", values_to = "cor") %>%
  pivot_longer(cols = starts_with("mse"), names_to = "rep", names_prefix = "mse.", values_to = "mse") %>%
  pivot_longer(c("acc","cor","mse"),
               names_to = 'metric',
               values_to = "value")
p = acc_data_long %>% filter (metric!="mse") %>% 
  mutate(method = factor(method, levels = method_list)) %>% 
  ggplot(aes(x = factor(r), y = value, fill = method)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  facet_wrap(metric~type, ncol = 3, scale = "free_y") +        # Create 2x2 layout
  #scale_x_discrete(labels = x_labels) +  # Custom x-axis labels
  labs(x = "r", y = "Metric", fill = "Sample size") +  # Custom legend title
  theme_minimal() +
  # ylim(0.2, 0.75) +  # Set y-axis limits
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Angle x-axis labels
    strip.text = element_text(size = 12),  # Increase facet label size
    legend.position = "bottom"                # Place the legend on top
  )
p

ggsave(fig_n,plot=p,width = 12, height = 6, dpi = 500)


