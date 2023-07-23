library(ggplot2)
library(data.table)
library(wesanderson)
library(forcats)
library(tidytext)
library(stringr)

### Functions
confidenceInterval <- function(pred, true){
  n   <- length(pred)
  ssd <- (pred - true)**2
  sd  <- sd(ssd)
  se  <- sd / sqrt(n)
  ci  <- 1.96*se
  return(ci)
  
}


getTheilCoefs <- function(pred, true){
  n <- length(pred)
  fit <- lm(pred ~ true)
  b1 <- coef(fit)[[2]]
  a <- coef(fit)[[1]]
  pred.lm <- a + b1*true
  
  u_d     <- mean(true**2)
  u_bias  <- n*mean((pred-true))**2
  u_slope <- (b1-1)**2 * sum((true - mean(true))**2)
  u_var   <- sum((pred - pred.lm)**2)
  u_bias  <- u_bias / (n*u_d)
  u_slope <- u_slope/ (n*u_d)
  u_var   <- u_var  / (n*u_d)
  ci      <- confidenceInterval(pred, true) / (u_d)
  
  # SSD_rec <- u_bias + u_slope + u_var
  SSD <- sum((pred-true)**2)/ (n*u_d)
  SSD_rec <- u_bias + u_slope + u_var
  
  theil_coef <- list("SSD" = SSD, "SSD_rec" = SSD_rec, "bias" = u_bias,
                     "slope" = u_slope, "var" = u_var, "ci"   = ci)
  
  return(theil_coef)
}


res <- readRDS("results/results_10000.RDS")
figure_path = "figures"

# res_mcmc <- readRDS("results/MCMC_1000.RDS")
# figure_path = "figures_mcmc"
# 
# mcmc_names <- c("distP", "b0_e", "b1_recr",
# "b2_recr", "b1_g", "b2_g")
# str(res_mcmc[[1]])
# res_mcmc <- lapply(res_mcmc, function(i)  { res_mcmc = as.data.table(i)})
# res_dt_mcmc = rbindlist(res_mcmc)
# res_dt_mcmc[["par"]] <- rep(mcmc_names, length(res_mcmc))
# res_dt_mcmc[["ID"]] <- rep(1:length(res_mcmc), each = length(mcmc_names))

if(!dir.exists(figure_path)) dir.create(figure_path)

# get data into one data.table
res_dt <- data.table()
i_comb = names(res)[3]
for(i_comb in names(res)){
  pred = res[[i_comb]]$results$pred
  obs = res[[i_comb]]$results$obs
  par_names = strsplit(i_comb,"\\+")[[1]]
  for(i in 1:length(par_names)){
    res_dt <- rbind(
      res_dt, 
      data.table(
        pred = pred[,i],
        obs = obs[,i],
        par = par_names[i],
        comb = i_comb,
        distP = grepl("distP", i_comb), 
        b0_e = grepl("b0_e", i_comb), 
        b1_recr = grepl("b1_recr", i_comb), 
        b2_recr = grepl("b2_recr", i_comb), 
        b1_g = grepl("b1_g", i_comb), 
        b2_g = grepl("b2_g", i_comb)
      ))
  }
}

res_dt[,":="(
  obs_z = as.numeric(scale(obs)),
  pred_z = as.numeric(scale(pred))
),by = .(comb,par)]

res_dt[,":="(
  cor = cor(obs,pred)
), by = .(comb, par)]

#res_dt[, c("pred_z", "obs_z") := split(as.numeric(scale(c(pred, obs))), rep(1:2, each = length(pred))) , by = .(par)]

res_dt[, c("SSD", "SSD_rec", "bias", "slope", "var", "ci") := getTheilCoefs(pred, obs) , by = .(comb, par)]
res_dt[, Npars := distP+b0_e+b1_recr+b2_recr+b1_g+b2_g,]
res_dt[, comb_group := fifelse(grepl("_g", comb),"growth", fifelse(grepl("_recr", comb), "regeneration", fifelse(grepl("_e", comb), "environment", fifelse(grepl("distP", comb), "disturbance", NA_character_)))),]


thails_coeffs <- unique(res_dt[,c("Npars","comb", "par","bias", "slope", "var", "SSD_rec")])
thails_coeffs = melt(thails_coeffs, id.vars = c("comb", "par","Npars", "SSD_rec"), variable.name = "error_name", value.name = "error")


thails_coeffs[, N_g_recr := str_count(comb, "_g")*grepl("_g", par)+str_count(comb, "_recr")*grepl("_recr", par), by = .(comb, par)]
res_dt[is.na(cor),cor := 0,]
res_dt[,r2:=cor^2]

label_fun <- function(x){
  x1 = strsplit(x, "___")[[1]]
  return(x1[-length(x1)])
}

pal = as.vector(wes_palette("FantasticFox1", n = 3, type = "discrete"))
Vectorize(label_fun)
names(pal) <- c("bias", "slope", "var")
for(i_par in c("_g", "_recr")){
  if(i_par == "_g") par_lab = "growth" else par_lab = "regeneration"
  p <- 
    ggplot(
    thails_coeffs[!(par %in% c("b0_e", "distP")) & grepl(i_par, par)],
    aes(x = abs(error), y = reorder_within(comb, SSD_rec, par), color = error_name, fill = error_name))+
    scale_y_discrete(labels = Vectorize(label_fun))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(name = element_text("error\ncomponent"), values = pal)+
    scale_color_manual(name = element_text("error\ncomponent"), values = pal)+
    facet_wrap(factor(N_g_recr, levels = 1:2, labels = c(paste0("one ",par_lab," parameter fixed"), "environment+density effects unknown"))~par, scales = "free_y")+
    theme_classic()+
    #theme(legend.title = element_text("Error\ncomponent"))+
    xlab("total error")+
    ylab("inferred parameter combination")+
    coord_cartesian(xlim = range(thails_coeffs$SSD_rec))
  png(paste0(figure_path, "/comb-errors_",par_lab,".png"), res = 300, units = "in", width = 12, height = 7)
  print(p)
  dev.off()
}

pal = as.vector(wes_palette("FantasticFox1", n = 3, type = "discrete"))[c(1,3)]

infer_label <- c("one driver\n(environment OR density)", "two drivers\n(environment AND density)") 
names(pal) <- infer_label
Vectorize(label_fun)

p_dat2 <- thails_coeffs[!(par %in% c("b0_e", "distP")) & grepl("b0_e", comb) & Npars == 4]
p_dat2 <- thails_coeffs[
  (grepl("_recr", par) & comb %in% c(
    "b0_e+b1_recr+b2_recr","b0_e+b2_recr","b0_e+b1_recr",
    "b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g"
    )) | 
    (grepl("_g", par) & comb %in% c(
      "b0_e+b1_g+b2_g","b0_e+b2_g","distP+b0_e+b1_g",
      "b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g"
      ))
  ]
p_dat2[grepl("_g", par), process := "growth",]
p_dat2[grepl("_recr", par), process := "regeneration",]
p_dat2[comb %in% c("b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g") & grepl("_g", par), setting := "two processes\n(growth AND regeneration)",]
p_dat2[comb %in% c("b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g") & grepl("_recr", par), setting := "two processes\n(growth AND regeneration)",]
p_dat2[!(comb %in% c("b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g")) & grepl("_g", par), setting := "one process\n(growth OR regeneration)",]
p_dat2[!(comb %in% c("b0_e+b1_recr+b2_recr+b1_g+b2_g","b0_e+b2_recr+b2_g","b0_e+b1_recr+b1_g")) & grepl("_recr", par), setting := "one process\n(growth OR regeneration)",]
p_dat2[grepl("b1", par), driver := "density",]
p_dat2[grepl("b2", par), driver := "environment",]
# for(i_par in c("_g", "_recr")){
#   if(i_par == "_g") par_lab = "growth" else par_lab = "regeneration"
  p <- 
    ggplot(
    p_dat2[error_name == "slope"],
    aes(y = abs(SSD_rec), x = driver, 
        color = factor(N_g_recr, levels = 1:2, labels = infer_label),
        fill = factor(N_g_recr, levels = 1:2, labels = infer_label)
        ))+
    # scale_y_discrete(labels = Vectorize(label_fun))+
    geom_bar(stat = "identity", position = position_dodge2(width = 0.4))+
    scale_fill_manual(name = element_text("inference effort"), values = pal)+
    scale_color_manual(name = element_text("inference effort"), values = pal)+
    # scale_color_manual(name = element_text("disturbances"), values = pal)+
    # facet_grid(process~factor(N_g_recr, levels = 1:2, labels = infer_label), scales = "free_y")+
    # facet_grid(process~factor(N_g_recr, levels = 1:2, labels = infer_label)
    #            , scales = "free_y")+
    facet_grid(process~setting, scales = "free_y")+
    theme_bw()+
    theme(legend.title = element_blank())+
    guides(
      color = guide_legend(direction = "horizontal"),
      fill = guide_legend(direction = "horizontal")
      )+
    xlab("driver")+
    ylab("total error")+
    coord_cartesian(ylim = c(0,1))+
    theme(legend.position = "bottom")
  png(paste0(figure_path, "/error-boxplots.png"), res = 300, units = "in", width = 5, height = 5)
  print(p)
  dev.off()
# }


# fit linear model for identifying effect of parameters for the r2
fm <- lm(var ~ (distP+b0_e + b1_recr + b2_recr + b2_g)*par -1, 
         data = res_dt)
coefs_dt <- 
  data.table(
    name = names(coefficients(fm)),
    coef = coefficients(fm)
    )

png(paste0(figure_path, "/pareffect.png"), res = 300, units = "in", width = 7, height = 7)
ggplot(coefs_dt, aes(y = name, x = coef))+
  geom_point()
dev.off()

pdf(paste0(figure_path, "/all_combs.pdf"), width = 8*2, height = 4.5*2)
i_par = "b2_recr"
for(i_par in c("distP", "b0_e", "b1_recr", "b2_recr", "b1_g", "b2_g")){
  cat(i_par, "\n")
  p_dat <- res_dt[par == i_par]
  axis_limits = range(c(p_dat$obs))
  p <- ggplot(p_dat, aes(x = obs, y = pred))+
    geom_point(alpha = 0.05)+
    coord_fixed(xlim = axis_limits, ylim = axis_limits)+
    facet_wrap(~fct_reorder(comb, cor), ncol = 8)+
    geom_abline(slope = 1, intercept = 0)+
    geom_smooth(method = "lm")+
    ggtitle(i_par)
  print(p)
}
dev.off()
