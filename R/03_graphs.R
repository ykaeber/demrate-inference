library(ggplot2)
library(data.table)
library(wesanderson)
library(forcats)

res <- readRDS("results/results_10000.RDS")
figures_path = "figures"

# res <- readRDS("results/MCMC_1000.RDS")
# figure_path = "figures_mcmc"

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

# calculate correlations
res_dt[,":="(
  cor = cor(obs,pred)
), by = .(comb, par)]
res_dt[is.na(cor),cor := 0,]
res_dt[,r2:=cor^2]

# fit linear model for identifying effect of parameters for the r2
fm <- lm(r2 ~ par*(distP + b0_e + b1_recr + b2_recr + b1_g + b2_g) -1, data = res_dt[par != comb])
#fm <- lm(r2 ~ par*comb -1, data = res_dt)

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
i_par = "b1_g"
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



# 
# par(mfrow = c(2,3),pty="s")
# one_one_allfixed_dt <- data.table()
# for(i_comb in c("distP", "b0_e", "b1_recr", "b2_recr", "b1_g", "b2_g")){
#   pred = res[[i_comb]]$results$pred
#   obs = res[[i_comb]]$results$obs
#   p_range = range(c(pred,obs))
#   plot(x = pred, y = obs,
#        xlab = paste("observed",i_comb),
#        ylab = paste("predicted",i_comb),
#        asp = 1, xlim = p_range, ylim = p_range
#   )
#   abline(0,1)
#   }
# 
# 
# i_comb2 = paste0("distP+b0_e+",i_comb)
# pred = res[[i_comb2]]$results$pred
# obs = res[[i_comb2]]$results$obs
# 
# 
# names(res)
# names(res)[names(res) == "b2_recr"]
# c("distP+b2_recr")
# data.table(res$`distP+b1_recr`$results$pred)
# 
# i_comb = "b2_recr"
# 
# par(mfrow = c(2,2),pty="s")
# for(i_comb in c("b1_recr", "b2_recr", "b1_g", "b2_g")){
#   pred = res[[i_comb]]$results$pred
#   obs = res[[i_comb]]$results$obs
#   p_range = range(c(pred,obs))
#   plot(x = pred, y = obs,
#        xlab = paste("observed",i_comb),
#        ylab = paste("predicted",i_comb),
#        asp = 1, xlim = p_range, ylim = p_range
#        )
#   abline(0,1)
# }
# 
# par(mfrow = c(2,2),pty="s")
# 
# for(i_comb in c("b1_recr", "b2_recr", "b1_g", "b2_g")){
#   #i_comb = "b1_g"
#   i_comb2 = paste0("distP+b0_e+",i_comb)
#   pred = res[[i_comb2]]$results$pred
#   obs = res[[i_comb2]]$results$obs
#   pal <- wes_palette("Zissou1", 100, type = "continuous")
#   p <- ggplot(
#     data = data.frame(x = obs[,1], y = obs[,3], z = log((pred[,3]-obs[,3])^2)), 
#     aes(cut(x, breaks = 20), cut(y, breaks = 20), fill = z))+
#     geom_tile()+
#     labs(
#       x =  paste("disturbance"),
#       y = paste("observed",i_comb)
#       )+
#     scale_fill_gradientn(colours = pal)
#   print(p)
# }
# 
# 
# 
# 
# ggplot(
#   data = data.frame(
#     pred = res[[i_comb]]$results$pred,
#     obs = res[[i_comb]]$results$obs
#   ), aes(x = obs, y = pred))+
#   geom_point()
# 
# 
# dt <- data.table(
#   obs = res$`distP+b1_recr`$results$pred[,1]
# )
# 
# plot(tail(res$distP$data$distP,3000), res$distP$results$obs)
# 
# ggplot(
#   data = data.frame(
#     distP_obs = res$distP$results$obs,
#     recr_pred = res$b2_recr$results$pred,
#     recr_obs = res$b2_recr$results$obs
#   ), aes(cut(recr_obs, breaks = 10), cut(distP_obs, breaks = 10), fill = (recr_pred-recr_obs)^2))+
#   geom_tile()
# 
# 
# ggplot(
#   data = data.frame(
#     distP_obs = res$b2_recr$results$obs,
#     recr_pred = res$b2_recr$results$pred,
#     recr_obs = res$b2_recr$results$obs
#   ), aes(cut(recr_obs, breaks = 10), cut(distP_obs, breaks = 10), fill = (recr_pred-recr_obs)^2))+
#   geom_tile()
# 



