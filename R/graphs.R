res <- readRDS("results/results_10000.RDS")
str(res)

library(ggplot2)
library(data.table)
library(wesanderson)

names(res)
c("distP+b2_recr")
data.table(res$`distP+b1_recr`$results$pred)

i_comb = "b2_recr"

par(mfrow = c(2,2),pty="s")
for(i_comb in c("b1_recr", "b2_recr", "b1_g", "b2_g")){
  pred = res[[i_comb]]$results$pred
  obs = res[[i_comb]]$results$obs
  p_range = range(c(pred,obs))
  plot(x = pred, y = obs,
       xlab = paste("observed",i_comb),
       ylab = paste("predicted",i_comb),
       asp = 1, xlim = p_range, ylim = p_range
       )
  abline(0,1)
}

par(mfrow = c(2,2),pty="s")

for(i_comb in c("b1_recr", "b2_recr", "b1_g", "b2_g")){
  i_comb = "b1_g"
  i_comb2 = paste0("distP+b0_e+",i_comb)
  pred = res[[i_comb2]]$results$pred
  obs = res[[i_comb2]]$results$obs
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  ggplot(
    data = data.frame(x = obs[,1], y = obs[,3], z = log((pred[,3]-obs[,3])^2)), 
    aes(cut(x, breaks = 20), cut(y, breaks = 20), fill = z))+
    geom_tile()+
    labs(
      x =  paste("disturbance"),
      y = paste("observed",i_comb)
      )+
    scale_fill_gradientn(colours = pal)
}




ggplot(
  data = data.frame(
    pred = res[[i_comb]]$results$pred,
    obs = res[[i_comb]]$results$obs
  ), aes(x = obs, y = pred))+
  geom_point()


dt <- data.table(
  obs = res$`distP+b1_recr`$results$pred[,1]
)

plot(tail(res$distP$data$distP,3000), res$distP$results$obs)

ggplot(
  data = data.frame(
    distP_obs = res$distP$results$obs,
    recr_pred = res$b2_recr$results$pred,
    recr_obs = res$b2_recr$results$obs
  ), aes(cut(recr_obs, breaks = 10), cut(distP_obs, breaks = 10), fill = (recr_pred-recr_obs)^2))+
  geom_tile()


ggplot(
  data = data.frame(
    distP_obs = res$b2_recr$results$obs,
    recr_pred = res$b2_recr$results$pred,
    recr_obs = res$b2_recr$results$obs
  ), aes(cut(recr_obs, breaks = 10), cut(distP_obs, breaks = 10), fill = (recr_pred-recr_obs)^2))+
  geom_tile()




