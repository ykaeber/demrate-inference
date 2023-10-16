library(sjSDM)
torch = sjSDM:::pkg.env$torch
dt = torch$float32
as_ft = function(x) torch$tensor(x, dtype = dt)
as_ftg = function(x) torch$tensor(x, dtype = dt, requires_grad = TRUE)



