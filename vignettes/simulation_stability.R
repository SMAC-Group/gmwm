library(GMWM)
N = 1000
model = AR1(phi = .99, sigma = 1) + WN(sigma2=20)
B = 100
phi = matrix(NA,B,10)
sigmaAR = phi
sigmaWN = phi

for (i in 1:B){
	set.seed(i)
	xt = gen.gts(model, N)
	mod1 = gmwm(AR1() + WN(), xt, model.type = "ssm")
	mod2 = gmwm(model, xt, model.type = "ssm")
	mod3 = gmwm(AR1() + WN(), xt, model.type = "ssm")
	mod4 = gmwm(model, xt, model.type = "ssm")
	mod5 = gmwm(AR1() + WN(), xt, model.type = "ssm", compute.v = "bootstrap")
	mod6 = gmwm(model, xt, model.type = "ssm", compute.v = "bootstrap")
	mod7 = gmwm(AR1() + WN(), xt, model.type = "ssm", compute.v = "bootstrap", K = 2)
	mod8 = gmwm(model, xt, model.type = "ssm", compute.v = "bootstrap", K = 2)
	mod9 = gmwm(AR1() + WN(), xt, model.type = "ssm", compute.v = "bootstrap", K = 3)
	mod10 = gmwm(model, xt, model.type = "ssm", compute.v = "bootstrap", K = 3)
	
	phi[i,] = c(mod1$estimate[1],mod2$estimate[1],mod3$estimate[1],mod4$estimate[1],mod5$estimate[1],
		mod6$estimate[1],mod7$estimate[1],mod8$estimate[1],mod9$estimate[1],mod10$estimate[1])
	sigmaAR[i,] = c(mod1$estimate[2],mod2$estimate[2],mod3$estimate[2],mod4$estimate[2],mod5$estimate[2],
		mod6$estimate[2],mod7$estimate[2],mod8$estimate[2],mod9$estimate[2],mod10$estimate[2])
	sigmaWN[i,] = c(mod1$estimate[3],mod2$estimate[3],mod3$estimate[3],mod4$estimate[3],mod5$estimate[3],
		mod6$estimate[3],mod7$estimate[3],mod8$estimate[3],mod9$estimate[3],mod10$estimate[3])
}


quartz()
par(mfrow = c(1,3))
boxplot(phi, col = "lightgrey")
abline(h = 0.99, lwd = 2, col = 2)

boxplot(sigmaAR, col = "lightgrey")
abline(h = 1, lwd = 2, col = 2)

boxplot(sigmaWN, col = "lightgrey")
abline(h = 20, lwd = 2, col = 2)

quartz()
par(mfrow = c(1,3))
boxplot(phi, col = "lightgrey", ylim = c(0.95,1))
abline(h = 0.99, lwd = 2, col = 2)

boxplot(sigmaAR, col = "lightgrey", ylim = c(0.5,1.7))
abline(h = 1, lwd = 2, col = 2)

boxplot(sigmaWN, col = "lightgrey", ylim = c(17,23))
abline(h = 20, lwd = 2, col = 2)
