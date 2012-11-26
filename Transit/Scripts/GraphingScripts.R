# Alex's graphing tools


# Loop over friction factor #
y <- vector(length=length(bval))
for(i in 1:length(bval)) {
	grav <- gravityModel(DATA$BOARD, DATA$ALIGHT , costMatrix(1,1), bval[i], 1E-9)
	y[i] <- output.cleaner(grav)[8] #the number in [] can be from 1:8#
}

# Code for plotting #

> plot(bval, y1, cex=0.5, xlab = "friction factor (b)", ylab = "transfer volume")
> points(bval, y5, cex=0.2, col="green")
> points(bval, y2, cex=0.2, col="red")
> points(bval, y8, cex=0.5, col="blue")
> points(bval, y7, cex=0.2, col="pink")
> points(bval, y3, cex=0.2, col="brown")
> points(bval, y6, cex=0.2, col="purple")
> points(bval, y4, cex=0.2, col="orange")
