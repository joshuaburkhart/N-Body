#USAGE:
#./nbody_plot.sh <data file>
#
#EXAMPLE:
#./nbody_plot.sh test_job.o74748
#
#NOTES:
#Use the .o file generated from the ACISS queuing system as the data file.

R --no-save << EOT
A = read.table("$1")
par(pty="s")
m = max(A)
plot(c(-m,m),c(-m,m),xlab="X",ylab="Y",main="N Body Project",pch=" ")
points(A[,2],A[,3],col="red")
points(A[,11],A[,12],col="red")
points(A[,14],A[,15],col="red")
points(A[,17],A[,18],col="red")
points(A[,20],A[,21],col="red")
quit("no")
EOT
evince "./Rplots.pdf"
