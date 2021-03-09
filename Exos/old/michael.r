postscript("speedup.ps",horizontal=F)
logx<-0:5
y<-c(1,1.9,3.9,6.6,3.9,3.3)
logy<-log(y,base=2)
plot(logx,logy,type="o",ylim=range(c(logx,logy)),
main="Speedup",xaxt="n",yaxt="n",xlab="OMP_NUM_THREADS",
ylab="speedup")
points(logx,logx,type="o",col="red")
axis(1,at=logx,lab=2^logx)
axis(2,at=logx,lab=2^logx)
for(i in 1:6)
{
 abline(h=logx[i],lty=2)
 abline(v=logx[i],lty=2)
}
dev.off()
