#PROGRAM 11.1
#Sample averages by treatment level
#Data from Figures 11.1 and 11.2


##Install packages "Hmisc" and "rms" if necessary.

library(Hmisc)
library(rms)

A<-c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
Y<-c(200,150,220,110,50,180,90,170,170,30,70,110,80,50,10,20)

plot(A,Y)
describe(Y[A==1])
describe(Y[A==0])

A2<-c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
Y2<-c(110,80,50,40,170,30,70,50,110,50,180,130,200,150,220,210)


plot(A2,Y2)
describe(Y2[A2==1])
describe(Y2[A2==2])
describe(Y2[A2==3])
describe(Y2[A2==4])

#PROGRAM 11.2
#2-parameter linear model
#Data from Figures 11.3 and 11.1

A3<-c(3,11,17,23,29,37,41,53,67,79,83,97,60,71,15,45)
Y3<-c(21,54,33,101,85,65,157,120,111,200,140,220,230,217,11,190)

ols(Y3~A3)

ols(Y~A)


#PROGRAM 11.3
#3-parameter linear model
#Data from Figure 11.3

Asq<-A3*A3

ols(Y3~A3+Asq)
