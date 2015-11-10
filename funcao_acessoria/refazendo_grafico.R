#### REFAZENDO GRAFICO PARA POSTER SIICUSP ####

#### FUNCAO CONTAGEM 
conta.sp=function(x)
{
  length(unique(x))
}
#### FUNCAO XI AO LONGO DO TEMPO + FUNCAO PARA PLOTAR GRAFICO
fert.t1 <- function(cod.sp,n.propag,fun=mean)
{ 
  especie <- unique(cod.sp[,1])
  sp.level<-factor(cod.sp,levels=cod.sp)
  t.a<-function(x){tapply(n.propag[,x],factor(cod.sp[,x],levels=especie),fun)}
  res<-sapply(1:ncol(n.propag), t.a)
  colnames(res) <- colnames(n.propag)
  rownames(res) <- paste("sp",especie, sep="")
  return(res)
}

matplot(t(tab5),type="l",xlab="Ciclo",ylab="",axes=F)
axis(1,at=c(0,2000,4000,6000,8000,10000),labels=c("0","200 mil","400 mil","600 mil","800 mil","1 milhão"))
axis(2,at=c(0,5,10,15,20,25,30),labels=c("","5","10","15","20","25","30"),las=2)
mtext("xi",2,las=1,line=2.7)

#### APLICANDO A FUNCAO XI AO LONGO DO TEMPO NAS SIMULACOES 
#### PREPARANDO PARA PLOTAR GRAFICO QUE QUERO

#### sd0
tab0 <- fert.t1(simula_sd_0$sp.list,simula_sd_0$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo0 <- apply(simula_sd_0$sp.list,2,conta.sp)
riqueza.tempo0[riqueza.tempo0==49]
tempo0<-as.integer(attributes(riqueza.tempo0[riqueza.tempo0==49][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre0 <- sd(as.vector(tab0[,tempo0]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos0 <- fert.t1(simula_sd_0$sp.list,simula_sd_0$sementes,sd)
sd.intra.todos0[,tempo0]
sd.intra0 <- mean(as.vector(sd.intra.todos0[,tempo0]),na.rm=T)
sd.intra.sd0 <- sd(as.vector(sd.intra.todos0[,tempo0]),na.rm=T)

#### sd0.5
tab0.5 <- fert.t1(simula_sd_0.5$sp.list,simula_sd_0.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo0.5 <- apply(simula_sd_0.5$sp.list,2,conta.sp)
riqueza.tempo0.5[riqueza.tempo0.5==50]
tempo0.5<-as.integer(attributes(riqueza.tempo0.5[riqueza.tempo0.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre0.5 <- sd(as.vector(tab0.5[,tempo0.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos0.5 <- fert.t1(simula_sd_0.5$sp.list,simula_sd_0.5$sementes,sd)
sd.intra.todos0.5[,tempo0.5]
sd.intra0.5 <- mean(as.vector(sd.intra.todos0.5[,tempo0.5]),na.rm=T)
sd.intra.sd0.5 <- sd(as.vector(sd.intra.todos0.5[,tempo0.5]),na.rm=T)

#### sd1
tab1 <- fert.t1(simula_sd_1$sp.list,simula_sd_1$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo1 <- apply(simula_sd_1$sp.list,2,conta.sp)
riqueza.tempo1[riqueza.tempo1==50]
tempo1<-as.integer(attributes(riqueza.tempo1[riqueza.tempo1==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre1 <- sd(as.vector(tab1[,tempo1]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos1 <- fert.t1(simula_sd_1$sp.list,simula_sd_1$sementes,sd)
sd.intra.todos1[,tempo1]
sd.intra1 <- mean(as.vector(sd.intra.todos1[,tempo1]),na.rm=T)
sd.intra.sd1 <- sd(as.vector(sd.intra.todos1[,tempo1]),na.rm=T)

#### sd1.5
tab1.5 <- fert.t1(simula_sd_1.5$sp.list,simula_sd_1.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo1.5 <- apply(simula_sd_1.5$sp.list,2,conta.sp)
riqueza.tempo1.5[riqueza.tempo1.5==50]
tempo1.5<-as.integer(attributes(riqueza.tempo1.5[riqueza.tempo1.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre1.5 <- sd(as.vector(tab1.5[,tempo1.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos1.5 <- fert.t1(simula_sd_1.5$sp.list,simula_sd_1.5$sementes,sd)
sd.intra.todos1.5[,tempo1.5]
sd.intra1.5 <- mean(as.vector(sd.intra.todos1.5[,tempo1.5]),na.rm=T)
sd.intra.sd1.5 <- sd(as.vector(sd.intra.todos1.5[,tempo1.5]),na.rm=T)

#### sd2
tab2 <- fert.t1(simula_sd_2$sp.list,simula_sd_2$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo2 <- apply(simula_sd_2$sp.list,2,conta.sp)
riqueza.tempo2[riqueza.tempo2==50]
tempo2<-as.integer(attributes(riqueza.tempo2[riqueza.tempo2==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre2 <- sd(as.vector(tab2[,tempo2]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos2 <- fert.t1(simula_sd_2$sp.list,simula_sd_2$sementes,sd)
sd.intra.todos2[,tempo2]
sd.intra2 <- mean(as.vector(sd.intra.todos2[,tempo2]),na.rm=T)
sd.intra.sd2 <- sd(as.vector(sd.intra.todos2[,tempo2]),na.rm=T)

#### sd2.5
tab2.5 <- fert.t1(simula_sd_2.5$sp.list,simula_sd_2.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo2.5 <- apply(simula_sd_2.5$sp.list,2,conta.sp)
riqueza.tempo2.5[riqueza.tempo2.5==50]
tempo2.5<-as.integer(attributes(riqueza.tempo2.5[riqueza.tempo2.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre2.5 <- sd(as.vector(tab2.5[,tempo2.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos2.5 <- fert.t1(simula_sd_2.5$sp.list,simula_sd_2.5$sementes,sd)
sd.intra.todos2.5[,tempo2.5]
sd.intra2.5 <- mean(as.vector(sd.intra.todos2.5[,tempo2.5]),na.rm=T)
sd.intra.sd2.5 <- sd(as.vector(sd.intra.todos2.5[,tempo2.5]),na.rm=T)

#### sd3
tab3 <- fert.t1(simula_sd_3$sp.list,simula_sd_3$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo3 <- apply(simula_sd_3$sp.list,2,conta.sp)
riqueza.tempo3[riqueza.tempo3==50]
tempo3<-as.integer(attributes(riqueza.tempo3[riqueza.tempo3==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre3 <- sd(as.vector(tab3[,tempo3]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos3 <- fert.t1(simula_sd_3$sp.list,simula_sd_3$sementes,sd)
sd.intra.todos3[,tempo3]
sd.intra3 <- mean(as.vector(sd.intra.todos3[,tempo3]),na.rm=T)
sd.intra.sd3 <- sd(as.vector(sd.intra.todos3[,tempo3]),na.rm=T)

#### sd3.5
tab3.5 <- fert.t1(simula_sd_3.5$sp.list,simula_sd_3.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo3.5 <- apply(simula_sd_3.5$sp.list,2,conta.sp)
riqueza.tempo3.5[riqueza.tempo3.5==50]
tempo3.5<-as.integer(attributes(riqueza.tempo3.5[riqueza.tempo3.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre3.5 <- sd(as.vector(tab3.5[,tempo3.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos3.5 <- fert.t1(simula_sd_3.5$sp.list,simula_sd_3.5$sementes,sd)
sd.intra.todos3.5[,tempo3.5]
sd.intra3.5 <- mean(as.vector(sd.intra.todos3.5[,tempo3.5]),na.rm=T)
sd.intra.sd3.5 <- sd(as.vector(sd.intra.todos3.5[,tempo3.5]),na.rm=T)

#### sd3.5
tab3.5 <- fert.t1(simula_sd_3.5$sp.list,simula_sd_3.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo3.5 <- apply(simula_sd_3.5$sp.list,2,conta.sp)
riqueza.tempo3.5[riqueza.tempo3.5==50]
tempo3.5<-as.integer(attributes(riqueza.tempo3.5[riqueza.tempo3.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre3.5 <- sd(as.vector(tab3.5[,tempo3.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos3.5 <- fert.t1(simula_sd_3.5$sp.list,simula_sd_3.5$sementes,sd)
sd.intra.todos3.5[,tempo3.5]
sd.intra3.5 <- mean(as.vector(sd.intra.todos3.5[,tempo3.5]),na.rm=T)
sd.intra.sd3.5 <- sd(as.vector(sd.intra.todos3.5[,tempo3.5]),na.rm=T)

#### sd4
tab4 <- fert.t1(simula_sd_4$sp.list,simula_sd_4$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo4 <- apply(simula_sd_4$sp.list,2,conta.sp)
riqueza.tempo4[riqueza.tempo4==50]
tempo4<-as.integer(attributes(riqueza.tempo4[riqueza.tempo4==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre4 <- sd(as.vector(tab4[,tempo4]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos4 <- fert.t1(simula_sd_4$sp.list,simula_sd_4$sementes,sd)
sd.intra.todos4[,tempo4]
sd.intra4 <- mean(as.vector(sd.intra.todos4[,tempo4]),na.rm=T)
sd.intra.sd4 <- sd(as.vector(sd.intra.todos4[,tempo4]),na.rm=T)

#### sd4.5
tab4.5 <- fert.t1(simula_sd_4.5$sp.list,simula_sd_4.5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo4.5 <- apply(simula_sd_4.5$sp.list,2,conta.sp)
riqueza.tempo4.5[riqueza.tempo4.5==50]
tempo4.5<-as.integer(attributes(riqueza.tempo4.5[riqueza.tempo4.5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre4.5 <- sd(as.vector(tab4.5[,tempo4.5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos4.5 <- fert.t1(simula_sd_4.5$sp.list,simula_sd_4.5$sementes,sd)
sd.intra.todos4.5[,tempo4.5]
sd.intra4.5 <- mean(as.vector(sd.intra.todos4.5[,tempo4.5]),na.rm=T)
sd.intra.sd4.5 <- sd(as.vector(sd.intra.todos4.5[,tempo4.5]),na.rm=T)

#### sd5
tab5 <- fert.t1(simula_sd_5$sp.list,simula_sd_5$sementes,fun=mean)
# Para saber tempo em que metade das especies foi extinta
riqueza.tempo5 <- apply(simula_sd_5$sp.list,2,conta.sp)
riqueza.tempo5[riqueza.tempo5==50]
tempo5<-as.integer(attributes(riqueza.tempo5[riqueza.tempo5==50][1])$names[1])/100+1
# Para calcular a variacao inter-especifica
sd.entre5 <- sd(as.vector(tab5[,tempo5]),na.rm=T)
# Para calcular a variacao intra-especifica
sd.intra.todos5 <- fert.t1(simula_sd_5$sp.list,simula_sd_5$sementes,sd)
sd.intra.todos5[,tempo5]
sd.intra5 <- mean(as.vector(sd.intra.todos5[,tempo5]),na.rm=T)
sd.intra.sd5 <- sd(as.vector(sd.intra.todos5[,tempo5]),na.rm=T)

# Montar o gráfico
plot(x=0,xlim=c(0,5),ylim=c(0,5),xlab="SD",ylab="Desvio padrão do xi",type="n",bty="l")
points(x=list(attributes(simula_sd_0)$start$sd,attributes(simula_sd_0.5)$start$sd,attributes(simula_sd_1)$start$sd,attributes(simula_sd_1.5)$start$sd,attributes(simula_sd_2)$start$sd,attributes(simula_sd_2.5)$start$sd,attributes(simula_sd_3)$start$sd,attributes(simula_sd_3.5)$start$sd,attributes(simula_sd_4)$start$sd,attributes(simula_sd_4.5)$start$sd,attributes(simula_sd_5)$start$sd),y=list(sd.entre0,sd.entre0.5,sd.entre1,sd.entre1.5,sd.entre2,sd.entre2.5,sd.entre3,sd.entre3.5,sd.entre4,sd.entre4.5,sd.entre5),col="blue",type="l")
points(x=list(attributes(simula_sd_0)$start$sd,attributes(simula_sd_0.5)$start$sd,attributes(simula_sd_1)$start$sd,attributes(simula_sd_1.5)$start$sd,attributes(simula_sd_2)$start$sd,attributes(simula_sd_2.5)$start$sd,attributes(simula_sd_3)$start$sd,attributes(simula_sd_3.5)$start$sd,attributes(simula_sd_4)$start$sd,attributes(simula_sd_4.5)$start$sd,attributes(simula_sd_5)$start$sd),y=list(sd.intra0,sd.intra0.5,sd.intra1,sd.intra1.5,sd.intra2,sd.intra2.5,sd.intra3,sd.intra3.5,sd.intra4,sd.intra4.5,sd.intra5),col="red",type="l")
xy <- locator(1)
legend(xy, legend=c("Intra-específico", "Interespecífico"), pch=c(16, 16), col=c("red","blue") , bty="n")

#### FUNCAO GRAFICO SAD AO LONGO DO TEMPO 
graf.abund=function(dados=resulta,...)
{
  info=attributes(dados)
  n.prop=info$start$nprop
  cv=info$start$cv
  tempo=dados$tempo
  ntempo=length(tempo)
  nspp=length(unique(dados$sp.list[,1]))
  #nmax=max(unlist(dados))
  nmax=max(table(dados$sp.list[,dim(dados$sp.list)[2]]))
  plot(x=c(1,nspp),y=c(1, nmax),log="y", type="n", ylab= "Abundance", xlab="Rank order", cex.lab=1.2, cex.axis=1.2, sub= paste("total seeds=", n.prop, "\t cv = ", cv),... )
  stempos=round(seq(2, (ntempo-100),length.out=100 ))
  colors=rainbow(length(stempos)*10)
  ncol=length(colors)
  ncol.meio=round(ncol*0.5)
  i=ncol.meio
  for (t in stempos )
  {
    i=i+1
    points(sort(table(factor(dados$sp.list[,t], levels=1:nspp)), decreasing=TRUE),type="l", lwd=0.5, col=colors[i])
  }
  points(sort(table(factor(dados$sp.list[,1], levels=1:nspp)), decreasing=TRUE), type="l", col="green", lwd=2)
  points(sort(table(factor(dados$sp.list[,ncol.meio], levels=1:nspp)), decreasing=TRUE), type="l", col="blue", lwd=2)
  points(sort(table(factor(dados$sp.list[,dim(dados$sp.list)[2]], levels=1:nspp)), decreasing=TRUE), type="l", col="red", lwd=2)
  legend("topright", lty=1, col=c("green", "blue", "red"), bty="n", legend=c("start", "middle", "end") )
}

graf.abund(simula_sd_1)


