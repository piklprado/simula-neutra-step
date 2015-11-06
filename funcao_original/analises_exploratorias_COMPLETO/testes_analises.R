#######################################
####### ANALISES EXPLORATORIAS ########
#######################################

###### OBSERVACAO 1: OLHAR O SCRIPT funcao_ajudaneutra.r PARA MAIS INFORMACOES SOBRE CADA GRAFICO DE ANALISE EXPLORATORIA
###### OBSERVACAO 2: AQUI, simulacao = RESULTADO (LISTA FINAL) DE DETERMINADA SIMULACAO COM simula.neutra.step

############################################################################
####### GRAFICO ESFORCO REPRODUTIVO INSTANTANEO POR ESPÉCIE X TEMPO  #######
############################################################################

##### VERSÃO ALE/PI RESUMIDA #####
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
matplot(t(TABELA RESULTANTE DE fert.t1),type="l",xlab="Tempo",ylab="Esforço reprodutivo instantâneo",main="Variação do ERI por espécie",labels=c())

#######################################
####### GRAFICO RIQUEZA X TEMPO #######
#######################################
conta.sp=function(x)
{
  length(unique(x))
}
riqueza.tempo <- apply(simulacao$sp.list,2,conta.sp)
plot(riqueza.tempo,type="l",xlab="Tempo",ylab="Riqueza",main="Variação da riqueza")

###################################
####### GRAFICO SAD X TEMPO #######
###################################
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
  points(sort(table(factor(dados$sp.list[,dim(dados$sp.list)[2]], levels=1:nspp)), decreasing=TRUE), type="p", col="red", lwd=2)
  legend("topright", lty=1, col=c("green", "blue", "red"), bty="n", legend=c("start", "middle", "end") )
}

graf.abund(simulacao$sp.list)

###################################
####### GRAFICO CV X VAR ERI ######
###################################

# Qual é o tempo em que metade das espécies são extintas?
# metade = riqueza inicial/2
riqueza.tempo <- apply(simulacao$sp.list,2,conta.sp)
riqueza.tempo[riqueza.tempo==metade]
# talvez: riqueza.tempo[riqueza.tempo==metade+1]
## Exemplo: conclusão: tempo 7000 -> usar tempo 8

# Para calcular a variacao inter-especifica
as.vector(tabela.final[,8])
mean(as.vector(tabela.final[,8]),na.rm=T)
sd.entre <- sd(as.vector(tabela.final[,8]),na.rm=T)

# Para calcular a variacao intra-especifica
sd.intra.todos <- fert.t1(simulacao$sp.list,simulacao$sementes,sd)
sd.intra.todos[,8]
sd.intra <- mean(as.vector(sd.intra.todos[,8]),na.rm=T)
sd.intra.sd <- sd(as.vector(sd.intra.todos[,8]),na.rm=T)

# Para chamar o cv 
attributes(simulacao)$start$cv

# Montar o gráfico
plot(x=0,xlim=c(0,10),ylim=c(0,5),xlab="CV",ylab="Desvio padrão do ERI",main="Variação intra e interespecífica do ERI em relação ao CV",type="n")
points(x=list(attributes(simulacao)$start$cv),y=list(sd.entre),col="blue",type="l")
points(x=list(attributes(simulacao)$start$cv),y=list(sd.intra),col="red",type="l")
