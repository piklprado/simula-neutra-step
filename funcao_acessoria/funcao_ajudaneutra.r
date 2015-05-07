##Funções acessorias e graficas do Modelo Neutro
## Conta espécie do arquivo de saída da simulação
# o objeto de entrada é o $sp.list do arquivo de saída
##AAO junho 2010

conta.sp=function(x)
 {
 length(unique(x))
 }


#################################################
#### Gráfico da variação na SAD ao longo do tempo
## entrada: dados é a tabela sp.list do objeto resultado da simulação e info são os attributos do resultado  
#################################################
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
##################################$$
#Rank-abundance Plot
############################
#Dado um vetor de abundăncia de espécies, faz um gráfico simples de Whitaker (log abundâncias x ranque).
########### possibilidade: as.vector(table(simulacao$sp.list[,101])), para o rank de abundancia do tempo final da simulacao
w.plot <- function(x,...)
{
  plot(1:length(x),sort(x,decreasing=T),
       ylab="Abundance",log="y",...)
}

################################
#Box-plot e rank-abundance plot
################################
#Funçăo feita para explorar graficamente a distribuiçăo da fertilidade. Plota box-plots da fertlidade 
#dos indivíduso de cada espécie, sobre um gr[áfico de rank-abundância da fertilidade de cada espécie, 
#num dado instante gravado da simulaçăo (argumento tempo).
## dados é o objeto (list) resultado da simulação

box.w <- function(dados, tempo=NULL )
 {
 if(is.null(tempo))
  {tempo=dim(dados$sp.list)[2]}
  #tempo,cod.sp,n.propag, info=attributes(dados)
  ast=attributes(dados)$start
  sp <- dados$sp.list[,tempo]
  prop <- dados$sementes[,tempo]
  ab <- tapply(sp,sp,length)
  par(mfrow=c(2,1), mar=c(2,4,4,3), mgp=c(2,1,0))
  w.plot(ab, xlab="", xaxt="n", main=paste("NEUTRALITY SIMULATION\n", paste(names(ast), ast, sep="=", collapse="; "), "\n(time = ",tempo, ")"), cex.main=0.8)
  par(mar=c(5,4,0,3))#, new=TRUE, oma=c(0,2,2,2))
  boxplot(prop~factor(sp,levels=names(ab)[order(ab,decreasing=T)]),xlab="Species Codes", ylab="Number of Seeds")
  par(mfrow=c(1,1))
  }
#################################
#Fertilidade Média por Espécie
##############################
#A simulaçăo grava uma matriz com as fertilidades de cada indivíduo a cada step ciclos. 
#Esta funçăo resume estes dados, calculando uma estatística descritiva (o default é a média)
# da fertilidade de cada espécie, em cada instante de tempo.

fert.t <- function(cod.sp,n.propag,fun=mean){
  especie <- unique(cod.sp[,1])
  medias <- matrix(ncol=ncol(n.propag),nrow=length(especie))
  colnames(medias) <- colnames(n.propag)
  rownames(medias) <- especie
  for(i in 1:ncol(n.propag)){
    medias[,i] <- tapply(n.propag[,i],factor(cod.sp[,i],levels=especie),fun)
  }
return(medias)
}

## uma versão resumida
# AAO 06 February 2012 
### parece não haver diferença de velocidade de processamento!
fert.t1 <- function(cod.sp,n.propag,fun=mean)
{ 
  especie <- unique(cod.sp[,1])
  #sp.level<-factor(cod.sp,levels=cod.sp)
  t.a<-function(x){tapply(n.propag[,x],factor(cod.sp[,x],levels=especie),fun)}
  res<-sapply(1:ncol(n.propag), t.a)
  colnames(res) <- colnames(n.propag)
  rownames(res) <- paste("sp",especie, sep="")
  return(res)
}

