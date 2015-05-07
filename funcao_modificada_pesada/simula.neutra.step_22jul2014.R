########################################
######### Simulando cenario neutro Hubbel
########################################
# versao modificada por Luísa Novara, em julho de 2014
# versao anterior criada por A. A de Oliveira, em outubro de 2009, e modificada por Paulo Inácio, em outubro de 2009
########################################
# funcao principal
########## Definidos ##########
#      S =  número de espécies da comunidade : inteiro positivo
#   *
#      j = número inicial de indivíduos por espécie: inteiro positivo
#   *
#      xi =  número de propágulos que o indivíduo i produz por intervalo (esforço reprodutivo instantâneo) : inteiro positivo
#   *
#      cv =  proporçăo de variaçăo no esforço reprodutivo instantâneo em relaçăo a sua esperança: real positivo
#    
########## Deduzidos ##########  
#      X =  total de propágulos que um indivíduo produz (esforço reprodutivo total): inteiro positivo
#      O X é dado por xi*ciclo ou por xi*J
#   *
#      Dado que o esforço reprodutivo total é o mesmo para todos os indivíduos (X), o tempo médio de vida de cada indivíduo é E[ti] = X/xi
#   *
#      A cada intervalo cada indivíduo tem uma probabilidade fixa pi de morrer. Portanto, a probabilidade de sobreviver t intervalos é dada por uma distribuiçăo geométrica com suporte 0,8. A esperança desta distribuiçăo é E[ti]=pi-1, portanto:
#            pi = xi/X
#   *
#      Dado s igual para todas as espécies, e uma distribuiçăo gaussiana de caracteres contínuos, o número de sementes por intervalo xi da prole de um indivíduo pode ser definido como uma variável normal, com parâmetros
#            µ = xi
#            s = cv·µ 

simula.neutra.step=function(S= 100, j=10, xi=10, cv=0.1, ciclo=1e6, step=100)
{
t0=proc.time()[[3]]
  ## Tamanho da comunidade
  J <- S*j
  ## Esforco reprodutivo total
  X <- xi*J 
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##Matriz de propagulos produzidos por individuo em cada ciclo
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step)
  ## Matriz de probabilidade de morte de cada individuo, por ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step)
  ##Vetor com numero de mortos por ciclo
  n.dead <- c()
  n.dead[1] <- 0
  ##CONDICOES INICIAIS##
  ##Deduzidas de acordo com o modelo de Hubbell:
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  ##A probabilidade inicial de morte, tomada de uma geométrica,
  ##dado que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  dead.mat[,1] <- 1/J
  p.death <- dead.mat[,1]
  ##O esforco reprodutivo instantaneo inicial é um dos definidos pelo modelo
  prop.mat[,1] <- xi
  n.propag <- prop.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step))
  {
    n.mortes <- 0
    for(j in 1:step) 
	{
      ## Sorteio dos que morrerao
      morte=rbinom(J, 1, prob=p.death)
      ##Total de mortos, que e armazenado em n.dead
      D=sum(morte)
      n.mortes <- n.mortes+D
      ## Se ha mortos comeca aqui a rotina de substituicao dos valores
      if(D>0) 
       {
           ##vetor de propagulos: cada propagulo tem o codigo numerico do individuo
          seed.bank <- rep(1:J,n.propag) 
          ##Indices dos individuos que morreram
          nascer= which(morte==1)
          ##Sorteio dos propagulos que irao repor os mortos
          mami=sample(seed.bank,D)
          ##Um vetor para armazenar o fenotipo do pai
          papi <- c()
          ##Um loop para sortear o pai entre os individuos da especie de cada
          ## propagulo-mae sorteado
          for(w in 1:D){
            papi[w] <- sample(n.propag[ cod.sp==cod.sp[mami[w]] ],1)
          }
          ##O valor esperado de propagulos dos filhotes eh
          ## a media do numero medio de propagulos produzidos pelo pai e pela mae
          medias.prop=(n.propag[mami]+papi)/2
          ##Codigos das especies dos mortos sao substituidos pelos nascidos
          ##dos propagulos sorteados
          cod.sp[nascer]<-cod.sp[mami]
          ## Numero medio de propagulos dos novos individuos nascidos e sorteada
          ## de uma normal discretizada e truncada entre 1 e J,
          ## com o coeficiente de variacao estabelecido
          n.propag[nascer] <- sapply(medias.prop,rnormt,cv=cv,min=1,max=X)
          ##A matriz de probabilidades de morrer eh atualizada para os novos individuos
          p.death[nascer] <- n.propag[nascer]/X
        }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
    dead.mat[,i] <- p.death
    prop.mat[,i] <- n.propag
    n.dead[i] <- n.mortes
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  colnames(dead.mat) <- tempo
  colnames(prop.mat) <- tempo
  names(n.dead) <- tempo
  ## Calculando riqueza final
  riq.final <- length(unique(ind.mat[,(1+(ciclo/step))]))
  ## Calculando o desvio padrão da prob. de morrer (que indica a variação entre espécies)
  prob.sp <- tapply(dead.mat[,(1+(ciclo/step))], list(ind.mat[,(1+(ciclo/step))]), mean)
  sd.dead <- sd(prob.sp)
  ## Calculando a correlação entre prob. de morrer e abundância por espécie
  ab.sp <- table(ind.mat[,(1+(ciclo/step))])
  cor.dead <- corr(cbind(ab.sp, prob.sp))
  ## Calculando o esforço reprodutivo final, por espécie
  prop.ind <- prop.mat[,2]
  prop.sp <- tapply(prop.mat[,(1+(ciclo/step))], list(ind.mat[,(1+(ciclo/step))]), mean)
  ## Chamando o resultado
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes=n.dead,riqueza.final=riq.final,sd.prob.morte=sd.dead,corr=cor.dead,esf.rep.inst.ind=prop.ind,esf.rep.inst.sp=prop.sp)
  t1=proc.time()[[3]]
  cat("\n\t tempo de processamento: ", round((t1-t0)/60,2),"\n") 
  ## incluindo atributos no arquivo resulta
  #attributes(resulta)$start <- list(espécies=S,indivíduos/sp=j,esforço reprodutivo instantâneo=xi,coeficente de variação=cv,ciclos=ciclo,passos=step)
  
return(resulta)
}
################################################
### funcăo modelo totalmente neutro - Hubell
#Esta é a funçăo simula.neutra.step editada para retirara o tradeoff reprodutivo e a variaçăo neste 
#caracter. O resulatdo é o modelo neutro do Hubbell, com a diferença que o número esperado de mortes
# a cada ciclo é um, mas há variaçăo estocástica nisto. no modelo original, há sempre uma morte por ciclo
################################################
simula.neutra.hub=function(S= 100, j=10, ciclo=1e4, step=100){ 
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##Vetor com numero de mortos por ciclo
  n.dead <- c()
  n.dead[1] <- 0
  ##CONDICOES INICIAIS##
  ##Deduzidas de acordo com o modelo de Hubbell:
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  ##A probabilidade inicial de morte, tomada de uma geométrica,
  ##dado que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  p.death <- 1/J
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    n.mortes <- 0
    for(j in 1:step){
      ## Sorteio dos que morrerao
      morte=rbinom(J, 1, prob=p.death)
      ##Total de mortos, que e armazenado em n.dead
      D=sum(morte)
      n.mortes <- n.mortes+D
      ## Se ha mortos comeca aqui a rotina de substituicao dos valores
      if(D>0) 
        {
          ##Indices dos individuos que morreram
          nascer= which(morte==1)
          ##Quem substitui os mortos
          novos <- sample(1:J,D,replace=T)
          ##Substituindo
          cod.sp[nascer]<-cod.sp[novos]
        }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
    n.dead[i] <- n.mortes
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  names(n.dead) <- tempo
  resulta=list(tempo=tempo,sp.list=ind.mat,n.mortes=n.dead)
  return(resulta)
}

############################################################
#########################################################
# Normal Discretizada Truncada

#Sorteia n números inteiros de uma aproximaçăo discreta da normal truncada em seus limites inferiores 
#e superiores dados sua média e coeficiente de variaçăo. 

##Funcao de sorteio de uma normal discretizada e truncada
rnormt <- function(mean,n=1,cv,min,max)
{
  vals <- (min-0.5):(max+0.5)
  p <- pnorm(vals,mean=mean,sd=mean*cv)
  p2 <- diff(p)
  sample(min:max,n,prob=p2,replace=T)
}

##########################################
###### Riqueza de espécies#######################
#Funçăo simples que retorna o número de elementos únicos de um vetor alfanumérico.
#Aplicado ao vetor de coódigo de espécie de cada indivíduo retorna a riqueza.

rich <- function(x)length(unique(x))

##################################$$
#Rannk-abundance Plot
############################
#Dado um vetor de abundăncia de espécies, faz um gráfico simples de Whitaker (log abundâncias x ranque).
###########
w.plot <- function(x,...)
{
  plot(1:length(x),sort(x,decreasing=T),
       xlab="Rank order",ylab="Abundance",log="y",...)
}

################################
#Box-plot e rank-abundance plot
################################
#Funçăo feita para explorar graficamente a distribuiçăo da fertilidade. Plota box-plots da fertlidade 
#dos indivíduso de cada espécie, sobre um gr[áfico de rank-abundância da fertilidade de cada espécie, 
#num dado instante gravado da simulaçăo (argumento tempo).

box.w <- function(tempo,cod.sp,n.propag)
 {
  sp <- cod.sp[,tempo]
  prop <- n.propag[,tempo]
  ab <- tapply(sp,sp,length)
  par(mfrow=c(2,1))
  boxplot(prop~factor(sp,levels=names(ab)[order(ab,decreasing=T)]))
  w.plot(ab)
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
  medias
}



