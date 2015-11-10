##########################################################
#################### HIPERCUBO LATINO ####################
##########################################################
#################### PRIMEIROS TESTES ####################
##########################################################

require(pse)

################# LISTANDO OS PARÂMETROS #################
factors <- c("S", "j", "X", "dp")
############### DEFININDO AS DISTRIBUIÇÕES ###############
############ DE PROBABILIDADE DOS PARÂMETROS #############
q <- c("qunif", "qunif", "qunif", "qlnorm")
################ DEFININDO OS PARÂMETROS #################
#################### DAS DISTRIBUIÇÕES ###################
q.arg <- list( list(min=2,max=100), list(min=10,max=1000),
               list(min=1, max=1000), list(meanlog=0.5,sdlog=1))
###################### MINHA FUNÇÃO: #####################
################### SIMULA.NEUTRA.STEP ###################
simula.neutra.step=function(S= 100, j=10, X=1000, dp=0.1, ciclo=1e6, step=100)
{
  t0=proc.time()[[3]]
  ## Tamanho da comunidade
  J <- S*j
  ##Verifica se X e multiplo de J
  if(abs(X/J - round(X/J)) > .Machine$double.eps^0.5) ## aqui verifica se o X/J e inteiro ou quase!
  {
    stop("\n\tO potencial reprodutivo (X) precisa ser multiplo do tamanho da comunidade (J). Tente novamente!\n\n")
  } 
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
  ##A probabilidade inicial de morte, tomada de uma geometrica,
  ##dado que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  dead.mat[,1] <- 1/J
  p.death <- dead.mat[,1]
  ##O esforco reprodutivo inicial e deduzido da esperanca de tempo de vida da geometrica E[t]=J
  ## o que portanto resulta em X/J propagulos produzidos a cada ciclo, por todos
  prop.mat[,1] <- X/J
  n.propag <- prop.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step))
  {
    n.mortes <- 0
    for(a in 1:step) 
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
        for(w in 1:D)
        {
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
        ## com o desvio padrao estabelecido
        n.propag[nascer] <- sapply(medias.prop,rnormt,dp=dp,min=1,max=X)
        ##A matriz de probabilidades de morrer eh atualizada para os novos individuos
        p.death[nascer] <- n.propag[nascer]/X
      }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
    dead.mat[,i] <- p.death
    prop.mat[,i] <- n.propag
    n.dead[i] <- n.mortes
    ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  colnames(dead.mat) <- tempo
  colnames(prop.mat) <- tempo
  names(n.dead) <- tempo
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes=n.dead)
  t1=proc.time()[[3]]
  cat("\n\t tempo de processamento: ", round((t1-t0)/60,2),"\n") 
  ## incluindo atributos no arquivo resulta
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=X, sd=dp, ciclos=ciclo, passos=step)
  return(resulta)
}
################### FUNÇÃO ACESSÓRIA À ###################
################### SIMULA.NEUTRA.STEP ###################
rnormt <- function(mean,n=1,dp,min,max)
{
  vals <- (min-0.5):(max+0.5)
  p <- pnorm(vals,mean=mean,sd=dp)
  p2 <- diff(p)
  sample(min:max,n,prob=p2,replace=T)
}
############ FUNÇÕES DE ANÁLISE EXPLORATÓRIA #############


############### SIMULA.NEUTRA.STEP.EXTERNA ############### ### modificar
simula.neutra.step.externa=function(S, j, X, dp){
  ciclo <- 100
  step <- 100
  if(abs(X/(S*j)) - round(X/(S*j)) > .Machine$double.eps^0.5)
  {stop}
  x <- simula.neutra.step(S, j, X, dp, ciclo, step)
  return(unlist(x[3]))
}
####################### "WRAPPER" ########################
modelRun <- function (dados) {
  mapply(simula.neutra.step.externa, dados[,1], dados[,2], dados[,3], dados[,4]) }
################## RODANDO O HIPERCUBO ###################
res.name <- c("propagulos/ciclo")
hipersuperincrivelcubolatino <- LHS(modelRun, factors, N=10, q, q.arg, res.name, nboot=50)
######### ACESSANDO OS VALORES USADOS COMO INPUT #########
get.data(hipersuperincrivelcubolatino)
################ ACESSANDO OS RESULTADOS #################
get.results(hipersuperincrivelcubolatino)
################# ANÁLISE DE INCERTEZA ###################
plotecdf(hipersuperincrivelcubolatino, stack=F)
############### ANÁLISE DE SENSIBILIDADE #################
plotscatter(hipersuperincrivelcubolatino, add.lm=FALSE)
plotprcc(hipersuperincrivelcubolatino)
############### SBMA: A AMOSTRA ESTÁ OK? #################
#targetLHS <- target.sbma (target=0.7, modelRun, factors, q, q.arg, res.names, FUN=min)

###
