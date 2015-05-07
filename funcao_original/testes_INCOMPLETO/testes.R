teste1<-simula.neutra.step(S=100, j=10,X=10000, cv=0.05, ciclo=1000, step=10)
teste2<-simula.neutra.step(S=100,j=10,X=10000,cv=0.05,ciclo=1e6,step=100)
teste3 <- simula.neutra.step(100,10,10000,0.5,1e6,1000)
teste4 <- simula.neutra.step(50,10,10000,0.5,1e6,1000)
teste5 <- simula.neutra.step(50,10,10000,2,1e6,1000) # rodado com a função em que o desvio padrão é o próprio cv
teste6 <- simula.neutra.step(40,10,10000,3,1e6,1000) # rodado com a função em que o desvio padrão é o próprio cv
teste7 <- simula.neutra.step(50,10,10000,3,1e6,1000) # rodado com a função em que o desvio padrão é o próprio cv
