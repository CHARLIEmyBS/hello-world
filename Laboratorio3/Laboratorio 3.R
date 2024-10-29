library("ape")
library ("phangorn")
library("phytools")

#Convertir archivo en compatible:
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                        format = "FASTA", type = "AA")
fraxatin

#crear una matriz de distancia para poder crear árboles de distancia o parsimonia a través de ella
#función dist.aa--> método o función de objeto o clase
#preparar los datos para que R pueda leer el algoritmo.
matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist

#ARBOLES DE DISTANCIA: estos árboles son árboles de similaridad, no representan inferencia evolutiva.

#MÉTODO UPGMA (grupo de pares no ponderados con media aritmética): 
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)
#Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.

#METODO NJ (unión de vecinos):
arbolNJ <- nj(matrizdist)
plot(arbolNJ)
#Modificaciones al árbol: cex (tamaño letra), edge.width (anchura de la rama), type (forma del esquema), font(cambia el tipo de letra)
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="brown", font=3) 
#node.pos (mueve la linea), edge.lty (cambia la continuidad de la linea), label.offset (para cambiar el espacio entre la palabra y la rama)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="lightblue", font=3)


#PERSONALIZACIÓN DEL ARBOL CON PYTHTOOLS:
plotTree(arbolNJ)
#diferentes modificaciones: 
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="pink", lwd=2)
plotTree(ladderize(arbolNJ))#cambiar el orden en que los grupos son visualizados


write.tree(arbolNJ, file = "file_name.nex")#Para guardar el árbol
read.tree(file="file_name.nex") #para verlo (me sale información sobre el árbol)

#ENRAIZAR:
#Hasta ahora, todos los árbolos  no están enraízados. Para enraizarlos podemos usar la función root del paquete ape. 
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
#Hacemos un árbol enraizado con el método UPGMA:
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)
#Visualizar los dos árboles a la vez con los siguientes comandos:
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)

#ÁRBOLES PARSIMONIA: usca disminuir el número de pasos que explican un árbol evolutivo contando el número de cambios de cada uno de los caracteres.
#Se pueden evaluar varios árboles.Se eliminan todos los caracteres onstantes en todos los taxones, o aquellos que son variables pero no informativos
#Tener en cuenta: caracter es informativo cuando si tiene al menos dos estados de caracter y por lo menos dos de estos estados ocurren con una frecuencia mínima de dos.
parsimony(arbolUPGMAraiz, fraxatin) #me da 313, que son los cambios evoliutivo. Lo que hace es compararme todos los árboles y me va a dar el que menos pasos evolutivos tenga.
parsimony(arbolUPGMA, fraxatin)#aunque esté con raíz o no el número de pasos debe ser el mismo (RAÍZ: se fija uno de los parámetros a partir del cual se genera el resto del árbol)

#Arboles con mejor parsimonia: menor número de pasos para generar ese árbol.
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin) #me da 307: y para generarlo ha usado solo dos operaiones (menor número)
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)#me da 307: para generarlo ha usado 1 operación, lo que significa que es mejor que UPGMA, ya que usa un menor numero de operaciones para generarlo)

#Otra metodología llamada pratched, que lo que hace es darte la misma solución que antes solo que te enseña lo que hace para generar dicho resultado.
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)

#Me da el árbol con mejor parsiominia, y el número de árboles que presentan dicho resultado.
fraxatin_parsimonia #Me da que 4 árboles cumplen la parsimonia de 307.

#Para compararlos es necesario enraizarlos:
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6) #me compara los cuatro árboles que me han salido antes

#A partir de los cuatro árboles que tenía antes, me los va a juntar en un único árbol que SOLO contienen las especies que tienen todos en común.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)

#Para un árbol menos estricto podemos cambiar el valor del parámetro p:
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)

#BOOTSTRAP:
#generamos seudoréplicas con remplazamiento de una matriz.Con cada réplica se hace un árbol consenso:
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
#La rutina anterior genera entonces 10 árboles pseudoréplicas.
plot(arbolesbootstrap, cex = .6)
#Generamos un consenso; en este caso con un consenso al 60%:
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

#MODELOS PROBABILÍSTICOS:
#Árboles de máxima verosimilitud:se calcula la probabilidad de obtener un dato según un modelo de un árbol, de acuerdo a un alineamiento de secuencias usando un modelo de sustitución de aminoácidos.
#Genera un árbol filogenético aleatorio. En este caso, se crea un árbol con 11 taxones aleatorias. establece las etiquetas para las hojas (especies) del árbol. Usa los nombres presentes en el objeto fraxatin para etiquetar las 11 especies.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)

#Lo enraizamos por la secuencias de ornitorrico. Además, los escalarizamos hacia la derecha, para que haya una caída y le agregamos escala.
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()

#A partir del arbol de arriba se puede iniciar la búsqueda del mejor árbol por máxima verosimilitud. Lo primero que se hace es carcular la verosimilitud del árbol dadas las secuencias.
ajustado <- pml(arbolazarR, fraxatin)
ajustado 
#La información que tiene el objeto ajustado nos reporta la verosimilitud del árbol al azar que habíamos creado, que es -4348.064. 

#
ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

#Para ver el árbol oculto usamos $tree. También lo enraizamos.
ajustadoconDay$tree

ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()

#El árbol anterior fue generado usando la matriz de sustitución de Dayhoff. Pero se pueden usar diferentes modelos.
justadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")

#Podemos comparar los modelos calculando el Criterio de información de Akaike AIC:
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)

#La primera columna corresponde a los grados libertad. Según el criterio anterior, el mejor modelo que se ajusta con los datos (con el AIC más bajo) es JTT modelo de Jones-Taylor-Thornton para evaluar la distancia entre secuencias de proteínas y optimiza la verosimilitud.
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")
mejorarbol
mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()






