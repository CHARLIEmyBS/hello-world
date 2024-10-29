library(tidyverse)
data(starwars)

#FILTRAR Y SELECCIONAR DATOS:
# Seleccionar todas las columnas menos el nombre
starwars %>% select(-name)
#Seleccionar sólo las columnas que tienen subraya (_)
starwars %>% select(contains("_"))
#Seleccionar sólo las columnas que empiezan con "s"
starwars %>% select(starts_with("s"))
#Crear un data frame con los nombres y planeta de origen (homeworld)
homeworld <- starwars %>% select(name, homeworld)
#Filtrar datos 
#Filtrar por especies: sólo humanos
human <- starwars %>% filter(species == "Human")
#Filtrar por especies: sólo humanos del planeta Tatooine
starwars %>% filter(species == "Human", homeworld == "Tatooine")
#Crear un nuevo datframe con todas las especies menos los Droides
starwars_nodroids <- starwars %>% filter(species != "Droid")

#¿Cuántos registros cumplen las condiciones finales?
#--->De los registros de homeworld aparecen 87. por otro lado, en los human aparecen 35 registros finales, mientras que en los no droids son 77 registros finales.

#SELECCIONAR Y AGRUPAR DATOS:
#Usamos group_by (agrupar por variables) y tally (cuenta las observaciones de cada grupo)
starwars %>% group_by(species) %>% tally()
#Añadiendo otra variable
starwars %>% group_by(species, gender) %>% tally()
#Si lo quieres guardar en el environment recuerda asignarle un nombre
table_gender <- starwars %>% group_by(species, gender) %>% tally()

#CALCULAR ALGUNOS ESTADÍSTICOS:
#na.rm=T quiere decir que elima los NA (valores No Asignados o sin datos)
starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))

#¿Cómo calcularías la desviación estándar (sd) de esos parámetros? Recuerda consultar con ? si no sabes como usar una función o comando. Por ejemplo: ?summarise() ,?sd().
#Paso1:
desv.estándar2<-desv.estándar1%>%select(contains("_")) #-->de la columna principal nos quedamos solo con aquellos valores numéricos
#Paso2: 
a<-desv.estándar2%>%pull(mean_height)#-->creamos un vector a partir de los datos de la columna mean_height
#Paso 3: 
sd (a)-->#calculamos la desviación estándar
#Resultado: sd=42.24324 (mean_height:a)
#Resultado:sd=229.3549 (mean_mass:b)
sd(b,na.rm=TRUE)

#CREAR GRÁFICOS Y MODIFICAR ELEMENTOS:
#Hacer un gráfico de la altura vs. la masa de los personajes
ggplot(starwars, aes(height, mass)) + geom_point()
#Puedes modificar el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "blue")
#Modificando el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)
#Modificando el color y el fondo 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()

#Actividad:
#StarwarssinJ<-starwars[-16,]-->hacemos la tabla sin Jaba 
#ggplot(StarwarssinJ, aes(height, mass)) + geom_point()

#EJERCICIO 1:
toy<-read_csv


#EJERCICIO 2:
getwd()
toy <- read_csv("toy.csv")
#Inspecciona el dataset, haz un resumen de la media (mean) de las variables (Peso, Altura,IMC, IAS, CCintura). Agrupando por sexo.
#Primer apartado (medias)
Mujeres<-toy%>%filter(Sex=="Women")
Hombres<-toy%>%filter(Sex=="Men")toy
toy %>% group_by(Sex) %>% summarise(mean_peso = mean(Weight_Kg, na.rm = TRUE),mean_altura = mean(Height_cm, na.rm = TRUE),mean_imc = mean(IMC, na.rm = TRUE),mean_ias = mean(IAS, na.rm = TRUE),mean_ccintura = mean(Ccintura, na.rm = TRUE))
#Segundo apartado (tabla)
Mujeres.nonormal<-Mujeres%>%filter(IMC_clas!="Normal")
Mujeres.obesas<-Mujeres.nonormal%>%filter(IMC_clas!="Obesity")
#Tercer apartado (gráfico)
ggplot(toy, aes(IMC, Weight_Kg)) + geom_point(colour = "red") + theme_light()
Obesity<-toy%>%filter(IMC_clas!="Normal")
ggplot(Obesity, aes(IMC, Weight_Kg)) + geom_point(colour = "blue") + theme_light()