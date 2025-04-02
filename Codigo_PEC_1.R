if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("metabolomicsWorkbenchR") 
#Necesitas el paquete BiocManager para acceder al repositorio
# donde se encuentra en paquete metabolomicsWorkbenchR

library("metabolomicsWorkbenchR") #Cargamos el paquete
df = do_query( # do_query necesita estos 4 argumentos siguientes:
  context = 'study',
  input_item = 'study_title',
  input_value = 'brain', #Obviamente, la busqueda se debe realizar en inglés
  output_item = 'summary')
View(df) #Para observar la lista de resultados

SE = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST003557',
  output_item = 'SummarizedExperiment' #Aquí indicamos en que formato 
  # queremos que se extraigan los datos
)
SE

if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")

se1 <- SE[[1]] #Separamos los elementos
se2 <- SE[[2]]

library(SummarizedExperiment)
library(Biobase)

exprs_matrix1 <- as.matrix(assays(se1)[[1]]) # Extraemos la matriz de datos con
exprs_matrix2 <- as.matrix(assays(se2)[[1]]) # assays()


expset1 <- ExpressionSet( #Creamos el objeto expressionSet rellenando cada
  assayData = exprs_matrix1, # atributo
  phenoData = AnnotatedDataFrame(as.data.frame(colData(se1))),
  featureData = AnnotatedDataFrame(as.data.frame(rowData(se1)))
)

expset2 <- ExpressionSet(
  assayData = exprs_matrix2,
  phenoData = AnnotatedDataFrame(as.data.frame(colData(se2))),
  featureData = AnnotatedDataFrame(as.data.frame(rowData(se2)))
)


expset_comb <- BiocGenerics::combine(expset1, expset2)

expset_comb <- BiocGenerics::combine(expset1, expset2)
expr_matrix<- rbind(exprs_matrix1, exprs_matrix2)
#Combinamos las matrices tambien

validObject(expset_comb) # Comprobación de que el objeto expressionSet es correcto

dim(expset_comb) # Se exploran las dimensiones

summary(pData(expset_comb)) #Se exploran variables

summary(fData(expset_comb))

metadata(se1)$description #Se muestra el título del estudio

expset_comb@featureData@data #Información de los metabolitos

metabolite_means <- apply(expr_matrix, 1, mean) # Calculamos la media
metabolite_sd <- apply(expr_matrix, 1, sd) # Calculamos la desviación

# Graficar media vs desviación estándar
plot(metabolite_means, metabolite_sd, xlab="Media", ylab="Desviación estándar", main="Dispersión de la media vs desviación estándar", pch=16, col=rgb(0,0,1,0.5))

hist(expr_matrix, breaks = 50, main = "Distribución de valores de expresión", xlab = "Expresión", col = "lightblue", border = "black")

boxplot(t(expr_matrix), main = "Boxplot de expresión por metabolito", col = rainbow(ncol(expr_matrix)), las = 2)

library(ggplot2)
library(ggfortify)

pca_result <- prcomp(t(expr_matrix), scale. = TRUE) # Aplicamos PCA

autoplot(pca_result, data = pData(expset_comb),
         colour = "Treatment",
         shape = "Sample_source"
) + ggtitle("Análisis de Componentes Principales") + theme_minimal()
#Graficamos diferenciando tratamiento y tejido
summary(pca_result)

