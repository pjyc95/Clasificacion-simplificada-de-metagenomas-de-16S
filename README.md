install.packages("readxl")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("purrr")
install.packages("readr")
install.packages("openxlsx")

# Cargar la librería necesaria
library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
##Garfico de dona
library(ggplot2)
library(purrr)
library(scales)
library(ggrepel)
library(openxlsx)

# Ruta del archivo de entrada para la lectura de los datos de Kraken 2
ruta_entrada <- "C:/Users/pjyc9/OneDrive/Escritorio/Master Bioinformatica/TFM/datos tesis/secuencias_filtradas_16S/PuertoVaras/QuebradaParque1_F_filt_report.reportk2"

# Leer los datos en un data.frame
datos <- read_delim(ruta_entrada, delim = "\t", col_names = FALSE)

# Definir nombres de columnas
nombres_columnas <- c("A", "B", "C", "D", "E", "Genero")
  # A es el Porcentaje de asignacion
  # B es Numero de lecturas que hay para el clado taxonomico
  # D es el Nivel taxonomica
  # Genero es la Asignación taxonómica

# Asignar nombres a las columnas
colnames(datos) <- nombres_columnas

# Filtrar los datos de Kraken en base a la identificacion de genero
head(datos)
filtrado1 <- datos  %>%
  filter (D == "G")

# Eliminar los espacios en blanco de la columna de genero
filtrado <- filtrado1 %>%
  mutate(Genero = trimws(Genero)) %>%
  arrange(desc(B))

# Ruta del archivo Excel de la Tabla maestra con la identificacion de Genero, Clasificacion y Clasificación Ligth
ruta_tabla_maestra <- "C:/Users/pjyc9/OneDrive/Escritorio/Master Bioinformatica/TFM/datos tesis/Tabla_Maestra.xlsx"

# Leer los datos del archivo Excel y crea un data frame
tabla_maestra <- read_excel(ruta_tabla_maestra)

# Visualizar las primeras filas de ambas tablas para entender su estructura
head(filtrado)
head(tabla_maestra)

# Fusionar las dos tablas en función del género de la especie
tabla_completa1 <- merge(filtrado, tabla_maestra, by = "Genero", all.x = TRUE)
tabla_completa <- tabla_completa1 %>%
  arrange(desc(B))
tabla_completa <- tabla_completa %>%
  slice(1:100)
# Guardar el resultado si es necesario
write.xlsx(tabla_completa, "C:/Users/pjyc9/OneDrive/Escritorio/Master Bioinformatica/TFM/datos tesis/tabla_completa_QuebradaParque1_F_filt_.xlsx", rowNames = FALSE)

# Visualizar las primeras filas del resultado
head(tabla_completa)

#Filtrado para las muestras que tienen clasificación en la tabla maestra y las que no.
# Filtrar los datos de Kraken en base a la identificacion de genero
tabla_completa2<- tabla_completa
head(tabla_completa2)
# Llenar los datos de ClasificacionLigth que no tienen valor con "NA"
tabla_completa2$ClasificacionLigth[is.na(tabla_completa2$ClasificacionLigth)] <- "NA"

# Sumar los valores de C agrupados por ClasificacionLigth
ident_tabla_maestra <- aggregate(B ~ ClasificacionLigth, data = tabla_completa2, FUN = function(x) sum(x, na.rm = TRUE))
ident_tabla_maestra$Porcentaje <- round(ident_tabla_maestra$B / sum(ident_tabla_maestra$B, na.rm = TRUE) * 100,2)
ident_tabla_maestra

## grafico de dona ident_tabla_maestra
ggplot(ident_tabla_maestra, aes(x = 2, y = Porcentaje, fill = ClasificacionLigth)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = percent(Porcentaje/100)), 
                  position = position_stack(vjust = 0.5), 
                  color = "black", size = 5.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("yellow", "red","blue", "green3")) +
  theme_void() +
  labs(title = "Gráfico de dona de la identificacion del analisis de secuencias inicial tabla maestra") +
  xlim(0.5, 2.5) 

#para clasificar entre NA y no NA en Clasificación Ligth 
filtrado_na <- tabla_completa %>%
  filter(is.na(ClasificacionLigth))
head(filtrado_na)
#para clasificar entre diferente a NA en Clasificación Ligth
filtrado_other1 <- tabla_completa %>%
  filter(!is.na(ClasificacionLigth))
head(filtrado_other1)
# Ordenar el data frame filtrado_other de mayor a menor según la columna B
filtrado_other1 <- filtrado_other1 %>%
  arrange(desc(B))
head(filtrado_other1)

# Seleccionar los primeros 100 datos y almacenarlos en filtrado_other1
filtrado_other <- filtrado_other1 %>%
  slice(1:100)

# Calcular el porcentaje por genero de manera individual
filtrado_other$porcentaje <- 0
filtrado_other$porcentaje <- round(filtrado_other$B / sum(filtrado_other$B) * 100, 2)

# Calcular el porcentaje final de cada clasificación
porcentaje_final_clasificacion <- aggregate(porcentaje ~ Clasificacion, data = filtrado_other, FUN = sum)

# Calcula el porcentaje total para normalizar los porcentajes finales
porcentaje_final_clasificacion$total <- sum(porcentaje_final_clasificacion$porcentaje)

porcentaje_final_clasificacion <- porcentaje_final_clasificacion %>% 
  mutate(csum = rev(cumsum(rev(porcentaje))), 
         pos = porcentaje/2 + lead(csum, 1),
         pos = if_else(is.na(pos), porcentaje/2, pos))

## grafico de pastel para la clasificacion basado en sus características biológicas y su relación con el medio.
ggplot(porcentaje_final_clasificacion, aes(x = 2, y = porcentaje, fill = Clasificacion)) +
  geom_bar(stat = "identity", color = "black") +
  geom_label_repel(data = porcentaje_final_clasificacion,
                   aes(y = pos, label = paste0(porcentaje, "%")),
                   size = 4, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("yellow", "blue", "#98F5FF","red")) +
  theme_void() +
  labs(title = "Gráfico de pastel para la clasificacion basado en sus características biológicas y su relación con el medio") 
  xlim(0.5, 2.5)

# Calcular el porcentaje final de cada clasificación
porcentaje_final_clasificacion_ligth <- aggregate(porcentaje ~ ClasificacionLigth, data = filtrado_other, FUN = sum)

# Calcula el porcentaje total para normalizar los porcentajes finales
porcentaje_final_clasificacion_ligth$total <- sum(porcentaje_final_clasificacion_ligth$porcentaje)

porcentaje_final_clasificacion_ligth <- porcentaje_final_clasificacion_ligth %>% 
  mutate(csum = rev(cumsum(rev(porcentaje))), 
         pos = porcentaje/2 + lead(csum, 1),
         pos = if_else(is.na(pos), porcentaje/2, pos))

# Gráfico de dona para la clasificacion basado en su clasificacion Ligth
ggplot(porcentaje_final_clasificacion_ligth, aes(x = 2, y = porcentaje, fill = ClasificacionLigth)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = percent(porcentaje/100)), 
                  position = position_stack(vjust = 0.5), 
                  color = "black", size = 4) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("yellow", "blue1")) +
  theme_void() +
  labs(title = "Gráfico de dona para la clasificacion basado en su clasificacion Ligth") +
  xlim(0.5, 2.5) 
