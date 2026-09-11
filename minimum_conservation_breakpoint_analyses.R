###############################################################
### UMBRAL ÓPTIMO DE CONSERVACIÓN - TRES ESCENARIOS ##########
###############################################################
rm(list = ls())

### 1. CARGAR DATOS ###########################################
DATOS <- read.csv("C:/Users/ASUS/Downloads/curve_performance_zonation_scenarios.csv",header = TRUE)

head(DATOS)
dim(DATOS)
names(DATOS)
str(DATOS)


### 2. PARÁMETROS GENERALES ###################################
PA_ACTUAL <- 13.7# Cobertura ACTUAL de áreas protegidas del NSDF
TARGET_30 <- 30# Meta internacional 30x30
VENTANA <- 10#Ventana para estimar pendientes locales


### 3. NOMBRES DE LOS TRES ESCENARIOS #########################
NOMBRES_ESCENARIOS <- c("Scenario 1", "Scenario 2", "Scenario 3")


### 4. COMPROBAR QUE EXISTEN LAS SEIS COLUMNAS ################
if(ncol(DATOS) < 6) {
  stop("El archivo necesita al menos seis columnas: ",
       "dos columnas para cada escenario.")
}


### 5. FUNCIÓN PARA ANALIZAR CADA ESCENARIO ###################
analizar_escenario <- function(X,
                               Y,
                               nombre,
                               pa_actual = 13.7,
                               target30 = 30,
                               ventana = 10,
                               hacer_plot = TRUE) {
  ### A. LIMPIAR DATOS #######
      temp <- data.frame(area = as.numeric(X),spp  = as.numeric(Y))
      temp <- temp[is.finite(temp$area) & is.finite(temp$spp),]  # Eliminar NA e infinitos
      temp <- temp[order(temp$area),]   # Ordenar por área protegida
      temp <- aggregate(spp ~ area, data = temp, FUN = mean)# Promediar valores repetidos de X, si existen
      temp <- temp[order(temp$area),]
      rownames(temp) <- NULL
  
  ### B. VERIFICAR RANGO #####
  if(min(temp$area) > pa_actual) {
    stop(nombre,": los datos no alcanzan el valor actual de ",
         pa_actual, "%.") }
  
  if(max(temp$area) < target30) {
    stop(nombre, ": los datos no alcanzan el target de ",
         target30,    "%.")}
  
  ### C. REPRESENTACIÓN EN 13.7% Y 30% #####
  SPP_ACTUAL <- approx(x = temp$area, y = temp$spp,  xout = pa_actual)$y
  SPP_30 <- approx(x = temp$area, y = temp$spp, xout = target30)$y
  
  ### D. CREAR CURVA ENTRE 13.7% Y 30% ######
  # Seleccionar puntos que se encuentran dentro del intervalo de interés para la expansión de las PAs
  tramo <- temp[temp$area >= pa_actual &  temp$area <= target30,]
  
  inicio <- data.frame(area = pa_actual,spp = SPP_ACTUAL)
  final <- data.frame(area = target30,spp = SPP_30)
  tramo <- rbind(inicio,tramo,final)
  
  # Eliminar posibles duplicados
  tramo <- aggregate(spp ~ area, data = tramo, FUN = mean)
  tramo <- tramo[order(tramo$area), ]
  rownames(tramo) <- NULL
  
  ### E. DETECCIÓN AUTOMÁTICA DEL KNEE / ELBOW POINT ##########
  tramo$x_norm <- (tramo$area - min(tramo$area)) / (max(tramo$area) - min(tramo$area))#Normalizar X entre 0 y 1
  tramo$y_norm <- (tramo$spp - min(tramo$spp)) / (max(tramo$spp) - min(tramo$spp))# Normalizar Y entre 0 y 1
  
  ### Distancia respecto a la línea entre inicio y final #######
  # Después de normalizar, la línea que conecta los extremos es: y = x. En una curva cóncava, el punto con mayor diferencia y_norm - x_norm corresponde al knee point.
  tramo$distance <- tramo$y_norm - tramo$x_norm
  
  ### Excluir los extremos #####
  candidatos <- 2:(nrow(tramo) - 1)
  POS_UMBRAL <- candidatos[which.max(tramo$distance[candidatos])]
  
  ### UMBRAL CALCULADO AUTOMÁTICAMENTE #########################
  UMBRAL <- tramo$area[POS_UMBRAL]
  SPP_UMBRAL <- tramo$spp[POS_UMBRAL]
  FUERZA_UMBRAL <- tramo$distance[POS_UMBRAL]
  
  ### F. PENDIENTES MÓVILES ####################################
  PENDIENTE <- numeric()
  INTERVAL <- numeric()

  if(nrow(temp) >= ventana) {
    for(k in ventana:nrow(temp)) {
      idx <- (k - ventana + 1):k
      modelo <- lm(
        spp ~ area,
        data = temp[idx, ])
      
      PENDIENTE <- c(PENDIENTE, unname(coef(modelo)[2]))
      INTERVAL <- c(INTERVAL,mean(temp$area[idx]))
    }
    
  }
  
  pendientes <- data.frame(area = INTERVAL, slope = PENDIENTE)
  
  
  ### G. GRÁFICA PRINCIPAL ########
  if(hacer_plot) {
    plot(temp$area,
         temp$spp,
         type = "l",
         col = "red",
         lwd = 2,
         xlim = c(0, 100),
         ylim = c(0, 100),
         las = 1,
         xlab = "Protected area of the NSDF (%)",
         ylab = "Protected distribution area of the species (%)",
         main = nombre,
         bty = "n")
    
    ### 13.7% - SITUACIÓN ACTUAL ##############################
    abline(a = 0, b = 1, lwd = 0.8, col = "black")### Línea 1:1 ######
    segments(x0 = pa_actual,
             y0 = 0,
             x1 = pa_actual,
             y1 = SPP_ACTUAL,
             lty = 3,
             lwd = 1.2,
             col = "grey40")    

    segments(x0 = 0,
             y0 = SPP_ACTUAL,
             x1 = pa_actual,
             y1 = SPP_ACTUAL,
             lty = 3,
             lwd = 1.2,
             col = "grey40")
    
    ### UMBRAL CALCULADO ########
    segments(x0 = UMBRAL,
             y0 = 0,
             x1 = UMBRAL,
             y1 = SPP_UMBRAL,
             lty = 2,
             lwd = 1.5,
             col = "blue")
    
    segments(x0 = 0,
             y0 = SPP_UMBRAL,
             x1 = UMBRAL,
             y1 = SPP_UMBRAL,
             lty = 2,
             lwd = 1.5,
             col = "blue")
    
    
    ### 30% #######
    segments(x0 = target30,
             y0 = 0,
             x1 = target30,
             y1 = SPP_30,
             lty = 3,
             lwd = 1.2,
             col = "blue")
    
    
    segments(x0 = 0,
             y0 = SPP_30,
             x1 = target30,
             y1 = SPP_30,
             lty = 3,
             lwd = 1.2,
             col = "blue")
    
    ### PUNTOS ########
    points(c(pa_actual, UMBRAL, target30),
           c(SPP_ACTUAL,SPP_UMBRAL, SPP_30),
           pch = 16,
           cex = 0.75)
    
    ### ETIQUETAS DEL EJE X #######
    mtext(round(pa_actual, 1), side = 1, at = pa_actual, line = -1, cex = 0.75)
    mtext(round(UMBRAL, 1), side = 1, at = UMBRAL,line = -1,  cex = 0.8,  font = 2)
    mtext(round(target30, 1), side = 1, at = target30,line = -1,cex = 0.75)
    
    ### ETIQUETAS DEL EJE Y #######
    mtext(round(SPP_ACTUAL, 1), side = 2, at = SPP_ACTUAL, line = -1.5, cex = 0.75, las = 1)
    mtext(round(SPP_UMBRAL, 1), side = 2, at = SPP_UMBRAL, line = -1.5, cex = 0.8, las = 1, font = 2)
    mtext(round(SPP_30, 1), side = 2, at = SPP_30, line = -1.5, cex = 0.75, las = 1)
  }
  
  ### H. RESULTADOS ##########
  resumen <- data.frame(scenario = nombre,
                        current_PA = pa_actual,
                        current_species_representation = SPP_ACTUAL,
                        estimated_threshold = UMBRAL,
                        threshold_species_representation = SPP_UMBRAL,
                        knee_strength = FUERZA_UMBRAL,
                        target_30 = target30,
                        representation_at_30 = SPP_30)
  
  
  return(
    list(datos = temp,
         tramo = tramo,
         pendientes = pendientes,
         resumen = resumen
    ))
}


### 6. CORRER LOS TRES ESCENARIOS ######
par(mfrow = c(1, 3),
    mar = c(5, 4.5, 3, 1))


### ESCENARIO 1 - COLUMNAS 1 Y 2 #######
RESULTADO_1 <- analizar_escenario(X = DATOS[[1]],
                                  Y = DATOS[[2]],
                                  nombre = NOMBRES_ESCENARIOS[1],
                                  pa_actual = PA_ACTUAL,
                                  target30 = TARGET_30,
                                  ventana = VENTANA)

### ESCENARIO 2 - COLUMNAS 3 Y 4 ##########
RESULTADO_2 <- analizar_escenario(X = DATOS[[3]],
                                  Y = DATOS[[4]],
                                  nombre = NOMBRES_ESCENARIOS[2],
                                  pa_actual = PA_ACTUAL,
                                  target30 = TARGET_30,
                                  ventana = VENTANA)

## ESCENARIO 3 - COLUMNAS 5 Y 6 ########
RESULTADO_3 <- analizar_escenario(X = DATOS[[5]],
                                  Y = DATOS[[6]],
                                  nombre = NOMBRES_ESCENARIOS[3],
                                  pa_actual = PA_ACTUAL,
                                  target30 = TARGET_30,
                                  ventana = VENTANA)

par(mfrow = c(1, 1))


### 7. TABLA FINAL CON LOS TRES UMBRALES #########
RESUMEN <- rbind(RESULTADO_1$resumen, RESULTADO_2$resumen, RESULTADO_3$resumen)
RESUMEN$current_PA <- round(RESUMEN$current_PA, 2)
RESUMEN$current_species_representation <- round(RESUMEN$current_species_representation, 2)
RESUMEN$estimated_threshold <- round(RESUMEN$estimated_threshold, 2)
RESUMEN$threshold_species_representation <- round(RESUMEN$threshold_species_representation, 2)
RESUMEN$knee_strength <- round(RESUMEN$knee_strength, 4)
RESUMEN$representation_at_30 <- round(RESUMEN$representation_at_30, 2)

print(RESUMEN)

### 8. GUARDAR TABLA Y FIGURA ####
write.csv(RESUMEN, "C:/Users/ASUS/Downloads/thresholds_three_scenarios.csv", row.names = FALSE)
png("C:/Users/ASUS/Downloads/thresholds_three_scenarios.png", width = 3300,  height = 1200,  res = 300)

par(mfrow = c(1, 3),  mar = c(5, 4.5, 3, 1))

analizar_escenario(DATOS[[1]], DATOS[[2]],  NOMBRES_ESCENARIOS[1],  PA_ACTUAL,  TARGET_30,  VENTANA)
analizar_escenario(DATOS[[3]], DATOS[[4]], NOMBRES_ESCENARIOS[2], PA_ACTUAL,  TARGET_30, VENTANA)
analizar_escenario(DATOS[[5]], DATOS[[6]], NOMBRES_ESCENARIOS[3], PA_ACTUAL,  TARGET_30,  VENTANA)
dev.off()


### 9. DEFINIR UMBRAL FINAL ####
UMBRALES <- RESUMEN$estimated_threshold
UMBRAL_MEDIANA <- median(UMBRALES, na.rm = TRUE)# Mediana de los tres knee points
UMBRAL_FINAL <- round(UMBRAL_MEDIANA * 2) / 2# Redondear al 0.5% más cercano

cat("\nKnee points:", paste(UMBRALES, collapse = ", "), "%\n")
cat("Median:", round(UMBRAL_MEDIANA, 2), "%\n")
cat("Final operational threshold:", UMBRAL_FINAL, "%\n")


### FIN ###