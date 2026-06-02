# Rose-transcriptomics-analysis
R scripts for filtering gene counts and running differential expression analysis using DESeq2
# ==============================================================================
# 1. SIMULACIÓN OPTIMIZADA DE CONTEOS DE RNA-SEQ (POSTCOSECHA DE ROSAS)
# ==============================================================================
```r
set.seed(123)
```
# Ampliamos la lista de genes para incluir controles internos y genes con ruido
```r
genes_rosa <- c(
  "Rhy_001254_Aquaporina_Hidratacion",
  "Rhy_003421_Receptor_Etileno_RhETR1",
  "Rhy_007845_Degradacion_Pared_Poligalacturonasa",
  "Rhy_012548_Auxinas_Apertura_RhARF8",
  "Rhy_009652_Proteina_Estres_RhSAUR72",
  "Rhy_002145_Sintesis_Antocianina_Color",
  "Rhy_005412_Respiracion_Celular_Hexokinasa",
  "Rhy_099999_Gen_Control_Actina",         # No debería cambiar (Significancia BAJA)
  "Rhy_088888_Factor_Transcripcion_Ruido"  # Mucha dispersión (No significativo)
)

# Creamos la matriz con una dispersión biológica simulada más realista
matriz_simulada <- data.frame(
  Control_Boton_R1   = as.integer(rpois(9, lambda = c(500, 100, 50,  80,  40,  600, 300, 1000, 150))),
  Control_Boton_R2   = as.integer(rpois(9, lambda = c(480, 110, 45,  85,  38,  580, 310, 1020, 80))),
  Control_Boton_R3   = as.integer(rpois(9, lambda = c(510, 95,  52,  78,  42,  610, 295,  980, 220))),
  Tratado_Florero_R1 = as.integer(rpois(9, lambda = c(150, 800, 400, 500, 300, 150, 90,  1010, 180))),
  Tratado_Florero_R2 = as.integer(rpois(9, lambda = c(140, 750, 420, 480, 320, 160, 85,  1005, 95))),
  Tratado_Florero_R3 = as.integer(rpois(9, lambda = c(160, 820, 390, 510, 290, 140, 95,   990, 290))),
  row.names = genes_rosa
)

metadatos_simulados <- data.frame(
  id_muestra = colnames(matriz_simulada),
  estado = factor(c("Boton", "Boton", "Boton", "Florero", "Florero", "Florero")),
  row.names = colnames(matriz_simulada)
)
```
# ==============================================================================
# 2. PIPELINE DE DESEQ2 Y VOLCANO PLOT CORREGIDO
# ==============================================================================
```r
library(DESeq2)
library(ggplot2)

dds <- DESeqDataSetFromMatrix(countData = matriz_simulada,
                              colData = metadatos_simulados,
                              design = ~ estado)
dds$estado <- relevel(dds$estado, ref = "Boton")
dds <- DESeq(dds)
resultados_postcosecha <- results(dds)
df_resultados <- as.data.frame(resultados_postcosecha)

# Definir umbrales estrictos para marcar los Verdaderos DEGs
df_resultados$Significativo <- df_resultados$pvalue < 0.01 & abs(df_resultados$log2FoldChange) > 1

# Gráfico de Volcán Mejorado
ggplot(df_resultados, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significativo), size = 4, alpha = 0.8) +
  scale_color_manual(values = c("darkgray", "firebrick3"), 
                     labels = c("No Significativo", "DEG Significativo")) +
  geom_text(aes(label = rownames(df_resultados)), 
            vjust = -0.8, size = 3, check_overlap = TRUE, fontface = "bold") +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot: Dinámica Transcriptómica en Postcosecha de Rosas",
       subtitle = "Comparación: Florero vs. Botón Floral",
       x = "Magnitud del Cambio (Log2 Fold Change)",
       y = "Significancia Estadística (-Log10 P-value)",
       color = "Filtro de Expresión") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))
```
<img width="811" height="511" alt="heatmap" src="https://github.com/user-attachments/assets/08362f57-f0a5-4528-a0af-230897012499" />
<img width="1000" height="625" alt="volcanoplot" src="https://github.com/user-attachments/assets/c7c4488e-817a-4293-ba0f-834dfd4c0cfc" />

# Bloque de Degradación e Inanición (Cuadrante Superior Derecho / Izquierdo):
Los genes como Aquaporina y Antocianina muestran un desplome masivo ($-\log_{10}(p\text{-value}) > 50$), confirmando molecularmente que el pétalo apaga su maquinaria de hidratación y color. Del mismo modo, el receptor de Etileno (RhETR1) y la Poligalacturonasa se disparan de forma estadísticamente ultra-significativa, marcando el inicio irreversible de la senescencia y el ablandamiento de la pared celular del pétalo.

Control Interno (Rh_099999_Actina): Se posiciona perfectamente en el centro del gráfico (cerca de cero en el eje X y abajo en el eje Y), demostrando que es un control de carga óptimo, ya que su expresión permanece inalterada entre el botón y el florero.Genes con Ruido Fisiológico (Rh_088888_Ruido): A pesar de mostrar cambios puntuales debido a la variabilidad entre las réplicas biológicas de florero, el modelo matemático de Wald en DESeq2 lo penaliza correctamente, situándolo por debajo de la línea de corte de significancia (color gris), tal como ocurre en muestras reales de campo.
