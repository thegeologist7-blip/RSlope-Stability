# ==============================================================================
# SCRIPT RSlope Stability PRO: PANEL DE CONTROL MULTICAPA Y EXPORTACIÓN TOTAL
# Proyecto: Código Abierto - Análisis de Confiabilidad Geotécnica
# ==============================================================================

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# ==============================================================================
# 📋 [PANEL DE CONTROL CENTRALIZADO] - CONFIGURA TU PROYECTO AQUÍ
# ==============================================================================

# 1. DEFINIR UBICACIÓN DE GUARDADO (Ruta de carpeta para tus resultados)
ruta_proyecto <- "W:/Articulos_escritos/SlopStability R/script/resultados_analisis"

# 2. COORDENADAS DE INTERFACES (Definidas de arriba hacia abajo)
capas_coords <- list(
  superficie = data.frame(
    x = c(0,  20.0, 30.0, 35.0, 50.0, 80.0),
    z = c(40.0, 40.0, 48.0, 48.0, 58.0, 58.0)  # Talud con Berma
  ),
  contacto_1 = data.frame(
    x = c(0, 80.0),
    z = c(46.0, 46.0)              # Interfaz horizontal superior
  ),
  contacto_2 = data.frame(
    x = c(0, 80.0),
    z = c(42.0, 42.0)              # Interfaz horizontal inferior
  )
)

# 3. PROPIEDADES MECÁNICAS DE LOS ESTRATOS
propiedades_estratos <- data.frame(
  nombre = c("Estrato Superior (Suelo)", "Estrato Medio (Tránsito)", "Sustrato Rocoso"),
  gamma  = c(17.5, 18.8, 20.5),
  c      = c(4.0,  15.0, 35.0),
  phi    = c(22.0,  18.0, 26.0),
  color  = c("burlywood1", "darkorange3", "gray70"),
  stringsAsFactors = FALSE
)

# 4. MALLA DE BÚSQUEDA CINEMÁTICA (GRID SEARCH)
grid <- expand.grid(
  xc = seq(20, 55, len=15), 
  zc = seq(60, 78, len=15), 
  R  = seq(20, 42, len=8)
)

# ==============================================================================
# ⚙️ GESTIÓN DE DIRECTORIOS Y LOGÍSTICA DE ARCHIVOS
# ==============================================================================
if(!dir.exists(ruta_proyecto)) dir.create(ruta_proyecto, recursive = TRUE)
archivo_reporte <- file.path(ruta_proyecto, "Reporte_Tecnico_Estabilidad.txt")
sink(archivo_reporte, split = TRUE) # Duplica la salida: Consola y Archivo .txt

# ==============================================================================
# ⚙️ MOTOR NUMÉRICO AUTOMATIZADO 
# ==============================================================================

limpiar_capa <- function(df) {
  df <- na.omit(df); df <- df[!duplicated(df$x), ]; df <- df[order(df$x), ]; return(df)
}
capas_coords <- lapply(capas_coords, limpiar_capa)
lineas_interfaz <- lapply(capas_coords, function(df) approxfun(df$x, df$z, rule = 2))

obtener_id_estrato <- function(px, pz) {
  n_interfaces <- length(lineas_interfaz)
  if(n_interfaces == 1) {
    if(pz <= lineas_interfaz$superficie(px)) return(1)
    return(1)
  }
  zs_niveles <- sapply(lineas_interfaz, function(f) f(px))
  if (pz > zs_niveles[1]) return(1)
  for(k in 1:(n_interfaces - 1)) {
    if(pz <= zs_niveles[k] && pz > zs_niveles[k+1]) return(k)
  }
  return(n_interfaces)
}

calc_fs_pro <- function(xc, zc, R, metodo = "Bishop") {
  x_s <- seq(xc - R + 0.01, xc + R - 0.01, length.out = 400)
  z_c <- zc - sqrt(R^2 - (x_s - xc)^2); z_t <- lineas_interfaz$superficie(x_s)
  idx <- which(z_c < z_t); if(length(idx) < 30) return(NULL)
  df <- data.frame(x = x_s[idx], z = z_c[idx], zs = z_t[idx])
  df$dx <- c(diff(df$x), mean(diff(df$x))); df$h  <- df$zs - df$z
  if(any(df$h <= 0)) return(NULL)
  dz <- c(diff(df$z), 0); df$alpha <- atan2(dz, df$dx)
  df$c_i <- 0; df$phi_i <- 0; df$gamma_i <- 0
  
  for(i in 1:nrow(df)) {
    z_medio <- (df$z[i] + df$zs[i]) / 2
    id <- obtener_id_estrato(df$x[i], z_medio)
    if(is.na(id) || id > nrow(propiedades_estratos)) return(NULL)
    df$c_i[i]     <- propiedades_estratos$c[id]
    df$phi_i[i]   <- propiedades_estratos$phi[id] * pi/180
    df$gamma_i[i] <- propiedades_estratos$gamma[id]
  }
  W <- df$gamma_i * df$h * df$dx; l <- sqrt(df$dx^2 + dz^2)
  
  if(metodo == "Fellenius") {
    resistencia <- df$c_i * l + (W * cos(df$alpha)) * tan(df$phi_i)
    fs_new <- sum(resistencia) / abs(sum(W * sin(df$alpha)))
  } else {
    fs_old <- 1.1
    for(iter in 1:60){
      m_alpha <- cos(df$alpha) + (sin(df$alpha) * tan(df$phi_i)) / fs_old
      m_alpha <- ifelse(m_alpha < 0.01, 0.01, m_alpha)
      if(metodo == "Bishop") {
        num <- (df$c_i * df$dx + W * tan(df$phi_i)) / m_alpha; den <- W * sin(df$alpha)
      } else if(metodo == "Janbu") {
        num <- (df$c_i * df$dx + W * tan(df$phi_i)) / (m_alpha * cos(df$alpha)); den <- W * tan(df$alpha)
      } else if(metodo == "Morgenstern") {
        ma_mp <- (cos(df$alpha) * (1 + tan(df$alpha) * tan(df$phi_i) / fs_old))
        num <- (df$c_i * df$dx + W * tan(df$phi_i)) / ma_mp; den <- W * sin(df$alpha)
      }
      fs_new <- sum(num) / abs(sum(den))
      if(is.nan(fs_new) || is.infinite(fs_new)) return(NULL)
      if(abs(fs_new - fs_old) < 1e-6) break
      fs_old <- fs_new
    }
  }
  ma_f <- cos(df$alpha) + (sin(df$alpha) * tan(df$phi_i)) / fs_new
  N_f  <- (W - (df$c_i * l * sin(df$alpha)) / fs_new) / ma_f
  S_f  <- df$c_i * l + N_f * tan(df$phi_i); T_f  <- W * sin(df$alpha)
  return(list(FS = fs_new, FS_force = sum(S_f)/abs(sum(T_f)), FS_moment = sum(S_f*R)/abs(sum(T_f*R)),
              T = sum(T_f), N = sum(N_f), S = sum(S_f), M_act = abs(sum(T_f*R)), M_res = sum(S_f*R)))
}

get_arc <- function(fs_val, id_val, i_idx) {
  if(is.null(fs_val) || is.infinite(fs_val) || is.nan(fs_val)) return(NULL)
  xc <- grid$xc[i_idx]; zc <- grid$zc[i_idx]; R <- grid$R[i_idx]
  x_arc_s <- seq(xc - R + 0.01, xc + R - 0.01, len=200); inside <- (R^2 - (x_arc_s - xc)^2) > 0; x_arc_s <- x_arc_s[inside]
  if(length(x_arc_s) < 20) return(NULL)
  z_arc_c <- zc - sqrt(R^2 - (x_arc_s - xc)^2); z_arc_t <- lineas_interfaz$superficie(x_arc_s)
  diff_z  <- z_arc_c - z_arc_t; cross <- which(diff_z[-1] * diff_z[-length(diff_z)] < 0)
  if(length(cross) < 2) return(NULL)
  interp_x <- function(idx) { x_arc_s[idx] - diff_z[idx] * (x_arc_s[idx+1] - x_arc_s[idx]) / (diff_z[idx+1] - diff_z[idx]) }
  x_in  <- interp_x(cross[1]); x_out <- interp_x(cross[length(cross)]); if(x_out <= x_in) return(NULL)
  x_final <- c(x_in, x_arc_s[(cross[1]+1):cross[length(cross)]], x_out)
  z_final <- c(lineas_interfaz$superficie(x_in), z_arc_c[(cross[1]+1):cross[length(cross)]], lineas_interfaz$superficie(x_out))
  return(data.frame(x = x_final, z = z_final, fs = fs_val, id = id_val))
}

res_b <- list(); res_f <- list(); res_j <- list(); res_m <- list()
fs_b_vec <- rep(Inf, nrow(grid)); fs_f_vec <- rep(Inf, nrow(grid)); fs_j_vec <- rep(Inf, nrow(grid)); fs_m_vec <- rep(Inf, nrow(grid))

for(i in 1:nrow(grid)) {
  if(grid$zc[i] <= lineas_interfaz$superficie(grid$xc[i])) next
  s_b <- calc_fs_pro(grid$xc[i], grid$zc[i], grid$R[i], metodo = "Bishop")
  s_f <- calc_fs_pro(grid$xc[i], grid$zc[i], grid$R[i], metodo = "Fellenius")
  s_j <- calc_fs_pro(grid$xc[i], grid$zc[i], grid$R[i], metodo = "Janbu")
  s_m <- calc_fs_pro(grid$xc[i], grid$zc[i], grid$R[i], metodo = "Morgenstern")
  fs_b_vec[i] <- if(is.null(s_b)) Inf else s_b$FS; fs_f_vec[i] <- if(is.null(s_f)) Inf else s_f$FS
  fs_j_vec[i] <- if(is.null(s_j)) Inf else s_j$FS; fs_m_vec[i] <- if(is.null(s_m)) Inf else s_m$FS
  if(fs_b_vec[i] < 2.5) res_b[[i]] <- get_arc(fs_b_vec[i], i, i)
  if(fs_f_vec[i] < 2.5) res_f[[i]] <- get_arc(fs_f_vec[i], i, i)
  if(fs_j_vec[i] < 2.5) res_j[[i]] <- get_arc(fs_j_vec[i], i, i)
  if(fs_m_vec[i] < 2.5) res_m[[i]] <- get_arc(fs_m_vec[i], i, i)
}

bind_safe <- function(lst) { lst <- lst[!sapply(lst, is.null)]; if(length(lst) == 0) return(data.frame()); return(do.call(rbind, lst)) }
df_b_all <- bind_safe(res_b); df_f_all <- bind_safe(res_f); df_j_all <- bind_safe(res_j); df_m_all <- bind_safe(res_m)
filtrar <- function(df) { if(nrow(df) == 0) return(df); subset(df, fs <= quantile(df$fs, 0.15, na.rm = TRUE)) }
df_b_f <- filtrar(df_b_all); df_f_f <- filtrar(df_f_all); df_j_f <- filtrar(df_j_all); df_m_f <- filtrar(df_m_all)

crear_plot <- function(df_nube, df_all_met, metodo_txt) {
  if(nrow(df_all_met) == 0) return(ggplot() + ggtitle(paste(metodo_txt, "- SIN SOLUCIONES VÁLIDAS")))
  fs_min <- min(df_all_met$fs, na.rm = TRUE)
  estado_msg <- if(fs_min < 1.0) "FALLA" else if(fs_min < 1.5) "CRÍTICO" else "ESTABLE"
  color_alerta <- if(fs_min < 1.0) "red3" else if(fs_min < 1.5) "orange3" else "darkgreen"
  id_min <- df_all_met$id[which.min(df_all_met$fs)]; f_crit <- subset(df_all_met, id == id_min)
  x_range <- range(capas_coords$superficie$x); x_seq <- seq(min(x_range), max(x_range), len = 400)
  df_plots_estratos <- data.frame(x = x_seq)
  n_total_capas <- length(lineas_interfaz)
  for(m in 1:n_total_capas) df_plots_estratos[[paste0("z_", m)]] <- lineas_interfaz[[m]](x_seq)
  df_plots_estratos$z_base_graf <- min(df_plots_estratos[[paste0("z_", n_total_capas)]]) - 10
  if(n_total_capas > 1) { for(m in 2:n_total_capas) df_plots_estratos[[paste0("z_", m)]] <- pmin(df_plots_estratos[[paste0("z_", m)]], df_plots_estratos[[paste0("z_", m-1)]]) }
  p <- ggplot()
  for(m in 1:(n_total_capas - 1)) p <- p + geom_ribbon(data = df_plots_estratos, aes_string(x="x", ymin=paste0("z_", m+1), ymax=paste0("z_", m)), fill=propiedades_estratos$color[m], alpha=0.6)
  p <- p + geom_ribbon(data = df_plots_estratos, aes_string(x="x", ymin="z_base_graf", ymax=paste0("z_", n_total_capas)), fill=propiedades_estratos$color[n_total_capas], alpha=0.6)
  if(nrow(df_nube) > 0) p <- p + geom_line(data = df_nube, aes(x=x, y=z, color=fs, group=id), alpha=0.15, linewidth=0.50)
  p <- p + geom_line(data = f_crit, aes(x=x, y=z), color="darkred", linewidth=1.5) + geom_line(data = df_plots_estratos, aes(x=x, y=z_1), color="black", linewidth=1) +
    scale_color_gradientn(colors = c("red", "orange", "yellow", "green"), name = "F.S.") + coord_fixed(ratio = 1, xlim = c(min(x_range), max(x_range)), ylim = c(min(df_plots_estratos$z_base_graf), max(df_plots_estratos$z_1) + 5)) +
    labs(title = metodo_txt, subtitle = paste0("FS Mínimo: ", round(fs_min, 3), " [", estado_msg, "]"), x = "Distancia (m)", y = "Elevación (m)") + theme_bw() + theme(plot.subtitle = element_text(color = color_alerta, face = "bold", size = 11))
  return(p)
}

imprimir_resumen_pro <- function(obj, nombre) {
  if(is.null(obj)) { cat(paste0("\n❌ MÉTODO: ", nombre, " - No convergió\n")); return() }
  alerta <- if(obj$FS < 1.0) " [FALLA]" else if(obj$FS < 1.5) " [CRÍTICO]" else " [ESTABLE]"
  cat(paste0("\n✅ MÉTODO: ", nombre, alerta, "\n"))
  cat("----------------------------------------------------\n")
  cat(sprintf("Factor de Seguridad Global:      %.3f\n", obj$FS))
  cat(sprintf("FS Equilibrio de Fuerzas:        %.3f\n", obj$FS_force))
  cat(sprintf("FS Equilibrio de Momentos:       %.3f\n", obj$FS_moment))
  cat(sprintf("Fuerza Normal Total (N):         %.1f kN\n", obj$N))
  cat(sprintf("Fuerza de Corte Actuante (T):    %.1f kN\n", obj$T))
  cat(sprintf("Fuerza Resistente Total (S):     %.1f kN\n", obj$S))
  cat(sprintf("Momento Actuante Total:          %.1f kN·m\n", obj$M_act))
  cat(sprintf("Momento Resistente Total:        %.1f kN·m\n", obj$M_res))
  cat("----------------------------------------------------\n")
}

# ==============================================================================
# 📊 REPORTE DE CONSOLA Y TEXTO AUTOMÁTICO
# ==============================================================================
valid_idx <- which(is.finite(fs_b_vec)); id_ref <- valid_idx[which.min(fs_b_vec[valid_idx])]
c_b <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Bishop")
c_f <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Fellenius")
c_j <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Janbu")
c_m <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Morgenstern")

cat("\n====================================================")
cat("\n   REPORTE TÉCNICO DE ESTABILIDAD - MOTOR INTEGRADO")
cat("\n====================================================\n")
imprimir_resumen_pro(c_b, "BISHOP SIMPLIFICADO")
imprimir_resumen_pro(c_f, "FELLENIUS (ORDINARIO)")
imprimir_resumen_pro(c_j, "JANBU SIMPLIFICADO")
imprimir_resumen_pro(c_m, "MORGENSTERN-PRICE")
sink() 

# ==============================================================================
# 📈 ANÁLISIS DE SENSIBILIDAD PARAMÉTRICA
# ==============================================================================
rango <- seq(0.7, 1.3, by = 0.05)
sens_df <- data.frame()
props_originales <- propiedades_estratos

for(f in rango) {
  propiedades_estratos$c <- props_originales$c * f
  res_c <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Bishop")
  
  propiedades_estratos <- props_originales
  propiedades_estratos$phi <- props_originales$phi * f
  res_phi <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Bishop")
  
  propiedades_estratos <- props_originales
  sens_df <- rbind(sens_df, 
                   data.frame(Variacion = (f-1)*100, FS = if(is.null(res_c)) NA else res_c$FS, Parametro = "Cohesión (c)"),
                   data.frame(Variacion = (f-1)*100, FS = if(is.null(res_phi)) NA else res_phi$FS, Parametro = "Fricción (phi)"))
}
sens_df <- sens_df[is.finite(sens_df$FS), ]

p_sens_pro <- ggplot(sens_df, aes(x = Variacion, y = FS, color = Parametro)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red3", linewidth = 0.8) +
  annotate("text", x = max(sens_df$Variacion), y = 1.0, label = "Falla Teórica (FS = 1.0) ", 
           color = "red3", fontface = "bold", hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = c("Cohesión (c)" = "royalblue", "Fricción (phi)" = "orange2")) +
  labs(title = "Análisis de Sensibilidad Geotécnica (Multi-Estrato)", subtitle = "Variación porcentual de parámetros en todas las capas simultáneamente", x = "% de Cambio en los Parámetros Originales", y = "Factor de Seguridad (FS)") +
  theme_minimal() + theme(legend.position = "bottom")

# ==============================================================================
# 📈 ANÁLISIS PROBABILÍSTICO (MONTE CARLO AVANZADO)
# ==============================================================================
set.seed(123); n_sim <- 1000
mc_data <- data.frame(sim = 1:n_sim, c_sim = numeric(n_sim), phi_sim = numeric(n_sim), fs = numeric(n_sim))
props_base <- propiedades_estratos

for(i in 1:n_sim) {
  f_c <- rlnorm(1, meanlog = 0, sdlog = 0.15)
  f_phi <- max(0.7, min(1.3, rnorm(1, mean = 1, sd = 0.10)))
  
  props_temp <- data.frame(nombre = props_base$nombre, gamma = props_base$gamma,
                           c = props_base$c * f_c, phi = pmax(5, pmin(45, props_base$phi * f_phi)), stringsAsFactors = FALSE)
  
  mc_data$c_sim[i]   <- props_temp$c[1]     
  mc_data$phi_sim[i] <- props_temp$phi[1]
  
  propiedades_estratos <- props_temp
  res <- calc_fs_pro(grid$xc[id_ref], grid$zc[id_ref], grid$R[id_ref], metodo = "Bishop")
  mc_data$fs[i] <- if(is.null(res)) NA else res$FS
}
propiedades_estratos <- props_base
mc_data <- mc_data[is.finite(mc_data$fs), ]
prob_falla <- mean(mc_data$fs < 1.0) * 100
mc_data$Estado <- ifelse(mc_data$fs < 1.0, "Falla (FS < 1.0)", "Estable (FS >= 1.0)")

p_prob_pro <- ggplot(mc_data, aes(x = fs)) +
  geom_histogram(aes(y = after_stat(density)), fill = "slategray3", color = "white", bins = 30, alpha=0.8) +
  geom_density(color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = 1.0, color = "red3", linewidth = 1.2, linetype = "longdash") +
  annotate("label", x = 1.0, y = 0.2, label = paste("P(Falla):", round(prob_falla, 2), "%"), color = "white", fill = "red3", fontface = "bold") +
  labs(title = "Distribución de Probabilidad del Factor de Seguridad", subtitle = paste("Simulación de Monte Carlo (N =", nrow(mc_data), ")"), x = "Factor de Seguridad (FS)", y = "Densidad de Probabilidad") +
  theme_minimal()

# ==============================================================================
# 📈 ANÁLISIS PROBABILÍSTICO (Dispersión Cruzada)
# ==============================================================================

p_scatter_pro <- ggplot(mc_data, aes(x = phi_sim, y = c_sim, color = Estado)) +
  geom_point(alpha = 0.7, size = 1.8) +
  # 🔹 AGREGAR ESTA LÍNEA PARA MARCAR EL VALOR DE DISEÑO ORIGINAL:
  geom_point(aes(x = props_base$phi[1], y = props_base$c[1]), color = "black", size = 4, shape = 18) +
  annotate("text", x = props_base$phi[1] + 1, y = props_base$c[1] + 0.5, label = "Línea Base (Diseño)", color = "black", fontface = "bold", hjust = 0) +
  scale_color_manual(values = c("Estable (FS >= 1.0)" = "seagreen3", "Falla (FS < 1.0)" = "red2")) +
  labs(title = "Nube de Dispersión Geotécnica y Envolvente de Falla", subtitle = "Interacción Cruzada de Parámetros Simulados (Estrato Superior)", x = "Ángulo de Fricción Interna \u03c6 (°)", y = "Cohesión c (kPa)") +
  theme_light() + theme(legend.position = "bottom")

# ==============================================================================
# 💾 EXPORTACIÓN AUTOMÁTICA EN ALTA RESOLUCIÓN (300 DPI)
# ==============================================================================
plot_bishop <- crear_plot(df_b_f, df_b_all, "BISHOP SIMPLIFICADO")
plot_fellen <- crear_plot(df_f_f, df_f_all, "FELLENIUS (ORDINARIO)")
plot_janbu  <- crear_plot(df_j_f, df_j_all, "JANBU SIMPLIFICADO")
plot_morgen <- crear_plot(df_m_f, df_m_all, "MORGENSTERN-PRICE")

g_cuadricula <- marrangeGrob(list(plot_bishop, plot_fellen, plot_janbu, plot_morgen), ncol = 2, nrow = 2, top="Comparativa Global LEM")

# Trazamos la lámina estadística extendida de 3 paneles (Sensibilidad + Histograma + Scatter)
g_estadistico <- marrangeGrob(list(p_sens_pro, p_prob_pro, p_scatter_pro), ncol = 1, nrow = 3, top="Análisis Estocástico y de Sensibilidad de Confiabilidad")

ggsave(file.path(ruta_proyecto, "Grafico_Comparativo_LEM.png"), g_cuadricula, width = 11, height = 8.5, dpi = 300)
ggsave(file.path(ruta_proyecto, "Grafico_Analisis_Probabilistico.png"), g_estadistico, width = 8.5, height = 14, dpi = 300)

cat(paste0("\n💾 ¡Proceso finalizado con éxito! Los reportes y gráficos han sido exportados a:\n👉 ", ruta_proyecto, "\n"))