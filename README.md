# RSlope Stability
### Motor de Equilibrio Límite (LEM) y Análisis Estocástico de Confiabilidad en Taludes Multicapa

Este repositorio contiene un algoritmo avanzado de código abierto desarrollado en **R** para la evaluación cuantitativa de la estabilidad de taludes y laderas. El sistema integra métodos de equilibrio límite tradicionales con un robusto módulo estocástico de **Simulación de Monte Carlo** y análisis de sensibilidad paramétrica no lineal, optimizado para perfiles litoestratigráficos complejos.

El software está diseñado bajo estándares académicos listos para soportar la reproducibilidad de datos en artículos científicos indexados y proyectos de ingeniería geotécnica de alta precisión.

## 🚀 Características Principales

* **Motor Multicapa Automatizado:** Algoritmo dinámico para la detección cinemática de interfaces y contactos estratigráficos mediante funciones de aproximación continua.
* **Malla de Búsqueda Flexible (Grid Search):** Evaluación sistemática de miles de combinaciones de centros $(X_c, Z_c)$ y radios $(R)$ de superficies de falla circulares.
* **Análisis de Sensibilidad Geotécnica:** Evaluación determinística del impacto individual del $\pm30\%$ de variación en la Cohesión ($c$) y el Ángulo de Fricción Interna ($\phi$) sobre el Factor de Seguridad base.
* **Análisis Probabilístico de Confiabilidad:** Simulación estocástica iterativa de Monte Carlo ($N = 1000$) para calcular la distribución real de densidad y la **Probabilidad de Falla ($P_f$)** exacta del sistema.
* **Exportación de Alta Resolución (300 DPI):** Rutas relativas automatizadas para compilar reportes estructurales `.txt` y láminas gráficas de calidad comercial o de publicación científica.

## 📁 Estructura del Repositorio

* `RSlope_Stability.R`: Código fuente principal en R con el panel de control y el motor de cálculo matemático.
* `README.md`: Documentación técnica y guía de usuario del repositorio.
* `LICENSE`: Licencia MIT de código abierto.

Una vez ejecutado el script, el sistema creará de forma dinámica la siguiente estructura en tu directorio local de trabajo:

```text
├── resultados_analisis/
│   ├── Reporte_Tecnico_Estabilidad.txt     <- Balances de momentos, fuerzas y factores de seguridad decimales.
│   ├── Grafico_Comparativo_LEM.png         <- Lámina HD con las cuñas y nubes de deslizamiento de los 4 métodos.
│   ├── Grafico_Analisis_Probabilistico.png <- Lámina HD de 3 niveles: Sensibilidad, Histograma y Scatter Plot.
