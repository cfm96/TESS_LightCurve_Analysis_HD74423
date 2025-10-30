# 🌌 Proyecto TESS — Análisis de Curva de Luz HD74423

Este proyecto implementa un análisis completo de la **curva de luz del sistema estelar HD74423**, utilizando datos del satélite **TESS (Transiting Exoplanet Survey Satellite)**.  
La implementación está basada en los tutoriales publicados en **U-Cursos** y desarrollada completamente en **Python 3**.

---

## 🚀 Ejecución

Para obtener los gráficos y resultados desde una consola **GNU/Linux**, ejecutar:

```bash
$ python3 tess_project.py
```

El programa analiza los datos contenidos en el archivo:
```
tess2019112060037-s0011-0000000355151781-0143-s_lc.fits
```

---

## 🧠 Funcionalidades principales

- **Lectura de archivos .FITS** descargados desde el catálogo MAST (NASA)
- **Limpieza y normalización de datos**
  - Eliminación de valores nulos y fuera de rango
  - Guardado de versiones intermedias (`light_curve.dat`, `normalized_data.dat`)
- **Análisis espectral mediante el método de Lomb–Scargle**
  - Identificación de frecuencias dominantes
  - Cálculo de periodos principales
- **Ajuste y graficación**
  - Curva de luz original y filtrada  
  - Periodograma Lomb–Scargle  
  - Ajuste de fase para los tres peaks principales

---

## 📊 Salidas generadas

El programa produce varios gráficos en formato **EPS**:

| Archivo | Descripción |
|----------|--------------|
| `exit_0.eps` | Curva de luz sin filtrar |
| `exit_1.eps` | Curva de luz normalizada |
| `exit_2.eps` | Periodograma Lomb–Scargle |
| `exit_3.eps` | Ajuste de fase para los tres peaks detectados |

También genera archivos `.dat` con los datos intermedios procesados.

---

## 🧩 Estructura del código

| Función | Descripción |
|----------|--------------|
| `data_handling()` | Controla la lectura, limpieza, análisis y graficación de la curva de luz. |
| `data_clean()` | Normaliza los datos y elimina valores fuera de rango. |
| `whitout_outliers()` | Filtra los valores anómalos según un umbral definido. |
| `ls_graph()` | Genera el periodograma de Lomb–Scargle y destaca los tres peaks principales. |
| `curve_graph()` | Grafica la curva de luz limpia y su posible ajuste. |
| `curve_graph_not_clean()` | Grafica la curva sin filtrar para comparación. |
| `data_phasing()` | Aplica ajuste de fase para los tres periodos dominantes. |
| `sin_fit()` | Ajusta una función seno a los datos para modelar la periodicidad. |

---

## 🧰 Requerimientos

Instala las dependencias necesarias ejecutando:

```bash
pip install numpy matplotlib astropy scipy
```

## 📁 Archivo principal

```
tess_project.py
```

Ejecuta todo el proceso de análisis, limpieza y graficación de resultados.
