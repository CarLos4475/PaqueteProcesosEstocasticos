# MEMORY.md — Contexto del Proyecto para Futuras Sesiones

## ¿Qué es este proyecto?

Paquete de software en **C estándar (C99)** que implementa de forma exacta la teoría y algoritmos de **Cadenas de Markov** y **Procesos Markovianos de Decisión (PMD)**. Todo el código y comentarios están en español.

## Ubicación

```
C:\Users\PC\Documents\PROYECTOS\Estocasticos\
├── MEMORY.md                    ← Este archivo
├── CASO_PUBLICIDAD.md           ← Caso de prueba: problema de publicidad
└── markov_lib/                  ← Paquete principal
    ├── matrix.h / matrix.c      ← TDA Matriz + Gauss
    ├── markov.h                 ← Estructuras y prototipos
    ├── markov_datos.c           ← Lectura (consola/archivo) + auto-gen C_ik
    ├── markov_cadenas.c         ← Módulo 1: CK, π, μ, absorción
    ├── markov_decision.c        ← Módulo 2: 5 algoritmos PMD + Simplex + batch
    ├── main.c                   ← Menú interactivo
    ├── Makefile / compilar.bat
    ├── ejemplo_pmd.txt          ← Ejemplo PMD básico (3 est, 2 dec, costos directos)
    ├── test_publicidad.txt      ← Caso publicidad (3 est, 3 dec, ingresos + max)
    └── README.md                ← Documentación completa
```

## Cómo compilar

```bash
# Windows (requiere MinGW/GCC)
gcc -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm

# Linux/Mac
make && ./markov
```

## Cómo ejecutar

```bash
./markov                              # Menú interactivo
./markov test_publicidad.txt          # Cargar problema de publicidad
./markov ejemplo_pmd.txt              # Cargar ejemplo PMD básico
```

## Estructura de datos principal (ACTUALIZADA)

```c
typedef struct {
    int     num_estados;        // m+1 (estados 0..m)
    int     num_decisiones;     // K
    int     tiene_decisiones;   // 0 = cadena simple, 1 = PMD
    double *prob_inicial;       // vector a
    Matriz *P;                  // matriz transición (cadena simple)
    Matriz **P_dec;             // K matrices de transición (PMD)
    Matriz  *C;                 // matriz de costos (m+1)×K
    double  alfa;               // factor de descuento
    double  tasa_interes;       // tasa i, alfa = 1/(1+i)

    // --- NUEVO: ingresos y maximización ---
    int     es_maximizacion;    // 1=MAX utilidad, 0=MIN costo
    int     usa_ingresos;       // 1=C generado desde ingresos
    Matriz **Ingreso;           // K matrices de ingreso (m+1)×(m+1)
    double  *costo_fijo;        // costo fijo por decisión (tam K)
} ModeloMarkov;
```

## Algoritmos implementados

### Módulo 1: Teoría de Cadenas de Markov
| # | Función | Archivo |
|---|---------|---------|
| 1 | Chapman-Kolmogorov (P^n) | `markov_cadenas.c` |
| 2 | Probabilidades incondicionales | `markov_cadenas.c` |
| 3 | Estado estable π | `markov_cadenas.c` |
| 4 | Tiempos de recurrencia μ_ii | `markov_cadenas.c` |
| 5 | Tiempos de primera pasada μ_ij | `markov_cadenas.c` |
| 6 | Probabilidades de absorción f_ik | `markov_cadenas.c` |

### Módulo 2: Procesos Markovianos de Decisión
| # | Algoritmo | Archivo |
|---|-----------|---------|
| 1 | Enumeración exhaustiva | `markov_decision.c` |
| 2 | Mejoramiento de políticas (sin desc) | `markov_decision.c` |
| 3 | Mejoramiento de políticas (con desc) | `markov_decision.c` |
| 4 | Aproximaciones sucesivas | `markov_decision.c` |
| 5 | Programación lineal (Simplex Big-M) | `markov_decision.c` |
| 6 | **Prueba completa batch** (ejecuta 1,2,3,5) | `markov_decision.c` |

## Novedades agregadas (sesión 2026-05-05)

### 1. Soporte para maximización / minimización
- El usuario elige si quiere MAXIMIZAR utilidad o MINIMIZAR costo
- Al maximizar, los costos se niegan internamente: min(-utilidad) = max(utilidad)
- El flag `es_maximizacion` en el modelo controla este comportamiento
- En los resultados se muestra "Utilidad esperada" en vez de "Costo esperado"

### 2. Auto-generación de C_ik desde matrices de ingreso
- Fórmula: `C_ik = sum_j P_ij(k) * Ingreso_ij(k) - Costo_fijo(k)`
- `modelo_calcular_costos_ingresos()` calcula todo automáticamente
- Soporte en consola (leer_matriz para cada Ingreso^(k)) y archivo (modo_costos=1 o 2)
- Imprime cada cálculo paso a paso para verificación

### 3. Modo batch — prueba completa
- `prueba_completa_pmd()` ejecuta enumeración, mejoramiento s/d, mejoramiento c/d, y PL en secuencia
- Imprime todos los resultados con encabezados claros (####)
- Menú Módulo 2, opción 6

### 4. Nuevo formato de archivo (línea 4: modo_costos)
```
1              ← tipo PMD
2              ← m
3              ← K
1              ← modo_costos: 0=directos, 1=ingresos+maxim, 2=ingresos+minim
0.4 0.3 0.3    ← a
...P^(1), P^(2), P^(3)...
...Ingreso^(1), Ingreso^(2), Ingreso^(3)...
200 900 300    ← costos fijos
0.95           ← alfa
```

## Verificaciones realizadas

- Compilación: exitosa con GCC (TDM-GCC-64), sin warnings ni errores
- Carga de modelo desde archivo: verificada (ejemplo_pmd.txt y test_publicidad.txt)
- Auto-generación C_ik: verificada, coincide con cálculo manual
- Estado estable: verificado (π = [0.5455, 0.2727, 0.1818])
- Mejoramiento con descuento: verificado (política [2,1,1])
- Aproximaciones sucesivas: verificado, converge a misma política
- **Los 4 métodos del caso publicidad convergen a la misma política óptima [1,3,3] con utilidad $282.00**
- Modo batch: ejecuta los 4 métodos secuencialmente sin errores

## Idioma y convenciones

- Todo el código, comentarios y mensajes en **español**
- Nombres de funciones/types en español: `Matriz`, `ModeloMarkov`, `Politica`, `estado_estable`, `mejoramiento_politicas`
- Comentarios de sección con `/* ====== */`
- Indentación: 4 espacios
- Sin dependencias externas (solo C estándar + math.h)
- Estándar: C99
