# MEMORY.md — Contexto del Proyecto para Futuras Sesiones

## ¿Qué es este proyecto?

Paquete de software en **C estándar (C99)** que implementa de forma exacta la teoría y algoritmos de **Cadenas de Markov** y **Procesos Markovianos de Decisión (PMD)**. Todo el código y comentarios están en español.

## Ubicación

```
C:\Users\PC\Documents\PROYECTOS\Estocasticos\
├── markov_lib/          ← Paquete principal
│   ├── matrix.h         ← TDA Matriz dinámica
│   ├── matrix.c         ← Operaciones matriciales + Gauss
│   ├── markov.h         ← Estructuras y prototipos
│   ├── markov_datos.c   ← Lectura de datos (consola + archivo)
│   ├── markov_cadenas.c ← Módulo 1: teoría básica de cadenas
│   ├── markov_decision.c← Módulo 2: 5 algoritmos PMD + Simplex
│   ├── main.c           ← Programa interactivo con menú
│   ├── Makefile         ← Compilación Linux/Mac/MinGW
│   ├── compilar.bat     ← Compilación Windows
│   ├── ejemplo_pmd.txt  ← Archivo de ejemplo (3 estados, 2 decisiones)
│   └── README.md        ← Documentación completa (fórmulas + manual)
└── .git/
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
./markov                    # Menú interactivo
./markov ejemplo_pmd.txt    # Cargar archivo directamente
```

## Estructura de datos principal

```c
typedef struct {
    int num_estados;        // m+1 (estados 0..m)
    int num_decisiones;     // K
    int tiene_decisiones;   // 0 = cadena simple, 1 = PMD
    double *prob_inicial;   // vector a
    Matriz *P;              // matriz transición (cadena simple)
    Matriz **P_dec;         // K matrices de transición (PMD)
    Matriz *C;              // matriz de costos (m+1)×K
    double alfa;            // factor de descuento
    double tasa_interes;    // tasa i, alfa = 1/(1+i)
} ModeloMarkov;
```

## Algoritmos implementados

### Módulo 1: Teoría de Cadenas de Markov
| # | Función | Archivo | Línea |
|---|---------|---------|-------|
| 1 | Chapman-Kolmogorov (P^n) | `markov_cadenas.c` | ~25 |
| 2 | Probabilidades incondicionales | `markov_cadenas.c` | ~45 |
| 3 | Estado estable π | `markov_cadenas.c` | ~75 |
| 4 | Tiempos de recurrencia μ_ii | `markov_cadenas.c` | ~125 |
| 5 | Tiempos de primera pasada μ_ij | `markov_cadenas.c` | ~140 |
| 6 | Probabilidades de absorción f_ik | `markov_cadenas.c` | ~190 |

### Módulo 2: Procesos Markovianos de Decisión
| # | Algoritmo | Archivo | Línea |
|---|-----------|---------|-------|
| 1 | Enumeración exhaustiva | `markov_decision.c` | ~120 |
| 2 | Mejoramiento de políticas | `markov_decision.c` | ~210 |
| 3 | Mejoramiento con descuento | `markov_decision.c` | ~380 |
| 4 | Aproximaciones sucesivas | `markov_decision.c` | ~530 |
| 5 | Programación lineal (Simplex Big-M) | `markov_decision.c` | ~670 |

## Dependencias internas

```
matrix.h → markov.h → {markov_datos.c, markov_cadenas.c, markov_decision.c}
                                  ↓
                              main.c
```

- `resolver_sistema_lineal()` (eliminación gaussiana con pivoteo) es la función más usada: aparece en estado estable, primera pasada, absorción, y los pasos de mejoramiento de políticas.
- `matriz_potencia()` usa exponenciación binaria O(log n).
- `markov_decision.c` incluye un solver Simplex completo con método de la Gran M para la PL.

## Formato de archivos de entrada

### Cadena simple
```
0              ← tipo
2              ← m (3 estados)
0.4 0.3 0.3    ← a
0.5 0.3 0.2    ← P fila 0
0.1 0.7 0.2    ← P fila 1
0.3 0.3 0.4    ← P fila 2
```

### PMD
```
1              ← tipo PMD
2              ← m
2              ← K
0.5 0.3 0.2    ← a
0.7 0.2 0.1    ← P^(1) fila 0
0.4 0.4 0.2    ← P^(1) fila 1
0.3 0.3 0.4    ← P^(1) fila 2
0.6 0.2 0.2    ← P^(2) fila 0
0.2 0.5 0.3    ← P^(2) fila 1
0.8 0.1 0.1    ← P^(2) fila 2
10.0 5.0       ← C fila 0
8.0 12.0       ← C fila 1
6.0 9.0        ← C fila 2
0.9            ← alfa
```

## Verificaciones realizadas

- Compilación: exitosa con GCC (TDM-GCC-64), sin warnings
- Carga de modelo desde archivo: verificada
- Estado estable: verificado (π = [0.5455, 0.2727, 0.1818] para P^(1))
- Mejoramiento con descuento: verificado, converge a política óptima [2, 1, 1]
- Aproximaciones sucesivas: verificado, converge a misma política
- Todos los métodos imprimen estado en cada iteración (requisito obligatorio)

## Idioma

Todo el código, comentarios, mensajes al usuario y documentación están en **español**, conforme a la solicitud del usuario.

## Convenciones del código

- Nombres de funciones en español (`matriz_crear`, `estado_estable`, `mejoramiento_politicas`)
- Nombres de tipos en español (`Matriz`, `ModeloMarkov`, `Politica`)
- Variables con nombres descriptivos en español cuando es natural
- Comentarios con `/* ====== */` para secciones grandes
- Indentación: 4 espacios
- Sin dependencias externas (solo C estándar + math.h)

## Posibles extensiones futuras

- Soporte para cadenas de Markov de tiempo continuo
- Exportación de resultados a CSV/JSON
- Interfaz gráfica (GTK/ncurses)
- Soporte para políticas aleatorias (no determinísticas)
- Paralelización de enumeración exhaustiva con OpenMP
- Más formatos de entrada (JSON, YAML)
