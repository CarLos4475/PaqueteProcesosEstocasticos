# Documentación y Manual de Usuario

## Paquete Estocástico — Cadenas de Markov y PMD en C

---

# Índice

1. [¿Qué es este paquete?](#1-qué-es-este-paquete)
2. [Fundamentos teóricos](#2-fundamentos-teóricos)
   - [Cadenas de Markov](#21-cadenas-de-markov)
   - [Procesos Markovianos de Decisión](#22-procesos-markovianos-de-decisión-pmd)
3. [Estructura del proyecto](#3-estructura-del-proyecto)
4. [Instalación y compilación](#4-instalación-y-compilación)
5. [Formato de archivos de entrada](#5-formato-de-archivos-de-entrada)
   - [Cadena simple](#51-cadena-de-markov-simple)
   - [PMD con costos directos](#52-pmd-con-costos-directos)
   - [PMD con ingresos y costos fijos](#53-pmd-con-ingresos-y-costos-fijos)
6. [TDA Matriz — El motor de cálculo](#6-tda-matriz--el-motor-de-cálculo)
7. [Módulo 1: Teoría de Cadenas de Markov](#7-módulo-1-teoría-de-cadenas-de-markov)
   - [1. Chapman-Kolmogorov](#71-chapman-kolmogorov-pn)
   - [2. Probabilidades Incondicionales](#72-probabilidades-incondicionales)
   - [3. Estado Estable π](#73-estado-estable-π)
   - [4. Tiempos de Recurrencia](#74-tiempos-de-recurrencia-μ_ii)
   - [5. Tiempos de Primera Pasada](#75-tiempos-de-primera-pasada-μ_ij)
   - [6. Probabilidades de Absorción](#76-probabilidades-de-absorción-f_ik)
8. [Módulo 2: Procesos Markovianos de Decisión](#8-módulo-2-procesos-markovianos-de-decisión)
   - [1. Enumeración Exhaustiva](#81-enumeración-exhaustiva)
   - [2. Mejoramiento de Políticas (sin descuento)](#82-mejoramiento-de-políticas-sin-descuento)
   - [3. Mejoramiento de Políticas (con descuento)](#83-mejoramiento-de-políticas-con-descuento)
   - [4. Aproximaciones Sucesivas](#84-aproximaciones-sucesivas)
   - [5. Programación Lineal](#85-programación-lineal)
   - [6. Prueba Completa (Batch)](#86-prueba-completa-batch)
9. [Manual de Usuario](#9-manual-de-usuario)
   - [Ejecución del programa](#91-ejecución-del-programa)
   - [Menú principal](#92-menú-principal)
   - [Módulo 1 paso a paso](#93-módulo-1-paso-a-paso)
   - [Módulo 2 paso a paso](#94-módulo-2-paso-a-paso)
   - [Carga por archivo vs consola](#95-carga-por-archivo-vs-consola)
10. [Ejemplos resueltos](#10-ejemplos-resueltos)
    - [Ejemplo 1: PMD básico (3 estados, 2 decisiones)](#101-ejemplo-1-pmd-básico)
    - [Ejemplo 2: Caso de publicidad (maximización)](#102-ejemplo-2-caso-de-publicidad)
    - [Ejemplo 3: Póker de los sábados (costo promedio)](#103-ejemplo-3-póker-de-los-sábados)
11. [Referencia rápida de API](#11-referencia-rápida-de-api)
12. [Limitaciones y recomendaciones](#12-limitaciones-y-recomendaciones)

---

# 1. ¿Qué es este paquete?

Es una biblioteca de software escrita en **C estándar (C99)** que implementa de forma exacta la teoría y los algoritmos de:

- **Cadenas de Markov** (6 algoritmos fundamentales)
- **Procesos Markovianos de Decisión (PMD)** (5 métodos de resolución)

El paquete no depende de ninguna biblioteca externa, solo usa `stdio.h`, `stdlib.h`, `math.h` y `string.h`. Todo el código, comentarios y mensajes están en **español**.

Puede resolver problemas de **minimización de costos** y **maximización de utilidades**, con o sin factor de descuento, y admite tanto costos directos como costos generados automáticamente a partir de matrices de ingreso.

---

# 2. Fundamentos teóricos

## 2.1 Cadenas de Markov

Una **Cadena de Markov** es un proceso estocástico `{X_n}` que cumple la **propiedad de Markov**: el futuro depende solo del presente, no del pasado.

```
P(X_{n+1} = j | X_n = i, X_{n-1} = i_{n-1}, ..., X_0 = i_0) = P(X_{n+1} = j | X_n = i) = p_{ij}
```

### Elementos de una cadena de Markov

| Elemento | Notación | Significado |
|----------|----------|-------------|
| **Estados** | `E = {0, 1, ..., m}` | Los `m+1` estados posibles del sistema |
| **Matriz de transición** | `P` de tamaño `(m+1)×(m+1)` | `p_ij = P(X_{n+1}=j | X_n=i)` — probabilidad de pasar de i a j en un paso |
| **Vector inicial** | `a` de tamaño `(m+1)` | `a_i = P(X_0 = i)` — probabilidad de empezar en cada estado |
| **Estado estable** | `π` de tamaño `(m+1)` | Distribución límite cuando `n → ∞`, cumple `π = π·P` |

### Propiedad fundamental de las filas de P

Cada fila de `P` debe sumar **exactamente 1**, porque desde un estado i, necesariamente se va a algún estado j:

```
p_{i0} + p_{i1} + ... + p_{im} = 1   para todo i
```

### ¿Qué calcula este módulo?

| # | Algoritmo | Fórmula | ¿Qué obtengo? |
|---|-----------|---------|---------------|
| 1 | Chapman-Kolmogorov | `P^(n) = P^n` | Probabilidades de transición en `n` pasos |
| 2 | Prob. Incondicionales | `a · P^n` | Probabilidad de estar en cada estado en el paso `n` |
| 3 | Estado estable | `π = π·P`, `Σπ_i = 1` | Distribución de largo plazo (estacionaria) |
| 4 | Tiempos de recurrencia | `μ_ii = 1/π_i` | Pasos esperados para regresar al estado `i` |
| 5 | Tiempos de 1ª pasada | Sistema lineal | Pasos esperados para llegar de `i` a `j` por 1ª vez |
| 6 | Prob. de absorción | Sistema lineal | Probabilidad de ser absorbido en `k` empezando en `i` |

---

## 2.2 Procesos Markovianos de Decisión (PMD)

Un **PMD** extiende una cadena de Markov añadiendo **decisiones** que el tomador de decisiones puede elegir en cada estado. Cada decisión tiene su propia matriz de transición y su propio costo.

### Elementos adicionales de un PMD

| Elemento | Notación | Significado |
|----------|----------|-------------|
| **Decisiones** | `K = {1, 2, ..., K}` | Acciones disponibles en cada estado |
| **Matrices de transición** | `P^(k)` para `k = 1..K` | Cada decisión `k` tiene su propia matriz |
| **Matriz de costos** | `C` de tamaño `(m+1)×K` | `C_{ik}` = costo inmediato de tomar decisión `k` en estado `i` |
| **Factor de descuento** | `α` (alfa), `0 < α ≤ 1` | Peso de costos futuros. `α = 1/(1+i)` donde `i` es tasa de interés |
| **Política** | `R = [d(0), d(1), ..., d(m)]` | Regla que asigna una decisión a cada estado |

### Objetivo de un PMD

Encontrar la **política óptima** `R*` que minimice (o maximice) el costo (o utilidad) esperado(a) a largo plazo.

```
min_R  E[ Σ_{n=0}^∞ α^n · C_{X_n, d(X_n)} ]
```

### ¿Qué calcula este módulo?

| # | Método | Tipo | ¿Cuándo usarlo? |
|---|--------|------|-----------------|
| 1 | Enumeración exhaustiva | Exacto | Pocos estados/decisiones (K^(m+1) ≤ ~10⁶) |
| 2 | Mejoramiento s/desc. | Exacto (Howard) | Problemas de **costo promedio** (α ≈ 1) |
| 3 | Mejoramiento c/desc. | Exacto (Howard) | Problemas con **descuento** (α < 1) |
| 4 | Aproximaciones sucesivas | Iterativo | Problemas grandes, converge por valor |
| 5 | Programación Lineal | Exacto (Simplex) | Verificación formal, formulación dual |
| 6 | Prueba completa | Batch | Ejecuta 1,2,3,5 en secuencia para verificar |

---

# 3. Estructura del proyecto

```
Estocasticos/
├── MEMORY.md                        ← Contexto para sesiones de IA
├── CASO_PUBLICIDAD.md               ← Planteamiento del problema de publicidad
├── PROBLEMA2.md                     ← Planteamiento del problema del póker
├── SOLUCION_PROBLEMA2.md            ← Solución detallada del problema 2
├── DOCUMENTACION/
│   └── MANUAL.md                    ← Este documento
└── markov_lib/                      ← Paquete principal
    ├── matrix.h / matrix.c          ← TDA Matriz + sistema de impresión + Gauss
    ├── markov.h                     ← Estructuras y prototipos (interfaz pública)
    ├── markov_datos.c               ← Lectura de datos (consola y archivo)
    ├── markov_cadenas.c             ← Módulo 1: 6 algoritmos de cadenas
    ├── markov_decision.c            ← Módulo 2: 5 algoritmos PMD + Simplex + batch
    ├── main.c                       ← Programa principal con menú interactivo
    ├── Makefile / compilar.bat      ← Scripts de compilación
    ├── ejemplo_pmd.txt              ← Ejemplo PMD (3 est, 2 dec, costos directos)
    ├── test_publicidad.txt          ← Caso publicidad (3 est, 3 dec, ingresos+max)
    └── problema2.txt                ← Caso póker (2 est, 2 dec, costos directos)
```

Cada archivo `.c` tiene una responsabilidad única:

| Archivo | Responsabilidad | Líneas |
|---------|----------------|--------|
| `matrix.c` | Álgebra lineal: crear, multiplicar, invertir, imprimir matrices | 326 |
| `markov_datos.c` | Parseo de entrada (consola interactiva + archivos de texto) | 468 |
| `markov_cadenas.c` | Algoritmos puros de cadenas de Markov | 332 |
| `markov_decision.c` | Algoritmos PMD + Simplex desde cero | 1221 |
| `main.c` | Interfaz de usuario (menús) | 345 |

---

# 4. Instalación y compilación

### Requisitos

- **Compilador GCC** (MinGW en Windows, GCC en Linux/Mac)
- Sin dependencias externas

### Compilación

**Windows (PowerShell / cmd):**
```bash
gcc -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm
```

**Linux / Mac:**
```bash
cd markov_lib
make && ./markov
```

El flag `-lm` enlaza la biblioteca matemática (`math.h`) necesaria para `pow()`, `fabs()`, etc.

### Verificar que compila correctamente

```bash
./markov ejemplo_pmd.txt
```

Si ves el mensaje `--- Datos leidos correctamente ---` y el menú principal, todo funciona.

---

# 5. Formato de archivos de entrada

El paquete lee problemas desde archivos `.txt`. Existen tres modalidades según el tipo de problema.

### Reglas generales

- Cada número se separa por espacios o saltos de línea.
- Las filas de matrices se escriben una por línea.
- Las probabilidades deben sumar 1 por fila (el programa normaliza automáticamente si hay pequeñas desviaciones).
- Los estados se indexan desde **0** hasta **m**.

---

## 5.1 Cadena de Markov simple

**Formato:**
```
0              ← tipo = 0 (cadena simple, sin decisiones)
m              ← último índice de estado (E = {0, 1, ..., m})
a_0 a_1 ... a_m   ← vector de probabilidad inicial
p_00 p_01 ... p_0m   ← fila 0 de P
p_10 p_11 ... p_1m   ← fila 1 de P
...
p_m0 p_m1 ... p_mm   ← fila m de P
```

**Ejemplo** (2 estados, E = {0, 1}):
```
0
1
0.6 0.4
0.7 0.3
0.2 0.8
```

---

## 5.2 PMD con costos directos

**Formato:**
```
1              ← tipo = 1 (PMD)
m              ← último índice de estado
K              ← número de decisiones
0              ← modo_costos = 0 (costos directos)
a_0 a_1 ... a_m   ← vector inicial
[fila 0 P^(1)]
[fila 1 P^(1)]
...
[fila m P^(1)]
[fila 0 P^(2)]
...
[fila m P^(K)]     ← K bloques de (m+1) filas cada uno
C_00 C_01 ... C_0,K-1   ← costos: estado 0
C_10 C_11 ... C_1,K-1   ← costos: estado 1
...
C_m0 C_m1 ... C_m,K-1   ← costos: estado m
alfa           ← factor de descuento (0 < α ≤ 1)
```

**Ejemplo real** (`ejemplo_pmd.txt`): 3 estados, 2 decisiones, costos directos, α = 0.9

```
1
2
2
0
0.5 0.3 0.2
0.7 0.2 0.1
0.4 0.4 0.2
0.3 0.3 0.4
0.6 0.2 0.2
0.2 0.5 0.3
0.8 0.1 0.1
10.0 5.0
8.0 12.0
6.0 9.0
0.9
```

**Explicación línea por línea:**

| Línea | Valor | Significado |
|-------|-------|-------------|
| 1 | `1` | Es un PMD |
| 2 | `2` | m=2 → 3 estados: {0, 1, 2} |
| 3 | `2` | K=2 decisiones |
| 4 | `0` | modo_costos=0 → costos directos |
| 5 | `0.5 0.3 0.2` | Vector a: P(X₀=0)=0.5, P(X₀=1)=0.3, P(X₀=2)=0.2 |
| 6-8 | `0.7 0.2 0.1` ... | P^(1): matriz de transición bajo decisión 1 |
| 9-11 | `0.6 0.2 0.2` ... | P^(2): matriz de transición bajo decisión 2 |
| 12-14 | `10.0 5.0` ... | Matriz C de costos (3 filas × 2 columnas) |
| 15 | `0.9` | α = 0.9 (factor de descuento) |

---

## 5.3 PMD con ingresos y costos fijos

Cuando el problema se plantea en términos de **ingresos** (ganancias) en vez de costos directos, se usa `modo_costos = 1` (maximizar) o `modo_costos = 2` (minimizar). El paquete **auto-genera** la matriz de costos `C` usando la fórmula:

```
C_ik = Σ_j P_ij(k) × Ingreso_ij(k)  −  CostoFijo(k)
```

Si `es_maximizacion = 1`, el valor se niega porque internamente el paquete siempre minimiza:

```
min (−utilidad) ≡ max (utilidad)
```

**Formato:**
```
1              ← tipo PMD
m              ← último índice de estado
K              ← número de decisiones
1              ← modo_costos: 1=maximizar, 2=minimizar
a_0 ... a_m    ← vector inicial
[K bloques de P^(k): (m+1)×(m+1) cada uno]
[K bloques de Ingreso^(k): (m+1)×(m+1) cada uno]
costo_fijo_1 costo_fijo_2 ... costo_fijo_K
alfa           ← factor de descuento
```

**Ejemplo real** (`test_publicidad.txt`): 3 estados, 3 decisiones, maximización con ingresos, α = 0.95

```
1
2
3
1
0.4 0.3 0.3
...P^(1), P^(2), P^(3)...   ← 3×3 filas = 9 líneas de transiciones
...Ingreso^(1), Ingreso^(2), Ingreso^(3)... ← 3×3 = 9 líneas de ingresos
200 900 300    ← costos fijos por decisión
0.95           ← alfa
```

---

# 6. TDA Matriz — El motor de cálculo

Todo el paquete se apoya en un **Tipo Abstracto de Datos (TDA)** llamado `Matriz` que encapsula las operaciones de álgebra lineal necesarias.

### Definición de la estructura

```c
typedef struct {
    int     filas;       /* número de filas     */
    int     columnas;    /* número de columnas  */
    double **datos;      /* arreglo 2D dinámico */
} Matriz;
```

`datos` es un **doble puntero**: un arreglo de punteros a filas, donde cada fila es un arreglo de `double`. Esto permite acceder a cualquier elemento con `m->datos[i][j]`.

### Funciones disponibles

```c
/* Crear y destruir */
Matriz* matriz_crear(int filas, int columnas);
void    matriz_destruir(Matriz *m);
Matriz* matriz_copiar(const Matriz *origen);
Matriz* matriz_identidad(int n);

/* Operaciones aritméticas */
Matriz* matriz_multiplicar(const Matriz *A, const Matriz *B);
Matriz* matriz_potencia(const Matriz *P, int n);
double* matriz_vector_multiplicar(const Matriz *A, const double *v);
Matriz* matriz_transponer(const Matriz *A);

/* Resolver Ax = b (eliminación gaussiana con pivoteo) */
int resolver_sistema_lineal(const Matriz *A, const double *b, double *x);
```

### Ejemplo: Crear, llenar y destruir una matriz

```c
#include "matrix.h"

int main() {
    /* Crear matriz 3×2 inicializada con ceros */
    Matriz *M = matriz_crear(3, 2);

    /* Llenar con valores */
    M->datos[0][0] = 0.7;  M->datos[0][1] = 0.3;
    M->datos[1][0] = 0.4;  M->datos[1][1] = 0.6;
    M->datos[2][0] = 1.0;  M->datos[2][1] = 0.0;

    /* Imprimir */
    matriz_imprimir(M, "Mi matriz");

    /* ¡Siempre liberar! */
    matriz_destruir(M);
    return 0;
}
```

### Eliminación gaussiana con pivoteo parcial

La función `resolver_sistema_lineal()` es el corazón de los métodos que requieren resolver sistemas de ecuaciones (estado estable, mejoramiento de políticas, absorción, primera pasada). Implementa:

1. **Pivoteo parcial**: en cada columna busca la fila con el mayor valor absoluto para minimizar errores numéricos.
2. **Eliminación hacia adelante**: convierte la matriz en triangular superior.
3. **Sustitución hacia atrás**: resuelve desde la última variable hasta la primera.

Si la matriz es singular (determinante ≈ 0), retorna 0 e informa del error.

```c
int resolver_sistema_lineal(const Matriz *A, const double *b, double *x);
// A: matriz n×n de coeficientes (NO se modifica — se trabaja sobre copia)
// b: vector de términos independientes (NO se modifica)
// x: vector solución (salida, debe estar preasignado con calloc)
// Retorna: 1 = éxito, 0 = matriz singular
```

---

# 7. Módulo 1: Teoría de Cadenas de Markov

> **Archivo:** `markov_cadenas.c` (332 líneas)
> **Archivo de cabecera:** `markov.h`

Cada función de este módulo toma una matriz de transición `P` y devuelve un resultado. Si el modelo actual es un PMD, se usa por defecto `P^(1)`.

---

## 7.1 Chapman-Kolmogorov (P^n)

**¿Qué hace?** Calcula las probabilidades de transición en **n pasos**.

**Fórmula:** `P^(n) = P^n` (la matriz P multiplicada por sí misma n veces)

**En el código:**

```c
Matriz* chapman_kolmogorov(const Matriz *P, int n) {
    return matriz_potencia(P, n);   /* usa exponenciación binaria */
}
```

La exponenciación binaria (`O(log n)` multiplicaciones en vez de `O(n)`) hace que calcular P^100 sea casi tan rápido como P^10.

**Uso en el menú:** Módulo 1 → Opción 1. Pide el número de pasos `n`.

**Salida esperada:** Una matriz `(m+1)×(m+1)` donde la celda `[i][j]` es `p_ij^(n)`.

---

## 7.2 Probabilidades Incondicionales

**¿Qué hace?** Dado un vector inicial `a`, calcula la probabilidad de estar en cada estado después de `n` pasos.

**Fórmula:** `P(X_n = j) = Σ_i  a_i · p_ij^(n)`  =  `a · P^n`

**En el código:**

```c
double* probabilidades_incondicionales(const double *a, const Matriz *P, int n) {
    // 1. Calcular P^n
    Matriz *Pn = matriz_potencia(P, n);
    // 2. Multiplicar a (fila) × Pn: prob[j] = Σ_i a[i] × Pn[i][j]
    int n_estados = P->filas;
    double *prob = (double*)calloc((size_t)n_estados, sizeof(double));
    for (int j = 0; j < n_estados; j++)
        for (int i = 0; i < n_estados; i++)
            prob[j] += a[i] * Pn->datos[i][j];
    matriz_destruir(Pn);
    return prob;
}
```

**Uso en el menú:** Módulo 1 → Opción 2. Pide `n`.

**Salida esperada:** Vector de tamaño `m+1` con la probabilidad de cada estado en el paso `n`. La suma debe dar 1.

---

## 7.3 Estado Estable (π)

**¿Qué hace?** Calcula la **distribución estacionaria** (o de equilibrio) de la cadena: a largo plazo, ¿en qué proporción del tiempo se visita cada estado?

**Fórmula:** `π = π·P`, con la restricción `Σ π_i = 1`.

**Método de resolución:**

1. Se construye el sistema `π·(P − I) = 0`.
2. Se reemplaza la última ecuación por `Σ π_i = 1` (para que tenga solución única).
3. Se resuelve con eliminación gaussiana.

**En el código:**

```c
double* estado_estable(const Matriz *P) {
    int n = P->filas;

    // 1. A = P^T - I
    Matriz *A = matriz_crear(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A->datos[i][j] = P->datos[j][i];   // transponer
    for (int i = 0; i < n; i++)
        A->datos[i][i] -= 1.0;                  // restar identidad

    // 2. Reemplazar última fila por 1's
    for (int j = 0; j < n; j++)
        A->datos[n-1][j] = 1.0;
    double *b = calloc(n, sizeof(double));
    b[n-1] = 1.0;

    // 3. Resolver y normalizar
    double *pi = calloc(n, sizeof(double));
    resolver_sistema_lineal(A, b, pi);
    // ...normalizar para que sume 1...
    return pi;
}
```

**Salida esperada:** Vector `π` con `π_0 + π_1 + ... + π_m = 1`.

---

## 7.4 Tiempos de Recurrencia (μ_ii)

**¿Qué hace?** Calcula el **número esperado de pasos** para regresar a un estado `i` partiendo de `i`.

**Fórmula:** `μ_ii = 1 / π_i`

Si `π_i` es muy pequeño, el tiempo de recurrencia es grande (el estado se visita poco). Si `π_i = 0`, `μ_ii = ∞`.

```c
double* tiempos_recurrencia(const double *pi, int n) {
    double *mu = calloc(n, sizeof(double));
    for (int i = 0; i < n; i++)
        mu[i] = (fabs(pi[i]) > 1e-15) ? 1.0 / pi[i] : INFINITY;
    return mu;
}
```

---

## 7.5 Tiempos de Primera Pasada (μ_ij)

**¿Qué hace?** Calcula el **número esperado de pasos** para llegar del estado `i` al estado `j` por **primera vez** (`i ≠ j`).

**Fórmula:** `μ_ij = 1 + Σ_{k ≠ j} p_ik · μ_kj`

Esto genera un sistema de `(m+1)` ecuaciones lineales para cada estado objetivo `j`. Se resuelve con Gauss.

**Salida esperada:** Matriz `M` de `(m+1)×(m+1)` donde `M[i][j] = μ_ij`.

---

## 7.6 Probabilidades de Absorción (f_ik)

**¿Qué hace?** Calcula la probabilidad de que la cadena sea **absorbida** en un estado `k`, partiendo de cada estado `i`.

**Requisito:** Hay que especificar cuáles estados son **absorbentes** (una vez que se entra en ellos, no se sale: `p_ii = 1`).

**Fórmula:** `f_ik = Σ_j p_ij · f_jk`, con condiciones de frontera `f_kk = 1` y `f_ik = 0` si `i` es absorbente y `i ≠ k`.

**Salida esperada:** Vector `f` de tamaño `m+1` con `f[i] = f_ik`.

---

# 8. Módulo 2: Procesos Markovianos de Decisión

> **Archivo:** `markov_decision.c` (1221 líneas)
> **Archivo de cabecera:** `markov.h`

Este módulo contiene los 5 algoritmos de resolución para PMD, más un sexto que ejecuta una prueba completa en lote.

---

### Estructura de una política

```c
typedef struct {
    int  num_estados;       /* m+1                         */
    int  num_decisiones;    /* K                           */
    int *decision;          /* decision[i] ∈ {0,...,K-1}   */
} Politica;
```

**Importante:** Internamente las decisiones se indexan desde **0** (C estándar), pero al usuario se muestran desde **1**. Por ejemplo, `decision[0] = 0` se imprime como `R = [1, ...]`.

---

### Estructura del modelo

```c
typedef struct {
    int     num_estados;        /* m+1                              */
    int     num_decisiones;     /* K                                */
    int     tiene_decisiones;   /* 0=cadena simple, 1=PMD           */
    double *prob_inicial;       /* vector a                         */
    Matriz **P_dec;             /* K matrices de transición         */
    Matriz  *C;                 /* matriz de costos (m+1)×K         */
    double  alfa;               /* factor de descuento              */
    int     es_maximizacion;    /* 1=MAX, 0=MIN                     */
    int     usa_ingresos;       /* 1=C auto-generado desde ingresos */
    Matriz **Ingreso;           /* K matrices de ingreso            */
    double  *costo_fijo;        /* costos fijos por decisión        */
} ModeloMarkov;
```

---

## 8.1 Enumeración Exhaustiva

**Tipo:** Exacto.  
**Complejidad:** `O(K^(m+1))` — evalúa **todas** las políticas posibles.

**¿Cómo funciona?**

1. Genera cada política posible (hay `K^(m+1)` en total).
2. Para cada política, construye la matriz combinada `P_R` y calcula su estado estable `π`.
3. Calcula el costo esperado: `E(C) = Σ_i C_{i, R(i)} · π_i`.
4. Se queda con la de menor costo.

**¿Cuándo usarlo?**
- Cuando `K^(m+1)` es manejable (ej. 2^4 = 16, 3^3 = 27).
- Como verificación de otros métodos.
- En problemas didácticos pequeños.

**Código clave:**

```c
// Generar políticas por índice en base K
for (long long idx = 0; idx < total; idx++) {
    long long temp = idx;
    for (int i = 0; i < n; i++) {
        actual.decision[i] = (int)(temp % K);
        temp /= K;
    }
    // evaluar política...
}
```

**En el menú:** Módulo 2 → Opción 1.

---

## 8.2 Mejoramiento de Políticas (sin descuento)

**Tipo:** Exacto.  
**Algoritmo:** Howard (1960) para **costo promedio** a largo plazo.  
**No usa** el factor de descuento `α` — es para problemas donde todas las semanas valen lo mismo.

**¿Cómo funciona?**

Es un algoritmo iterativo de dos pasos:

#### Paso 1 — Determinación del valor (Policy Evaluation)

Dada una política `R`, se resuelve el sistema:

```
g + V_i = C_{i, R(i)} + Σ_j P_{ij}(R(i)) · V_j      para i = 0..m
con V_m = 0
```

Donde:
- `g` es el **costo promedio** (ganancia) de la política.
- `V_i` son los **valores relativos** (sesgos) de cada estado.

Es un sistema de `(m+1)` ecuaciones con `(m+1)` incógnitas (`V_0, ..., V_{m-1}, g`).

**En el código:**

```c
// Construir sistema para la política actual
for (int i = 0; i < n; i++) {
    int k = R_actual->decision[i];
    for (int j = 0; j < n - 1; j++) {
        if (i == j)
            A->datos[i][j] = 1.0 - modelo->P_dec[k]->datos[i][j];
        else
            A->datos[i][j] = -modelo->P_dec[k]->datos[i][j];
    }
    A->datos[i][n-1] = 1.0;          // coeficiente de g
    b[i] = modelo->C->datos[i][k];   // costo inmediato
}
resolver_sistema_lineal(A, b, x);
// x = [V_0, V_1, ..., V_{m-1}, g]
```

#### Paso 2 — Mejoramiento (Policy Improvement)

Para cada estado `i`, se busca la decisión `k` que minimiza:

```
test = C_{ik} + Σ_j P_{ij}(k) · V_j  −  V_i
```

Si para algún estado se encuentra una decisión mejor, se actualiza la política y se repite desde el Paso 1.

**Criterio de parada:** Cuando `R_{n+1} = R_n` (la política no cambia), se encontró el óptimo.

**En el código:**

```c
for (int i = 0; i < n; i++) {
    double mejor_valor = INFINITY;
    int mejor_k = R_actual->decision[i];
    for (int k = 0; k < K; k++) {
        double suma = 0.0;
        for (int j = 0; j < n; j++)
            suma += modelo->P_dec[k]->datos[i][j] * V[j];
        double test = modelo->C->datos[i][k] + suma - V[i];
        if (test < mejor_valor - 1e-10) {
            mejor_valor = test;
            mejor_k = k;
        }
    }
    R_nueva->decision[i] = mejor_k;
}
```

**En el menú:** Módulo 2 → Opción 2. Pregunta si usar política inicial por defecto `[1,1,...,1]` o ingresar una manualmente.

---

## 8.3 Mejoramiento de Políticas (con descuento)

**Tipo:** Exacto.  
**Algoritmo:** Howard con factor de descuento `α`.  
**Usa:** `α < 1` (costo futuro vale menos que el presente).

**Diferencia clave con la versión sin descuento:**

En el Paso 1, el sistema a resolver es:

```
V_i = C_{i, R(i)} + α · Σ_j P_{ij}(R(i)) · V_j      para i = 0..m
```

No hay variable `g` (el costo está "incluido" en `V_i`). Son `(m+1)` ecuaciones con `(m+1)` incógnitas.

En el Paso 2, el test es:

```
test = C_{ik} + α · Σ_j P_{ij}(k) · V_j
```

(No se resta `V_i` porque `V_i` ya incorpora el valor de estar en `i`.)

**En el código:**

```c
// Paso 1: Sistema con descuento
for (int i = 0; i < n; i++) {
    int k = R_actual->decision[i];
    for (int j = 0; j < n; j++) {
        if (i == j)
            A->datos[i][j] = 1.0 - alfa * modelo->P_dec[k]->datos[i][j];
        else
            A->datos[i][j] = -alfa * modelo->P_dec[k]->datos[i][j];
    }
    b[i] = modelo->C->datos[i][k];
}

// Paso 2: Mejoramiento
for (int k = 0; k < K; k++) {
    double suma = 0.0;
    for (int j = 0; j < n; j++)
        suma += modelo->P_dec[k]->datos[i][j] * V[j];
    vals[k] = modelo->C->datos[i][k] + alfa * suma;
}
```

**En el menú:** Módulo 2 → Opción 3.

---

## 8.4 Aproximaciones Sucesivas

**Tipo:** Iterativo (aproximado).  
**Alias:** Iteración de Valor (Value Iteration).  
**Usa:** Factor de descuento `α`.

**¿Cómo funciona?**

1. **Inicialización:** `V_i^1 = min_k C_{ik}` (el costo mínimo inmediato en cada estado).

2. **Iteración n:** Para cada estado, calcular:

   ```
   V_i^n = min_k [ C_{ik} + α · Σ_j P_{ij}(k) · V_j^{n-1} ]
   ```

3. **Convergencia:** Cuando `|V_i^n − V_i^{n-1}| < ε` para todo `i`, o se alcanza el máximo de iteraciones `N`.

**Tabla de iteraciones:** El programa imprime una tabla con cada iteración, mostrando `V_i`, la diferencia máxima (`delta_max`), y la política implícita en esa iteración.

```
+------+-------------+-------------+---------------+---------------+---------------+
| iter |    V_0      |    V_1      |  delta_max    |   politica    |
+------+-------------+-------------+---------------+---------------+---------------+
|   1  |   0.000000  |  14.000000  |   (inicial)   |  [2,1]        |
|   2  |   6.993000  |  20.979000  |    6.9930e+00 |  [2,1]        |
...
```

**Parámetros configurables:**

| Parámetro | Típico | Significado |
|-----------|--------|-------------|
| `max_iter` | 100-1000 | Máximo de iteraciones |
| `epsilon` | 1e-6 | Tolerancia para convergencia |

**¿Cuándo usarlo?**
- Problemas grandes donde la enumeración es inviable.
- Cuando se necesita una aproximación rápida.
- La política suele converger antes que los valores V_i.

**En el menú:** Módulo 2 → Opción 4. Pide `N` y `ε`.

---

## 8.5 Programación Lineal

**Tipo:** Exacto.  
**Método:** Simplex con técnica de la **Gran M** (Big-M).  
**Implementación propia:** ~400 líneas de C, sin dependencias externas.

### Variables del PL

La formulación usa variables `Y_{ik}` = probabilidad estacionaria conjunta de estar en el estado `i` y tomar la decisión `k`.

Total de variables: `(m+1) × K`

### Formulación matemática

**Minimizar:**
```
Z = Σ_i Σ_k  C_{ik} · Y_{ik}
```

**Sujeto a:**

1. **Normalización:** `Σ_i Σ_k Y_{ik} = 1`

2. **Balance de flujo** (para cada estado j):
   ```
   Σ_k Y_{jk} − Σ_i Σ_k Y_{ik} · P_{ij}(k) = 0
   ```
   Esto asegura que el flujo que entra al estado j es igual al que sale.

3. **No negatividad:** `Y_{ik} ≥ 0`

### Resolución con Simplex

El programa construye la tabla simplex internamente:

```c
// Variables: Y_{ik}  +  variable artificial a  +  holguras s_j
int num_var = n_est * K;
int col_a = num_var;                          // artificial
int col_s_inicio = num_var + 1;               // holguras
int total_columnas = num_var + 1 + n_est;

// Fila 0: Σ Y_{ik} + a = 1
for (int i = 0; i < n_est; i++)
    for (int k = 0; k < K; k++)
        tabla->tabla[0][IDX(i,k)] = 1.0;
tabla->tabla[0][col_a] = 1.0;
tabla->rhs[0] = 1.0;

// Filas de balance: balance_j + s_j = 0
// ...coeficientes según P_{ij}(k)...
```

El costo de la variable artificial es `M = 1e10` (Gran M) para forzar que salga de la base.

### Extracción de política determinística

Al obtener los `Y_{ik}` óptimos, se extrae la política:

```
D_{ik} = Y_{ik} / Σ_k Y_{ik}        (probabilidad de elegir k en estado i)
R(i) = argmax_k D_{ik}              (política determinística)
```

**En el menú:** Módulo 2 → Opción 5.

---

## 8.6 Prueba Completa (Batch)

Ejecuta **4 métodos en secuencia** sin intervención del usuario:

1. Enumeración Exhaustiva
2. Mejoramiento sin descuento (desde política `[1,1,...,1]`)
3. Mejoramiento con descuento (desde política `[1,1,...,1]`)
4. Programación Lineal

**Objetivo:** Verificar que todos los métodos convergen a la misma política óptima (consistencia).

**En el menú:** Módulo 2 → Opción 6.

---

# 9. Manual de Usuario

## 9.1 Ejecución del programa

### Desde archivo (recomendado)

```bash
# Windows
.\markov.exe ejemplo_pmd.txt

# Linux/Mac
./markov ejemplo_pmd.txt
```

Al pasar un archivo como argumento, el modelo se carga automáticamente y se muestra el menú principal listo para usar.

### Desde consola interactiva

```bash
./markov
```

Luego usar la opción 1 del menú principal para ingresar todos los datos a mano.

---

## 9.2 Menú principal

```
+===============================================================+
|         PAQUETE ESTOCASTICO - CADENAS DE MARKOV Y PMD         |
+===============================================================+

  MENU PRINCIPAL
  +----+------------------------------------------------+
  | 1  | Cargar modelo desde consola                    |
  | 2  | Cargar modelo desde archivo                    |
  | 3  | Mostrar modelo cargado                         |
  | 4  | Modulo 1: Teoria basica de Cadenas de Markov   |
  | 5  | Modulo 2: Procesos Markovianos de Decision     |
  | 0  | Salir                                          |
  +----+------------------------------------------------+
  Opcion:
```

| Opción | ¿Qué hace? |
|--------|-----------|
| **1** | Ingresar el problema paso a paso por teclado |
| **2** | Cargar desde un archivo `.txt` |
| **3** | Ver un resumen de los datos cargados (matrices, costos, α) |
| **4** | Entrar al Módulo 1 (cadenas de Markov) |
| **5** | Entrar al Módulo 2 (PMD) |
| **0** | Salir |

---

## 9.3 Módulo 1 paso a paso

```
  MODULO 1 - TEORIA BASICA DE CADENAS DE MARKOV
  +----+------------------------------------------------+
  | 1  | Ecuaciones de Chapman-Kolmogorov (P^n)         |
  | 2  | Probabilidades Incondicionales                 |
  | 3  | Vector de Estado Estable                       |
  | 4  | Tiempos de Recurrencia                         |
  | 5  | Tiempos de Primera Pasada                      |
  | 6  | Probabilidades de Absorcion                    |
  | 0  | Volver al menu principal                       |
  +----+------------------------------------------------+
```

**Flujo típico para analizar una cadena:**

1. **Opción 3** → Obtener el estado estable `π`. ¿La cadena tiene distribución límite?
2. **Opción 4** → ¿Cada cuánto se visita cada estado?
3. **Opción 1** → Ver cómo evolucionan las probabilidades en `n` pasos.
4. **Opción 2** → Partiendo de una distribución inicial `a`, ¿dónde estaremos en `n` pasos?
5. **Opción 5** → ¿Cuánto tarda en promedio ir de un estado a otro?
6. **Opción 6** → Si hay estados absorbentes, ¿con qué probabilidad se cae en cada uno?

**Nota:** Si el modelo es un PMD, el Módulo 1 usa `P^(1)` (la primera decisión) como matriz de transición.

---

## 9.4 Módulo 2 paso a paso

```
  MODULO 2 - PROCESOS MARKOVIANOS DE DECISION
  +----+------------------------------------------------+
  | 1  | Enumeracion Exhaustiva de Politicas            |
  | 2  | Mejoramiento de Politicas (sin descuento)      |
  | 3  | Mejoramiento de Politicas con Descuento        |
  | 4  | Metodo de Aproximaciones Sucesivas             |
  | 5  | Solucion por Programacion Lineal               |
  | 6  | PRUEBA COMPLETA (ejecuta 1,2,3,5 en secuencia) |
  | 0  | Volver al menu principal                       |
  +----+------------------------------------------------+
```

**Flujo recomendado para resolver un PMD:**

1. **Opción 6** (Prueba completa) → Ejecuta todo de una vez.Compara resultados.
2. Si los resultados coinciden → **confianza total** en la política óptima.
3. Si hay discrepancias → Investigar: ¿matriz singular?, ¿problema de convergencia?

**Opciones 2 y 3** preguntan:

```
  Elegir politica inicial para Mejoramiento de Politicas:
    1. Ingresar politicas manualmente
    2. Usar politica por defecto (decision 1 para todo)
```

Para el problema del póker, por ejemplo, si quieres empezar desde `(Ofrecer, Ofrecer)`, eliges la opción 2 (coincide porque la decisión 1 es Ofrecer).

**Opción 4** pide:

```
  Maximo de iteraciones N: 1000
  Tolerancia epsilon: 1e-6
```

Valores típicos: `N = 1000`, `ε = 0.000001`. Si `α` está cerca de 1, la convergencia es más lenta; aumentar `N`.

---

## 9.5 Carga por archivo vs consola

### Ventajas del archivo

- **Reproducible:** El mismo archivo siempre da el mismo resultado.
- **Compartible:** Se puede enviar el `.txt` a otra persona.
- **Rápido:** No hay que re-ingresar 20+ números.
- **Versionable:** Se puede poner en git.

### Ventajas de la consola

- **Exploratoria:** Para probar ideas rápidamente.
- **Didáctica:** Para entender qué dato se pide en cada paso.

### Formato del archivo (resumen visual)

```
1        ← ¿PMD? (0=no, 1=sí)
m        ← estados 0..m
K        ← decisiones (solo si PMD)
modo     ← 0=costos directos, 1=ingresos+maximizar, 2=ingresos+minimizar
a0..am   ← probabilidad inicial
[P^(1)]  ← matriz (m+1)×(m+1) ...
[P^(K)]  ← ... K matrices
[C]      ← matriz (m+1)×K de costos  (si modo=0)
[Ingreso^(1..K)] [costos_fijos]      (si modo=1/2)
alfa     ← factor de descuento
```

---

# 10. Ejemplos resueltos

## 10.1 Ejemplo 1: PMD básico

**Archivo:** `ejemplo_pmd.txt`

```
1
2
2
0
0.5 0.3 0.2
0.7 0.2 0.1
0.4 0.4 0.2
0.3 0.3 0.4
0.6 0.2 0.2
0.2 0.5 0.3
0.8 0.1 0.1
10.0 5.0
8.0 12.0
6.0 9.0
0.9
```

**Interpretación:**
- 3 estados, 2 decisiones, costos directos, α = 0.9.
- Vector inicial: `a = (0.5, 0.3, 0.2)`.
- P^(1): la decisión 1 favorece transiciones al estado 0.
- P^(2): la decisión 2 favorece transiciones al estado 2.
- Costos: decisión 1 es más cara en estado 0 ($10 vs $5), pero más barata en estado 1 ($8 vs $12) y estado 2 ($6 vs $9).

---

## 10.2 Ejemplo 2: Caso de publicidad

**Archivo:** `test_publicidad.txt`

**Problema:** Una empresa puede elegir entre 3 niveles de publicidad (decisión 1, 2, 3) en 3 estados del mercado (0, 1, 2). Se busca **maximizar** la utilidad esperada.

**Particularidades de este ejemplo:**

- Usa `modo_costos = 1` (ingresos + maximización).
- Tiene 3 matrices de ingreso (una por decisión) y costos fijos por decisión.
- El paquete auto-genera `C_ik` a partir de los ingresos.
- α = 0.95.

**Resultado:** Los 4 métodos convergen a la misma política óptima `[1, 3, 3]` con utilidad esperada de **$282.00**.

---

## 10.3 Ejemplo 3: Póker de los sábados

**Archivo:** `problema2.txt`  
**Solución:** `SOLUCION_PROBLEMA2.md`

```
1
1
2
0
0.5 0.5
0.875 0.125
0.875 0.125
0.125 0.875
0.125 0.875
14 0
14 75
0.999
```

**Problema resumido:**

- 2 estados (buen humor=0, mal humor=1), 2 decisiones (ofrecer=1, no ofrecer=2).
- P^(1): al ofrecer, 7/8 de probabilidad de buen humor la próxima semana.
- P^(2): al no ofrecer, solo 1/8 de probabilidad de buen humor.
- Costos: ofrecer cuesta $14; no ofrecer en mal humor cuesta $75 (lo molestan).
- α = 0.999 ≈ 1 (costo promedio, sin descuento real).

**Política óptima:** `[2, 1]` = No ofrecer si están de buen humor, ofrecer si están de mal humor.  
**Costo esperado:** $7.00 / semana.

---

# 11. Referencia rápida de API

### Estructuras principales

```c
// Matriz dinámica (matrix.h)
typedef struct { int filas, columnas; double **datos; } Matriz;

// Política determinística (markov.h)
typedef struct { int num_estados, num_decisiones; int *decision; } Politica;

// Modelo completo (markov.h)
typedef struct {
    int num_estados, num_decisiones, tiene_decisiones;
    double *prob_inicial;
    Matriz *P;           // cadena simple
    Matriz **P_dec;      // PMD: K matrices
    Matriz *C;           // costos (m+1)×K
    double alfa;
    int es_maximizacion, usa_ingresos;
    Matriz **Ingreso;
    double *costo_fijo;
} ModeloMarkov;
```

### Funciones públicas (markov.h)

```c
// --- Carga de datos ---
ModeloMarkov* modelo_crear(void);
void          modelo_destruir(ModeloMarkov *m);
ModeloMarkov* modelo_leer_consola(void);
ModeloMarkov* modelo_leer_archivo(const char *archivo);
void          modelo_imprimir(const ModeloMarkov *m);
void          modelo_calcular_costos_ingresos(ModeloMarkov *m);

// --- Módulo 1: Cadenas ---
Matriz* chapman_kolmogorov(const Matriz *P, int n);
double* probabilidades_incondicionales(const double *a, const Matriz *P, int n);
double* estado_estable(const Matriz *P);
double* tiempos_recurrencia(const double *pi, int n);
Matriz* tiempos_primera_pasada(const Matriz *P);
double* probabilidades_absorcion(const Matriz *P, const int *absorbentes, int k);

// --- Módulo 2: PMD ---
Politica* politica_crear(int num_estados, int num_decisiones);
void      politica_destruir(Politica *p);
void      politica_imprimir(const Politica *p);
Matriz*   politica_matriz_transicion(const ModeloMarkov *m, const Politica *p);
void      enumeracion_exhaustiva(ModeloMarkov *m);
Politica* mejoramiento_politicas(const ModeloMarkov *m, const Politica *ini);
Politica* mejoramiento_politicas_descuento(const ModeloMarkov *m, const Politica *ini);
void      aproximaciones_sucesivas(const ModeloMarkov *m, int N, double eps);
void      programacion_lineal(const ModeloMarkov *m);
void      prueba_completa_pmd(ModeloMarkov *m);
```

### Funciones de matriz (matrix.h)

```c
Matriz* matriz_crear(int filas, int cols);
void    matriz_destruir(Matriz *m);
Matriz* matriz_copiar(const Matriz *orig);
Matriz* matriz_identidad(int n);
Matriz* matriz_multiplicar(const Matriz *A, const Matriz *B);
Matriz* matriz_potencia(const Matriz *P, int n);
double* matriz_vector_multiplicar(const Matriz *A, const double *v);
Matriz* matriz_transponer(const Matriz *A);
int     resolver_sistema_lineal(const Matriz *A, const double *b, double *x);
void    matriz_imprimir(const Matriz *m, const char *titulo);
void    vector_imprimir(const double *v, int n, const char *titulo);
```

---

# 12. Limitaciones y recomendaciones

### Limitaciones

| Aspecto | Límite | Observación |
|---------|--------|-------------|
| Estados (`m+1`) | 500 | Definido por `leer_entero(..., 1, 500)` |
| Decisiones (`K`) | 100 | Definido por `leer_entero(..., 1, 100)` |
| Políticas (exhaustiva) | ~10⁶ | Tiempos prácticos; más allá es muy lento |
| Precisión numérica | `double` (64 bits) | ~15 dígitos decimales |
| Simplex | 10000 iteraciones | Normalmente converge en mucho menos |
| Memoria | Dinámica | Crece con `(m+1)² × K` |

### Recomendaciones

1. **Usar `α = 0.999` para problemas de costo promedio.** El paquete no tiene un modo específico para costo promedio sin descuento, pero `α ≈ 1` da resultados equivalentes. El método de mejoramiento sin descuento (opción 2 del Módulo 2) **no usa α** y es el apropiado para estos casos.

2. **Verificar con prueba completa.** Antes de confiar en un resultado, ejecutar la opción 6 del Módulo 2 para confirmar que los 4 métodos coinciden.

3. **Normalizar probabilidades.** Aunque el paquete normaliza automáticamente desviaciones menores a 0.001, es buena práctica que cada fila de P sume exactamente 1.

4. **Preferir archivos sobre consola.** Para problemas con más de 2 estados, es mucho más práctico y menos propenso a errores.

5. **Cuidado con la indexación.** Las decisiones se muestran al usuario como 1..K, pero internamente son 0..K-1. La política `[2, 1]` significa: estado 0 → decisión 2, estado 1 → decisión 1.

6. **No modificar los archivos del paquete.** Si se necesita extender funcionalidad, crear nuevos archivos que incluyan `markov.h`.

---

<div align="center">

**— Fin de la documentación —**

*Paquete Estocástico v1.0 · C99 · Sin dependencias externas*

</div>
