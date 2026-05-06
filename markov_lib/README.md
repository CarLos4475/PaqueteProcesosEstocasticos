# Cadenas de Markov y Procesos Markovianos de Decisión (PMD)

## Paquete de Software Estocástico en C

---

## Índice

1. [Introducción](#1-introducción)
2. [Fundamentos Matemáticos](#2-fundamentos-matemáticos)
   - [2.1 Cadenas de Markov](#21-cadenas-de-markov)
   - [2.2 Procesos Markovianos de Decisión](#22-procesos-markovianos-de-decisión)
3. [Estructura del Paquete](#3-estructura-del-paquete)
4. [Arquitectura y Lógica de los Algoritmos](#4-arquitectura-y-lógica-de-los-algoritmos)
   - [4.1 TDA Matriz y Eliminación Gaussiana](#41-tda-matriz-y-eliminación-gaussiana)
   - [4.2 Algoritmo de Lectura de Datos](#42-algoritmo-de-lectura-de-datos)
   - [4.3 Módulo 1: Teoría Básica de Cadenas](#43-módulo-1-teoría-básica-de-cadenas)
   - [4.4 Módulo 2: Resolución de PMD](#44-módulo-2-resolución-de-pmd)
5. [Manual de Usuario](#5-manual-de-usuario)
   - [5.1 Requisitos](#51-requisitos)
   - [5.2 Compilación](#52-compilación)
   - [5.3 Ejecución](#53-ejecución)
   - [5.4 Formato de Archivos de Entrada](#54-formato-de-archivos-de-entrada)
   - [5.5 Guía Paso a Paso](#55-guía-paso-a-paso)
6. [Ejemplos](#6-ejemplos)
7. [Referencias](#7-referencias)

---

## 1. Introducción

Este paquete implementa de forma exacta la teoría, métodos y algoritmos sobre **Cadenas de Markov** y **Procesos Markovianos de Decisión (PMD)**. Está escrito completamente en **C estándar (C99)** con módulos independientes, estructuras de datos dinámicas y una interfaz interactiva de consola.

Cada algoritmo sigue estrictamente las definiciones matemáticas del documento de especificación, sin recurrir a fórmulas inventadas ni métodos externos. Los métodos iterativos imprimen obligatoriamente el estado del sistema en cada iteración para fines de depuración y análisis.

---

## 2. Fundamentos Matemáticos

### 2.1 Cadenas de Markov

Una **cadena de Markov** es un proceso estocástico $\{X_n\}_{n \geq 0}$ con espacio de estados $E = \{0, 1, \ldots, m\}$ que satisface la **propiedad de Markov**:

$$P(X_{n+1} = j \mid X_n = i, X_{n-1} = i_{n-1}, \ldots, X_0 = i_0) = P(X_{n+1} = j \mid X_n = i) = p_{ij}$$

#### Variables Fundamentales

| Variable | Significado |
|----------|-------------|
| $E = \{0, 1, \ldots, m\}$ | Espacio de estados |
| $P = [p_{ij}]$ | Matriz de transición $(m+1) \times (m+1)$, estocástica por filas |
| $a = (a_0, a_1, \ldots, a_m)$ | Vector de probabilidad inicial |
| $\pi = (\pi_0, \pi_1, \ldots, \pi_m)$ | Distribución estacionaria |
| $\mu_{ij}$ | Tiempo esperado de primera pasada de $i$ a $j$ |
| $f_{ik}$ | Probabilidad de absorción desde $i$ al estado absorbente $k$ |

#### Fórmulas Implementadas

**1. Ecuaciones de Chapman-Kolmogorov (transiciones en $n$ pasos):**

$$p_{ij}^{(n)} = \sum_{k=0}^{m} p_{ik}^{(l)} \cdot p_{kj}^{(n-l)}$$

Equivalentemente: $P^{(n)} = P^n$ (potencia $n$-ésima de la matriz de transición).

**2. Probabilidades Incondicionales en el paso $n$:**

$$P(X_n = j) = \sum_{i=0}^{m} a_i \cdot p_{ij}^{(n)} = (a \cdot P^n)_j$$

**3. Vector de Estado Estable $\pi$:**

$$\pi = \pi \cdot P, \quad \sum_{i=0}^{m} \pi_i = 1$$

Se resuelve el sistema lineal $(P^T - I)\pi^T = 0$ sustituyendo la última ecuación por la restricción de normalización $\sum \pi_i = 1$.

**4. Tiempos de Recurrencia:**

$$\mu_{ii} = \frac{1}{\pi_i}$$

**5. Tiempos de Primera Pasada:**

$$\mu_{ij} = 1 + \sum_{k \neq j} p_{ik} \cdot \mu_{kj}$$

Para cada estado objetivo $j$ se resuelve un sistema lineal de tamaño $(m) \times (m)$.

**6. Probabilidades de Absorción:**

$$f_{ik} = \sum_{j=0}^{m} p_{ij} \cdot f_{jk}$$

Con condiciones de frontera: $f_{kk} = 1$ y $f_{ik} = 0$ si $i$ es absorbente con $i \neq k$.

---

### 2.2 Procesos Markovianos de Decisión

Un **PMD** extiende la cadena de Markov añadiendo un conjunto de decisiones $K = \{1, 2, \ldots, K\}$. En cada estado $i$, se elige una decisión $k$ que determina:

- Una **matriz de transición** $P_{ij}(k)$ específica para esa decisión
- Un **costo esperado** $C_{ik}$ por estar en el estado $i$ y tomar la decisión $k$

Una **política determinística** $R$ asigna a cada estado $i$ una decisión $R(i) \in K$. El objetivo es encontrar la política que minimiza el costo esperado a largo plazo.

#### Variables Adicionales

| Variable | Significado |
|----------|-------------|
| $K = \{1, 2, \ldots, K\}$ | Conjunto de decisiones posibles |
| $P_{ij}(k)$ | Probabilidad de transición de $i$ a $j$ bajo decisión $k$ |
| $C_{ik}$ | Costo esperado en estado $i$ con decisión $k$ |
| $\alpha$ | Factor de descuento, $\alpha = (1 + i)^{-1}$ |
| $V_i(R)$ | Función de valor en estado $i$ bajo política $R$ |
| $g(R)$ | Ganancia promedio a largo plazo bajo política $R$ |
| $Y_{ik}$ | Probabilidad estacionaria conjunta de estado $i$ y decisión $k$ |

#### Los 5 Algoritmos de Resolución

---

##### Algoritmo 1: Enumeración Exhaustiva de Políticas

**Lógica:** Generar todas las $K^{m+1}$ políticas posibles, evaluar cada una y seleccionar la óptima.

**Pasos:**
1. Para cada política $R$, construir $P_R[i][j] = P_{ij}(R(i))$
2. Calcular el estado estacionario $\pi$ de $P_R$
3. Calcular costo esperado: $E(C) = \sum_{i=0}^{m} C_{i, R(i)} \cdot \pi_i$
4. Seleccionar la política con mínimo $E(C)$

**Complejidad:** $O(K^{m+1} \cdot (m+1)^3)$ — factible solo para $m$ y $K$ pequeños.

---

##### Algoritmo 2: Mejoramiento de Políticas (Sin Descuento) — Howard

**Lógica:** Alternar entre evaluar una política y mejorarla, garantizando convergencia al óptimo en un número finito de pasos.

$$
\boxed{\text{Paso 1 (Valor): } g(R_n) = C_{ik} + \sum_{j=0}^{m} P_{ij}(k) \cdot V_j(R_n) - V_i(R_n), \quad V_m(R_n) = 0}
$$

Se resuelve un sistema lineal de $(m+1) \times (m+1)$ para las incógnitas $[V_0, V_1, \ldots, V_{m-1}, g]$.

$$
\boxed{\text{Paso 2 (Mejora): } \argmin_k \left[ C_{ik} + \sum_{j=0}^{m} P_{ij}(k) \cdot V_j(R_n) - V_i(R_n) \right]}
$$

**Criterio de parada:** $R_{n+1} = R_n$ (política óptima alcanzada).

---

##### Algoritmo 3: Mejoramiento de Políticas con Descuento

**Lógica:** Misma estructura que el algoritmo 2, pero incorporando el factor de descuento $\alpha$.

$$
\boxed{\text{Paso 1 (Valor): } V_i(R_n) = C_{ik} + \alpha \cdot \sum_{j=0}^{m} P_{ij}(k) \cdot V_j(R_n)}
$$

Sistema lineal de $(m+1) \times (m+1)$: $V_i - \alpha \sum_j P_{ij}(k) V_j = C_{ik}$.

$$
\boxed{\text{Paso 2 (Mejora): } \argmin_k \left[ C_{ik} + \alpha \cdot \sum_{j=0}^{m} P_{ij}(k) \cdot V_j(R_n) \right]}
$$

---

##### Algoritmo 4: Método de Aproximaciones Sucesivas

**Lógica:** Iteración de valor. Comienza con una estimación inicial y la refina iterativamente usando programación dinámica.

$$
\boxed{V_i^1 = \min_k \; C_{ik}}
$$

$$
\boxed{V_i^n = \min_k \left( C_{ik} + \alpha \cdot \sum_{j=0}^{m} P_{ij}(k) \cdot V_j^{n-1} \right)}
$$

**Criterio de parada:** $|V_i^n - V_i^{n-1}| < \epsilon$ para todo $i$, o $n > N_{\max}$.

**Nota:** Este método imprime obligatoriamente $V_i^n$ en cada iteración.

---

##### Algoritmo 5: Programación Lineal

**Lógica:** Reformular el PMD como un problema de programación lineal con variables $Y_{ik}$ (probabilidad estacionaria conjunta) y resolverlo mediante el método Simplex con la **técnica de la Gran M**.

**Función objetivo** (minimizar):

$$Z = \sum_{i=0}^{m} \sum_{k=1}^{K} C_{ik} \cdot Y_{ik}$$

**Restricciones:**

1. Normalización: $\displaystyle\sum_{i=0}^{m} \sum_{k=1}^{K} Y_{ik} = 1$

2. Balance de flujo ($\forall j \in \{0, \ldots, m\}$):
   $$\sum_{k=1}^{K} Y_{jk} - \sum_{i=0}^{m} \sum_{k=1}^{K} Y_{ik} \cdot P_{ij}(k) = 0$$

3. No negatividad: $Y_{ik} \geq 0$

**Transformación a política determinística:**
$$D_{ik} = \frac{Y_{ik}}{\sum_{k=1}^{K} Y_{ik}}$$

La decisión para el estado $i$ es $\argmax_k D_{ik}$.

**Método de la Gran M:** La restricción de igualdad $\sum Y = 1$ se maneja con una variable artificial $a \geq 0$ y un costo penalizador $M$ (muy grande) en la función objetivo: $\min \sum C_{ik}Y_{ik} + M \cdot a$. El Simplex fuerza $a \to 0$, garantizando que $\sum Y = 1$.

---

## 3. Estructura del Paquete

```
markov_lib/
├── matrix.h              # TDA Matriz: operaciones + eliminación gaussiana
├── matrix.c              # Implementación del TDA Matriz
├── markov.h              # Estructuras de datos y prototipos de funciones
├── markov_datos.c        # Algoritmo de Lectura de Datos (núcleo del paquete)
├── markov_cadenas.c      # Módulo 1: Teoría Básica de Cadenas de Markov
├── markov_decision.c     # Módulo 2: Algoritmos de Resolución de PMD
├── main.c                # Programa principal con menú interactivo
├── Makefile              # Script de compilación (Linux/Mac/MinGW)
├── compilar.bat          # Script de compilación (Windows)
├── ejemplo_pmd.txt       # Archivo de ejemplo para PMD
└── README.md             # Este documento
```

### Dependencias entre módulos

```
matrix.h ────────────────► markov.h
                               │
          ┌────────────────────┼────────────────────┐
          ▼                    ▼                    ▼
   markov_datos.c       markov_cadenas.c      markov_decision.c
          │                    │                    │
          └────────────────────┼────────────────────┘
                               ▼
                           main.c
```

- `matrix.h/c`: Capa más baja. Todas las operaciones matriciales y el solver de sistemas lineales.
- `markov.h`: Define `ModeloMarkov`, `Politica`, `Matriz` y declara todas las funciones públicas.
- `markov_datos.c`: Depende de `matrix.h` para crear matrices. Lee datos de consola o archivo.
- `markov_cadenas.c`: Usa `matrix.h` para operaciones matriciales y `markov.h` para el modelo.
- `markov_decision.c`: Usa `matrix.h` (sistemas lineales) + `markov_cadenas.c` (estado estable). Contiene también el Simplex.
- `main.c`: Orquesta todos los módulos mediante el menú interactivo.

---

## 4. Arquitectura y Lógica de los Algoritmos

### 4.1 TDA Matriz y Eliminación Gaussiana

**Estructura `Matriz`:**
```c
typedef struct {
    int filas;
    int columnas;
    double **datos;  // arreglo dinámico de punteros a filas
} Matriz;
```

**Operaciones implementadas:**
- `matriz_crear(filas, cols)` — asigna memoria dinámica con `calloc`
- `matriz_destruir(m)` — libera toda la memoria
- `matriz_multiplicar(A, B)` — producto matricial $C_{ij} = \sum_k A_{ik} B_{kj}$
- `matriz_potencia(P, n)` — $P^n$ usando exponenciación binaria $O(\log n)$
- `matriz_vector_multiplicar(A, v)` — $u_i = \sum_j A_{ij} v_j$
- `resolver_sistema_lineal(A, b, x)` — eliminación gaussiana con pivoteo parcial

**Eliminación Gaussiana con pivoteo parcial** es el corazón del paquete. Se usa en:
- Cálculo del estado estacionario $\pi$
- Tiempos de primera pasada $\mu_{ij}$
- Probabilidades de absorción $f_{ik}$
- Determinación del valor en el mejoramiento de políticas
- Determinación del valor en el mejoramiento con descuento

**Lógica del pivoteo:** En cada columna, se busca la fila con el máximo valor absoluto y se intercambia. Esto minimiza errores numéricos. Si el pivote es menor que $10^{-12}$, el sistema es singular.

---

### 4.2 Algoritmo de Lectura de Datos

Es el módulo más importante. Soporta dos modos:

**Modo consola (`modelo_leer_consola`):**
1. Pregunta si es cadena simple o PMD
2. Lee $m$ (último estado → $m+1$ estados)
3. Si es cadena simple: vector $a$, matriz $P$
4. Si es PMD: vector $a$, $K$ matrices $P^{(k)}$, matriz $C$, factor $\alpha$
5. Valida que las matrices estocásticas sumen 1 por fila (normaliza si no)
6. Pregunta si ingresa tasa de interés $i$ o factor $\alpha$ directamente

**Modo archivo (`modelo_leer_archivo`):**
- Lee el mismo formato desde un archivo de texto
- Útil para modelos grandes o para automatizar pruebas

---

### 4.3 Módulo 1: Teoría Básica de Cadenas

**Chapman-Kolmogorov (`chapman_kolmogorov`):**
- Simplemente llama a `matriz_potencia(P, n)` para obtener $P^n$
- La exponenciación binaria evita $n-1$ multiplicaciones

**Estado Estable (`estado_estable`):**
- Construye el sistema $(P^T - I)\pi^T = 0$ con sustitución de última fila
- Resuelve con eliminación gaussiana
- La restricción $\sum \pi_i = 1$ reemplaza la última ecuación

**Tiempos de Primera Pasada (`tiempos_primera_pasada`):**
- Para cada estado objetivo $j$, construye un sistema de $m$ ecuaciones
- Las incógnitas son $\mu_{ij}$ para cada $i \neq j$
- La ecuación para cada $i$ es: $\mu_{ij} - \sum_{k \neq j} p_{ik} \mu_{kj} = 1$

**Probabilidades de Absorción (`probabilidades_absorcion`):**
- Identifica estados transientes (no absorbentes) excluyendo $k$
- Construye sistema: $f_{ik} - \sum_{j \in T} p_{ij} f_{jk} = p_{ik}$
- Los estados absorbentes tienen $f$ fijo (1 para $k$, 0 para otros)

---

### 4.4 Módulo 2: Resolución de PMD

**Enumeración Exhaustiva:**
- Itera desde $0$ hasta $K^{m+1}-1$, interpretando cada número en base $K$ como una política
- Para cada política, ejecuta el pipeline: construir $P_R$ → calcular $\pi$ → calcular $E(C)$
- Mantiene registro de la mejor política encontrada

**Mejoramiento de Políticas (sin y con descuento):**
- Bucle infinito que alterna Paso 1 (resolver sistema) y Paso 2 (elegir mejores decisiones)
- Sin descuento: $m+1$ ecuaciones con $V_m = 0$ fijo, incógnitas $V_0\ldots V_{m-1}, g$
- Con descuento: $m+1$ ecuaciones, incógnitas $V_0\ldots V_m$
- Mejora greedy: para cada estado, evaluar las $K$ decisiones y elegir la mínima
- Convergencia en número finito de iteraciones

**Aproximaciones Sucesivas:**
- Inicializa $V_i^1 = \min_k C_{ik}$
- Itera: $V_i^n = \min_k (C_{ik} + \alpha \sum_j P_{ij}(k) V_j^{n-1})$
- Verifica tolerancia $\epsilon$ componente a componente
- Imprime $V^n$ en cada iteración obligatoriamente

**Programación Lineal (Simplex con Gran M):**
- Variables $Y_{ik}$ ($n \times K$ variables), variables de holgura, variable artificial
- Restricciones: 1 igualdad (normalización) + $m+1$ desigualdades (balance)
- La igualdad se maneja con variable artificial $a$ y costo $M$ en el objetivo
- Inicialización: $a=1$ en primera fila, holguras en las demás
- El Simplex pivotea para minimizar el costo, forzando $a \to 0$
- Al terminar, transforma $Y_{ik}$ a política determinística vía $D_{ik}$

---

## 5. Manual de Usuario

### 5.1 Requisitos

- **Compilador GCC** (MinGW en Windows, GCC en Linux/Mac)
- **Sistema operativo:** Windows, Linux o macOS
- **Sin dependencias externas:** solo biblioteca estándar de C y `math.h`

### 5.2 Compilación

#### Windows (con MinGW)

**Opción A — Usar el script:**
```
doble clic en compilar.bat
```

**Opción B — Línea de comandos:**
```bash
gcc -Wall -O2 -std=c99 -I. -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm
```

#### Linux / Mac

```bash
make          # compila
make clean    # limpia archivos objeto
./markov      # ejecuta
```

O manualmente:
```bash
gcc -Wall -O2 -std=c99 -I. -o markov matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm
```

### 5.3 Ejecución

```bash
# Iniciar el programa (menú interactivo)
./markov

# Cargar directamente un archivo de datos
./markov ejemplo_pmd.txt

# Cargar otro archivo
./markov mi_modelo.txt
```

Al iniciar, verás el menú principal:

```
=============================================================
  CADENAS DE MARKOV Y PROCESOS MARKOVIANOS DE DECISION
  Paquete de Software Estocastico
=============================================================

  MENU PRINCIPAL:
    1. Cargar modelo desde consola
    2. Cargar modelo desde archivo
    3. Mostrar modelo cargado
    4. Modulo 1: Teoria basica de Cadenas de Markov
    5. Modulo 2: Procesos Markovianos de Decision (PMD)
    0. Salir
  Opcion:
```

### 5.4 Formato de Archivos de Entrada

#### Cadena de Markov Simple

```
0              ← tipo (0 = cadena simple)
2              ← m (último índice de estado → 3 estados: 0,1,2)
0.4 0.3 0.3    ← vector de probabilidad inicial a
0.5 0.3 0.2    ← matriz P, fila 0
0.1 0.7 0.2    ← matriz P, fila 1
0.3 0.3 0.4    ← matriz P, fila 2
```

#### Proceso Markoviano de Decisión (PMD)

```
1              ← tipo (1 = PMD)
2              ← m (3 estados: 0,1,2)
2              ← K (2 decisiones: 1,2)
0.5 0.3 0.2    ← vector a
0.7 0.2 0.1    ← P^(1) fila 0
0.4 0.4 0.2    ← P^(1) fila 1
0.3 0.3 0.4    ← P^(1) fila 2
0.6 0.2 0.2    ← P^(2) fila 0
0.2 0.5 0.3    ← P^(2) fila 1
0.8 0.1 0.1    ← P^(2) fila 2
10.0 5.0       ← C, fila 0 (costo estado 0 con decisión 1 y 2)
8.0 12.0       ← C, fila 1
6.0 9.0        ← C, fila 2
0.9            ← factor de descuento alfa
```

**Reglas del formato:**
- Cada número separado por espacios o saltos de línea
- Las matrices se leen por filas, en orden
- Las matrices de transición deben ser estocásticas (cada fila suma 1)
- Los costos pueden ser cualquier número real
- $\alpha$ debe estar en $(0, 1]$

### 5.5 Guía Paso a Paso

#### Para alguien que nunca ha usado C:

1. **Instala MinGW** (solo Windows):
   - Descarga de https://sourceforge.net/projects/mingw-w64/
   - O usa TDM-GCC: https://jmeubank.github.io/tdm-gcc/

2. **Abre una terminal** (PowerShell, CMD, o Git Bash)

3. **Navega a la carpeta del proyecto:**
   ```bash
   cd C:\Users\PC\Documents\PROYECTOS\Estocasticos\markov_lib
   ```

4. **Compila el programa** (elige una opción):
   - Windows: `gcc -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm`
   - O simplemente: `.\compilar.bat`

5. **Ejecuta el programa:**
   ```bash
   .\markov.exe
   ```

6. **Flujo de trabajo típico:**

   **Paso 1 — Cargar datos:**
   - Opción 1 del menú: ingresar datos manualmente por consola
   - O opción 2: cargar desde archivo (escribe el nombre, ej. `ejemplo_pmd.txt`)

   **Paso 2 — Verificar datos:**
   - Opción 3: muestra un resumen completo del modelo cargado

   **Paso 3 — Ejecutar Módulo 1 (cadenas simples) o Módulo 2 (PMD):**
   - Opción 4: accede a las 6 funciones de teoría de cadenas
   - Opción 5: accede a los 5 algoritmos de PMD

   **Ejemplo con Módulo 1:**
   ```
   Opción 4 → Opción 1 (Chapman-Kolmogorov) → Ingresa n=5 → Muestra P^5
   Opción 4 → Opción 3 (Estado Estable) → Muestra π
   Opción 4 → Opción 5 (Primera Pasada) → Muestra matriz μ
   ```

   **Ejemplo con Módulo 2:**
   ```
   Opción 5 → Opción 1 (Enumeración Exhaustiva) → Evalúa todas las políticas
   Opción 5 → Opción 3 (Mejoramiento con Descuento) → Elige política inicial → Converge
   Opción 5 → Opción 4 (Aproximaciones Sucesivas) → Ingresa N=100, ε=0.0001
   Opción 5 → Opción 5 (Programación Lineal) → Resuelve con Simplex
   ```

   **Paso 4 — Salir:**
   - Opción 0 en cualquier menú para retroceder o salir

---

## 6. Ejemplos

### Ejemplo 1: Cadena de Markov Simple

**Archivo `cadena.txt`:**
```
0
2
0.4 0.3 0.3
0.5 0.3 0.2
0.1 0.7 0.2
0.3 0.3 0.4
```

**Ejecución:**
```bash
./markov cadena.txt
# Menú: 3 (mostrar modelo) → 4 (Módulo 1) → 3 (estado estable)
```

**Resultado esperado:** $\pi = [0.2692, 0.4615, 0.2692]$

---

### Ejemplo 2: PMD con 3 estados y 2 decisiones

Usando el archivo incluido `ejemplo_pmd.txt`:

```
Tipo: PMD
Estados: 3 (0, 1, 2)
Decisiones: 2
α = 0.9

P^(1):
  0.7  0.2  0.1
  0.4  0.4  0.2
  0.3  0.3  0.4

P^(2):
  0.6  0.2  0.2
  0.2  0.5  0.3
  0.8  0.1  0.1

Costos C:
  [10,  5]
  [ 8, 12]
  [ 6,  9]
```

**Algoritmo de Mejoramiento con Descuento:**
- **Política inicial:** $R = [1, 1, 1]$ (decisión 1 para todos)
- **Iteración 1:** Valores $V = [89.18, 86.10, 83.32]$, mejora a $R = [2, 1, 1]$
- **Iteración 2:** Valores $V = [59.50, 63.16, 61.12]$, $R$ se mantiene → **ÓPTIMO**
- **Política óptima:** $[2, 1, 1]$ (decisión 2 en estado 0, decisión 1 en estados 1 y 2)
- **Costo:** $6.09$

**Aproximaciones Sucesivas:** Converge a la misma política $[2, 1, 1]$ tras ~150 iteraciones con $\epsilon = 10^{-6}$ y $\alpha = 0.9$.

---

## 7. Referencias

Las fórmulas y algoritmos implementados siguen las siguientes referencias del documento de especificación:

| # | Contenido |
|---|-----------|
| [1] | Espacio de estados $E$ y conjunto de decisiones $K$ |
| [2] | Vector de probabilidad inicial $a$ |
| [3] | Ecuaciones de Chapman-Kolmogorov para transiciones de $n$ pasos |
| [4] | Matrices de transición por decisión $P_{ij}(k)$ |
| [5] | Matriz de costos esperados $C_{ik}$ |
| [6] | Factor de descuento $\alpha = (1+i)^{-1}$ |
| [7] | Mejoramiento de políticas con descuento |
| [8] | Probabilidades incondicionales usando vector inicial |
| [9] | Vector de estado estable y tiempos de recurrencia |
| [10] | Tiempos de primera pasada |
| [11] | Probabilidades de absorción |
| [12] | Políticas viables e ingreso interactivo de políticas |
| [13-14] | Enumeración exhaustiva de políticas |
| [15-17] | Mejoramiento de políticas (sin descuento, Howard) |
| [18-19] | Mejoramiento con descuento |
| [20-21] | Método de aproximaciones sucesivas |
| [22-24] | Programación lineal con variables $Y_{ik}$ y transformación $D_{ik}$ |

---

*Paquete desarrollado conforme a la especificación de Cadenas de Markov y Procesos Markovianos de Decisión. Todos los comentarios en el código fuente están en español.*
