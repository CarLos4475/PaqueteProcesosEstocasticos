# Caso de Prueba: Optimización de Publicidad con PMD

## Resumen

Una empresa debe elegir semanalmente entre tres medios publicitarios (Radio, TV, Periódico) para maximizar sus utilidades. El volumen de ventas se clasifica en tres estados (Regular, Bueno, Excelente). Este documento muestra cómo el paquete de Cadenas de Markov y PMD resolvió el problema usando 4 métodos distintos, todos convergiendo a la misma política óptima.

---

## 1. Definición del Problema

### 1.1 Espacio de Estados

| Índice | Estado | Descripción |
|--------|--------|-------------|
| 0 | Regular | Ventas semanales bajas |
| 1 | Bueno | Ventas semanales medias |
| 2 | Excelente | Ventas semanales altas |

$$E = \{0, 1, 2\} \quad (m = 2)$$

### 1.2 Conjunto de Decisiones y Costos Fijos

| Índice | Decisión | Medio | Costo Fijo |
|--------|----------|-------|-------------|
| 0 | $k=1$ | Radio | \$200 |
| 1 | $k=2$ | TV | \$900 |
| 2 | $k=3$ | Periódico | \$300 |

$$K = \{1, 2, 3\} \quad (K = 3)$$

### 1.3 Matrices de Transición $P_{ij}(k)$

**Radio ($k=1$):**
$$ P^{(1)} = \begin{bmatrix} 0.4 & 0.5 & 0.1 \\ 0.1 & 0.7 & 0.2 \\ 0.1 & 0.2 & 0.7 \end{bmatrix}$$

**TV ($k=2$):**
$$ P^{(2)} = \begin{bmatrix} 0.7 & 0.2 & 0.1 \\ 0.3 & 0.6 & 0.1 \\ 0.1 & 0.7 & 0.2 \end{bmatrix}$$

**Periódico ($k=3$):**
$$ P^{(3)} = \begin{bmatrix} 0.2 & 0.5 & 0.3 \\ 0.0 & 0.7 & 0.3 \\ 0.0 & 0.2 & 0.8 \end{bmatrix}$$

> **Observación clave:** Bajo Periódico ($k=3$), los estados 1 y 2 **nunca transicionan al estado 0**. Esto hace que el estado Regular sea transitorio bajo ciertas políticas.

### 1.4 Matrices de Ingreso $Ingreso_{ij}(k)$

**Radio ($k=1$):**
$$ Ingreso^{(1)} = \begin{bmatrix} 400 & 520 & 600 \\ 300 & 400 & 700 \\ 200 & 250 & 500 \end{bmatrix}$$

**TV ($k=2$):**
$$ Ingreso^{(2)} = \begin{bmatrix} 1000 & 1300 & 1600 \\ 800 & 1000 & 1700 \\ 600 & 700 & 1100 \end{bmatrix}$$

**Periódico ($k=3$):**
$$ Ingreso^{(3)} = \begin{bmatrix} 400 & 530 & 710 \\ 350 & 450 & 800 \\ 250 & 400 & 650 \end{bmatrix}$$

---

## 2. Cálculo de la Matriz de Costos $C_{ik}$

El paquete **auto-genera** $C_{ik}$ usando la fórmula:

$$C_{ik} = \sum_{j=0}^{2} P_{ij}(k) \cdot Ingreso_{ij}(k) - CostoFijo(k)$$

Como el objetivo es **maximizar**, los valores se niegan ($-C_{ik}$) para que el algoritmo de minimización encuentre el máximo.

### Cálculo detallado:

**Radio ($k=0$):**

| Estado $i$ | $\sum_j P_{ij} \cdot Ing_{ij}$ | $-$ Costo | $C_{ik}$ real | $C_{ik}$ interno |
|------------|-------------------------------|-----------|---------------|-------------------|
| 0 (Regular) | $0.4\times400+0.5\times520+0.1\times600 = 480$ | $-200$ | $280$ | $-280$ |
| 1 (Bueno) | $0.1\times300+0.7\times400+0.2\times700 = 450$ | $-200$ | $250$ | $-250$ |
| 2 (Excelente) | $0.1\times200+0.2\times250+0.7\times500 = 420$ | $-200$ | $220$ | $-220$ |

**TV ($k=1$):**

| Estado $i$ | $\sum_j P_{ij} \cdot Ing_{ij}$ | $-$ Costo | $C_{ik}$ real | $C_{ik}$ interno |
|------------|-------------------------------|-----------|---------------|-------------------|
| 0 (Regular) | $0.7\times1000+0.2\times1300+0.1\times1600 = 1120$ | $-900$ | $220$ | $-220$ |
| 1 (Bueno) | $0.3\times800+0.6\times1000+0.1\times1700 = 1010$ | $-900$ | $110$ | $-110$ |
| 2 (Excelente) | $0.1\times600+0.7\times700+0.2\times1100 = 770$ | $-900$ | $-130$ | $+130$ |

> **Nota:** TV en estado Excelente genera **pérdida** ($-130). El alto costo de TV (\$900) no se justifica en ese estado.

**Periódico ($k=2$):**

| Estado $i$ | $\sum_j P_{ij} \cdot Ing_{ij}$ | $-$ Costo | $C_{ik}$ real | $C_{ik}$ interno |
|------------|-------------------------------|-----------|---------------|-------------------|
| 0 (Regular) | $0.2\times400+0.5\times530+0.3\times710 = 558$ | $-300$ | $258$ | $-258$ |
| 1 (Bueno) | $0.0\times350+0.7\times450+0.3\times800 = 555$ | $-300$ | $255$ | $-255$ |
| 2 (Excelente) | $0.0\times250+0.2\times400+0.8\times650 = 600$ | $-300$ | $300$ | $-300$ |

### Matriz $C_{ik}$ final (interna, para minimizar):

$$C = \begin{bmatrix} -280 & -220 & -258 \\ -250 & -110 & -255 \\ -220 & +130 & -300 \end{bmatrix}$$

---

## 3. Resultados de los 4 Métodos

### 3.1 Enumeración Exhaustiva

Se evaluaron las $3^3 = 27$ políticas posibles. La matriz de transición bajo cada política $R$ es $P_R[i][j] = P_{ij}(R(i))$, y el costo esperado es:

$$E(C) = \sum_{i=0}^{2} C_{i,R(i)} \cdot \pi_i \quad \text{donde } \pi = \pi \cdot P_R$$

**Mejor política encontrada:** $R = [1, 3, 3]$
- Estado 0 (Regular) → Decisión 1 (Radio, \$200)
- Estado 1 (Bueno) → Decisión 3 (Periódico, \$300)
- Estado 2 (Excelente) → Decisión 3 (Periódico, \$300)

**Costo mínimo:** $E(C) = -282.00$ → **Utilidad máxima: \$282.00 / semana**

---

### 3.2 Mejoramiento de Políticas (Sin Descuento)

Algoritmo de Howard. Partiendo de $R_0 = [1, 1, 1]$ (Radio para todos):

| Iteración | Política | Ganancia $g$ | Valores $V$ |
|-----------|----------|-------------|-------------|
| 1 | $[1, 1, 1]$ | $-243.14$ | $V = [-111.4, -60.0, 0.0]$ |
| 2 | $[1, 1, 3]$ | $-276.30$ | $V = [92.6, 118.5, 0.0]$ |
| 3 | $[1, 3, 3]$ | $\mathbf{-282.00}$ | $V = [78.3, 90.0, 0.0]$ |
| 4 | $[1, 3, 3]$ | — | **ÓPTIMO** |

**Sistema resuelto en cada iteración (ej. iteración 3):**
$$g + V_i - \sum_j P_{ij}(R(i)) \cdot V_j = C_{i,R(i)}, \quad V_2 = 0$$

$$g = -282.00, \quad V_0 = 78.33, \quad V_1 = 90.00$$

**Paso de mejoramiento:** Para cada estado $i$, se evalúan las 3 decisiones y se elige la que minimiza:
$$C_{ik} + \sum_j P_{ij}(k) \cdot V_j - V_i$$

En la iteración 3, para todo estado $i$, la decisión actual ya es óptima → convergencia.

**Resultado:** $R = [1, 3, 3]$ con utilidad **\$282.00**

---

### 3.3 Mejoramiento de Políticas con Descuento ($\alpha = 0.95$)

Mismo algoritmo pero con factor de descuento:

$$V_i = C_{ik} + \alpha \cdot \sum_j P_{ij}(k) \cdot V_j$$

| Iteración | Política | $V_0$ |
|-----------|----------|-------|
| 1 | $[1, 1, 1]$ | $-4926.55$ |
| 2 | $[1, 1, 3]$ | $-5498.13$ |
| 3 | $[1, 3, 3]$ | $\mathbf{-5602.63}$ |
| 4 | $[1, 3, 3]$ | **ÓPTIMO** |

**Resultado:** Misma política $[1, 3, 3]$. El valor $V_0 = -5602.63$ es el valor presente descontado de la utilidad futura con $\alpha = 0.95$. La utilidad por período es $\$282.00$.

---

### 3.4 Programación Lineal

Formulación como PL con variables $Y_{ik}$ (probabilidad estacionaria conjunta):

**Minimizar:**

$$ Z = -280Y_{00} -220Y_{01} -258Y_{02} -250Y_{10} -110Y_{11} -255Y_{12} -220Y_{20} +130Y_{21} -300Y_{22}$$

**Sujeto a:**

$$ \sum_i \sum_k Y_{ik} = 1$$

$$ \sum_k Y_{jk} - \sum_i \sum_k Y_{ik} \cdot P_{ij}(k) = 0 \quad (j = 0,1,2)$$

$$ Y_{ik} \geq 0$$

**Resuelto con Simplex (Método de la Gran M):**

| Variable | Valor |
|----------|-------|
| $Y_{0,0}$ | $0$ |
| $Y_{1,2}$ | $0.4$ |
| $Y_{2,2}$ | $0.6$ |
| Resto | $0$ |

**Valor óptimo:** $Z = -282.00$

**Transformación a política determinística:**
$$D_{ik} = \frac{Y_{ik}}{\sum_k Y_{ik}}$$

| Estado $i$ | $D_{i,1}$ | $D_{i,2}$ | $D_{i,3}$ | Decisión |
|------------|-----------|-----------|-----------|----------|
| 0 (Regular) | $0$ | $0$ | $0$ | 1 (por defecto)* |
| 1 (Bueno) | $0$ | $0$ | $1.0$ | 3 |
| 2 (Excelente) | $0$ | $0$ | $1.0$ | 3 |

> \* El estado 0 tiene $Y_{0k} = 0$ para todo $k$ porque es **transitorio** bajo la política óptima: una vez que las ventas mejoran a Bueno/Excelente, nunca regresan a Regular. Las 3 decisiones son equivalentes en ese estado.

---

## 4. Interpretación de la Política Óptima

$$\boxed{R^* = [\text{Radio}, \text{Periódico}, \text{Periódico}]}$$

| Estado actual | Acción recomendada | Costo | Lógica |
|---------------|-------------------|-------|--------|
| **Regular** | 📻 Radio | \$200 | Radio es barato y tiene buena probabilidad de sacar de Regular (60% a Bueno/Excelente) |
| **Bueno** | 📰 Periódico | \$300 | Periódico maximiza retención en Bueno/Excelente (nunca cae a Regular) |
| **Excelente** | 📰 Periódico | \$300 | Periódico da el mayor ingreso esperado neto en Excelente (\$300) |

**¿Por qué no TV?**
- TV cuesta \$900, demasiado caro
- En estado Excelente, TV genera **pérdida** neta de -\$130
- Solo en estado Regular TV compite (utilidad \$220 vs Radio \$280), pero Radio gana

**¿Por qué no Radio en Bueno/Excelente?**
- Radio cuesta \$200 pero en estados buenos el Periódico (\$300) genera más ingreso neto:
  - Bueno: Radio = \$250, Periódico = \$255
  - Excelente: Radio = \$220, Periódico = \$300

---

## 5. Verificación de Consistencia

Los 4 métodos convergen al mismo resultado:

| Método | Política | Utilidad | ¿Coincide? |
|--------|----------|----------|------------|
| Enumeración Exhaustiva | $[1, 3, 3]$ | \$282.00 | ✓ |
| Mejoramiento s/desc | $[1, 3, 3]$ | \$282.00 | ✓ |
| Mejoramiento c/desc ($\alpha=0.95$) | $[1, 3, 3]$ | \$282.00 | ✓ |
| Programación Lineal | $[1, 3, 3]$ | \$282.00 | ✓ |

---

## 6. Propiedades de la Política Óptima

**Distribución estacionaria:**
$$\pi = [0, \ 0.4, \ 0.6]$$

El estado Regular tiene $\pi_0 = 0$ porque la cadena es **absorbente** hacia los estados {Bueno, Excelente} bajo esta política. Una vez que se sale de Regular, no se regresa.

**Matriz de transición bajo $R^*$:**
$$P_{R^*} = \begin{bmatrix} 0.4 & 0.5 & 0.1 \\ 0.0 & 0.7 & 0.3 \\ 0.0 & 0.2 & 0.8 \end{bmatrix}$$

**Comportamiento a largo plazo:**
- 40% del tiempo en estado Bueno
- 60% del tiempo en estado Excelente
- 0% en estado Regular (transitorio)

**Utilidad desglosada:**
$$\text{Utilidad} = 0.4 \times 255 + 0.6 \times 300 = 102 + 180 = \$282$$

---

## 7. Cómo Reproducir este Análisis

```bash
# Compilar
gcc -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm

# Ejecutar con el archivo de prueba
./markov.exe test_publicidad.txt

# En el menú:
#   5 (Módulo 2: PMD)
#   6 (Prueba Completa)
#   → Ejecuta los 4 métodos y muestra todos los resultados
```

O ejecutar métodos individuales desde el menú interactivo del Módulo 2.

El archivo `test_publicidad.txt` contiene todos los datos del problema en el formato que el paquete lee automáticamente, incluyendo las matrices de ingreso para que $C_{ik}$ se auto-genere.

---

*Documento generado a partir de la ejecución real del paquete sobre los datos del problema de publicidad. Todos los valores numéricos fueron calculados por el programa y verificados manualmente.*
