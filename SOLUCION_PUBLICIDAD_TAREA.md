# Tarea: Optimización de Estrategia Publicitaria con PMD

---

## Enunciado

Una empresa debe elegir semanalmente entre tres medios publicitarios para maximizar sus utilidades. El volumen de ventas se clasifica en tres estados: **Regular** (0), **Bueno** (1) y **Excelente** (2). Los medios disponibles y sus costos fijos semanales son:

| Decisión | Medio | Costo fijo |
|:---:|-------|:-----------:|
| 1 | Radio | \$200 |
| 2 | TV | \$900 |
| 3 | Periódico | \$300 |

Las probabilidades de transición entre estados de ventas dependen del medio elegido y están dadas por las siguientes matrices. Además, cada transición genera un ingreso según las matrices correspondientes.

Se ha estimado una tasa de interés del 5.26% semanal, equivalente a un factor de descuento α = 0.95. La distribución inicial de ventas es P(X₀ = 0) = 0.4, P(X₀ = 1) = 0.3, P(X₀ = 2) = 0.3.

Se pide:

> **A)** Definir formalmente el espacio de estados E, el conjunto de decisiones K, y los costos fijos por decisión.
>
> **B)** Escribir las matrices de transición P^(k) y las matrices de ingreso Ingreso^(k) para cada medio publicitario.
>
> **C)** Calcular *a mano* la matriz de costos netos C_ik para cada estado y decisión, usando la fórmula:
> $$C_{ik} = \left(\sum_j P_{ij}(k) \cdot Ingreso_{ij}(k)\right) - CostoFijo(k)$$
> Mostrar el desarrollo completo de los 9 cálculos. Explicar por qué los valores se niegan internamente.
>
> **D)** ¿Cuántas políticas determinísticas existen? Enumerarlas conceptualmente (no una por una).
>
> **E)** Realizar *a mano* **una iteración** del algoritmo de mejoramiento de políticas sin descuento, partiendo de la política inicial R₀ = [Radio, Radio, Radio] = [1, 1, 1]. Resolver el sistema de determinación del valor y aplicar el paso de mejoramiento, indicando la nueva política obtenida.
>
> **F)** Plantear el problema de programación lineal asociado (variables, función objetivo, restricciones) y, *con ayuda del paquete de programas*, reportar los valores óptimos de Y_ik y la política determinística resultante.
>
> **G)** *Con ayuda del paquete de programas*, ejecutar los 4 métodos y reportar la política óptima. Interpretar el resultado: ¿por qué es óptima cada decisión? ¿Qué sucede con el estado Regular a largo plazo?

---

## Solución

---

### A) Espacio de estados, decisiones y costos fijos

**Espacio de estados:**

$$E = \{0,\ 1,\ 2\}, \qquad m = 2 \quad \text{(3 estados)}$$

| Índice | Estado | Descripción |
|:---:|--------|-------------|
| 0 | Regular | Ventas semanales bajas |
| 1 | Bueno | Ventas semanales medias |
| 2 | Excelente | Ventas semanales altas |

**Conjunto de decisiones:** $K = 3$

| Índice | Decisión | Medio | Costo fijo semanal |
|:---:|:---:|-------|:---:|
| 1 | $k = 1$ | Radio | \$200 |
| 2 | $k = 2$ | TV | \$900 |
| 3 | $k = 3$ | Periódico | \$300 |

**Factor de descuento:** $\alpha = 0.95$, equivalente a una tasa de interés $i = \frac{1}{\alpha} - 1 \approx 5.26\%$ semanal.

**Vector de probabilidad inicial:** $a = [0.4,\ 0.3,\ 0.3]$.

**Objetivo:** *Maximizar* la utilidad esperada descontada a largo plazo. El paquete internamente minimiza, por lo que los costos se niegan para convertir el problema: $\min(-\text{utilidad}) \equiv \max(\text{utilidad})$.

---

### B) Matrices de transición y de ingreso

#### Matrices de transición $P^{(k)}$

**Radio ($k = 1$):**
$$P^{(1)} = \begin{bmatrix} 0.4 & 0.5 & 0.1 \\ 0.1 & 0.7 & 0.2 \\ 0.1 & 0.2 & 0.7 \end{bmatrix}$$

**TV ($k = 2$):**
$$P^{(2)} = \begin{bmatrix} 0.7 & 0.2 & 0.1 \\ 0.3 & 0.6 & 0.1 \\ 0.1 & 0.7 & 0.2 \end{bmatrix}$$

**Periódico ($k = 3$):**
$$P^{(3)} = \begin{bmatrix} 0.2 & 0.5 & 0.3 \\ 0.0 & 0.7 & 0.3 \\ 0.0 & 0.2 & 0.8 \end{bmatrix}$$

> **Observación importante:** Bajo Periódico, desde los estados Bueno y Excelente **nunca se regresa al estado Regular** (filas 1 y 2 tienen $P_{i0}(3) = 0$). Esto vuelve transitorio al estado 0 bajo ciertas políticas.

#### Matrices de ingreso $Ingreso^{(k)}$

**Radio ($k = 1$):**
$$Ingreso^{(1)} = \begin{bmatrix} 400 & 520 & 600 \\ 300 & 400 & 700 \\ 200 & 250 & 500 \end{bmatrix}$$

**TV ($k = 2$):**
$$Ingreso^{(2)} = \begin{bmatrix} 1000 & 1300 & 1600 \\ 800 & 1000 & 1700 \\ 600 & 700 & 1100 \end{bmatrix}$$

**Periódico ($k = 3$):**
$$Ingreso^{(3)} = \begin{bmatrix} 400 & 530 & 710 \\ 350 & 450 & 800 \\ 250 & 400 & 650 \end{bmatrix}$$

---

### C) Cálculo manual de la matriz de costos netos $C_{ik}$

La fórmula usada por el paquete es:

$$C_{ik}^{\text{interno}} = -\left[ \sum_{j=0}^{2} P_{ij}(k) \cdot Ingreso_{ij}(k) \;-\; CostoFijo(k) \right]$$

El signo negativo convierte la maximización de utilidad en minimización de costo (el paquete siempre minimiza). La **utilidad neta real** es el valor sin negar:

$$UtilidadNeta(i,k) = \underbrace{\sum_j P_{ij}(k) \cdot Ingreso_{ij}(k)}_{\text{Ingreso esperado}} \;-\; CostoFijo(k)$$

---

#### Radio ($k = 1$), costo fijo = \$200

**Estado 0 (Regular):**
$$\begin{aligned} E[Ingreso] &= 0.4(400) + 0.5(520) + 0.1(600) \\ &= 160 + 260 + 60 = 480 \\ UtilidadNeta &= 480 - 200 = \mathbf{280} \\ C_{0,1}^{\text{interno}} &= -280 \end{aligned}$$

**Estado 1 (Bueno):**
$$\begin{aligned} E[Ingreso] &= 0.1(300) + 0.7(400) + 0.2(700) \\ &= 30 + 280 + 140 = 450 \\ UtilidadNeta &= 450 - 200 = \mathbf{250} \\ C_{1,1}^{\text{interno}} &= -250 \end{aligned}$$

**Estado 2 (Excelente):**
$$\begin{aligned} E[Ingreso] &= 0.1(200) + 0.2(250) + 0.7(500) \\ &= 20 + 50 + 350 = 420 \\ UtilidadNeta &= 420 - 200 = \mathbf{220} \\ C_{2,1}^{\text{interno}} &= -220 \end{aligned}$$

---

#### TV ($k = 2$), costo fijo = \$900

**Estado 0 (Regular):**
$$\begin{aligned} E[Ingreso] &= 0.7(1000) + 0.2(1300) + 0.1(1600) \\ &= 700 + 260 + 160 = 1120 \\ UtilidadNeta &= 1120 - 900 = \mathbf{220} \\ C_{0,2}^{\text{interno}} &= -220 \end{aligned}$$

**Estado 1 (Bueno):**
$$\begin{aligned} E[Ingreso] &= 0.3(800) + 0.6(1000) + 0.1(1700) \\ &= 240 + 600 + 170 = 1010 \\ UtilidadNeta &= 1010 - 900 = \mathbf{110} \\ C_{1,2}^{\text{interno}} &= -110 \end{aligned}$$

**Estado 2 (Excelente):**
$$\begin{aligned} E[Ingreso] &= 0.1(600) + 0.7(700) + 0.2(1100) \\ &= 60 + 490 + 220 = 770 \\ UtilidadNeta &= 770 - 900 = \mathbf{-130} \quad \text{(¡pérdida!)} \\ C_{2,2}^{\text{interno}} &= -(-130) = \mathbf{+130} \end{aligned}$$

> **Nota:** En estado Excelente, TV genera **pérdida** neta de \$130. Su alto costo fijo (\$900) no se justifica. Esto se refleja en que $C_{2,2} = +130$ es el **único valor positivo** de toda la matriz: el algoritmo de minimización lo evitará.

---

#### Periódico ($k = 3$), costo fijo = \$300

**Estado 0 (Regular):**
$$\begin{aligned} E[Ingreso] &= 0.2(400) + 0.5(530) + 0.3(710) \\ &= 80 + 265 + 213 = 558 \\ UtilidadNeta &= 558 - 300 = \mathbf{258} \\ C_{0,3}^{\text{interno}} &= -258 \end{aligned}$$

**Estado 1 (Bueno):**
$$\begin{aligned} E[Ingreso] &= 0.0(350) + 0.7(450) + 0.3(800) \\ &= 0 + 315 + 240 = 555 \\ UtilidadNeta &= 555 - 300 = \mathbf{255} \\ C_{1,3}^{\text{interno}} &= -255 \end{aligned}$$

**Estado 2 (Excelente):**
$$\begin{aligned} E[Ingreso] &= 0.0(250) + 0.2(400) + 0.8(650) \\ &= 0 + 80 + 520 = 600 \\ UtilidadNeta &= 600 - 300 = \mathbf{300} \\ C_{2,3}^{\text{interno}} &= -300 \end{aligned}$$

---

#### Matriz de costos resultante

**Utilidad neta real $U_{ik}$ (lo que gana la empresa):**

$$U = \begin{bmatrix} 280 & 220 & 258 \\ 250 & 110 & 255 \\ 220 & -130 & 300 \end{bmatrix}$$

**Matriz $C_{ik}$ interna (negada, para minimizar):**

$$C = \begin{bmatrix} -280 & -220 & -258 \\ -250 & -110 & -255 \\ -220 & \mathbf{+130} & -300 \end{bmatrix}$$

---

### D) Número de políticas determinísticas

Cada política asigna una decisión $d(i) \in \{1, 2, 3\}$ a cada estado $i \in \{0, 1, 2\}$:

$$R = [\,d(0),\ d(1),\ d(2)\,]$$

Número total de políticas: $K^{\,m+1} = 3^3 = \mathbf{27}$ políticas.

---

### E) Una iteración de mejoramiento de políticas (sin descuento) desde $R_0 = [1, 1, 1]$

**Política inicial:** $R_0 = [1, 1, 1]$ — Radio en los tres estados.

#### Paso 1 — Determinación del valor

Bajo $R_0$, la matriz de transición es $P_{R_0} = P^{(1)}$:

$$P_{R_0} = \begin{bmatrix} 0.4 & 0.5 & 0.1 \\ 0.1 & 0.7 & 0.2 \\ 0.1 & 0.2 & 0.7 \end{bmatrix}$$

El sistema a resolver es ($V_2 = 0$):

$$\begin{cases} g + V_0 - (0.4V_0 + 0.5V_1 + 0.1V_2) = C_{0,1} = -280 \\ g + V_1 - (0.1V_0 + 0.7V_1 + 0.2V_2) = C_{1,1} = -250 \\ g + V_2 - (0.1V_0 + 0.2V_1 + 0.7V_2) = C_{2,1} = -220 \end{cases}$$

Con $V_2 = 0$:

$$\begin{cases} g + 0.6V_0 - 0.5V_1 = -280 & (1) \\ g - 0.1V_0 + 0.3V_1 = -250 & (2) \\ g - 0.1V_0 - 0.2V_1 = -220 & (3) \end{cases}$$

De (3): $\quad g = -220 + 0.1V_0 + 0.2V_1$

Sustituyendo en (2):
$$(-220 + 0.1V_0 + 0.2V_1) - 0.1V_0 + 0.3V_1 = -250$$
$$-220 + 0.5V_1 = -250$$
$$0.5V_1 = -30$$
$$\boxed{V_1 = -60}$$

Sustituyendo $V_1$ en la expresión de $g$:
$$g = -220 + 0.1V_0 + 0.2(-60) = -232 + 0.1V_0$$

Sustituyendo en (1):
$$(-232 + 0.1V_0) + 0.6V_0 - 0.5(-60) = -280$$
$$-232 + 0.1V_0 + 0.6V_0 + 30 = -280$$
$$-202 + 0.7V_0 = -280$$
$$0.7V_0 = -78$$
$$\boxed{V_0 = -\frac{78}{0.7} = -111.4286}$$

$$\boxed{g = -232 + 0.1(-111.4286) = -243.1429}$$

**Resultado de la evaluación:**

$$\begin{aligned} g &= -243.1429 \quad \text{(costo promedio semanal)} \\ V_0 &= -111.4286 \\ V_1 &= -60.0000 \\ V_2 &= 0 \end{aligned}$$

La utilidad promedio actual es $\$243.14$ / semana.

#### Paso 2 — Mejoramiento

Para cada estado $i$, se evalúa cada decisión $k$ con:

$$\text{test} = C_{ik} + \sum_{j} P_{ij}(k) \cdot V_j - V_i$$

Se elige la que da el **menor** valor de test.

---

**Estado 0** ($V_0 = -111.4286$):

| $k$ | Medio | Cálculo de $\sum P_{0j}V_j$ | test |
|:---:|-------|------------------------------|:---:|
| 1 | Radio | $0.4(-111.43) + 0.5(-60) + 0.1(0) = -44.57 - 30 = \mathbf{-74.57}$ | $-280 + (-74.57) - (-111.43) = \mathbf{-243.14}$ |
| 2 | TV | $0.7(-111.43) + 0.2(-60) + 0.1(0) = -78.00 - 12 = \mathbf{-90.00}$ | $-220 + (-90.00) - (-111.43) = \mathbf{-198.57}$ |
| 3 | Periódico | $0.2(-111.43) + 0.5(-60) + 0.3(0) = -22.29 - 30 = \mathbf{-52.29}$ | $-258 + (-52.29) - (-111.43) = \mathbf{-198.86}$ |

→ **Mejor: Radio ($k = 1$)** con test = $-243.14$.

---

**Estado 1** ($V_1 = -60.0000$):

| $k$ | Medio | Cálculo de $\sum P_{1j}V_j$ | test |
|:---:|-------|------------------------------|:---:|
| 1 | Radio | $0.1(-111.43) + 0.7(-60) + 0.2(0) = -11.14 - 42 = \mathbf{-53.14}$ | $-250 + (-53.14) - (-60) = \mathbf{-243.14}$ |
| 2 | TV | $0.3(-111.43) + 0.6(-60) + 0.1(0) = -33.43 - 36 = \mathbf{-69.43}$ | $-110 + (-69.43) - (-60) = \mathbf{-119.43}$ |
| 3 | Periódico | $0.0(-111.43) + 0.7(-60) + 0.3(0) = 0 - 42 = \mathbf{-42.00}$ | $-255 + (-42.00) - (-60) = \mathbf{-237.00}$ |

→ **Mejor: Radio ($k = 1$)** con test = $-243.14$.

---

**Estado 2** ($V_2 = 0$):

| $k$ | Medio | Cálculo de $\sum P_{2j}V_j$ | test |
|:---:|-------|------------------------------|:---:|
| 1 | Radio | $0.1(-111.43) + 0.2(-60) + 0.7(0) = -11.14 - 12 = \mathbf{-23.14}$ | $-220 + (-23.14) - 0 = \mathbf{-243.14}$ |
| 2 | TV | $0.1(-111.43) + 0.7(-60) + 0.2(0) = -11.14 - 42 = \mathbf{-53.14}$ | $+130 + (-53.14) - 0 = \mathbf{+76.86}$ |
| 3 | Periódico | $0.0(-111.43) + 0.2(-60) + 0.8(0) = 0 - 12 = \mathbf{-12.00}$ | $-300 + (-12.00) - 0 = \mathbf{-312.00}$ |

→ **Mejor: Periódico ($k = 3$)** con test = $-312.00$.

---

#### Resultado de la iteración

| Estado | Decisión anterior | Nueva decisión | ¿Cambió? |
|:---:|:---:|:---:|:---:|
| 0 (Regular) | 1 (Radio) | 1 (Radio) | No |
| 1 (Bueno) | 1 (Radio) | 1 (Radio) | No |
| 2 (Excelente) | 1 (Radio) | **3 (Periódico)** | **Sí** |

$$\boxed{R_1 = [1,\ 1,\ 3] = [\text{Radio},\ \text{Radio},\ \text{Periódico}]}$$

La utilidad promedio estimada mejora de \$243.14 a aproximadamente \$276.30 / semana (calculado por el paquete en la iteración siguiente).

---

### F) Planteamiento del PPL y solución con el paquete

#### Variables de decisión

$Y_{ik}$ = probabilidad estacionaria conjunta de estar en el estado $i$ y tomar la decisión $k$.

Total: $3 \times 3 = 9$ variables: $Y_{0,1}, Y_{0,2}, Y_{0,3}, Y_{1,1}, Y_{1,2}, Y_{1,3}, Y_{2,1}, Y_{2,2}, Y_{2,3}$.

#### Función objetivo (minimizar)

$$\begin{aligned} \min Z = &-280\,Y_{0,1} - 220\,Y_{0,2} - 258\,Y_{0,3} \\ &-250\,Y_{1,1} - 110\,Y_{1,2} - 255\,Y_{1,3} \\ &-220\,Y_{2,1} + 130\,Y_{2,2} - 300\,Y_{2,3} \end{aligned}$$

#### Restricciones

**1. Normalización:**
$$\sum_{i=0}^{2}\sum_{k=1}^{3} Y_{ik} = 1$$

**2. Balances de flujo** (para cada $j = 0, 1, 2$):
$$\sum_{k=1}^{3} Y_{jk} \;-\; \sum_{i=0}^{2}\sum_{k=1}^{3} Y_{ik} \cdot P_{ij}(k) = 0$$

**3. No negatividad:**
$$Y_{ik} \geq 0 \quad \forall i,k$$

#### Solución del PPL (vía paquete — Simplex Big-M)

El paquete resolvió el sistema en **5 iteraciones** del simplex. Variables óptimas no nulas:

| Variable | Valor |
|----------|:-----:|
| $Y_{1,3}$ | $0.4$ |
| $Y_{2,3}$ | $0.6$ |
| Resto | $0$ |

**Valor óptimo:** $Z^* = -282.00$

**Transformación a política determinística:**

$$D_{ik} = \frac{Y_{ik}}{\sum_k Y_{ik}}$$

| Estado $i$ | $D_{i,1}$ | $D_{i,2}$ | $D_{i,3}$ | Decisión |
|:---:|:---:|:---:|:---:|:---:|
| 0 (Regular) | — | — | — | 1 (por defecto)* |
| 1 (Bueno) | 0 | 0 | **1.0** | **3** |
| 2 (Excelente) | 0 | 0 | **1.0** | **3** |

> \* El estado 0 tiene $Y_{0k} = 0$ para todo $k$ porque es **transitorio** bajo la política óptima. Cualquier decisión es equivalente en él; el paquete asigna la 1 por defecto.

**Política resultante del PPL:** $\boxed{R = [1,\ 3,\ 3]}$

---

### G) Resultados completos y análisis

#### Ejecución de los 4 métodos (vía `prueba_completa_pmd`)

| Método | Iteraciones | Política óptima | Utilidad esperada |
|--------|:---:|:---:|:---:|
| **Enumeración exhaustiva** | 27 políticas evaluadas | $[1, 3, 3]$ | \$282.00 |
| **Mejoramiento s/desc.** | 3 iteraciones | $[1, 3, 3]$ | \$282.00 |
| **Mejoramiento c/desc.** ($\alpha=0.95$) | 3 iteraciones | $[1, 3, 3]$ | \$282.00 |
| **Programación Lineal** | 5 iter. simplex | $[1, 3, 3]$ | \$282.00 |

**Los 4 métodos convergen exactamente al mismo resultado.**

---

#### Política óptima

$$\boxed{R^* = [\text{Radio},\ \text{Periódico},\ \text{Periódico}]}$$

| Estado actual | Acción | Costo fijo | Utilidad neta | Razón |
|---------------|:---:|:---:|:---:|---|
| **Regular** | 📻 Radio | \$200 | \$280 | Mayor prob. de sacar de Regular (60% a Bueno/Excelente). Es la opción más barata con buena transición. |
| **Bueno** | 📰 Periódico | \$300 | \$255 | Periódico da más utilidad neta que Radio (\$255 vs \$250) y nunca degrada a Regular. |
| **Excelente** | 📰 Periódico | \$300 | \$300 | Máxima utilidad neta entre las 3 opciones. TV genera **pérdida** de \$130 en este estado. |

---

#### ¿Por qué no TV?

- **Costo fijo altísimo** (\$900/semana).
- En estado Excelente genera **pérdida neta** de \$130 (es el único $C_{ik}$ positivo en la matriz interna: $+130$).
- En estado Regular, su utilidad neta (\$220) es inferior a Radio (\$280).
- En estado Bueno, su utilidad neta (\$110) es la peor de las tres opciones.

TV queda descartada en **todos** los estados por la política óptima.

---

#### ¿Por qué no Radio en Bueno/Excelente?

| Estado | Radio | Periódico | Ganador |
|--------|:---:|:---:|:---:|
| Bueno | \$250 | \$255 | Periódico (+$5) |
| Excelente | \$220 | \$300 | Periódico (+$80) |

Aunque Radio es más barato (\$200 vs \$300), en los estados altos el Periódico genera suficiente ingreso extra para compensar su mayor costo fijo.

---

#### Comportamiento a largo plazo

**Distribución estacionaria bajo $R^*$:**

$$\pi = [\,0,\ 0.4,\ 0.6\,]$$

El estado Regular tiene $\pi_0 = 0$: una vez que las ventas salen de Regular hacia Bueno o Excelente, **nunca regresan**. La cadena es absorbente hacia el subconjunto {Bueno, Excelente}.

**Matriz de transición bajo $R^*$:**

$$P_{R^*} = \begin{bmatrix} 0.4 & 0.5 & 0.1 \\ 0.0 & 0.7 & 0.3 \\ 0.0 & 0.2 & 0.8 \end{bmatrix}$$

- Fila 0: Radio en Regular (puede transicionar a cualquier estado).
- Filas 1-2: Periódico en Bueno y Excelente ($P_{i0} = 0$: nunca a Regular).

**Utilidad desglosada:**

$$\begin{aligned} \text{Utilidad} &= \pi_1 \cdot U_{1,3} + \pi_2 \cdot U_{2,3} \\ &= 0.4 \times 255 + 0.6 \times 300 \\ &= 102 + 180 \\ &= \mathbf{\$282.00\ / \text{semana}} \end{aligned}$$

A largo plazo, la empresa pasa el **40%** de las semanas en estado Bueno y el **60%** en Excelente, con una utilidad estable de **\$282 semanales**.

---

<div align="center">

**— Fin de la tarea —**

*Cálculos manuales verificados contra el paquete `markov_lib`. Resultados del Simplex, enumeración exhaustiva (27 políticas) y convergencia completa del mejoramiento de políticas obtenidos con ayuda del programa.*

</div>
