# Solución — Problema 2 (Póker de los Sábados)

**Problema:** PMD de costo promedio a largo plazo, sin descuento.  
**Objetivo:** Minimizar el costo esperado semanal promedio a largo plazo.

---

## A) Estados, Decisiones y Costos

### Estados (E = {0, 1})

| Estado | Significado |
|--------|-------------|
| 0 | El grupo está de **buen humor** |
| 1 | El grupo está de **mal humor** |

m = 1 → 2 estados (E = {0, 1})

### Decisiones (K = 2)

| Decisión | Significado |
|----------|-------------|
| k = 1 | **Ofrecer refrescos** |
| k = 2 | **No ofrecer refrescos** |

### Matriz de Costos C_ik

| Estado i | k = 1 (Ofrecer) | k = 2 (No ofrecer) |
|----------|:---------------:|:------------------:|
| 0 (Buen humor) | $14 | $0 |
| 1 (Mal humor)  | $14 | $75 |

**Justificación de costos:**

- **Ofrecer refrescos (k = 1):** Cuesta $14 fijos, sin importar el humor. El juego queda tablas ($0 de ganancia/pérdida en el póker).
- **No ofrecer, buen humor (k = 2, i = 0):** No hay costo de refrescos, el juego queda tablas → $0.
- **No ofrecer, mal humor (k = 2, i = 1):** Lo molestan y tiene una pérdida esperada de $75 en el póker.

---

## B) Matrices de Transición para Cada Decisión

Las transiciones **no dependen del estado actual** (solo de la decisión tomada).

### P^(1) — Ofrecer refrescos (k = 1)

Al ofrecer refrescos, el siguiente sábado el grupo estará de buen humor con probabilidad 7/8.

```
         → Bueno(0)  → Malo(1)
Desde 0 [   7/8         1/8    ]
Desde 1 [   7/8         1/8    ]
```

En decimal:

```
         → Bueno(0)   → Malo(1)
Desde 0 [  0.875       0.125  ]
Desde 1 [  0.875       0.125  ]
```

### P^(2) — No ofrecer refrescos (k = 2)

Al no ofrecer refrescos, el siguiente sábado el grupo estará de buen humor con probabilidad 1/8, sin importar el humor actual.

```
         → Bueno(0)  → Malo(1)
Desde 0 [   1/8         7/8    ]
Desde 1 [   1/8         7/8    ]
```

En decimal:

```
         → Bueno(0)   → Malo(1)
Desde 0 [  0.125       0.875  ]
Desde 1 [  0.125       0.875  ]
```

---

## C) Políticas Determinísticas

Cada política asigna una decisión a cada estado.

Total de políticas: K^(m+1) = 2² = **4 políticas**.

| # | Política R = [d(0), d(1)] | Significado | π estacionario | Costo esperado E(C) |
|---|--------------------------|-------------|:---:|:---:|
| 1 | [1, 1] | Siempre ofrecer    | π = (0.875, 0.125) | $14.00 |
| 2 | [1, 2] | Ofrecer si buen humor, no si mal humor | π = (0.5, 0.5) | $44.50 |
| 3 | [2, 1] | **No ofrecer si buen humor, ofrecer si mal humor** | π = (0.5, 0.5) | **$7.00** ✓ |
| 4 | [2, 2] | Nunca ofrecer      | π = (0.125, 0.875) | $65.625 |

**Conclusión:** La **política óptima es R = [2, 1]**:
- Estado 0 (buen humor) → **No ofrecer** refrescos (no hace falta).
- Estado 1 (mal humor) → **Ofrecer** refrescos (cuesta $14 pero evita perder $75).

Costo promedio esperado a largo plazo: **$7.00 / semana**.

---

## D) Una Iteración de Mejoramiento de Políticas desde R₀ = (Ofrecer, Ofrecer)

Se aplica el **algoritmo de Howard para costo promedio** (sin descuento).

### Paso 0 — Política inicial

R₀ = [1, 1] → Ofrecer refrescos en ambos estados.

### Paso 1 — Determinación del valor

Bajo R₀, la matriz de transición es P^(1):

```
        [0.875   0.125]
P_R₀ =  [0.875   0.125]
```

Se resuelve el sistema:

> g + V_i = C_{i, k(i)} + Σⱼ P_{ij}(k(i)) · V_j,   con V_1 = 0

Para i = 1:
```
g + 0 = 14 + (0.875)V_0 + (0.125)(0)
g = 14 + 0.875 V_0          ... (1)
```

Para i = 0:
```
g + V_0 = 14 + (0.875)V_0 + (0.125)(0)
g + V_0 = 14 + 0.875 V_0
```

Sustituyendo (1) en la ecuación de i = 0:
```
(14 + 0.875 V_0) + V_0 = 14 + 0.875 V_0
14 + 1.875 V_0 = 14 + 0.875 V_0
V_0 = 0
```

De (1): **g = 14**, V_0 = 0, V_1 = 0.

**Costo promedio actual:** g = $14/semana.

### Paso 2 — Mejoramiento

Para cada estado i, se evalúa:  **test = C_{ik} + Σⱼ P_{ij}(k) · Vⱼ − V_i**

| Estado | Decisión | Cálculo | test | Marca |
|--------|----------|---------|:----:|:-----:|
| 0 | k = 1 (Ofrecer) | 14 + 0.875(0) + 0.125(0) − 0 | **14** | |
| 0 | k = 2 (No ofrecer) | 0 + 0.125(0) + 0.875(0) − 0 | **0** | **← MEJOR** |
| 1 | k = 1 (Ofrecer) | 14 + 0.875(0) + 0.125(0) − 0 | **14** | **← MEJOR** |
| 1 | k = 2 (No ofrecer) | 75 + 0.125(0) + 0.875(0) − 0 | **75** | |

### Nueva política

```
R₁ = [2, 1] = (No ofrecer en buen humor, Ofrecer en mal humor)
```

**La política cambió en una iteración.** La siguiente iteración confirmaría que R₁ es óptima (R₂ = R₁ con g = $7).

---

## E) Planteamiento del Problema de Programación Lineal

### Variables de decisión

Y_{ik} = probabilidad estacionaria conjunta de estar en el estado i y tomar la decisión k.

| Variable | Significado |
|----------|-------------|
| Y₀₁ | Estado 0, decisión 1 (Ofrecer en buen humor) |
| Y₀₂ | Estado 0, decisión 2 (No ofrecer en buen humor) |
| Y₁₁ | Estado 1, decisión 1 (Ofrecer en mal humor) |
| Y₁₂ | Estado 1, decisión 2 (No ofrecer en mal humor) |

### Función objetivo (Minimizar)

```
Min Z = 14·Y₀₁  +  0·Y₀₂  +  14·Y₁₁  +  75·Y₁₂
```

### Sujeto a

**1. Normalización:**

```
Y₀₁ + Y₀₂ + Y₁₁ + Y₁₂ = 1
```

**2. Balance para el estado 0 (buen humor):**

```
(Y₀₁ + Y₀₂) − [0.875·Y₀₁ + 0.125·Y₀₂ + 0.875·Y₁₁ + 0.125·Y₁₂] = 0
```

En fracciones:
```
(Y₀₁ + Y₀₂) − [(7/8)Y₀₁ + (1/8)Y₀₂ + (7/8)Y₁₁ + (1/8)Y₁₂] = 0
```

**3. Balance para el estado 1 (mal humor):**

```
(Y₁₁ + Y₁₂) − [0.125·Y₀₁ + 0.875·Y₀₂ + 0.125·Y₁₁ + 0.875·Y₁₂] = 0
```

En fracciones:
```
(Y₁₁ + Y₁₂) − [(1/8)Y₀₁ + (7/8)Y₀₂ + (1/8)Y₁₁ + (7/8)Y₁₂] = 0
```

**4. No negatividad:**

```
Y₀₁, Y₀₂, Y₁₁, Y₁₂ ≥ 0
```

### Solución del PPL (obtenida del paquete)

| Variable | Valor óptimo | Significado |
|----------|:---:|------|
| Y₀₁ | 0 | (No conviene ofrecer en buen humor) |
| Y₀₂ | **0.5** | En buen humor → no ofrecer |
| Y₁₁ | **0.5** | En mal humor → ofrecer |
| Y₁₂ | 0 | (No conviene no ofrecer en mal humor) |

**Valor óptimo:** Z* = **$7.00**  
**Política determinística:** R = [2, 1]

---

## Verificación con el paquete `markov_lib`

Se creó el archivo de entrada [`problema2.txt`](markov_lib/problema2.txt) con α = 0.999 (aproximando el caso sin descuento) y se ejecutó la **prueba completa** (batch):

| Método | Política óptima | Costo/Uso |
|--------|:---:|:---:|
| Enumeración exhaustiva | [2, 1] | $7.00 |
| Mejoramiento s/desc. | [2, 1] | $7.00 |
| Mejoramiento c/desc. (α=0.999) | [2, 1] | $7.00 |
| Programación Lineal | [2, 1] | Z* = $7.00 |

**Los 4 métodos convergen a la misma política óptima R = [2, 1] con costo $7.00/semana.**

### Resultado de la iteración D del paquete

La salida del paquete para la primera iteración del mejoramiento sin descuento partiendo de R₀ = [1, 1] coincide exactamente con el cálculo manual:

```
V₀ = -0.000000,  V₁ = 0.000000,  g = 14.000000

Estado 0:  k=1 → test = 14.00    k=2 → test = 0.00  ← MEJOR
Estado 1:  k=1 → test = 14.00 ← MEJOR    k=2 → test = 75.00

R₁ = [2, 1]
```

---

**Respuesta final:** La política óptima es **no ofrecer refrescos cuando el grupo está de buen humor y ofrecerlos cuando está de mal humor**, con un costo promedio esperado de **$7.00 por semana**.
