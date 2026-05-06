#include "markov.h"
#include <string.h>

/* ===================================================================
 * MODULO 2: PROCESOS MARKOVIANOS DE DECISION (PMD)
 *
 * Implementa los 5 algoritmos de resolucion para PMD:
 *   1. Enumeracion Exhaustiva de Politicas
 *   2. Mejoramiento de Politicas (sin descuento)
 *   3. Mejoramiento de Politicas con Descuento
 *   4. Metodo de Aproximaciones Sucesivas
 *   5. Solucion por Programacion Lineal
 * =================================================================== */

/* ===================================================================
 * FUNCIONES AUXILIARES PARA POLITICAS
 * =================================================================== */

Politica* politica_crear(int num_estados, int num_decisiones) {
    Politica *p = (Politica*)malloc(sizeof(Politica));
    if (!p) return NULL;
    p->num_estados    = num_estados;
    p->num_decisiones = num_decisiones;
    p->decision       = (int*)calloc((size_t)num_estados, sizeof(int));
    return p;
}

void politica_destruir(Politica *p) {
    if (!p) return;
    free(p->decision);
    free(p);
}

void politica_imprimir(const Politica *p) {
    printf("\n  Politica R: [");
    for (int i = 0; i < p->num_estados; i++) {
        printf("%d", p->decision[i] + 1); /* +1 para mostrar decision 1..K */
        if (i < p->num_estados - 1) printf(", ");
    }
    printf("]  (decision para cada estado 0..%d)\n", p->num_estados - 1);
}

Politica* politica_copiar(const Politica *origen) {
    Politica *copia = politica_crear(origen->num_estados, origen->num_decisiones);
    if (!copia) return NULL;
    memcpy(copia->decision, origen->decision,
           (size_t)origen->num_estados * sizeof(int));
    return copia;
}

int politica_igual(const Politica *a, const Politica *b) {
    if (a->num_estados != b->num_estados) return 0;
    for (int i = 0; i < a->num_estados; i++)
        if (a->decision[i] != b->decision[i]) return 0;
    return 1;
}

/* ===================================================================
 * Construye la matriz de transicion bajo una politica especifica.
 * P_pol[i][j] = P_dec[decision[i]][i][j]
 * =================================================================== */
Matriz* politica_matriz_transicion(const ModeloMarkov *modelo, const Politica *p) {
    Matriz *Ppol = matriz_crear(modelo->num_estados, modelo->num_estados);
    if (!Ppol) return NULL;
    for (int i = 0; i < modelo->num_estados; i++) {
        int k = p->decision[i];
        for (int j = 0; j < modelo->num_estados; j++)
            Ppol->datos[i][j] = modelo->P_dec[k]->datos[i][j];
    }
    return Ppol;
}

/* ===================================================================
 * politica_leer_consola
 *
 * Permite al usuario ingresar politicas una por una para decidir
 * cuales son viables y cuales absurdas/irrelevantes.
 *
 * El usuario ingresa la decision para cada estado.
 * Se le pregunta si la politica es valida.
 * =================================================================== */
void politica_leer_consola(ModeloMarkov *modelo, Politica **politicas,
                            int *num_politicas) {
    int capacidad = 10;
    *politicas = (Politica*)malloc((size_t)capacidad * sizeof(Politica));
    *num_politicas = 0;

    printf("\n--- Ingreso de Politicas ---\n");
    printf("Estados: E = {0, 1, ..., %d}\n", modelo->num_estados - 1);
    printf("Decisiones: K = {1, 2, ..., %d}\n\n", modelo->num_decisiones);

    int continuar = 1;
    while (continuar) {
        Politica nueva;
        nueva.num_estados    = modelo->num_estados;
        nueva.num_decisiones = modelo->num_decisiones;
        nueva.decision       = (int*)calloc((size_t)modelo->num_estados,
                                            sizeof(int));

        printf("Politica #%d:\n", *num_politicas + 1);
        for (int i = 0; i < modelo->num_estados; i++) {
            char msg[64];
            snprintf(msg, sizeof(msg), "  Decision para estado %d (1-%d): ",
                     i, modelo->num_decisiones);
            printf("%s", msg);
            int d;
            if (scanf("%d", &d) == 1 && d >= 1 && d <= modelo->num_decisiones) {
                nueva.decision[i] = d - 1;
            } else {
                printf("  Invalido, usando decision 1.\n");
                nueva.decision[i] = 0;
            }
            int c;
            while ((c = getchar()) != '\n' && c != EOF);
        }

        politica_imprimir(&nueva);
        printf("  Esta politica es valida? (S/N): ");
        char resp;
        scanf(" %c", &resp);
        int c;
        while ((c = getchar()) != '\n' && c != EOF);

        if (resp == 'S' || resp == 's') {
            if (*num_politicas >= capacidad) {
                capacidad *= 2;
                *politicas = (Politica*)realloc(*politicas,
                    (size_t)capacidad * sizeof(Politica));
            }
            memcpy(&(*politicas)[*num_politicas], &nueva, sizeof(Politica));
            (*num_politicas)++;
        } else {
            free(nueva.decision);
        }

        printf("  Ingresar otra politica? (S/N): ");
        scanf(" %c", &resp);
        while ((c = getchar()) != '\n' && c != EOF);
        if (resp != 'S' && resp != 's') continuar = 0;
    }

    printf("\nSe aceptaron %d politicas.\n", *num_politicas);
    for (int i = 0; i < *num_politicas; i++)
        politica_imprimir(&(*politicas)[i]);
}

/* ===================================================================
 * 1. ENUMERACION EXHAUSTIVA DE POLITICAS
 *
 * Genera todas las K^(m+1) politicas posibles, calcula el costo
 * esperado a largo plazo E(C) = sum_i C_{ik} * pi(i) para cada una,
 * y selecciona la de minimo costo.
 *
 * Cada iteracion del algoritmo se imprime en pantalla.
 * =================================================================== */
void enumeracion_exhaustiva(ModeloMarkov *modelo) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;

    printf("\n=============================================================\n");
    printf("  1. ENUMERACION EXHAUSTIVA DE POLITICAS\n");
    printf("=============================================================\n");
    printf("  Numero de estados: %d\n", n);
    printf("  Numero de decisiones: %d\n", K);
    printf("  Total de politicas a evaluar: %.0f\n\n", pow((double)K, (double)n));

    /* Calcular total de politicas */
    long long total = 1;
    for (int i = 0; i < n; i++) total *= K;

    Politica *mejor_pol = politica_crear(n, K);
    double mejor_costo  = INFINITY;

    Politica actual;
    actual.num_estados    = n;
    actual.num_decisiones = K;
    actual.decision       = (int*)calloc((size_t)n, sizeof(int));

    int iteracion = 0;
    for (long long idx = 0; idx < total; idx++) {
        /* Convertir indice a politica (base K) */
        long long temp = idx;
        for (int i = 0; i < n; i++) {
            actual.decision[i] = (int)(temp % K);
            temp /= K;
        }

        iteracion++;
        printf("\n--- Iteracion %d / %lld ---\n", iteracion, total);
        politica_imprimir(&actual);

        /* Construir matriz de transicion bajo esta politica */
        Matriz *Ppol = politica_matriz_transicion(modelo, &actual);

        /* Calcular estado estacionario pi */
        double *pi = estado_estable(Ppol);
        matriz_destruir(Ppol);

        if (!pi) {
            printf("  (Politica con matriz singular, se omite)\n");
            continue;
        }

        printf("  Estado estacionario pi: [");
        for (int i = 0; i < n; i++) {
            printf("%.6f", pi[i]);
            if (i < n - 1) printf(", ");
        }
        printf("]\n");

        /* Calcular costo esperado: E(C) = sum_i C_{i,k(i)} * pi(i) */
        double costo = 0.0;
        for (int i = 0; i < n; i++) {
            int k = actual.decision[i];
            costo += modelo->C->datos[i][k] * pi[i];
        }

        printf("  Costo esperado E(C) = %.8f\n", costo);

        if (costo < mejor_costo) {
            mejor_costo = costo;
            memcpy(mejor_pol->decision, actual.decision,
                   (size_t)n * sizeof(int));
            printf("  *** NUEVA MEJOR POLITICA! Costo = %.8f ***\n", mejor_costo);
        }

        free(pi);
    }

    printf("\n=============================================================\n");
    printf("  RESULTADO FINAL - Enumeracion Exhaustiva\n");
    printf("=============================================================\n");
    printf("  Mejor politica encontrada:\n");
    politica_imprimir(mejor_pol);
    printf("  Costo minimo esperado: %.8f\n", mejor_costo);

    free(actual.decision);
    politica_destruir(mejor_pol);
}

/* ===================================================================
 * 2. MEJORAMIENTO DE POLITICAS (SIN DESCUENTO)
 *
 * Algoritmo iterativo de Howard:
 *
 * Paso 0: Elegir politica inicial arbitraria R_n, hacer n=1.
 * Paso 1 (Determinacion del valor): Resolver el sistema:
 *          g(R_n) = c_{ik} + sum_j p_{ij}(k) * V_j(R_n) - V_i(R_n)
 *          para i = 0, 1, ..., m, con V_m(R_n) = 0.
 * Paso 2 (Mejoramiento): Para cada estado i, encontrar k que minimiza:
 *          test = c_{ik} + sum_j p_{ij}(k) * V_j(R_n) - V_i(R_n)
 *          Si mejora, actualizar.
 * Prueba de optimalidad: Si R_{n+1} = R_n, detener.
 * =================================================================== */
Politica* mejoramiento_politicas(const ModeloMarkov *modelo,
                                 const Politica *inicial) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;

    printf("\n=============================================================\n");
    printf("  2. MEJORAMIENTO DE POLITICAS (Sin Descuento)\n");
    printf("=============================================================\n");

    Politica *R_actual = politica_copiar(inicial);
    int iter = 0;

    while (1) {
        iter++;
        printf("\n--- ITERACION %d ---\n", iter);
        printf("Politica actual:");
        politica_imprimir(R_actual);

        /* ========================================================
         * Paso 1: Determinacion del valor
         * Resolver: g + V_i - sum_j p_{ij}(k_i) * V_j = c_{i,k_i}
         * para i = 0..m, con V_m = 0.
         *
         * Incognitas: V_0, ..., V_{m-1}, g   (total m+1 = n)
         * Ecuaciones: una para cada i = 0..m (total m+1 = n)
         * ======================================================== */

        Matriz *A = matriz_crear(n, n);
        double *b = (double*)calloc((size_t)n, sizeof(double));
        double *x = (double*)calloc((size_t)n, sizeof(double));
        if (!A || !b || !x) {
            matriz_destruir(A); free(b); free(x);
            politica_destruir(R_actual); return NULL;
        }

        /* Construir el sistema para la politica actual */
        for (int i = 0; i < n; i++) {
            int k = R_actual->decision[i];

            /* Coeficientes de V_0..V_{m-1} */
            for (int j = 0; j < n - 1; j++) {
                if (i == j)
                    A->datos[i][j] = 1.0 - modelo->P_dec[k]->datos[i][j];
                else
                    A->datos[i][j] = -modelo->P_dec[k]->datos[i][j];
            }
            /* Coeficiente de g (ultima columna) */
            A->datos[i][n-1] = 1.0;

            /* Termino independiente (V_m=0 asi que p_im * 0 = 0) */
            b[i] = modelo->C->datos[i][k];
        }

        printf("  Resolviendo sistema lineal de %d x %d...\n", n, n);

        if (!resolver_sistema_lineal(A, b, x)) {
            fprintf(stderr, "  Error: sistema singular en iteracion %d.\n", iter);
            matriz_destruir(A); free(b); free(x);
            politica_destruir(R_actual); return NULL;
        }

        matriz_destruir(A);
        free(b);

        /* Extraer V y g de la solucion */
        /* x = [V_0, V_1, ..., V_{m-1}, g] */
        double *V = (double*)calloc((size_t)n, sizeof(double));
        for (int i = 0; i < n - 1; i++)
            V[i] = x[i];
        V[n-1] = 0.0; /* V_m = 0 */
        double g = x[n-1];

        printf("  Valores V_i:\n");
        for (int i = 0; i < n; i++)
            printf("    V_%d = %12.6f\n", i, V[i]);
        printf("  Ganancia g = %12.6f\n", g);
        free(x);

        /* ========================================================
         * Paso 2: Mejoramiento
         * Para cada estado i, encontrar k que minimiza:
         *   test = c_{ik} + sum_j p_{ij}(k) * V_j  -  V_i
         * ======================================================== */

        Politica *R_nueva = politica_crear(n, K);
        printf("\n  Paso de mejoramiento:\n");

        for (int i = 0; i < n; i++) {
            double mejor_valor  = INFINITY;
            int    mejor_decision = R_actual->decision[i];

            /* Evaluar cada decision posible para el estado i */
            for (int k = 0; k < K; k++) {
                double suma = 0.0;
                for (int j = 0; j < n; j++)
                    suma += modelo->P_dec[k]->datos[i][j] * V[j];

                double valor = modelo->C->datos[i][k] + suma - V[i];

                printf("    Estado %d, decision %d: test = %12.6f",
                       i, k + 1, valor);

                if (valor < mejor_valor - 1e-10) {
                    mejor_valor    = valor;
                    mejor_decision = k;
                    printf("  <-- mejor");
                }
                printf("\n");
            }

            R_nueva->decision[i] = mejor_decision;
        }

        /* ========================================================
         * Prueba de optimalidad
         * ======================================================== */
        if (politica_igual(R_actual, R_nueva)) {
            printf("\n  *** POLITICA OPTIMA ENCONTRADA! (R_%d = R_%d) ***\n",
                   iter + 1, iter);
            printf("  Politica optima:");
            politica_imprimir(R_actual);
            printf("  Ganancia g = %.8f\n", g);
            free(V);
            politica_destruir(R_nueva);
            break;
        }

        /* Actualizar para siguiente iteracion */
        printf("\n  R_%d != R_%d, continuando...\n", iter + 1, iter);
        free(V);
        politica_destruir(R_actual);
        R_actual = R_nueva;
    }

    return R_actual;
}

/* ===================================================================
 * 3. MEJORAMIENTO DE POLITICAS CON DESCUENTO
 *
 * Algoritmo de Howard con factor de descuento alfa.
 *
 * Paso 0: Politica inicial R_n, n=1.
 * Paso 1: Resolver V_i = c_{ik} + alfa * sum_j p_{ij}(k) * V_j
 * Paso 2: Encontrar k que minimiza c_{ik} + alfa * sum_j p_{ij}(k) * V_j
 * Optimalidad: Si R_{n+1} = R_n, terminar.
 * =================================================================== */
Politica* mejoramiento_politicas_descuento(const ModeloMarkov *modelo,
                                           const Politica *inicial) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;
    double alfa = modelo->alfa;

    printf("\n=============================================================\n");
    printf("  3. MEJORAMIENTO DE POLITICAS CON DESCUENTO\n");
    printf("=============================================================\n");
    printf("  Factor de descuento alfa = %.6f\n\n", alfa);

    Politica *R_actual = politica_copiar(inicial);
    int iter = 0;

    while (1) {
        iter++;
        printf("\n--- ITERACION %d ---\n", iter);
        printf("Politica actual:");
        politica_imprimir(R_actual);

        /* ========================================================
         * Paso 1: Determinacion del valor con descuento
         * Sistema: V_i - alfa * sum_j p_{ij}(k) * V_j = c_{ik}
         * para i = 0..m. (m+1 ecuaciones, m+1 incognitas)
         * ======================================================== */

        Matriz *A = matriz_crear(n, n);
        double *b = (double*)calloc((size_t)n, sizeof(double));
        double *V = (double*)calloc((size_t)n, sizeof(double));

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

        printf("  Resolviendo sistema lineal con descuento...\n");
        if (!resolver_sistema_lineal(A, b, V)) {
            fprintf(stderr, "  Error: sistema singular.\n");
            matriz_destruir(A); free(b); free(V);
            politica_destruir(R_actual); return NULL;
        }

        matriz_destruir(A);
        free(b);

        printf("  Valores V_i:\n");
        for (int i = 0; i < n; i++)
            printf("    V_%d = %12.6f\n", i, V[i]);

        /* ========================================================
         * Paso 2: Mejoramiento con descuento
         * Encontrar k que minimiza: c_{ik} + alfa * sum_j p_{ij}(k) * V_j
         * ======================================================== */

        Politica *R_nueva = politica_crear(n, K);
        printf("\n  Paso de mejoramiento (con descuento):\n");

        for (int i = 0; i < n; i++) {
            double mejor_valor    = INFINITY;
            int    mejor_decision = R_actual->decision[i];

            for (int k = 0; k < K; k++) {
                double suma = 0.0;
                for (int j = 0; j < n; j++)
                    suma += modelo->P_dec[k]->datos[i][j] * V[j];

                double valor = modelo->C->datos[i][k] + alfa * suma;

                printf("    Estado %d, decision %d: valor = %12.6f",
                       i, k + 1, valor);

                if (valor < mejor_valor - 1e-10) {
                    mejor_valor    = valor;
                    mejor_decision = k;
                    printf("  <-- mejor");
                }
                printf("\n");
            }

            R_nueva->decision[i] = mejor_decision;
        }

        /* Prueba de optimalidad */
        if (politica_igual(R_actual, R_nueva)) {
            printf("\n  *** POLITICA OPTIMA ENCONTRADA! (R_%d = R_%d) ***\n",
                   iter + 1, iter);
            printf("  Politica optima:");
            politica_imprimir(R_actual);
            printf("  Valor optimo V_0 = %.8f\n", V[0]);

            /* Calcular costo esperado bajo la politica optima */
            Matriz *Popt = politica_matriz_transicion(modelo, R_actual);
            double *pi_opt = estado_estable(Popt);
            if (pi_opt) {
                double costo_opt = 0.0;
                for (int i = 0; i < n; i++)
                    costo_opt += modelo->C->datos[i][R_actual->decision[i]] * pi_opt[i];
                printf("  Costo esperado a largo plazo = %.8f\n", costo_opt);
                free(pi_opt);
            }
            matriz_destruir(Popt);

            free(V);
            politica_destruir(R_nueva);
            break;
        }

        printf("\n  R_%d != R_%d, continuando...\n", iter + 1, iter);
        free(V);
        politica_destruir(R_actual);
        R_actual = R_nueva;
    }

    return R_actual;
}

/* ===================================================================
 * 4. METODO DE APROXIMACIONES SUCESIVAS
 *
 * Algoritmo de iteracion de valor con descuento.
 * Imprime obligatoriamente el valor de V_i^n en cada iteracion.
 *
 * Paso 1: n=1. V_i^1 = min_k c_{ik}
 * Paso 2: n = n+1. V_i^n = min_k (c_{ik} + alfa * sum_j p_{ij}(k) * V_j^{n-1})
 * Tolerancia: |V_i^n - V_i^{n-1}| < epsilon para todo i, o n > N.
 * =================================================================== */
void aproximaciones_sucesivas(const ModeloMarkov *modelo,
                              int max_iter, double epsilon) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;
    double alfa = modelo->alfa;

    printf("\n=============================================================\n");
    printf("  4. METODO DE APROXIMACIONES SUCESIVAS\n");
    printf("=============================================================\n");
    printf("  Maximo de iteraciones N = %d\n", max_iter);
    printf("  Tolerancia epsilon = %.10f\n", epsilon);
    printf("  Factor de descuento alfa = %.6f\n\n", alfa);

    double *V_anterior = (double*)calloc((size_t)n, sizeof(double));
    double *V_actual   = (double*)calloc((size_t)n, sizeof(double));
    int    *politica   = (int*)calloc((size_t)n, sizeof(int));

    /* ========================================================
     * Paso 1: Inicializacion, n=1.
     * V_i^1 = min_k c_{ik}
     * ======================================================== */
    printf("--- ITERACION 1 (inicializacion) ---\n");
    for (int i = 0; i < n; i++) {
        double min_costo = INFINITY;
        for (int k = 0; k < K; k++) {
            if (modelo->C->datos[i][k] < min_costo) {
                min_costo     = modelo->C->datos[i][k];
                politica[i]   = k;
            }
        }
        V_anterior[i] = min_costo;
    }
    printf("  V^1: [");
    for (int i = 0; i < n; i++) {
        printf("%.6f", V_anterior[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");

    int iter = 1;
    int convergio = 0;

    while (iter < max_iter && !convergio) {
        iter++;

        /* ========================================================
         * Paso 2: Calcular V^n
         * V_i^n = min_k (c_{ik} + alfa * sum_j p_{ij}(k) * V_j^{n-1})
         * ======================================================== */
        for (int i = 0; i < n; i++) {
            double min_valor = INFINITY;
            int    mejor_k   = 0;

            for (int k = 0; k < K; k++) {
                double suma = 0.0;
                for (int j = 0; j < n; j++)
                    suma += modelo->P_dec[k]->datos[i][j] * V_anterior[j];

                double valor = modelo->C->datos[i][k] + alfa * suma;

                if (valor < min_valor - 1e-12) {
                    min_valor = valor;
                    mejor_k   = k;
                }
            }

            V_actual[i] = min_valor;
            politica[i] = mejor_k;
        }

        printf("\n--- ITERACION %d ---\n", iter);
        printf("  V^%d: [", iter);
        for (int i = 0; i < n; i++) {
            printf("%.6f", V_actual[i]);
            if (i < n - 1) printf(", ");
        }
        printf("]\n");

        /* Verificar tolerancia */
        convergio = 1;
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(V_actual[i] - V_anterior[i]);
            if (diff > max_diff) max_diff = diff;
            if (diff >= epsilon) convergio = 0;
        }
        printf("  Diferencia maxima: %.10f\n", max_diff);

        if (convergio) {
            printf("  *** CONVERGENCIA ALCANZADA en iteracion %d ***\n", iter);
        }

        /* Intercambiar V */
        double *temp = V_anterior;
        V_anterior = V_actual;
        V_actual = temp;
    }

    if (!convergio && iter >= max_iter) {
        printf("\n  *** Se alcanzo el maximo de iteraciones (%d) ***\n", max_iter);
    }

    printf("\n=============================================================\n");
    printf("  RESULTADO FINAL - Aproximaciones Sucesivas\n");
    printf("=============================================================\n");
    printf("  Iteraciones realizadas: %d\n", iter);
    printf("  Valores V_i:\n");
    for (int i = 0; i < n; i++)
        printf("    V_%d = %12.6f\n", i, V_anterior[i]);
    printf("\n  Politica resultante:\n");
    printf("  [");
    for (int i = 0; i < n; i++) {
        printf("%d", politica[i] + 1);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");

    free(V_anterior);
    free(V_actual);
    free(politica);
}

/* ===================================================================
 * 5. SOLUCION POR PROGRAMACION LINEAL
 *
 * Formula el PMD como un problema de programacion lineal y lo resuelve
 * usando el metodo Simplex.
 *
 * Variables de decision: Y_{ik} (prob. estacionaria conjunta de
 *                         estado i y decision k)
 *
 * Minimizar: Z = sum_i sum_k C_{ik} * Y_{ik}
 *
 * Sujeto a:
 *   sum_i sum_k Y_{ik} = 1
 *   sum_k Y_{jk} - sum_i sum_k Y_{ik} * p_{ij}(k) = 0  (para j=0..m)
 *   Y_{ik} >= 0
 *
 * Al encontrar el optimo, transformar a politica deterministica:
 *   D_{ik} = Y_{ik} / sum_k Y_{ik}
 * =================================================================== */

/* Tabla simplex */
typedef struct {
    int     filas;
    int     columnas;
    double *costos;    /* funcion objetivo (costos reducidos) */
    double *rhs;       /* lado derecho */
    double **tabla;    /* matriz de restricciones */
    int    *base;      /* variables basicas */
    int    *no_base;   /* variables no basicas */
} TablaSimplex;

static TablaSimplex* simplex_crear(int filas, int columnas) {
    TablaSimplex *t = (TablaSimplex*)malloc(sizeof(TablaSimplex));
    t->filas    = filas;
    t->columnas = columnas;
    t->costos   = (double*)calloc((size_t)columnas, sizeof(double));
    t->rhs      = (double*)calloc((size_t)filas, sizeof(double));
    t->tabla    = (double**)malloc((size_t)filas * sizeof(double*));
    for (int i = 0; i < filas; i++)
        t->tabla[i] = (double*)calloc((size_t)columnas, sizeof(double));
    t->base    = (int*)malloc((size_t)filas * sizeof(int));
    t->no_base = (int*)malloc((size_t)columnas * sizeof(int));
    return t;
}

static void simplex_destruir(TablaSimplex *t) {
    if (!t) return;
    for (int i = 0; i < t->filas; i++) free(t->tabla[i]);
    free(t->tabla);
    free(t->costos);
    free(t->rhs);
    free(t->base);
    free(t->no_base);
    free(t);
}

static void simplex_resolver(TablaSimplex *t, double *solucion) {
    int m = t->filas;
    int n = t->columnas;

    /* Inicializar variables basicas y no basicas */
    for (int i = 0; i < n; i++)
        t->no_base[i] = i;
    for (int i = 0; i < m; i++)
        t->base[i] = n - m + i; /* ultimas m variables como basicas */
    for (int i = 0; i < m; i++)
        t->no_base[n - m + i] = -1;

    int iter = 0;
    int max_iter_simplex = 10000;

    while (iter < max_iter_simplex) {
        iter++;

        /* Encontrar variable entrante (costo reducido mas negativo) */
        int    entra = -1;
        double min_costo = 0.0;
        for (int j = 0; j < n; j++) {
            if (t->no_base[j] >= 0 && t->costos[j] < min_costo - 1e-12) {
                min_costo = t->costos[j];
                entra = j;
            }
        }

        if (entra < 0) break; /* optimo encontrado */

        /* Encontrar variable saliente (prueba de razon minima) */
        int    sale  = -1;
        double razon_min = INFINITY;
        for (int i = 0; i < m; i++) {
            if (t->tabla[i][entra] > 1e-12) {
                double razon = t->rhs[i] / t->tabla[i][entra];
                if (razon < razon_min - 1e-12) {
                    razon_min = razon;
                    sale = i;
                }
            }
        }

        if (sale < 0) {
            /* Problema no acotado */
            printf("  Advertencia: PL no acotada.\n");
            break;
        }

        /* Pivotear */
        double pivote = t->tabla[sale][entra];

        /* Normalizar fila pivote */
        for (int j = 0; j < n; j++)
            t->tabla[sale][j] /= pivote;
        t->rhs[sale] /= pivote;

        /* Actualizar otras filas */
        for (int i = 0; i < m; i++) {
            if (i == sale) continue;
            double factor = t->tabla[i][entra];
            if (fabs(factor) < 1e-15) continue;
            for (int j = 0; j < n; j++)
                t->tabla[i][j] -= factor * t->tabla[sale][j];
            t->rhs[i] -= factor * t->rhs[sale];
        }

        /* Actualizar costos reducidos */
        double factor_costo = t->costos[entra];
        if (fabs(factor_costo) > 1e-15) {
            for (int j = 0; j < n; j++)
                t->costos[j] -= factor_costo * t->tabla[sale][j];
        }

        /* Actualizar conjuntos basicos/no basicos */
        int var_sale = t->base[sale];
        t->base[sale]     = entra;
        t->no_base[entra] = -1; /* ahora es basica */
        if (var_sale >= 0)
            t->no_base[var_sale] = var_sale; /* ahora es no basica */
    }

    /* Extraer solucion */
    for (int j = 0; j < n; j++)
        solucion[j] = 0.0;
    for (int i = 0; i < m; i++)
        if (t->base[i] < n)
            solucion[t->base[i]] = t->rhs[i];
}

void programacion_lineal(const ModeloMarkov *modelo) {
    int n_est = modelo->num_estados;
    int K     = modelo->num_decisiones;

    printf("\n=============================================================\n");
    printf("  5. SOLUCION POR PROGRAMACION LINEAL\n");
    printf("=============================================================\n");

    /* Numero de variables Y_{ik}: n_est * K */
    int num_var = n_est * K;

    /* Numero de restricciones:
       - 1 restriccion de normalizacion: sum_i sum_k Y_{ik} = 1
       - n_est restricciones de balance: para cada j
       Total: 1 + n_est restricciones de igualdad.
       Pero el simplex necesita restricciones de desigualdad.
       Convertimos cada igualdad en dos desigualdades, o usamos variables
       de holgura. Para simplificar, convertimos en restricciones de
       desigualdad usando la formulacion estandar.

       Como las restricciones son igualdades, las reemplazamos por:
       sum_k Y_{jk} - sum_i sum_k Y_{ik} p_{ij}(k) >= 0  y  <= 0
       Es decir, >= 0 y -(...) >= 0.
    */

    int num_rest = 1 + 2 * n_est; /* dos desigualdades por cada balance */

    printf("\n  Formulacion del problema lineal:\n");
    printf("  Variables: Y_{ik} para i=0..%d, k=1..%d  (total: %d)\n",
           n_est - 1, K, num_var);
    printf("  Restricciones: %d\n\n", num_rest);

    /* Imprimir funcion objetivo */
    printf("  Funcion Objetivo (Minimizar):\n");
    printf("  Z = ");
    int primero = 1;
    for (int i = 0; i < n_est; i++) {
        for (int k = 0; k < K; k++) {
            if (!primero) printf(" + ");
            printf("%.4f * Y_%d,%d",
                   modelo->C->datos[i][k], i, k + 1);
            primero = 0;
        }
    }
    printf("\n\n");

    /* Imprimir restricciones */
    printf("  Restricciones:\n");
    printf("  (1) sum_i sum_k Y_{ik} = 1\n");
    for (int j = 0; j < n_est; j++) {
        printf("  (%d) sum_k Y_%d,k - sum_i sum_k Y_{ik} * P_i,%d(k) = 0\n",
               j + 2, j, j);
    }

    printf("\n  Resolviendo con metodo Simplex (Big M)...\n\n");

    /* ---------------------------------------------------------------
     * Construccion de la tabla simplex con metodo de la Gran M.
     *
     * Variables:
     *   Y_{ik} para i=0..n_est-1, k=0..K-1   (total: num_var)
     *   s_0: holgura de restriccion sum Y <= 1 (se usa la igualdad
     *        con variable artificial en su lugar)
     *   a  : variable artificial para sum Y = 1
     *   s_j: holgura para balance_j <= 0 (para j=0..n_est-1)
     *
     * Restricciones:
     *   Fila 0: sum Y_{ik} + a = 1            (igualdad, artificial)
     *   Fila 1+j: balance_j + s_j = 0         (desigualdad <= 0)
     *
     * Nota: solo una desigualdad por balance porque sum_j balance_j=0
     * implica que si cada una es <= 0, entonces cada una es = 0.
     * --------------------------------------------------------------- */

    /* Numero de restricciones: 1 (normalizacion) + n_est (balances) */
    int num_rest_simplex = 1 + n_est;

    /* Variables: num_var (Y) + 1 (artificial a) + n_est (holguras s_j) */
    int col_a = num_var;                       /* indice de variable artificial */
    int col_s_inicio = num_var + 1;              /* indices de holguras s_j */
    int total_columnas = num_var + 1 + n_est;

    #define IDX(i, k) ((i) * K + (k))
    #define M_GRANDE 1e10

    TablaSimplex *tabla = simplex_crear(num_rest_simplex, total_columnas);

    /* Funcion objetivo: min sum C_{ik}*Y_{ik} + M*a */
    for (int i = 0; i < n_est; i++)
        for (int k = 0; k < K; k++)
            tabla->costos[IDX(i, k)] = modelo->C->datos[i][k];
    tabla->costos[col_a] = M_GRANDE;  /* costo artificial */
    /* holguras tienen costo 0 */

    /* Fila 0: sum Y_{ik} + a = 1 */
    int fila = 0;
    for (int i = 0; i < n_est; i++)
        for (int k = 0; k < K; k++)
            tabla->tabla[fila][IDX(i, k)] = 1.0;
    tabla->tabla[fila][col_a] = 1.0;
    tabla->rhs[fila] = 1.0;
    fila++;

    /* Filas 1..n_est: balance_j <= 0  =>  balance_j + s_j = 0 */
    for (int j = 0; j < n_est; j++) {
        /*
         * balance_j = sum_k Y_{jk} - sum_i sum_k Y_{ik} * p_{ij}(k)
         * Coeficiente de Y_{ik}: (i==j ? 1 : 0) - p_{ij}(k)
         */
        for (int i = 0; i < n_est; i++) {
            for (int k = 0; k < K; k++) {
                double coef = -modelo->P_dec[k]->datos[i][j];
                if (i == j) coef += 1.0;
                tabla->tabla[fila][IDX(i, k)] = coef;
            }
        }
        /* Coeficiente de a (artificial) es 0 (no aparece en balance) */
        tabla->tabla[fila][col_a] = 0.0;
        /* Holgura s_j */
        tabla->tabla[fila][col_s_inicio + j] = 1.0;
        tabla->rhs[fila] = 0.0;
        fila++;
    }

    /* ============================================================
     * Inicializacion del simplex: variables basicas iniciales
     * La variable artificial a es basica en fila 0 (valor = 1)
     * Las holguras s_j son basicas en filas 1+j (valor = 0)
     *
     * Se deben ajustar los costos reducidos porque "a" tiene
     * costo M en la funcion objetivo y es basica.
     * Para cada columna j: costo[j] -= costo_basico[fila_a] * coef[fila_a][j]
     * =================================================================== */

    /* Fila 0 es donde esta la artificial "a" como basica */
    for (int j = 0; j < total_columnas; j++) {
        tabla->costos[j] -= M_GRANDE * tabla->tabla[0][j];
    }

    /* Inicializar conjuntos base/no_base */
    for (int i = 0; i < total_columnas; i++)
        tabla->no_base[i] = i;
    tabla->base[0] = col_a;       /* artificial en fila 0 */
    tabla->no_base[col_a] = -1;
    for (int j = 0; j < n_est; j++) {
        tabla->base[1 + j] = col_s_inicio + j;  /* holguras en filas 1+j */
        tabla->no_base[col_s_inicio + j] = -1;
    }

    double *solucion = (double*)calloc((size_t)total_columnas, sizeof(double));
    simplex_resolver(tabla, solucion);

    printf("=============================================================\n");
    printf("  RESULTADO - Programacion Lineal\n");
    printf("=============================================================\n");

    double valor_optimo = 0.0;
    printf("\n  Variables Y_{ik} optimas:\n");
    for (int i = 0; i < n_est; i++) {
        for (int k = 0; k < K; k++) {
            double val = solucion[IDX(i, k)];
            if (fabs(val) > 1e-8) {
                printf("    Y_%d,%d = %.8f\n", i, k + 1, val);
                valor_optimo += modelo->C->datos[i][k] * val;
            }
        }
    }

    printf("\n  Valor optimo Z = %.8f\n", valor_optimo);

    /* Transformar a politica deterministica:
       D_{ik} = Y_{ik} / sum_k Y_{ik}
       decision[i] = argmax_k D_{ik} */
    printf("\n  Politica deterministica resultante:\n");
    printf("  [");
    for (int i = 0; i < n_est; i++) {
        double suma_fila = 0.0;
        for (int k = 0; k < K; k++)
            suma_fila += solucion[IDX(i, k)];

        int mejor_k = 0;
        double mejor_valor = -1.0;
        for (int k = 0; k < K; k++) {
            double dik = (suma_fila > 1e-12) ?
                solucion[IDX(i, k)] / suma_fila : 0.0;
            if (dik > mejor_valor) {
                mejor_valor = dik;
                mejor_k = k;
            }
        }

        printf("%d", mejor_k + 1);
        if (i < n_est - 1) printf(", ");
    }
    printf("]\n");

    printf("\n  Probabilidades D_{ik}:\n");
    for (int i = 0; i < n_est; i++) {
        double suma_fila = 0.0;
        for (int k = 0; k < K; k++)
            suma_fila += solucion[IDX(i, k)];
        printf("    Estado %d: ", i);
        for (int k = 0; k < K; k++) {
            double dik = (suma_fila > 1e-12) ?
                solucion[IDX(i, k)] / suma_fila : 0.0;
            printf("D_%d,%d=%.4f  ", i, k + 1, dik);
        }
        printf("\n");
    }

    #undef IDX

    free(solucion);
    simplex_destruir(tabla);
}

/* ===================================================================
 * PRUEBA COMPLETA DE PMD (BATCH TEST)
 *
 * Ejecuta los 4 metodos de resolucion en secuencia, imprimiendo
 * todos los resultados de una sola vez. Util para verificar
 * consistencia entre metodos.
 *
 * Metodos ejecutados:
 *   1. Enumeracion Exhaustiva
 *   2. Mejoramiento de Politicas (sin descuento)
 *   3. Mejoramiento de Politicas con Descuento
 *   4. Programacion Lineal
 * =================================================================== */
void prueba_completa_pmd(ModeloMarkov *modelo) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;

    printf("\n");
    printf("################################################################\n");
    printf("#                                                              #\n");
    printf("#   PRUEBA COMPLETA DE PROCESOS MARKOVIANOS DE DECISION        #\n");
    printf("#                                                              #\n");
    printf("################################################################\n");
    printf("# Estados: %d   Decisiones: %d   alfa=%.4f\n",
           n, K, modelo->alfa);
    if (modelo->es_maximizacion)
        printf("# Objetivo: MAXIMIZAR (costos negados internamente)\n");
    else
        printf("# Objetivo: MINIMIZAR\n");
    printf("################################################################\n");

    /* ================================================================
     * 1. ENUMERACION EXHAUSTIVA
     * ================================================================ */
    enumeracion_exhaustiva(modelo);

    /* ================================================================
     * 2. MEJORAMIENTO DE POLITICAS (SIN DESCUENTO)
     * ================================================================ */
    printf("\n");
    printf("################################################################\n");
    printf("#   MEJORAMIENTO DE POLITICAS (Sin Descuento)                  #\n");
    printf("################################################################\n");
    Politica *R_ini = politica_crear(n, K);
    for (int i = 0; i < n; i++) R_ini->decision[i] = 0;
    printf("Politica inicial: ");
    politica_imprimir(R_ini);
    Politica *opt_sin_desc = mejoramiento_politicas(modelo, R_ini);
    if (opt_sin_desc) {
        printf("\n>>> Politica optima (sin descuento): ");
        politica_imprimir(opt_sin_desc);

        Matriz *Popt = politica_matriz_transicion(modelo, opt_sin_desc);
        double *pi_opt = estado_estable(Popt);
        if (pi_opt) {
            double valor_obj = 0.0;
            for (int i = 0; i < n; i++)
                valor_obj += modelo->C->datos[i][opt_sin_desc->decision[i]]
                           * pi_opt[i];
            if (modelo->es_maximizacion)
                printf("  Utilidad esperada = %.4f\n", -valor_obj);
            else
                printf("  Costo esperado = %.4f\n", valor_obj);
            free(pi_opt);
        }
        matriz_destruir(Popt);
        politica_destruir(opt_sin_desc);
    }
    politica_destruir(R_ini);

    /* ================================================================
     * 3. MEJORAMIENTO DE POLITICAS CON DESCUENTO
     * ================================================================ */
    printf("\n");
    printf("################################################################\n");
    printf("#   MEJORAMIENTO DE POLITICAS CON DESCUENTO (alfa=%.4f)        #\n",
           modelo->alfa);
    printf("################################################################\n");
    R_ini = politica_crear(n, K);
    for (int i = 0; i < n; i++) R_ini->decision[i] = 0;
    printf("Politica inicial: ");
    politica_imprimir(R_ini);
    Politica *opt_con_desc = mejoramiento_politicas_descuento(modelo, R_ini);
    if (opt_con_desc) {
        printf("\n>>> Politica optima (con descuento): ");
        politica_imprimir(opt_con_desc);

        Matriz *Popt = politica_matriz_transicion(modelo, opt_con_desc);
        double *pi_opt = estado_estable(Popt);
        if (pi_opt) {
            double valor_obj = 0.0;
            for (int i = 0; i < n; i++)
                valor_obj += modelo->C->datos[i][opt_con_desc->decision[i]]
                           * pi_opt[i];
            if (modelo->es_maximizacion)
                printf("  Utilidad esperada a LP = %.4f\n", -valor_obj);
            else
                printf("  Costo esperado a LP = %.4f\n", valor_obj);
            free(pi_opt);
        }
        matriz_destruir(Popt);
        politica_destruir(opt_con_desc);
    }
    politica_destruir(R_ini);

    /* ================================================================
     * 4. PROGRAMACION LINEAL
     * ================================================================ */
    printf("\n");
    printf("################################################################\n");
    printf("#   PROGRAMACION LINEAL                                        #\n");
    printf("################################################################\n");
    programacion_lineal(modelo);

    /* ================================================================
     * RESUMEN FINAL
     * ================================================================ */
    printf("\n");
    printf("################################################################\n");
    printf("#                                                              #\n");
    printf("#   PRUEBA COMPLETA FINALIZADA                                 #\n");
    printf("#   Todos los metodos ejecutados exitosamente.                 #\n");
    printf("#   Compare los resultados para verificar consistencia.        #\n");
    printf("#                                                              #\n");
    printf("################################################################\n\n");
}
