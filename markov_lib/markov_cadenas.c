#include "markov.h"
#include <string.h>

/* ===================================================================
 * MODULO 1: TEORIA BASICA DE CADENAS DE MARKOV
 *
 * Implementa los calculos teoricos fundamentales de cadenas de Markov
 * que serviran como soporte para los metodos de PMD.
 * =================================================================== */

/* ===================================================================
 * 1. ECUACIONES DE CHAPMAN-KOLMOGOROV
 *
 * p_{ij}^{(n)} = sum_{k=0}^{m} p_{ik}^{(l)} * p_{kj}^{(n-l)}
 *
 * Equivalentemente: P^(n) = P^n  (potencia n-esima de P)
 *
 * Parametros:
 *   P : matriz de transicion de un paso
 *   n : numero de pasos
 *
 * Retorna: matriz P^(n) (las transiciones en n pasos)
 * =================================================================== */
Matriz* chapman_kolmogorov(const Matriz *P, int n) {
    return matriz_potencia(P, n);
}

/* ===================================================================
 * 2. PROBABILIDADES INCONDICIONALES
 *
 * P(X_n = j) = a_0 * p_{0j}^{(n)} + a_1 * p_{1j}^{(n)} + ... + a_m * p_{mj}^{(n)}
 *            = (a * P^n)[j]
 *
 * Es decir, multiplicamos el vector de probabilidad inicial a
 * por la matriz de transicion en n pasos P^n.
 *
 * Parametros:
 *   a : vector de probabilidad inicial (tamano num_estados)
 *   P : matriz de transicion de un paso
 *   n : numero de pasos
 *
 * Retorna: vector de probabilidades incondicionales en el paso n
 * =================================================================== */
double* probabilidades_incondicionales(const double *a, const Matriz *P, int n) {
    Matriz *Pn = matriz_potencia(P, n);
    if (!Pn) return NULL;

    double *prob = matriz_vector_multiplicar(Pn, a);
    /* Nota: matriz_vector_multiplicar hace A * v (como vector columna)
       pero a es vector fila. Para a * P^n necesitariamos transponer.
       Corregimos: (a * M)[j] = sum_i a_i * M[i][j] = (M^T * a)[j]
       Por lo tanto, hacemos M^T * a donde M = P^n */

    /* Reemplazamos con calculo correcto: */
    free(prob);
    int n_estados = P->filas;
    prob = (double*)calloc((size_t)n_estados, sizeof(double));
    for (int j = 0; j < n_estados; j++)
        for (int i = 0; i < n_estados; i++)
            prob[j] += a[i] * Pn->datos[i][j];

    matriz_destruir(Pn);
    return prob;
}

/* ===================================================================
 * 3. VECTOR DE ESTADO ESTABLE (pi) Y TIEMPOS DE RECURRENCIA
 *
 * El estado estable satisface: pi = pi * P,  sum(pi_i) = 1
 * Es decir: pi * (P - I) = 0  con la restriccion sum(pi) = 1.
 *
 * Esto equivale a resolver el sistema:
 *   (P^T - I) * pi^T = 0
 * sustituyendo la ultima ecuacion por sum(pi_i) = 1.
 *
 * Tiempos de recurrencia: mu_{ii} = 1 / pi_i
 *
 * Parametros:
 *   P : matriz de transicion
 *
 * Retorna: vector pi de tamano num_estados con las prob. estacionarias
 * =================================================================== */
double* estado_estable(const Matriz *P) {
    int n = P->filas;

    /* Construir sistema: (P^T - I) pero reemplazar ultima fila por 1's */
    Matriz *A = matriz_crear(n, n);
    double *b = (double*)calloc((size_t)n, sizeof(double));
    if (!A || !b) { matriz_destruir(A); free(b); return NULL; }

    /* A = P^T, luego A = A - I, y reemplazar ultima fila */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A->datos[i][j] = P->datos[j][i];  /* A = P^T */
        }
    }

    /* A = P^T - I (restar identidad) */
    for (int i = 0; i < n; i++)
        A->datos[i][i] -= 1.0;

    /* Reemplazar ultima fila por la restriccion sum(pi) = 1 */
    for (int j = 0; j < n; j++)
        A->datos[n-1][j] = 1.0;
    b[n-1] = 1.0;

    double *pi = (double*)calloc((size_t)n, sizeof(double));
    if (!pi) { matriz_destruir(A); free(b); return NULL; }

    if (!resolver_sistema_lineal(A, b, pi)) {
        fprintf(stderr, "Error: No se pudo calcular el estado estable "
                "(matriz singular).\n");
        free(pi);
        pi = NULL;
    }

    /* Normalizar para asegurar que sumen 1 */
    if (pi) {
        double suma = 0.0;
        for (int i = 0; i < n; i++) suma += pi[i];
        if (fabs(suma) > 1e-12)
            for (int i = 0; i < n; i++) pi[i] /= suma;
    }

    matriz_destruir(A);
    free(b);
    return pi;
}

/* ===================================================================
 * 4. TIEMPOS DE RECURRENCIA
 *
 * mu_{ii} = 1 / pi_i   para cada estado i.
 *
 * Parametros:
 *   pi : vector de estado estable
 *   n  : numero de estados
 *
 * Retorna: vector mu donde mu[i] = 1 / pi[i]
 * =================================================================== */
double* tiempos_recurrencia(const double *pi, int n) {
    double *mu = (double*)calloc((size_t)n, sizeof(double));
    if (!mu) return NULL;
    for (int i = 0; i < n; i++) {
        if (fabs(pi[i]) > 1e-15)
            mu[i] = 1.0 / pi[i];
        else
            mu[i] = INFINITY; /* estado no recurrente o con prob. cero */
    }
    return mu;
}

/* ===================================================================
 * 5. TIEMPOS DE PRIMERA PASADA (TRANSICION)
 *
 * mu_{ij} = 1 + sum_{k != j} p_{ik} * mu_{kj}
 *
 * Para cada estado objetivo j, resolvemos un sistema lineal de tamano
 * (n-1) para las incognitas mu_{ij} con i != j.
 *
 * Parametros:
 *   P : matriz de transicion
 *
 * Retorna: matriz M de (num_estados) x (num_estados) donde
 *          M[i][j] = mu_{ij} (tiempo esperado para ir de i a j)
 *          M[i][i] = mu_{ii} (tiempo de recurrencia)
 * =================================================================== */
Matriz* tiempos_primera_pasada(const Matriz *P) {
    int n = P->filas;
    Matriz *M = matriz_crear(n, n);
    if (!M) return NULL;

    /* Para cada estado objetivo j, calcular mu_{ij} para todo i != j */
    for (int obj = 0; obj < n; obj++) {
        int tam = n - 1; /* numero de incognitas (todas las i != obj) */

        Matriz *A = matriz_crear(tam, tam);
        double *b = (double*)calloc((size_t)tam, sizeof(double));
        double *u = (double*)calloc((size_t)tam, sizeof(double));
        if (!A || !b || !u) {
            matriz_destruir(A); free(b); free(u); matriz_destruir(M);
            return NULL;
        }

        /* Construir el sistema para el objetivo obj */
        int fila = 0;
        for (int i = 0; i < n; i++) {
            if (i == obj) continue; /* saltar el estado objetivo */

            /* Ecuacion: mu_{i,obj} - sum_{k!=obj} p_{ik} * mu_{k,obj} = 1 */
            int col = 0;
            for (int k = 0; k < n; k++) {
                if (k == obj) continue; /* saltar en la suma */
                if (i == k)
                    A->datos[fila][col] = 1.0 - P->datos[i][k];
                else
                    A->datos[fila][col] = -P->datos[i][k];
                col++;
            }
            b[fila] = 1.0;
            fila++;
        }

        if (!resolver_sistema_lineal(A, b, u)) {
            /* Si no se puede resolver, rellenar con INFINITY */
            fila = 0;
            for (int i = 0; i < n; i++) {
                if (i == obj) {
                    M->datos[i][obj] = 0.0; /* tiempo de i a si mismo es 0 */
                } else {
                    M->datos[i][obj] = INFINITY;
                    fila++;
                }
            }
        } else {
            /* Mapear solucion de vuelta a M */
            fila = 0;
            for (int i = 0; i < n; i++) {
                if (i == obj) {
                    /* Calcular tiempo de recurrencia: mu_{jj} */
                    double suma = 0.0;
                    for (int k = 0; k < n; k++)
                        if (k != obj)
                            suma += P->datos[obj][k] * u[
                                (k > obj) ? k-1 : k];
                    M->datos[obj][obj] = 1.0 + suma;
                } else {
                    M->datos[i][obj] = u[fila];
                    fila++;
                }
            }
        }

        matriz_destruir(A);
        free(b);
        free(u);
    }

    return M;
}

/* ===================================================================
 * 6. PROBABILIDADES DE ABSORCION
 *
 * f_{ik} = sum_{j=0}^{m} p_{ij} * f_{jk}
 *
 * Condiciones de frontera:
 *   f_{kk} = 1  (si empezamos en k, ya estamos absorbidos en k)
 *   f_{ik} = 0  si i es absorbente y distinto de k
 *
 * Para estados no absorbentes i != k:
 *   f_{ik} = p_{ik} * 1 + sum_{j no absorbente, j != k} p_{ij} * f_{jk}
 *   f_{ik} - sum_{j en T} p_{ij} * f_{jk} = p_{ik}
 *
 * donde T = {j : j no es absorbente, j != k}
 *
 * Parametros:
 *   P           : matriz de transicion
 *   absorbentes : arreglo booleano (1 si el estado es absorbente)
 *   k           : indice del estado absorbente objetivo
 *
 * Retorna: vector f de tamano num_estados con las probabilidades f_{ik}
 * =================================================================== */
double* probabilidades_absorcion(const Matriz *P, const int *absorbentes, int k) {
    int n = P->filas;
    double *f = (double*)calloc((size_t)n, sizeof(double));
    if (!f) return NULL;

    /* Contar estados transientes (no absorbentes) excluyendo k */
    int *indice = (int*)malloc((size_t)n * sizeof(int));
    int *a_nuevo = (int*)malloc((size_t)n * sizeof(int)); /* mapeo inverso */
    int num_trans = 0;
    for (int i = 0; i < n; i++) {
        if (i == k) {
            f[i] = 1.0;
            indice[i] = -1;
        } else if (absorbentes[i]) {
            f[i] = 0.0;
            indice[i] = -1;
        } else {
            indice[i] = num_trans;
            a_nuevo[num_trans] = i;
            num_trans++;
        }
    }

    if (num_trans == 0) {
        /* Todos los estados son absorbentes o k ya esta determinado */
        free(indice);
        free(a_nuevo);
        return f;
    }

    /* Construir sistema para estados transientes */
    Matriz *A = matriz_crear(num_trans, num_trans);
    double *b = (double*)calloc((size_t)num_trans, sizeof(double));
    double *u = (double*)calloc((size_t)num_trans, sizeof(double));
    if (!A || !b || !u) {
        matriz_destruir(A); free(b); free(u);
        free(indice); free(a_nuevo); free(f);
        return NULL;
    }

    for (int r = 0; r < num_trans; r++) {
        int i = a_nuevo[r];
        for (int c = 0; c < num_trans; c++) {
            int j = a_nuevo[c];
            if (i == j)
                A->datos[r][c] = 1.0 - P->datos[i][j];
            else
                A->datos[r][c] = -P->datos[i][j];
        }
        b[r] = P->datos[i][k]; /* p_{ik} * 1 */
    }

    if (!resolver_sistema_lineal(A, b, u)) {
        fprintf(stderr, "Error: No se pudo calcular las probabilidades "
                "de absorcion.\n");
        free(u); u = NULL;
    }

    /* Asignar resultados */
    for (int r = 0; r < num_trans; r++)
        f[a_nuevo[r]] = u ? u[r] : INFINITY;

    matriz_destruir(A);
    free(b);
    free(u);
    free(indice);
    free(a_nuevo);
    return f;
}
