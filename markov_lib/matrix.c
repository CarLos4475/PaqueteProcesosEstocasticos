#include "matrix.h"
#include <string.h>

/* ===================================================================
 * Crear una matriz de filas x columnas inicializada con ceros.
 * =================================================================== */
Matriz* matriz_crear(int filas, int columnas) {
    Matriz *m = (Matriz*)malloc(sizeof(Matriz));
    if (!m) return NULL;

    m->filas    = filas;
    m->columnas = columnas;

    m->datos = (double**)malloc((size_t)filas * sizeof(double*));
    if (!m->datos) { free(m); return NULL; }

    for (int i = 0; i < filas; i++) {
        m->datos[i] = (double*)calloc((size_t)columnas, sizeof(double));
        if (!m->datos[i]) {
            for (int k = 0; k < i; k++) free(m->datos[k]);
            free(m->datos);
            free(m);
            return NULL;
        }
    }
    return m;
}

/* ===================================================================
 * Liberar la memoria de una matriz.
 * =================================================================== */
void matriz_destruir(Matriz *m) {
    if (!m) return;
    for (int i = 0; i < m->filas; i++) free(m->datos[i]);
    free(m->datos);
    free(m);
}

/* ===================================================================
 * Crear una copia exacta de una matriz.
 * =================================================================== */
Matriz* matriz_copiar(const Matriz *origen) {
    Matriz *copia = matriz_crear(origen->filas, origen->columnas);
    if (!copia) return NULL;
    for (int i = 0; i < origen->filas; i++)
        for (int j = 0; j < origen->columnas; j++)
            copia->datos[i][j] = origen->datos[i][j];
    return copia;
}

/* ===================================================================
 * Crear una matriz identidad de tamano n x n.
 * =================================================================== */
Matriz* matriz_identidad(int n) {
    Matriz *I = matriz_crear(n, n);
    if (!I) return NULL;
    for (int i = 0; i < n; i++)
        I->datos[i][i] = 1.0;
    return I;
}

/* ===================================================================
 * Multiplicar dos matrices: A (m x n) * B (n x p) = C (m x p).
 * Retorna NULL si las dimensiones no coinciden.
 * =================================================================== */
Matriz* matriz_multiplicar(const Matriz *A, const Matriz *B) {
    if (A->columnas != B->filas) return NULL;

    Matriz *C = matriz_crear(A->filas, B->columnas);
    if (!C) return NULL;

    for (int i = 0; i < A->filas; i++)
        for (int k = 0; k < A->columnas; k++)
            if (fabs(A->datos[i][k]) > 1e-15)
                for (int j = 0; j < B->columnas; j++)
                    C->datos[i][j] += A->datos[i][k] * B->datos[k][j];

    return C;
}

/* ===================================================================
 * Calcular la potencia n-esima de una matriz P: P^n.
 * Usa exponenciacion binaria para eficiencia.
 * =================================================================== */
Matriz* matriz_potencia(const Matriz *P, int n) {
    if (P->filas != P->columnas) return NULL;

    Matriz *resultado = matriz_identidad(P->filas);  /* P^0 = I */
    Matriz *base      = matriz_copiar(P);

    int exp = n;
    while (exp > 0) {
        if (exp & 1) {
            Matriz *temp = matriz_multiplicar(resultado, base);
            matriz_destruir(resultado);
            resultado = temp;
        }
        exp >>= 1;
        if (exp > 0) {
            Matriz *temp = matriz_multiplicar(base, base);
            matriz_destruir(base);
            base = temp;
        }
    }
    matriz_destruir(base);
    return resultado;
}

/* ===================================================================
 * Multiplicar matriz A (m x n) por vector v (tamano n).
 * Retorna vector resultado de tamano m.
 * La memoria del vector devuelto debe ser liberada con free().
 * =================================================================== */
double* matriz_vector_multiplicar(const Matriz *A, const double *v) {
    double *resultado = (double*)calloc((size_t)A->filas, sizeof(double));
    if (!resultado) return NULL;

    for (int i = 0; i < A->filas; i++)
        for (int j = 0; j < A->columnas; j++)
            resultado[i] += A->datos[i][j] * v[j];

    return resultado;
}

/* ===================================================================
 * Transponer una matriz: A^T.
 * =================================================================== */
Matriz* matriz_transponer(const Matriz *A) {
    Matriz *T = matriz_crear(A->columnas, A->filas);
    if (!T) return NULL;
    for (int i = 0; i < A->filas; i++)
        for (int j = 0; j < A->columnas; j++)
            T->datos[j][i] = A->datos[i][j];
    return T;
}

/* ===================================================================
 * Producto punto de dos vectores de tamano n.
 * =================================================================== */
double vector_producto_punto(const double *a, const double *b, int n) {
    double suma = 0.0;
    for (int i = 0; i < n; i++)
        suma += a[i] * b[i];
    return suma;
}

/* ===================================================================
 * Imprimir una matriz en consola con formato.
 * =================================================================== */
void matriz_imprimir(const Matriz *m, const char *titulo) {
    if (titulo) printf("\n%s [%d x %d]:\n", titulo, m->filas, m->columnas);
    for (int i = 0; i < m->filas; i++) {
        for (int j = 0; j < m->columnas; j++)
            printf("%10.6f ", m->datos[i][j]);
        printf("\n");
    }
}

/* ===================================================================
 * Imprimir un vector en consola.
 * =================================================================== */
void vector_imprimir(const double *v, int n, const char *titulo) {
    if (titulo) printf("\n%s [%d]:\n", titulo, n);
    for (int i = 0; i < n; i++)
        printf("  [%d] = %12.8f\n", i, v[i]);
}

/* ===================================================================
 * Resolver sistema lineal Ax = b mediante eliminacion gaussiana
 * con pivoteo parcial.
 *
 * Parametros:
 *   A : matriz de coeficientes n x n (se modifica internamente)
 *   b : vector de terminos independientes (se modifica internamente)
 *   x : vector solucion de tamano n (resultado de salida)
 *
 * Retorna:
 *   1 si el sistema fue resuelto exitosamente.
 *   0 si la matriz es singular.
 * =================================================================== */
int resolver_sistema_lineal(const Matriz *A_orig, const double *b_orig, double *x) {
    int n = A_orig->filas;
    if (n != A_orig->columnas) return 0;

    /* Crear copia local de A y b para no modificar los originales */
    Matriz *A = matriz_copiar(A_orig);
    double *b = (double*)malloc((size_t)n * sizeof(double));
    if (!A || !b) { matriz_destruir(A); free(b); return 0; }
    memcpy(b, b_orig, (size_t)n * sizeof(double));

    /* Eliminacion hacia adelante con pivoteo parcial */
    for (int col = 0; col < n; col++) {
        /* Buscar pivote (maximo absoluto en esta columna) */
        int    fila_pivote = col;
        double val_max     = fabs(A->datos[col][col]);
        for (int f = col + 1; f < n; f++) {
            if (fabs(A->datos[f][col]) > val_max) {
                val_max     = fabs(A->datos[f][col]);
                fila_pivote = f;
            }
        }

        /* Si el pivote es casi cero, el sistema es singular */
        if (val_max < 1e-12) {
            matriz_destruir(A); free(b);
            return 0;
        }

        /* Intercambiar filas si es necesario */
        if (fila_pivote != col) {
            double *temp_fila    = A->datos[col];
            A->datos[col]        = A->datos[fila_pivote];
            A->datos[fila_pivote] = temp_fila;

            double temp_b = b[col];
            b[col]        = b[fila_pivote];
            b[fila_pivote] = temp_b;
        }

        /* Eliminar filas inferiores */
        for (int f = col + 1; f < n; f++) {
            double factor = A->datos[f][col] / A->datos[col][col];
            A->datos[f][col] = 0.0;
            for (int c = col + 1; c < n; c++)
                A->datos[f][c] -= factor * A->datos[col][c];
            b[f] -= factor * b[col];
        }
    }

    /* Sustitucion hacia atras */
    for (int i = n - 1; i >= 0; i--) {
        double suma = b[i];
        for (int j = i + 1; j < n; j++)
            suma -= A->datos[i][j] * x[j];
        x[i] = suma / A->datos[i][i];
    }

    matriz_destruir(A);
    free(b);
    return 1;
}
