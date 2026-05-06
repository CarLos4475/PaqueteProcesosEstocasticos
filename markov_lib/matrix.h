#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ===================================================================
 * TDA Matriz Dinamica
 * Representa una matriz de doubles con dimensiones filas x columnas.
 * =================================================================== */
typedef struct {
    int filas;
    int columnas;
    double **datos;
} Matriz;

/* Crear y destruir */
Matriz* matriz_crear(int filas, int columnas);
void    matriz_destruir(Matriz *m);
Matriz* matriz_copiar(const Matriz *origen);
Matriz* matriz_identidad(int n);

/* Operaciones aritmeticas */
Matriz* matriz_multiplicar(const Matriz *A, const Matriz *B);
Matriz* matriz_potencia(const Matriz *P, int n);
double* matriz_vector_multiplicar(const Matriz *A, const double *v);
Matriz* matriz_transponer(const Matriz *A);
double  vector_producto_punto(const double *a, const double *b, int n);

/* Impresion */
void matriz_imprimir(const Matriz *m, const char *titulo);
void vector_imprimir(const double *v, int n, const char *titulo);

/* Resolver sistema lineal Ax = b mediante eliminacion gaussiana
   con pivoteo parcial. A es n x n, b es vector de tamano n,
   x es vector resultado de tamano n (debe estar preasignado).
   Retorna 1 si tuvo exito, 0 si el sistema es singular. */
int resolver_sistema_lineal(const Matriz *A, const double *b, double *x);

#endif
