#ifndef MARKOV_H
#define MARKOV_H

#include "matrix.h"

/* ===================================================================
 * Estructuras de Datos Principales
 * =================================================================== */

/* Representa una politica deterministica:
   decision[i] = indice de la decision (0..K-1) para el estado i */
typedef struct {
    int  num_estados;
    int  num_decisiones;
    int *decision;       /* decision[i] en {0, 1, ..., K-1} */
} Politica;

/* Modelo principal que almacena toda la informacion del sistema */
typedef struct {
    int     num_estados;        /* m + 1  (estados: 0, 1, ..., m) */
    int     num_decisiones;     /* K      (decisiones: 1, 2, ..., K) */
    int     tiene_decisiones;   /* 0 = cadena simple, 1 = PMD */

    double *prob_inicial;       /* vector a de tamano m+1 */

    /* Para cadenas de Markov simples (sin decisiones) */
    Matriz *P;                  /* matriz de transicion (m+1) x (m+1) */

    /* Para Procesos Markovianos de Decision (PMD) */
    Matriz **P_dec;             /* arreglo de K matrices de transicion,
                                   cada una (m+1) x (m+1) */
    Matriz  *C;                 /* matriz de costos (m+1) x K */
    double   alfa;              /* factor de descuento alfa */
    double   tasa_interes;      /* tasa de interes i (para calcular alfa) */

    /* ---- Nuevo: soporte para ingresos y maximizacion ---- */
    int      es_maximizacion;   /* 1 = maximizar utilidad, 0 = minimizar costo */
    int      usa_ingresos;      /* 1 = C se genero desde matrices de ingreso */
    Matriz **Ingreso;           /* arreglo de K matrices de ingreso (m+1)x(m+1) */
    double  *costo_fijo;        /* costo fijo por decision (tamano K) */
    char    *nombres_estados;   /* nombres opcionales de estados (arreglo de char*) */
    char    *nombres_decisiones;/* nombres opcionales de decisiones (arreglo de char*) */
} ModeloMarkov;

/* ===================================================================
 * Funciones de Lectura de Datos (markov_datos.c)
 * =================================================================== */
ModeloMarkov* modelo_crear(void);
void          modelo_destruir(ModeloMarkov *modelo);
ModeloMarkov* modelo_leer_consola(void);
ModeloMarkov* modelo_leer_archivo(const char *nombre_archivo);
void          modelo_imprimir(const ModeloMarkov *modelo);

/* Calcula C_ik a partir de matrices de ingreso y costos fijos:
   C_ik = sum_j P_ij(k) * Ingreso_ij(k) - costo_fijo[k]
   Si es_maximizacion=1, los valores se niegan (min -utilidad = max utilidad) */
void modelo_calcular_costos_ingresos(ModeloMarkov *modelo);

/* ===================================================================
 * Funciones del Modulo 1: Teoria Basica de Cadenas de Markov
 * (markov_cadenas.c)
 * =================================================================== */

/* 1. Ecuaciones de Chapman-Kolmogorov: P^(n) = P^n */
Matriz* chapman_kolmogorov(const Matriz *P, int n);

/* 2. Probabilidades Incondicionales: P(X_n = j) = a * P^n */
double* probabilidades_incondicionales(const double *a, const Matriz *P, int n);

/* 3. Vector de Estado Estable: pi = pi * P, sum(pi) = 1 */
double* estado_estable(const Matriz *P);

/* 4. Tiempos de Recurrencia: mu_ii = 1/pi_i */
double* tiempos_recurrencia(const double *pi, int n);

/* 5. Tiempos de Primera Pasada: mu_ij para cada i != j objetivo.
   Retorna una matriz de (num_estados) x (num_estados) con mu_ij. */
Matriz* tiempos_primera_pasada(const Matriz *P);

/* 6. Probabilidades de Absorcion hacia el estado absorbente k.
   absorbentes es un arreglo booleano de tamano num_estados
   (1 = absorbente, 0 = no absorbente). k es el indice del estado objetivo.
   Retorna vector f de tamano num_estados con f_ik. */
double* probabilidades_absorcion(const Matriz *P, const int *absorbentes, int k);

/* ===================================================================
 * Funciones del Modulo 2: Procesos Markovianos de Decision
 * (markov_decision.c)
 * =================================================================== */

/* Politica */
Politica* politica_crear(int num_estados, int num_decisiones);
void      politica_destruir(Politica *p);
void      politica_imprimir(const Politica *p);
Politica* politica_copiar(const Politica *origen);
int       politica_igual(const Politica *a, const Politica *b);
void      politica_leer_consola(ModeloMarkov *modelo, Politica **politicas,
                                int *num_politicas);

/* Construye la matriz de transicion bajo una politica */
Matriz* politica_matriz_transicion(const ModeloMarkov *modelo, const Politica *p);

/* 1. Enumeracion Exhaustiva de Politicas */
void enumeracion_exhaustiva(ModeloMarkov *modelo);

/* 2. Mejoramiento de Politicas (sin descuento) */
Politica* mejoramiento_politicas(const ModeloMarkov *modelo, const Politica *inicial);

/* 3. Mejoramiento de Politicas con Descuento */
Politica* mejoramiento_politicas_descuento(const ModeloMarkov *modelo,
                                           const Politica *inicial);

/* 4. Metodo de Aproximaciones Sucesivas */
void aproximaciones_sucesivas(const ModeloMarkov *modelo,
                              int max_iter, double epsilon);

/* 5. Solucion por Programacion Lineal */
void programacion_lineal(const ModeloMarkov *modelo);

/* ===================================================================
 * Prueba completa (batch test) — Ejecuta los 4 metodos de PMD
 * de una sola vez e imprime todos los resultados.
 * (markov_decision.c)
 * =================================================================== */
void prueba_completa_pmd(ModeloMarkov *modelo);

#endif
