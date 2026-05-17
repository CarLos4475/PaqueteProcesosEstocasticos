#include "markov.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>

/* ===================================================================
 * Programa principal: Menu interactivo para Cadenas de Markov y PMD.
 *
 * Permite al usuario cargar un modelo y ejecutar cualquiera de los
 * metodos implementados.
 * =================================================================== */

/* Limpiar pantalla con portabilidad Windows/Linux */
static void limpiar_pantalla(void) {
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

static void menu_principal(void) {
    imprimir_titulo_box("PAQUETE ESTOCASTICO - CADENAS DE MARKOV Y PMD");
    printf("\n  MENU PRINCIPAL\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  | 1  | Cargar modelo desde consola                    |\n");
    printf("  | 2  | Cargar modelo desde archivo                    |\n");
    printf("  | 3  | Mostrar modelo cargado                         |\n");
    printf("  | 4  | Modulo 1: Teoria basica de Cadenas de Markov   |\n");
    printf("  | 5  | Modulo 2: Procesos Markovianos de Decision     |\n");
    printf("  | 0  | Salir                                          |\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  Opcion: ");
}

static void menu_modulo1(void) {
    printf("\n  MODULO 1 - TEORIA BASICA DE CADENAS DE MARKOV\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  | 1  | Ecuaciones de Chapman-Kolmogorov (P^n)         |\n");
    printf("  | 2  | Probabilidades Incondicionales                 |\n");
    printf("  | 3  | Vector de Estado Estable                       |\n");
    printf("  | 4  | Tiempos de Recurrencia                         |\n");
    printf("  | 5  | Tiempos de Primera Pasada                      |\n");
    printf("  | 6  | Probabilidades de Absorcion                    |\n");
    printf("  | 0  | Volver al menu principal                       |\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  Opcion: ");
}

static void menu_modulo2(void) {
    printf("\n  MODULO 2 - PROCESOS MARKOVIANOS DE DECISION\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  | 1  | Enumeracion Exhaustiva de Politicas            |\n");
    printf("  | 2  | Mejoramiento de Politicas (sin descuento)      |\n");
    printf("  | 3  | Mejoramiento de Politicas con Descuento        |\n");
    printf("  | 4  | Metodo de Aproximaciones Sucesivas             |\n");
    printf("  | 5  | Solucion por Programacion Lineal               |\n");
    printf("  | 6  | PRUEBA COMPLETA (ejecuta 1,2,3,5 en secuencia) |\n");
    printf("  | 0  | Volver al menu principal                       |\n");
    printf("  +----+------------------------------------------------+\n");
    printf("  Opcion: ");
}

static void ejecutar_modulo1(ModeloMarkov *modelo) {
    int n = modelo->num_estados;
    if (!modelo->P && !modelo->tiene_decisiones) {
        printf("  Error: El modelo no tiene matriz de transicion.\n");
        return;
    }

    /* Para modulo 1, si es PMD usamos la primera matriz de decision */
    Matriz *P_usar = modelo->P;
    if (modelo->tiene_decisiones) {
        P_usar = modelo->P_dec[0];
        printf("  Nota: Usando P^(1) como matriz de transicion.\n");
    }

    int opcion;
    do {
        limpiar_pantalla();
        menu_modulo1();
        if (scanf("%d", &opcion) != 1) opcion = -1;
        while (getchar() != '\n');

        switch (opcion) {
        case 1: {
            int pasos;
            printf("  Ingrese el numero de pasos n: ");
            if (scanf("%d", &pasos) == 1 && pasos >= 0) {
                while (getchar() != '\n');
                Matriz *Pn = chapman_kolmogorov(P_usar, pasos);
                if (Pn) {
                    char titulo[64];
                    snprintf(titulo, sizeof(titulo),
                             "P^(%d) - Transiciones en %d pasos", pasos, pasos);
                    matriz_imprimir(Pn, titulo);
                    matriz_destruir(Pn);
                }
            } else { while (getchar() != '\n'); }
            break;
        }
        case 2: {
            int pasos;
            printf("  Ingrese el numero de pasos n: ");
            if (scanf("%d", &pasos) == 1 && pasos >= 0) {
                while (getchar() != '\n');
                double *prob = probabilidades_incondicionales(
                    modelo->prob_inicial, P_usar, pasos);
                if (prob) {
                    char titulo[64];
                    snprintf(titulo, sizeof(titulo),
                             "P(X_%d = j) - Probabilidades Incondicionales", pasos);
                    vector_imprimir(prob, n, titulo);
                    free(prob);
                }
            } else { while (getchar() != '\n'); }
            break;
        }
        case 3: {
            double *pi = estado_estable(P_usar);
            if (pi) {
                vector_imprimir(pi, n, "pi - Vector de Estado Estable");
                free(pi);
            }
            break;
        }
        case 4: {
            double *pi = estado_estable(P_usar);
            if (pi) {
                double *mu = tiempos_recurrencia(pi, n);
                if (mu) {
                    vector_imprimir(mu, n, "mu_ii - Tiempos de Recurrencia");
                    free(mu);
                }
                free(pi);
            }
            break;
        }
        case 5: {
            Matriz *M = tiempos_primera_pasada(P_usar);
            if (M) {
                matriz_imprimir(M, "M - Matriz de Tiempos de Primera Pasada");
                matriz_destruir(M);
            }
            break;
        }
        case 6: {
            /* Preguntar cuales estados son absorbentes */
            int *absorbentes = (int*)calloc((size_t)n, sizeof(int));
            printf("  Estados absorbentes (E = {0..%d}):\n", n - 1);
            for (int i = 0; i < n; i++) {
                char resp;
                printf("    El estado %d es absorbente? (S/N): ", i);
                if (scanf(" %c", &resp) == 1 && (resp == 'S' || resp == 's'))
                    absorbentes[i] = 1;
                else
                    absorbentes[i] = 0;
                while (getchar() != '\n');
            }
            int k;
            printf("  Ingrese el indice del estado absorbente objetivo: ");
            if (scanf("%d", &k) == 1 && k >= 0 && k < n) {
                while (getchar() != '\n');
                double *f = probabilidades_absorcion(P_usar, absorbentes, k);
                if (f) {
                    char titulo[64];
                    snprintf(titulo, sizeof(titulo),
                             "f_{i,%d} - Probabilidades de Absorcion", k);
                    vector_imprimir(f, n, titulo);
                    free(f);
                }
            } else { while (getchar() != '\n'); }
            free(absorbentes);
            break;
        }
        case 0: break;
        default: printf("  Opcion invalida.\n");
        }
    } while (opcion != 0);
}

static void ejecutar_modulo2(ModeloMarkov *modelo) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;

    if (!modelo->tiene_decisiones) {
        printf("  Error: El modelo actual es una cadena de Markov simple.\n");
        printf("  Se requiere un PMD con decisiones.\n");
        return;
    }

    int opcion;
    do {
        limpiar_pantalla();
        menu_modulo2();
        if (scanf("%d", &opcion) != 1) opcion = -1;
        while (getchar() != '\n');

        switch (opcion) {
        case 1:
            enumeracion_exhaustiva(modelo);
            break;
        case 2: {
            /* Preguntar al usuario si quiere ingresar politicas o usar
               una politica inicial por defecto */
            printf("\n  Elegir politica inicial para Mejoramiento de Politicas:\n");
            printf("    1. Ingresar politicas manualmente\n");
            printf("    2. Usar politica por defecto (decision 1 para todo)\n");
            printf("    Opcion: ");
            int sub;
            if (scanf("%d", &sub) == 1) {
                while (getchar() != '\n');
                Politica *R_inicial = NULL;
                if (sub == 1) {
                    Politica *politicas = NULL;
                    int num_pol = 0;
                    politica_leer_consola(modelo, &politicas, &num_pol);
                    if (num_pol > 0) {
                        R_inicial = politica_copiar(&politicas[0]);
                        for (int i = 0; i < num_pol; i++)
                            free(politicas[i].decision);
                        free(politicas);
                    }
                } else {
                    R_inicial = politica_crear(n, K);
                    for (int i = 0; i < n; i++)
                        R_inicial->decision[i] = 0;
                }

                if (R_inicial) {
                    Politica *optima = mejoramiento_politicas(modelo, R_inicial);
                    politica_destruir(optima);
                    politica_destruir(R_inicial);
                }
            } else { while (getchar() != '\n'); }
            break;
        }
        case 3: {
            printf("\n  Elegir politica inicial para Mejoramiento con Descuento:\n");
            printf("    1. Ingresar politicas manualmente\n");
            printf("    2. Usar politica por defecto (decision 1 para todo)\n");
            printf("    Opcion: ");
            int sub;
            if (scanf("%d", &sub) == 1) {
                while (getchar() != '\n');
                Politica *R_inicial = NULL;
                if (sub == 1) {
                    Politica *politicas = NULL;
                    int num_pol = 0;
                    politica_leer_consola(modelo, &politicas, &num_pol);
                    if (num_pol > 0) {
                        R_inicial = politica_copiar(&politicas[0]);
                        for (int i = 0; i < num_pol; i++)
                            free(politicas[i].decision);
                        free(politicas);
                    }
                } else {
                    R_inicial = politica_crear(n, K);
                    for (int i = 0; i < n; i++)
                        R_inicial->decision[i] = 0;
                }

                if (R_inicial) {
                    Politica *optima = mejoramiento_politicas_descuento(
                        modelo, R_inicial);
                    politica_destruir(optima);
                    politica_destruir(R_inicial);
                }
            } else { while (getchar() != '\n'); }
            break;
        }
        case 4: {
            int max_iter;
            double epsilon;
            printf("  Maximo de iteraciones N: ");
            if (scanf("%d", &max_iter) != 1 || max_iter <= 0) max_iter = 100;
            while (getchar() != '\n');
            printf("  Tolerancia epsilon: ");
            if (scanf("%lf", &epsilon) != 1 || epsilon <= 0.0) epsilon = 1e-6;
            while (getchar() != '\n');
            aproximaciones_sucesivas(modelo, max_iter, epsilon);
            break;
        }
        case 5:
            programacion_lineal(modelo);
            break;
        case 6:
            printf("\n  Ejecutando PRUEBA COMPLETA de PMD...\n");
            prueba_completa_pmd(modelo);
            break;
        case 0: break;
        default: printf("  Opcion invalida.\n");
        }
    } while (opcion != 0);
}

int main(int argc, char *argv[]) {
    setlocale(LC_ALL, "");  /* Soporte para caracteres en espanol */

    ModeloMarkov *modelo = NULL;

    /* Si se pasa archivo como argumento, cargarlo directamente */
    if (argc > 1) {
        modelo = modelo_leer_archivo(argv[1]);
    }

    int opcion;
    do {
        limpiar_pantalla();
        if (modelo) {
            printf("\n  [Modelo cargado: %d estados, %s]\n",
                   modelo->num_estados,
                   modelo->tiene_decisiones ? "PMD" : "Cadena simple");
        }
        menu_principal();
        if (scanf("%d", &opcion) != 1) opcion = -1;
        while (getchar() != '\n');

        switch (opcion) {
        case 1:
            if (modelo) modelo_destruir(modelo);
            modelo = modelo_leer_consola();
            break;
        case 2: {
            char archivo[256];
            printf("  Nombre del archivo: ");
            if (fgets(archivo, sizeof(archivo), stdin)) {
                archivo[strcspn(archivo, "\n")] = '\0';
                if (modelo) modelo_destruir(modelo);
                modelo = modelo_leer_archivo(archivo);
            }
            break;
        }
        case 3:
            if (modelo) modelo_imprimir(modelo);
            else printf("  No hay modelo cargado.\n");
            break;
        case 4:
            limpiar_pantalla();
            if (modelo) ejecutar_modulo1(modelo);
            else printf("  No hay modelo cargado. Cargue uno primero.\n");
            break;
        case 5:
            limpiar_pantalla();
            if (modelo) ejecutar_modulo2(modelo);
            else printf("  No hay modelo cargado. Cargue uno primero.\n");
            break;
        case 0:
            printf("\n  Saliendo del programa.\n");
            break;
        default:
            printf("  Opcion invalida.\n");
        }
    } while (opcion != 0);

    if (modelo) modelo_destruir(modelo);
    return 0;
}
