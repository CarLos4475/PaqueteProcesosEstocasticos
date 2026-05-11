#include "markov.h"
#include <string.h>
#include <ctype.h>

#define TAM_BUFFER 1024

/* ===================================================================
 * ALGORITMO DE LECTURA DE DATOS
 * Nucleo del paquete. Lee y almacena toda la configuracion del sistema.
 * Soporta cadenas de Markov simples, PMD con costos directos, y PMD
 * con matrices de ingreso + costos fijos (auto-generacion de C_ik).
 * =================================================================== */

/* Crear un modelo vacio */
ModeloMarkov* modelo_crear(void) {
    ModeloMarkov *m = (ModeloMarkov*)calloc(1, sizeof(ModeloMarkov));
    return m;
}

/* Destruir un modelo y liberar toda la memoria asociada */
void modelo_destruir(ModeloMarkov *modelo) {
    if (!modelo) return;
    free(modelo->prob_inicial);
    matriz_destruir(modelo->P);
    if (modelo->P_dec) {
        for (int k = 0; k < modelo->num_decisiones; k++)
            matriz_destruir(modelo->P_dec[k]);
        free(modelo->P_dec);
    }
    matriz_destruir(modelo->C);
    if (modelo->Ingreso) {
        for (int k = 0; k < modelo->num_decisiones; k++)
            matriz_destruir(modelo->Ingreso[k]);
        free(modelo->Ingreso);
    }
    free(modelo->costo_fijo);
    if (modelo->nombres_estados) {
        for (int i = 0; i < modelo->num_estados; i++)
            free(((char**)modelo->nombres_estados)[i]);
        free(modelo->nombres_estados);
    }
    if (modelo->nombres_decisiones) {
        for (int k = 0; k < modelo->num_decisiones; k++)
            free(((char**)modelo->nombres_decisiones)[k]);
        free(modelo->nombres_decisiones);
    }
    free(modelo);
}

/* ===================================================================
 * Funciones auxiliares de entrada
 * =================================================================== */
static void limpiar_buffer(void) {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

static int leer_entero(const char *mensaje, int min, int max) {
    int valor;
    while (1) {
        printf("%s", mensaje);
        if (scanf("%d", &valor) == 1 && valor >= min && valor <= max) {
            limpiar_buffer();
            return valor;
        }
        printf("  Valor invalido. Ingrese un entero entre %d y %d.\n", min, max);
        limpiar_buffer();
    }
}

static double leer_double(const char *mensaje, double min, double max) {
    double valor;
    while (1) {
        printf("%s", mensaje);
        if (scanf("%lf", &valor) == 1 && valor >= min && valor <= max) {
            limpiar_buffer();
            return valor;
        }
        printf("  Valor invalido. Ingrese un numero entre %.4f y %.4f.\n", min, max);
        limpiar_buffer();
    }
}

static int leer_si_no(const char *mensaje) {
    char buffer[TAM_BUFFER];
    while (1) {
        printf("%s", mensaje);
        if (fgets(buffer, TAM_BUFFER, stdin)) {
            if (buffer[0] == 'S' || buffer[0] == 's') return 1;
            if (buffer[0] == 'N' || buffer[0] == 'n') return 0;
        }
        printf("  Ingrese S (Si) o N (No).\n");
    }
}

static double* leer_vector(int n, const char *nombre) {
    double *v = (double*)calloc((size_t)n, sizeof(double));
    if (!v) return NULL;
    printf("\n--- Ingrese el vector %s (%d componentes) ---\n", nombre, n);
    for (int i = 0; i < n; i++) {
        char msg[128];
        snprintf(msg, sizeof(msg), "  Componente %d: ", i);
        v[i] = leer_double(msg, -1e9, 1e9);
    }
    return v;
}

static Matriz* leer_matriz(int filas, int columnas, const char *nombre, int estocastica) {
    Matriz *M = matriz_crear(filas, columnas);
    if (!M) return NULL;
    printf("\n--- Ingrese la matriz %s (%d x %d) ---\n", nombre, filas, columnas);
    for (int i = 0; i < filas; i++) {
        printf("  Fila %d:\n", i);
        double suma_fila = 0.0;
        for (int j = 0; j < columnas; j++) {
            char msg[128];
            snprintf(msg, sizeof(msg), "    [%d][%d]: ", i, j);
            double val = leer_double(msg, estocastica ? 0.0 : -1e9,
                                          estocastica ? 1.0 : 1e9);
            M->datos[i][j] = val;
            if (estocastica) suma_fila += val;
        }
        if (estocastica && fabs(suma_fila - 1.0) > 0.001) {
            printf("  Advertencia: la fila %d suma %.6f (se esperaba 1.0). "
                   "Se normalizara.\n", i, suma_fila);
            for (int j = 0; j < columnas; j++)
                M->datos[i][j] /= suma_fila;
        }
    }
    return M;
}

/* ===================================================================
 * modelo_calcular_costos_ingresos
 *
 * Auto-genera la matriz C_ik a partir de las matrices de ingreso
 * y los costos fijos:
 *    C_ik = sum_j P_ij(k) * Ingreso_ij(k) - costo_fijo[k]
 *
 * Si es_maximizacion == 1, los valores se niegan porque el paquete
 * minimiza: min (-utilidad) = max utilidad.
 * =================================================================== */
void modelo_calcular_costos_ingresos(ModeloMarkov *modelo) {
    int n = modelo->num_estados;
    int K = modelo->num_decisiones;

    /* Liberar C anterior si existe */
    matriz_destruir(modelo->C);
    modelo->C = matriz_crear(n, K);
    if (!modelo->C) return;

    imprimir_titulo_seccion("Calculo de la matriz de costos netos C_ik");
    printf("  Formula:  C_ik = ( sum_j P_ij(k) * R_ij(k) )  -  Costo_fijo(k)\n");
    if (modelo->es_maximizacion)
        printf("  Nota   :  valor negado (paquete minimiza)  =>  C = -costo_neto\n");

    printf("\n  +-------+-----+-------------+-------------+-------------+-------------+\n");
    printf("  |  i,k  |  k  |    E[R_ik]  |   C_fijo(k) |  Neto (i,k) |    C_ik     |\n");
    printf("  +-------+-----+-------------+-------------+-------------+-------------+\n");

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < K; k++) {
            double ingreso_esperado = 0.0;
            for (int j = 0; j < n; j++) {
                ingreso_esperado += modelo->P_dec[k]->datos[i][j]
                                  * modelo->Ingreso[k]->datos[i][j];
            }
            double costo_neto = ingreso_esperado - modelo->costo_fijo[k];
            double c_final = modelo->es_maximizacion ? -costo_neto : costo_neto;
            modelo->C->datos[i][k] = c_final;

            printf("  | (%d,%d) | %3d | %11.6f | %11.6f | %11.6f | %11.6f |\n",
                   i, k+1, k+1, ingreso_esperado, modelo->costo_fijo[k],
                   costo_neto, c_final);
        }
    }
    printf("  +-------+-----+-------------+-------------+-------------+-------------+\n");
    modelo->usa_ingresos = 1;
}

/* ===================================================================
 * modelo_leer_consola: lee todo el modelo desde la entrada estandar.
 * Pregunta interactivamente por cada componente.
 * =================================================================== */
ModeloMarkov* modelo_leer_consola(void) {
    ModeloMarkov *modelo = modelo_crear();
    if (!modelo) return NULL;

    printf("\n");
    printf("=============================================================\n");
    printf("  ALGORITMO DE LECTURA DE DATOS\n");
    printf("  Cadenas de Markov y Procesos Markovianos de Decision\n");
    printf("=============================================================\n\n");

    /* Preguntar si es cadena simple o PMD */
    int es_pmd = leer_si_no("Es un Proceso Markoviano de Decision (PMD)? (S/N): ");

    /* Espacio de estados: E = {0, 1, ..., m} */
    int m = leer_entero("Ingrese el valor de m (ultimo estado, E = {0..m}): ", 1, 500);
    modelo->num_estados = m + 1;
    modelo->tiene_decisiones = es_pmd;

    if (!es_pmd) {
        /* ----- CADENA DE MARKOV SIMPLE ----- */
        modelo->num_decisiones = 0;

        modelo->prob_inicial = leer_vector(modelo->num_estados, "a (probabilidad inicial)");
        double suma = 0.0;
        for (int i = 0; i < modelo->num_estados; i++)
            suma += modelo->prob_inicial[i];
        if (fabs(suma) > 1e-12)
            for (int i = 0; i < modelo->num_estados; i++)
                modelo->prob_inicial[i] /= suma;

        modelo->P = leer_matriz(modelo->num_estados, modelo->num_estados,
                                "P (matriz de transicion)", 1);
    } else {
        /* ----- PROCESO MARKOVIANO DE DECISION ----- */
        modelo->num_decisiones = leer_entero(
            "Ingrese el numero de decisiones K: ", 1, 100);

        /* Vector de probabilidad inicial */
        modelo->prob_inicial = leer_vector(modelo->num_estados,
                                           "a (probabilidad inicial)");
        double suma = 0.0;
        for (int i = 0; i < modelo->num_estados; i++)
            suma += modelo->prob_inicial[i];
        if (fabs(suma) > 1e-12)
            for (int i = 0; i < modelo->num_estados; i++)
                modelo->prob_inicial[i] /= suma;

        /* Matrices de transicion por decision */
        modelo->P_dec = (Matriz**)malloc(
            (size_t)modelo->num_decisiones * sizeof(Matriz*));
        for (int k = 0; k < modelo->num_decisiones; k++) {
            char nombre[64];
            snprintf(nombre, sizeof(nombre),
                     "P^(%d) (matriz de transicion decision %d)", k+1, k+1);
            modelo->P_dec[k] = leer_matriz(modelo->num_estados, modelo->num_estados,
                                           nombre, 1);
        }

        /* Preguntar si usar ingresos o costos directos */
        int usar_ingresos = leer_si_no(
            "Desea ingresar matrices de ingreso (en lugar de costos directos)? (S/N): ");

        if (usar_ingresos) {
            /* Preguntar maximizacion o minimizacion */
            modelo->es_maximizacion = leer_si_no(
                "Es un problema de MAXIMIZACION de utilidad? (S/N, N=minimizacion): ");

            /* Matrices de ingreso */
            modelo->Ingreso = (Matriz**)malloc(
                (size_t)modelo->num_decisiones * sizeof(Matriz*));
            for (int k = 0; k < modelo->num_decisiones; k++) {
                char nombre[64];
                snprintf(nombre, sizeof(nombre),
                         "Ingreso^(%d) (matriz de ingresos, decision %d)", k+1, k+1);
                modelo->Ingreso[k] = leer_matriz(modelo->num_estados,
                    modelo->num_estados, nombre, 0);
            }

            /* Costos fijos por decision */
            modelo->costo_fijo = (double*)malloc(
                (size_t)modelo->num_decisiones * sizeof(double));
            if (modelo->costo_fijo)
                memset(modelo->costo_fijo, 0,
                       (size_t)modelo->num_decisiones * sizeof(double));
            printf("\n--- Ingrese los costos fijos por decision ---\n");
            for (int k = 0; k < modelo->num_decisiones; k++) {
                char msg[64];
                snprintf(msg, sizeof(msg), "  Costo fijo decision %d: ", k+1);
                modelo->costo_fijo[k] = leer_double(msg, -1e9, 1e9);
            }

            /* Auto-generar C_ik desde ingresos */
            modelo_calcular_costos_ingresos(modelo);

            /* Factor de descuento */
            int usar_tasa = leer_si_no(
                "Desea ingresar tasa de interes i en lugar de alfa? (S/N): ");
            if (usar_tasa) {
                modelo->tasa_interes = leer_double(
                    "  Ingrese i (ej. 0.10 para 10%%): ", 0.0, 10.0);
                modelo->alfa = 1.0 / (1.0 + modelo->tasa_interes);
            } else {
                modelo->alfa = leer_double(
                    "  Ingrese alfa (0 < alfa <= 1): ", 0.0001, 1.0);
                modelo->tasa_interes = (1.0 / modelo->alfa) - 1.0;
            }
        } else {
            /* Costos directos */
            modelo->es_maximizacion = 0;
            modelo->C = leer_matriz(modelo->num_estados, modelo->num_decisiones,
                                    "C (matriz de costos directos)", 0);

            int usar_tasa = leer_si_no(
                "Desea ingresar tasa de interes i en lugar de alfa? (S/N): ");
            if (usar_tasa) {
                modelo->tasa_interes = leer_double(
                    "  Ingrese i (ej. 0.10 para 10%%): ", 0.0, 10.0);
                modelo->alfa = 1.0 / (1.0 + modelo->tasa_interes);
            } else {
                modelo->alfa = leer_double(
                    "  Ingrese alfa (0 < alfa <= 1): ", 0.0001, 1.0);
                modelo->tasa_interes = (1.0 / modelo->alfa) - 1.0;
            }
        }
    }

    printf("\n--- Datos leidos correctamente ---\n\n");
    return modelo;
}

/* ===================================================================
 * modelo_leer_archivo: lee el modelo desde un archivo de texto.
 *
 * Formato del archivo (nuevo, extendido):
 *   Linea 1: tipo (0=cadena simple, 1=PMD)
 *   Linea 2: m (ultimo indice de estado)
 *   Si cadena simple:
 *     Linea 3: a_0 a_1 ... a_m
 *     Lineas 4..(4+m): matriz P (m+1 filas de m+1 numeros)
 *   Si PMD:
 *     Linea 3: K (numero de decisiones)
 *     Linea 4: modo_costos (0=directos, 1=ingresos+maxim, 2=ingresos+minim)
 *     Linea 5: a_0 a_1 ... a_m
 *     Luego K matrices P (una por decision)
 *     Si modo_costos=0: luego matriz C (m+1 filas de K numeros)
 *     Si modo_costos=1/2: luego K matrices Ingreso, luego K costos fijos
 *     Luego alfa (factor de descuento)
 * =================================================================== */
ModeloMarkov* modelo_leer_archivo(const char *nombre_archivo) {
    FILE *f = fopen(nombre_archivo, "r");
    if (!f) {
        fprintf(stderr, "Error: No se pudo abrir el archivo '%s'.\n", nombre_archivo);
        return NULL;
    }

    ModeloMarkov *modelo = modelo_crear();
    if (!modelo) { fclose(f); return NULL; }

    int tipo, m;
    if (fscanf(f, "%d", &tipo) != 1) { fclose(f); modelo_destruir(modelo); return NULL; }
    if (fscanf(f, "%d", &m) != 1)    { fclose(f); modelo_destruir(modelo); return NULL; }

    modelo->num_estados      = m + 1;
    modelo->tiene_decisiones = tipo;
    modelo->prob_inicial = (double*)calloc((size_t)modelo->num_estados, sizeof(double));

    if (tipo == 0) {
        modelo->num_decisiones = 0;
        for (int i = 0; i < modelo->num_estados; i++)
            fscanf(f, "%lf", &modelo->prob_inicial[i]);
        modelo->P = matriz_crear(modelo->num_estados, modelo->num_estados);
        for (int i = 0; i < modelo->num_estados; i++)
            for (int j = 0; j < modelo->num_estados; j++)
                fscanf(f, "%lf", &modelo->P->datos[i][j]);
    } else {
        int K, modo_costos;
        if (fscanf(f, "%d", &K) != 1 || fscanf(f, "%d", &modo_costos) != 1) {
            fclose(f); modelo_destruir(modelo); return NULL;
        }
        modelo->num_decisiones = K;

        for (int i = 0; i < modelo->num_estados; i++)
            fscanf(f, "%lf", &modelo->prob_inicial[i]);

        modelo->P_dec = (Matriz**)malloc((size_t)K * sizeof(Matriz*));
        for (int k = 0; k < K; k++) {
            modelo->P_dec[k] = matriz_crear(modelo->num_estados, modelo->num_estados);
            for (int i = 0; i < modelo->num_estados; i++)
                for (int j = 0; j < modelo->num_estados; j++)
                    fscanf(f, "%lf", &modelo->P_dec[k]->datos[i][j]);
        }

        if (modo_costos == 0) {
            modelo->es_maximizacion = 0;
            modelo->usa_ingresos = 0;
            modelo->C = matriz_crear(modelo->num_estados, modelo->num_decisiones);
            for (int i = 0; i < modelo->num_estados; i++)
                for (int k = 0; k < modelo->num_decisiones; k++)
                    fscanf(f, "%lf", &modelo->C->datos[i][k]);
        } else {
            modelo->es_maximizacion = (modo_costos == 1);
            modelo->usa_ingresos = 1;

            modelo->Ingreso = (Matriz**)malloc((size_t)K * sizeof(Matriz*));
            for (int k = 0; k < K; k++) {
                modelo->Ingreso[k] = matriz_crear(modelo->num_estados, modelo->num_estados);
                for (int i = 0; i < modelo->num_estados; i++)
                    for (int j = 0; j < modelo->num_estados; j++)
                        fscanf(f, "%lf", &modelo->Ingreso[k]->datos[i][j]);
            }

            modelo->costo_fijo = (double*)malloc((size_t)K * sizeof(double));
            if (modelo->costo_fijo)
                memset(modelo->costo_fijo, 0, (size_t)K * sizeof(double));
            for (int k = 0; k < K; k++)
                fscanf(f, "%lf", &modelo->costo_fijo[k]);

            modelo_calcular_costos_ingresos(modelo);
        }

        if (fscanf(f, "%lf", &modelo->alfa) != 1)
            modelo->alfa = 0.95;
        modelo->tasa_interes = (fabs(modelo->alfa) > 1e-12)
            ? (1.0 / modelo->alfa) - 1.0 : 0.0;
    }

    fclose(f);
    printf("\n--- Datos leidos correctamente desde '%s' ---\n\n", nombre_archivo);
    return modelo;
}

/* ===================================================================
 * modelo_imprimir: muestra en pantalla un resumen de los datos cargados.
 * =================================================================== */
void modelo_imprimir(const ModeloMarkov *modelo) {
    imprimir_titulo_box("RESUMEN DEL MODELO");
    printf("  Tipo                : %s\n",
           modelo->tiene_decisiones ? "PMD" : "Cadena de Markov simple");
    printf("  Numero de estados   : %d  (E = {0, 1, ..., %d})\n",
           modelo->num_estados, modelo->num_estados - 1);
    printf("  Numero de decisiones: %d\n", modelo->num_decisiones);

    if (modelo->tiene_decisiones) {
        printf("  Objetivo            : %s\n",
               modelo->es_maximizacion ? "MAXIMIZAR utilidad" : "MINIMIZAR costo");
        if (modelo->usa_ingresos)
            printf("  Matriz C            : auto-generada desde ingresos\n");
    }

    vector_imprimir(modelo->prob_inicial, modelo->num_estados,
                    "Vector de probabilidad inicial a");

    if (!modelo->tiene_decisiones) {
        matriz_imprimir(modelo->P, "Matriz de transicion P");
    } else {
        for (int k = 0; k < modelo->num_decisiones; k++) {
            char titulo[64];
            snprintf(titulo, sizeof(titulo), "Matriz de transicion P^(%d)", k+1);
            matriz_imprimir(modelo->P_dec[k], titulo);
        }

        if (modelo->usa_ingresos && modelo->Ingreso) {
            for (int k = 0; k < modelo->num_decisiones; k++) {
                char titulo[64];
                snprintf(titulo, sizeof(titulo),
                         "Matriz de Ingreso^(%d)", k+1);
                matriz_imprimir(modelo->Ingreso[k], titulo);
            }
            imprimir_titulo_sub("Costos fijos por decision");
            printf("  +----------+----------+\n");
            printf("  | decision |   costo  |\n");
            printf("  +----------+----------+\n");
            for (int k = 0; k < modelo->num_decisiones; k++)
                printf("  |    %2d    | %8.2f |\n", k+1, modelo->costo_fijo[k]);
            printf("  +----------+----------+\n");
        }

        matriz_imprimir(modelo->C, "Matriz de costos C_ik (usada en algoritmos)");
        printf("\n  Factor de descuento alfa : %.6f\n", modelo->alfa);
        printf("  Tasa de interes i        : %.6f\n", modelo->tasa_interes);
    }
    imprimir_separador('=', 65);
    printf("\n");
}
