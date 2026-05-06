#include "markov.h"
#include <string.h>
#include <ctype.h>

#define TAM_BUFFER 1024

/* ===================================================================
 * ALGORITMO DE LECTURA DE DATOS
 * Nucleo del paquete. Lee y almacena toda la configuracion del sistema.
 * Soporta cadenas de Markov simples y Procesos Markovianos de Decision.
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
    free(modelo);
}

/* ===================================================================
 * limpiar_buffer: descarta caracteres sobrantes en stdin.
 * =================================================================== */
static void limpiar_buffer(void) {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

/* ===================================================================
 * leer_entero: lee un entero desde stdin con validacion.
 * =================================================================== */
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

/* ===================================================================
 * leer_double: lee un double desde stdin con validacion.
 * =================================================================== */
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

/* ===================================================================
 * leer_si_no: lee respuesta S/N.
 * =================================================================== */
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

/* ===================================================================
 * leer_vector: lee un vector de n componentes desde consola.
 * =================================================================== */
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

/* ===================================================================
 * leer_matriz: lee una matriz de filas x columnas desde consola.
 *              si estocastica == 1, valida que las filas sumen 1.0 y
 *              cada elemento este en [0, 1].
 * =================================================================== */
static Matriz* leer_matriz(int filas, int columnas, const char *nombre, int estocastica) {
    Matriz *M = matriz_crear(filas, columnas);
    if (!M) return NULL;
    printf("\n--- Ingrese la matriz %s (%d x %d) ---\n", nombre, filas, columnas);
    for (int i = 0; i < filas; i++) {
        printf("  Fila %d:\n", i);
        double suma_fila = 0.0;
        for (int j = 0; j < columnas; j++) {
            char msg[128];
            if (estocastica)
                snprintf(msg, sizeof(msg), "    P[%d][%d]: ", i, j);
            else
                snprintf(msg, sizeof(msg), "    [%d][%d]: ", i, j);
            double val = leer_double(msg, estocastica ? 0.0 : -1e9,
                                          estocastica ? 1.0 : 1e9);
            M->datos[i][j] = val;
            if (estocastica) suma_fila += val;
        }
        if (estocastica && fabs(suma_fila - 1.0) > 0.001) {
            printf("  Advertencia: la fila %d suma %.6f (se esperaba 1.0)\n",
                   i, suma_fila);
            printf("  Los valores seran normalizados automaticamente.\n");
            for (int j = 0; j < columnas; j++)
                M->datos[i][j] /= suma_fila;
        }
    }
    return M;
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

        /* Vector de probabilidad inicial a */
        modelo->prob_inicial = leer_vector(modelo->num_estados, "a (probabilidad inicial)");

        /* Normalizar si no suma 1 */
        double suma = 0.0;
        for (int i = 0; i < modelo->num_estados; i++)
            suma += modelo->prob_inicial[i];
        if (fabs(suma) > 1e-12)
            for (int i = 0; i < modelo->num_estados; i++)
                modelo->prob_inicial[i] /= suma;

        /* Matriz de transicion P */
        modelo->P = leer_matriz(modelo->num_estados, modelo->num_estados,
                                "P (matriz de transicion)", 1);

    } else {
        /* ----- PROCESO MARKOVIANO DE DECISION ----- */

        /* Conjunto de decisiones K = {1, 2, ..., K} */
        modelo->num_decisiones = leer_entero(
            "Ingrese el numero de decisiones K: ", 1, 100);

        /* Vector de probabilidad inicial */
        modelo->prob_inicial = leer_vector(modelo->num_estados,
                                           "a (probabilidad inicial)");
        /* Normalizar */
        double suma = 0.0;
        for (int i = 0; i < modelo->num_estados; i++)
            suma += modelo->prob_inicial[i];
        if (fabs(suma) > 1e-12)
            for (int i = 0; i < modelo->num_estados; i++)
                modelo->prob_inicial[i] /= suma;

        /* Matrices de transicion por decision: P^k para k=1..K */
        modelo->P_dec = (Matriz**)malloc(
            (size_t)modelo->num_decisiones * sizeof(Matriz*));
        for (int k = 0; k < modelo->num_decisiones; k++) {
            char nombre[64];
            snprintf(nombre, sizeof(nombre), "P^(%d) (matriz de transicion decision %d)", k+1, k+1);
            modelo->P_dec[k] = leer_matriz(modelo->num_estados, modelo->num_estados,
                                           nombre, 1);
        }

        /* Matriz de costos esperados C_ik */
        modelo->C = leer_matriz(modelo->num_estados, modelo->num_decisiones,
                                "C (matriz de costos)", 0);

        /* Factor de descuento alfa o tasa de interes i */
        int usar_tasa = leer_si_no(
            "Desea ingresar la tasa de interes i en lugar del factor de descuento? (S/N): ");
        if (usar_tasa) {
            modelo->tasa_interes = leer_double(
                "  Ingrese la tasa de interes i (ej. 0.10 para 10%%): ", 0.0, 10.0);
            modelo->alfa = 1.0 / (1.0 + modelo->tasa_interes);
        } else {
            modelo->alfa = leer_double(
                "  Ingrese el factor de descuento alfa (0 < alfa <= 1): ", 0.0001, 1.0);
            modelo->tasa_interes = (1.0 / modelo->alfa) - 1.0;
        }
    }

    printf("\n--- Datos leidos correctamente ---\n\n");
    return modelo;
}

/* ===================================================================
 * modelo_leer_archivo: lee el modelo desde un archivo de texto.
 *
 * Formato del archivo:
 *   Linea 1: tipo (0=cadena simple, 1=PMD)
 *   Linea 2: m (ultimo indice de estado)
 *   Si cadena simple:
 *     Linea 3: a_0 a_1 ... a_m
 *     Lineas 4..(4+m): matriz P (m+1 filas de m+1 numeros)
 *   Si PMD:
 *     Linea 3: K (numero de decisiones)
 *     Linea 4: a_0 a_1 ... a_m
 *     Luego K matrices P (una por decision)
 *     Luego matriz C (m+1 filas de K numeros)
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
        /* Cadena simple */
        modelo->num_decisiones = 0;

        for (int i = 0; i < modelo->num_estados; i++)
            fscanf(f, "%lf", &modelo->prob_inicial[i]);

        modelo->P = matriz_crear(modelo->num_estados, modelo->num_estados);
        for (int i = 0; i < modelo->num_estados; i++)
            for (int j = 0; j < modelo->num_estados; j++)
                fscanf(f, "%lf", &modelo->P->datos[i][j]);

    } else {
        /* PMD */
        int K;
        fscanf(f, "%d", &K);
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

        modelo->C = matriz_crear(modelo->num_estados, modelo->num_decisiones);
        for (int i = 0; i < modelo->num_estados; i++)
            for (int k = 0; k < modelo->num_decisiones; k++)
                fscanf(f, "%lf", &modelo->C->datos[i][k]);

        fscanf(f, "%lf", &modelo->alfa);
        modelo->tasa_interes = (1.0 / modelo->alfa) - 1.0;
    }

    fclose(f);
    printf("\n--- Datos leidos correctamente desde '%s' ---\n\n", nombre_archivo);
    return modelo;
}

/* ===================================================================
 * modelo_imprimir: muestra en pantalla un resumen de los datos cargados.
 * =================================================================== */
void modelo_imprimir(const ModeloMarkov *modelo) {
    printf("\n=============================================================\n");
    printf("  RESUMEN DEL MODELO\n");
    printf("=============================================================\n");
    printf("  Tipo: %s\n", modelo->tiene_decisiones ? "PMD" : "Cadena de Markov simple");
    printf("  Numero de estados: %d (E = {0, 1, ..., %d})\n",
           modelo->num_estados, modelo->num_estados - 1);
    printf("  Numero de decisiones: %d\n", modelo->num_decisiones);

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
        matriz_imprimir(modelo->C, "Matriz de costos esperados C_ik");
        printf("\n  Factor de descuento alfa: %.6f\n", modelo->alfa);
        printf("  Tasa de interes i: %.6f\n", modelo->tasa_interes);
    }
    printf("=============================================================\n\n");
}
