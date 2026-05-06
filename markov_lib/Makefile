# Makefile para el paquete de Cadenas de Markov y PMD
# Compilador y banderas
CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -I.
LDFLAGS = -lm

# Archivos fuente
SRC = matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c
OBJ = $(SRC:.c=.o)

# Nombre del ejecutable
TARGET = markov

# Regla principal
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.c markov.h matrix.h
	$(CC) $(CFLAGS) -c $< -o $@

# Limpiar archivos objeto y ejecutable
clean:
	rm -f $(OBJ) $(TARGET)

# Limpiar solo objetos
cleanobj:
	rm -f $(OBJ)

# Compilar y ejecutar con ejemplo de prueba
test: $(TARGET)
	./$(TARGET)

.PHONY: all clean cleanobj test
