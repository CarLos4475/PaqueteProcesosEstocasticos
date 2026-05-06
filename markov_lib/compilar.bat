@echo off
REM Script de compilacion para Windows (usando GCC/MinGW)
REM Alternativa: usar "gcc -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm"

echo Compilando paquete de Cadenas de Markov y PMD...
echo.

gcc -Wall -O2 -std=c99 -I. -o markov.exe matrix.c markov_datos.c markov_cadenas.c markov_decision.c main.c -lm

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilacion exitosa. Ejecutable: markov.exe
    echo Uso: markov.exe [archivo_entrada.txt]
) else (
    echo.
    echo Error en la compilacion. Verifique que GCC/MinGW este instalado.
)

pause
