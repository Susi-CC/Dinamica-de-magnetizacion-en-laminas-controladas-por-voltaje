# Dynamics of Magnetization in Thin Films Controlled by Voltage

This repository contains the codes developed for the thesis: **"Dynamics of Magnetization in Thin Films Controlled by Voltage"**. The programs are implemented in **C** and **Python**, focusing on the numerical integration of the **Landau-Lifshitz equation** to study the dynamics of isolated and coupled magnetic oscillators.

## Repository Structure

* **`Lyapunov-oscilador-desacoplado.py`**: Calculates the **Lyapunov exponent** for an isolated magnetic oscillator under an oscillatory voltage. The code is designed to integrate over a β1 interval from 0.12 to 0.23. It outputs a `.csv` file with the λLLE for each β1 in the interval, with data separated by commas.

* **`Lyapunov-osciladores-acoplados.c`**: Calculates the **Lyapunov exponent** for multiple magnetic oscillators coupled via dipolar fields. The code is designed to couple 6 oscillators with a specific initial condition. At the end, the code outputs the maximum Lyapunov exponent to the terminal and generates `.csv` files for **mx**, **my**, and **mz** at a given β1. These files contain the trajectory difference vector information for the six oscillators before normalization, organized in comma-separated tables.

  **Compilation and execution:**
  ```bash
  gcc -o Lyapunov-osciladores-acoplados Lyapunov-osciladores-acoplados.c -lm
  ./Lyapunov-osciladores-acoplados


# Dinámica de la Magnetización en Láminas Delgadas Controladas por Voltaje

Este repositorio contiene los códigos desarrollados para la tesis: **"Dinámica de la magnetización en láminas delgadas controladas por voltaje"**. Los programas están implementados en C y Python, enfocados en la integración numérica de la ecuación de Landau-Lifshitz para estudiar la dinámica de osciladores magnéticos aislados y acoplados.

## Estructura del Repositorio

- **`Lyapunov-oscilador-desacoplado.py`**: Calcula el exponente de Lyapunov para un oscilador magnético aislado bajo un voltaje oscilatorio. El código está diseñado para integrar dando un intervalo de β1 de 0.12 a 0.23. El código devuelve un archivo .csv con el λLLE para cada β1 de intervalo, con los datos separados por comas.

- **`Lyapunov-osciladores-acoplados.c`**: Calcula el exponente de Lyapunov para múltiples osciladores magnéticos acoplados mediante campos dipolares. El código está diseñado para acoplar 6 osciladores con una condición inicial específica. Al finalizar, el código devuelve por el terminal el máximo exponente de Lyapunov, además una vez transcurrido el tiempo transiente genera documentos .csv para mx, my y mz en determinado β1, con la información del vector de diferencia de trayectoria para los seis osciladores antes de normalizar, organizados en tablas separadas por comas.

  **Compilación y ejecución:**
  ```bash
  gcc -o Lyapunov-osciladores-acoplados Lyapunov-osciladores-acoplados.c -lm
  ./Lyapunov-osciladores-acoplados
  ```

- **`Oscilador-desacoplado.py`**: Simula la dinámica de un oscilador magnético aislado bajo un voltaje oscilatorio. El código está diseñado para integrar dando un intervalo de β1 de 0.0 a 0.123 y con una condición aleatoria originada por una semilla. El código devuelve tres archivos .csv con los datos de la magnetización en el eje x, y y z para cada β1 del intervalo, los datos están separados por comas.

  **Ejecución:**
  ```bash
  python Oscilador-desacoplado.py
  ```

- **`Osciladores-acoplados.c`**: Integra las ecuaciones de múltiples osciladores magnéticos acoplados mediante campos dipolares. El código está diseñado para acoplar 6 osciladores. Para acoplar una cantidad distinta de osciladores, se cambia el valor de N en la línea 7 del código. Al finalizar, el código devuelve documentos .csv para mx, my y mz en determinado β1, con la información de la magnetización para los seis osciladores posterior al tiempo transiente, esta información está organizada en tablas separadas por comas.

  **Compilación y ejecución:**
  ```bash
  gcc -o Osciladores-acoplados Osciladores-acoplados.c -lm
  ./Osciladores-acoplados
  ```

## Publicaciones Relacionadas

Los resultados de este trabajo han sido publicados en los siguientes artículos científicos:

- "Magnetic chimeras in voltage-driven nano-oscillators", *Communications in Nonlinear Science and Numerical Simulation*, 2025. [https://doi.org/10.1016/j.cnsns.2024.108420](https://doi.org/10.1016/j.cnsns.2024.108420). (Autor correspondiente: S. Contreras-Celada)

- "Voltage-driven multistability and chaos in magnetic films", *Journal of Magnetism and Magnetic Materials*, 2022. [https://doi.org/10.1016/j.jmmm.2022.169793](https://doi.org/10.1016/j.jmmm.2022.169793). (Autor correspondiente: S. Contreras-Celada)

## Dependencias

Para compilar y ejecutar los códigos en C:
- Compilador de C (e.g., `gcc`).

Para los scripts en Python:
- Python 3.8 o superior.
- Librerías:
  - `numpy`
  - `matplotlib`

Instalación de librerías de Python:
```bash
pip install numpy matplotlib
```

## Resultados Esperados

- **Lyapunov-oscilador-desacoplado.py:** Genera un archivo .csv con el λLLE para cada β1 en el intervalo definido.
- **Lyapunov-osciladores-acoplados.c:** Imprime el máximo exponente de Lyapunov en el terminal y genera archivos .csv con datos organizados por comas para las magnetizaciones mx, my, mz en determinado β1.
- **Oscilador-desacoplado.py:** Produce tres archivos .csv con los datos de la magnetización en los ejes x, y y z para cada β1 del intervalo.
- **Osciladores-acoplados.c:** Crea documentos .csv para mx, my y mz, mostrando las magnetizaciones para múltiples osciladores después del tiempo transiente.

## Autor

**Susana Alejandra Contreras Celada**

## Licencia

Este proyecto está bajo la licencia MIT. Más detalles en el archivo `LICENSE`.

---

¡Espero que este repositorio sea útil para investigaciones futuras sobre dinámicas de magnetización y osciladores acoplados! Para consultas, no dudes en contactarme.

