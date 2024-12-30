# Dinámica de la Magnetización en Láminas Delgadas Controladas por Voltaje

Este repositorio contiene los códigos desarrollados para la tesis: **"Dinámica de la magnetización en láminas delgadas controladas por voltaje"**. Los programas están implementados en C y Python, enfocados en la integración numérica de la ecuación de Landau-Lifshitz para estudiar la dinámica de osciladores magnéticos aislados y acoplados.

## Estructura del Repositorio

- **`Lyapunov-oscilador-desacoplado.py`**: Calcula el exponente de Lyapunov para un oscilador magnético aislado bajo un voltaje oscilatorio.

- **`Lyapunov-osciladores-acoplados.c`**: Calcula el exponente de Lyapunov para múltiples osciladores magnéticos acoplados mediante campos dipolares.

- **`Oscilador-desacoplado.py`**: Simula la dinámica de un oscilador magnético aislado bajo un voltaje oscilatorio.

- **`Osciladores-acoplados.c`**: Integra las ecuaciones de múltiples osciladores magnéticos acoplados mediante campos dipolares.

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

