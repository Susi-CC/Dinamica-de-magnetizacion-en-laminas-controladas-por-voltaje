# Dinámica de la Magnetización en Láminas Delgadas Controladas por Voltaje

Este repositorio contiene los códigos utilizados en el desarrollo de mi tesis de magíster: **"Dinámica de la magnetización en láminas delgadas controladas por voltaje"**. Los códigos están implementados principalmente en lenguaje C y Python y se enfocan en la integración numérica de la ecuación de Landau-Lifshitz para estudiar las dinámicas de osciladores magnéticos aislados y acoplados.

## Contenido

- **`oscilador_aislado.c`**: Código en C para simular la dinámica de un oscilador magnético desacoplado bajo un voltaje oscilatorio.
- **`osciladores_acoplados.c`**: Código en C que integra las ecuaciones de múltiples osciladores acoplados mediante campos dipolares.
- **`calculo_lyapunov.py`**: Script en Python para el cálculo del máximo exponente de Lyapunov.
- **`graficar_resultados.py`**: Script en Python para la visualización de resultados.

## Publicaciones relacionadas

Los resultados de este trabajo han sido publicados en los siguientes artículos científicos:

- [Voltage-controlled magnetic anisotropy: Dynamics and chaos](https://doi.org/10.1016/j.cnsns.2024.106051)
- [Magnetic oscillators and their synchronization: Theory and applications](https://doi.org/10.1016/j.jmmm.2022.169023)

## Dependencias

Para compilar y ejecutar los códigos en C:
- Un compilador de C, como `gcc`.

Para ejecutar los scripts de Python:
- Python 3.8 o superior.
- Librerías necesarias:
  - `numpy`
  - `matplotlib`

Puedes instalar las librerías de Python usando:
```bash
pip install numpy matplotlib
```

## Uso

### 1. Simulación de un oscilador aislado

1. Compila el código:
   ```bash
   gcc -o oscilador_aislado oscilador_aislado.c -lm
   ```
2. Ejecuta el programa:
   ```bash
   ./oscilador_aislado
   ```
3. Los resultados se guardarán en archivos CSV en el directorio de ejecución.

### 2. Simulación de múltiples osciladores acoplados

1. Compila el código:
   ```bash
   gcc -o osciladores_acoplados osciladores_acoplados.c -lm
   ```
2. Ejecuta el programa:
   ```bash
   ./osciladores_acoplados
   ```
3. Los resultados se guardarán en archivos CSV en el directorio de ejecución.

### 3. Análisis y visualización de resultados

1. Edita los scripts de Python para apuntar a los archivos CSV generados por las simulaciones.
2. Ejecuta los scripts:
   ```bash
   python calculo_lyapunov.py
   python graficar_resultados.py
   ```

## Autor

**Susana Alejandra Contreras Celada**

## Licencia

Este proyecto está licenciado bajo la licencia MIT. Puedes encontrar más detalles en el archivo `LICENSE`.

---

¡Espero que este repositorio sea útil para futuros trabajos relacionados con dinámicas de magnetización y osciladores acoplados! Si tienes preguntas, no dudes en contactarme.
