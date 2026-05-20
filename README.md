# 2D Steady-State Heat Conduction Solver (Fixed Boundaries)

## Aim
To solve the 2D steady-state heat conduction equation (Laplace's equation) on a square plate with fixed boundary temperatures (Dirichlet conditions) using the Finite Difference Method (FDM).

## Theory
The governing equation for 2D steady-state heat conduction without heat generation is:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

Assuming a square grid ($\Delta x = \Delta y$), we apply central difference discretization to the spatial derivatives at interior nodes:

$$\frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{\Delta x^2} + \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{\Delta y^2} = 0 \implies -4T_{i,j} + T_{i+1,j} + T_{i-1,j} + T_{i,j+1} + T_{i,j-1} = 0$$

At the corners, the boundaries are averaged. The script solves the large linear system ($A T_{flat} = B$) and plots a 2D temperature contour.

## File Structure
- `conduction.py` - Core solver initializing the boundary conditions, setting up the $A$ matrix, solving the equations, and displaying the annotated temperature contour.
- `output values.txt` - Generated equation tables and tabular internal nodal temperatures.
- `Graph conduction.png` - Visual plot containing the 2D heat contour grid mapped with temperature readings.

## How to Run
Ensure you have the required dependencies:
```bash
pip install numpy matplotlib
python conduction.py
```
