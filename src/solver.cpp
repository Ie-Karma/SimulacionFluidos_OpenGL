#include "solver.h"
#include <stdlib.h>
#include <stdio.h>

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
}
/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
	free(u);
	free(v);
	free(dens);

	free(u_prev);
	free(v_prev);
	free(dens_prev);
}

void Solver::ClearData(void)
{

	for (int i = 0; i < (N + 2) * (N + 2); i++) {
		u[i] = v[i] = dens[i] = 0.0f;
	}
	ClearPrevData();
	
}

bool Solver::AllocateData(void)
{

	u = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
	v = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
	dens = (float *)malloc((N + 2) * (N + 2) * sizeof(float));

	u_prev = (float *)malloc((N + 2) * (N + 2) * sizeof(float));
	v_prev = (float *)malloc((N + 2) * (N + 2) * sizeof(float));

	dens_prev = (float *)malloc((N + 2) * (N + 2) * sizeof(float));

	if (u != NULL || v != NULL || dens != NULL || u_prev != NULL || v_prev != NULL || dens_prev != NULL) {
		ClearData();
		return true;
	}
	
    return false;
}
void Solver::ClearPrevData()
{

	for (int i = 0; i < (N + 2) * (N + 2); i++)
	{
		u_prev[i] = v_prev[i] = dens_prev[i] = 0.0f;
	}
	
	
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	dens[XY_TO_ARRAY(x, y)] += source;
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	u[XY_TO_ARRAY(x, y)] += forceX;
	v[XY_TO_ARRAY(x, y)] += forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

void Solver::DensStep()
{
	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	//SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	//Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	//SWAP (u_prev,u)			
	//SWAP (v_prev,v)
	//Advect(1, u, u_prev, u_prev, v_prev);
	//Advect(2, v, v_prev, u_prev, v_prev);
	//Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
		base[i] += dt * source[i];
}


void Solver::SetBounds(int b, float * x)
{
/*TODO:
Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
*/

	int i;
	for (i = 1; i <= N; i++) {
		x[XY_TO_ARRAY(0, i)] = b == 1 ? -x[XY_TO_ARRAY(1, i)] : x[XY_TO_ARRAY(1, i)];
		x[XY_TO_ARRAY(N + 1, i)] = b == 1 ? -x[XY_TO_ARRAY(N, i)] : x[XY_TO_ARRAY(N, i)];
		x[XY_TO_ARRAY(i, 0)] = b == 2 ? -x[XY_TO_ARRAY(i, 1)] : x[XY_TO_ARRAY(i, 1)];
		x[XY_TO_ARRAY(i, N + 1)] = b == 2 ? -x[XY_TO_ARRAY(i, N)] : x[XY_TO_ARRAY(i, N)];
	}
	x[XY_TO_ARRAY(0, 0)] = 0.5f*(x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]);
	x[XY_TO_ARRAY(0, N + 1)] = 0.5f*(x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]);
	x[XY_TO_ARRAY(N + 1, 0)] = 0.5f*(x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]);
	x[XY_TO_ARRAY(N + 1, N + 1)] = 0.5f*(x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]);
	
}

/*
https://www.youtube.com/watch?v=62_RUX_hrT4
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii)
{
//TODO: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
	int i, j, k;
	for (k = 0; k < 20; k++) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (x0[XY_TO_ARRAY(i, j)] + aij * (x[XY_TO_ARRAY(i - 1, j)] + x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i, j - 1)] + x[XY_TO_ARRAY(i, j + 1)])) / aii;
		END_FOR
		SetBounds(b, x);
	}
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
//TODO: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.

	float a = dt * diff * N * N;
	LinSolve(b, x, x0, a, 1 + 4 * a);
	
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v)
{
//TODO: Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.


	

}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{
	int i, j;

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}