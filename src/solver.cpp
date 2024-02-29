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
//TODO: Libera los buffers de memoria.
}

void Solver::ClearData(void)
{
//TODO: Borra todo el contenido de los buffers
}

bool Solver::AllocateData(void)
{
//TODO:
//Reservamos memoria, en caso de fallo devlvemos false.
//Antes de devolver true, hay que limpiar la memoria reservada con un ClearData().
	
	return false;
}

void Solver::ClearPrevData() 
{
//TODO: Borra el contenido de los buffers _prev
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
//TODO: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
//TODO: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
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
//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
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
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
//TODO: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
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