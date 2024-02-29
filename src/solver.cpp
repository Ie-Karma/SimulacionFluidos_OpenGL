#include "solver.h"
#include <stdlib.h>
#include <stdio.h>

void solver::init(const unsigned n, const float dt, const float diff, const float visc, SolverMethod method)
{
    // Inicializa las variables de la simulación con los valores proporcionados.
	this->dt_ = dt;  // Intervalo de tiempo para cada paso de la simulación.
	this->diff_ = diff;  // Coeficiente de difusión.
	this->visc_ = visc;  // Viscosidad del fluido.
	this->N = n;  // Tamaño de la cuadrícula.
	this->method = method;  // Método de resolución a utilizar.
}

bool solver::allocate_data(void)
{
    // Reserva memoria para las matrices de velocidad y densidad.
	u = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	v = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	dens = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));

    // Reserva memoria para las matrices de velocidad y densidad del paso anterior.
	u_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	v_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	dens_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));

    // Si la asignación de memoria fue exitosa, inicializa las matrices a cero.
	if (CORRECT_DATA) {
		clear_data();
		return true;
	}
	
    return false;
}

void solver::free_data() const
{
    // Libera la memoria reservada para las matrices de velocidad y densidad.
	free(u);
	free(v);
	free(dens);
	free(u_prev);
	free(v_prev);
	free(dens_prev);
}

void solver::clear_data(void) const
{
    // Inicializa las matrices de velocidad y densidad a cero.
	for (int i = 0; i < (N + 2) * (N + 2); i++) {
		u[i] = v[i] = dens[i] = 0.0f;
	}
	clear_prev_data();
}

void solver::clear_prev_data() const
{
    // Inicializa las matrices de velocidad y densidad del paso anterior a cero.
	for (int i = 0; i < (N + 2) * (N + 2); i++)
	{
		u_prev[i] = v_prev[i] = dens_prev[i] = 0.0f;
	}
}

void solver::add_density(const unsigned i, const unsigned j, const float source) const
{
    // Añade densidad en la celda (i, j).
	dens[XY_TO_ARRAY(i, j)] += source;
}

void solver::add_velocity(const unsigned i, const unsigned j, const float force_x, const float force_y) const
{
    // Añade velocidad en la celda (i, j).
	u[XY_TO_ARRAY(i, j)] += force_x;
	v[XY_TO_ARRAY(i, j)] += force_y;
}

void solver::add_gravity() const
{
    // Añade la fuerza de gravedad a cada celda.
	int i, j;
	FOR_EACH_CELL
		add_velocity(i, j, -grav_,0);
	END_FOR
}

void solver::solve()
{
    // Realiza un paso de la simulación.
	vel_step();  // Actualiza las velocidades.
	dens_step();  // Actualiza las densidades.
	add_gravity();  // Añade la gravedad.
}

void solver::dens_step()
{
    // Realiza un paso de la simulación para la densidad.
	add_source(dens, dens_prev);  // Añade la densidad de la fuente.
	SWAP(dens_prev, dens)  // Intercambia las matrices de densidad.
	diffuse(0, dens, dens_prev);  // Difunde la densidad.
	SWAP(dens_prev, dens)  // Intercambia las matrices de densidad.
	advect(0, dens, dens_prev, u, v);  // Advierte la densidad.
}

void solver::vel_step()
{
    // Realiza un paso de la simulación para la velocidad.
	add_source(u, u_prev);  // Añade la velocidad de la fuente en la dirección x.
	add_source(v, v_prev);  // Añade la velocidad de la fuente en la dirección y.
	SWAP (u_prev,u)  // Intercambia las matrices de velocidad en la dirección x.
	SWAP (v_prev,v)  // Intercambia las matrices de velocidad en la dirección y.
	diffuse(1, u, u_prev);  // Difunde la velocidad en la dirección x.
	diffuse(2, v, v_prev);  // Difunde la velocidad en la dirección y.
	project(u, v, u_prev, v_prev);  // Proyecta las velocidades para hacerlas masivamente conservadoras.
	SWAP (u_prev,u)  // Intercambia las matrices de velocidad en la dirección x.
	SWAP (v_prev,v)  // Intercambia las matrices de velocidad en la dirección y.
	advect(1, u, u_prev, u_prev, v_prev);  // Advierte la velocidad en la dirección x.
	advect(2, v, v_prev, u_prev, v_prev);  // Advierte la velocidad en la dirección y.
	project(u, v, u_prev, v_prev);  // Proyecta las velocidades para hacerlas masivamente conservadoras.
}

void solver::add_source(float * x, const float * s) const
{
    // Añade la fuente a la matriz x.
	const int size = (N + 2) * (N + 2);
	for (int i = 0; i < size; i++)
		x[i] += dt_ * s[i];
}

auto solver::set_bounds(const int b, float* x) const -> void
{
    // Establece las condiciones de contorno para la matriz x.
	for (int i = 1; i <= N; i++) {
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

void solver::lin_solve(const int b, float * x, const float * x0, const float a, const float c) const
{
    // Resuelve el sistema lineal utilizando el método de Jacobi.
	int i, j;
	for (int k = 0; k < 20; k++) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (x0[XY_TO_ARRAY(i, j)] + a * (x[XY_TO_ARRAY(i - 1, j)] + x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i, j - 1)] + x[XY_TO_ARRAY(i, j + 1)])) / c;
		END_FOR
		set_bounds(b, x);
	}
}

void solver::diffuse(const int b, float * x, const float * x0) const
{
    // Difunde la matriz x.
	const float a = dt_ * diff_ * N * N;
	lin_solve(b, x, x0, a, 1 + 4 * a);
}

void solver::advect(const int b, float * d, const float * d0, const float * u, const float * v) const
{
    // Calculamos el paso de tiempo multiplicado por el tamaño de la cuadrícula.
	const float dt0 = dt_ * N;

    // Recorremos cada celda de la cuadrícula.
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {

            // Calculamos las coordenadas de la celda de origen de la advección.
            float x = i - dt0 * u[XY_TO_ARRAY(i, j)];
            float y = j - dt0 * v[XY_TO_ARRAY(i, j)];

            // Nos aseguramos de que las coordenadas de la celda de origen estén dentro de los límites de la cuadrícula.
            if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; int i0 = (int)x; int i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; int j0 = (int)y; int j1 = j0 + 1;

            // Calculamos los pesos para la interpolación lineal.
            const float s1 = x - i0; float s0 = 1 - s1; float t1 = y - j0; float t0 = 1 - t1;

            // Realizamos la interpolación lineal para obtener el valor de la celda en el paso de tiempo siguiente.
            d[XY_TO_ARRAY(i, j)] = s0 * (t0 * d0[XY_TO_ARRAY(i0, j0)] + t1 * d0[XY_TO_ARRAY(i0, j1)]) +
                                   s1 * (t0 * d0[XY_TO_ARRAY(i1, j0)] + t1 * d0[XY_TO_ARRAY(i1, j1)]);
        }
    }

    // Aplicamos las condiciones de contorno a la matriz d.
    set_bounds(b, d);
}
void solver::project(float * u, float * v, float * p, float * div) const
{
	switch (method)
	{
	case JACOBI:
		jacobi_project(u, v, p, div);
		break;
	case GAUSS_SEIDEL:
		gauss_seidel_project(u, v, p, div);
		break;
	case SOR:
		sor_project(u, v, p, div);
		break;
	}
}

void solver::jacobi_project(float* u, float* v, float* p, float* div) const
{
    int i, j;

    // Para cada celda, calculamos la divergencia del campo de velocidad
    FOR_EACH_CELL
        div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
        p[XY_TO_ARRAY(i, j)] = 0;
    END_FOR
    set_bounds(0, div);
    set_bounds(0, p);

    // Realizamos 20 iteraciones del método de Jacobi para resolver el sistema lineal
    for (int k = 0; k < 20; k++) {
        FOR_EACH_CELL
            p[XY_TO_ARRAY(i, j)] = (div[XY_TO_ARRAY(i, j)] + p[XY_TO_ARRAY(i - 1, j)] + p[XY_TO_ARRAY(i + 1, j)] + p[XY_TO_ARRAY(i, j - 1)] + p[XY_TO_ARRAY(i, j + 1)]) / 4;
        END_FOR
        set_bounds(0, p);
    }

    // Ajustamos el campo de velocidad para que sea libre de divergencia (masa conservada)
    FOR_EACH_CELL
        u[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
        v[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
    END_FOR
    set_bounds(1, u);
    set_bounds(2, v);
}

void solver::gauss_seidel_project(float * u, float * v, float * p, float * div) const
{

	int i, j;

	// Para cada celda, calculamos la divergencia del campo de velocidad
	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
	p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	set_bounds(0, div);
	set_bounds(0, p);

	// Resolvemos el sistema lineal utilizando el método de Gauss-Seidel
	lin_solve(0, p, div, 1, 4);

	// Ajustamos el campo de velocidad para que sea libre de divergencia (masa conservada)
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
	v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	set_bounds(1, u);
	set_bounds(2, v);
	
}

void solver::sor_project(float* u, float* v, float* p, float* div) const
{
    int i, j;

    // Para cada celda, calculamos la divergencia del campo de velocidad
    FOR_EACH_CELL
        div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
        p[XY_TO_ARRAY(i, j)] = 0;
    END_FOR
    set_bounds(0, div);
    set_bounds(0, p);

    // Realizamos 20 iteraciones del método SOR (Successive Over-Relaxation) para resolver el sistema lineal
    for (int k = 0; k < 20; k++) {
	    constexpr float w = 1.7f;
	    FOR_EACH_CELL
		    const float p_new = (div[XY_TO_ARRAY(i, j)] + p[XY_TO_ARRAY(i - 1, j)] + p[XY_TO_ARRAY(i + 1, j)] + p[XY_TO_ARRAY(i, j - 1)] + p[XY_TO_ARRAY(i, j + 1)]) / 4;
            p[XY_TO_ARRAY(i, j)] = (1 - w) * p[XY_TO_ARRAY(i, j)] + w * p_new;
        END_FOR
        set_bounds(0, p);
    }

    // Ajustamos el campo de velocidad para que sea libre de divergencia (masa conservada)
    FOR_EACH_CELL
        u[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
        v[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
    END_FOR
    set_bounds(1, u);
    set_bounds(2, v);
}