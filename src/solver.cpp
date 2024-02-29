#include "solver.h"
#include <stdlib.h>
#include <stdio.h>

void solver::init(const unsigned n, const float dt, const float diff, const float visc, SolverMethod method)
{
	this->dt_ = dt;
	this->diff_ = diff;
	this->visc_ = visc;
	this->N = n;
	this->method = method;
}

bool solver::allocate_data(void)
{
	u = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	v = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	dens = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));

	u_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	v_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));
	dens_prev = static_cast<float*>(malloc((N + 2) * (N + 2) * sizeof(float)));

	if (CORRECT_DATA) {
		clear_data();
		return true;
	}
	
    return false;
}

void solver::free_data() const
{
	free(u);
	free(v);
	free(dens);
	free(u_prev);
	free(v_prev);
	free(dens_prev);
}

void solver::clear_data(void)
{
	for (int i = 0; i < (N + 2) * (N + 2); i++) {
		u[i] = v[i] = dens[i] = 0.0f;
	}
	clear_prev_data();
}

void solver::clear_prev_data() const
{
	for (int i = 0; i < (N + 2) * (N + 2); i++)
	{
		u_prev[i] = v_prev[i] = dens_prev[i] = 0.0f;
	}
}

void solver::add_density(const unsigned i, const unsigned j, const float source) const
{
	dens[XY_TO_ARRAY(i, j)] += source;
}

void solver::add_velocity(const unsigned i, const unsigned j, const float force_x, const float force_y) const
{
	u[XY_TO_ARRAY(i, j)] += force_x;
	v[XY_TO_ARRAY(i, j)] += force_y;
}

void solver::add_gravity() const
{
	int i, j;
	FOR_EACH_CELL
		add_velocity(i, j, -grav_,0);
	END_FOR
}

void solver::solve()
{
	vel_step();
	dens_step();
	add_gravity();
}

void solver::dens_step()
{
	add_source(dens, dens_prev);
	SWAP(dens_prev, dens)
	diffuse(0, dens, dens_prev);
	SWAP(dens_prev, dens)
	advect(0, dens, dens_prev, u, v);
}

void solver::vel_step()
{
	add_source(u, u_prev);
	add_source(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev,v)
	diffuse(1, u, u_prev);  
	diffuse(2, v, v_prev); 
	project(u, v, u_prev, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev,v)
	advect(1, u, u_prev, u_prev, v_prev);
	advect(2, v, v_prev, u_prev, v_prev);
	project(u, v, u_prev, v_prev);
}

void solver::add_source(float * x, const float * s) const
{
	const int size = (N + 2) * (N + 2);
	for (int i = 0; i < size; i++)
		x[i] += dt_ * s[i];
}

auto solver::set_bounds(const int b, float* x) const -> void
{
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
	const float a = dt_ * diff_ * N * N;
	lin_solve(b, x, x0, a, 1 + 4 * a);
}

void solver::advect(const int b, float * d, const float * d0, const float * u, const float * v) const
{
	const float dt0 = dt_ * N;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            float x = i - dt0 * u[XY_TO_ARRAY(i, j)];
            float y = j - dt0 * v[XY_TO_ARRAY(i, j)];
            if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; int i0 = (int)x; int i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; int j0 = (int)y; int j1 = j0 + 1;
            const float s1 = x - i0; float s0 = 1 - s1; float t1 = y - j0; float t0 = 1 - t1;
            d[XY_TO_ARRAY(i, j)] = s0 * (t0 * d0[XY_TO_ARRAY(i0, j0)] + t1 * d0[XY_TO_ARRAY(i0, j1)]) +
                                   s1 * (t0 * d0[XY_TO_ARRAY(i1, j0)] + t1 * d0[XY_TO_ARRAY(i1, j1)]);
        }
    }
    set_bounds(b, d);
}

void solver::project(float * u, float * v, float * p, float * div)
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

    // Calculate divergence
    FOR_EACH_CELL
        div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
        p[XY_TO_ARRAY(i, j)] = 0;
    END_FOR
    set_bounds(0, div);
    set_bounds(0, p);

    // Perform Jacobi iteration
    for (int k = 0; k < 20; k++) {
        FOR_EACH_CELL
            p[XY_TO_ARRAY(i, j)] = (div[XY_TO_ARRAY(i, j)] + p[XY_TO_ARRAY(i - 1, j)] + p[XY_TO_ARRAY(i + 1, j)] + p[XY_TO_ARRAY(i, j - 1)] + p[XY_TO_ARRAY(i, j + 1)]) / 4;
        END_FOR
        set_bounds(0, p);
    }

    // Subtract gradient of pressure from velocity to make it mass-conserving
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

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
	p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	set_bounds(0, div);
	set_bounds(0, p);

	lin_solve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
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

    // Calculate divergence
    FOR_EACH_CELL
        div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
        p[XY_TO_ARRAY(i, j)] = 0;
    END_FOR
    set_bounds(0, div);
    set_bounds(0, p);

    // Perform SOR iteration
    for (int k = 0; k < 20; k++) {
	    constexpr float w = 1.7f;
	    FOR_EACH_CELL
		    const float p_new = (div[XY_TO_ARRAY(i, j)] + p[XY_TO_ARRAY(i - 1, j)] + p[XY_TO_ARRAY(i + 1, j)] + p[XY_TO_ARRAY(i, j - 1)] + p[XY_TO_ARRAY(i, j + 1)]) / 4;
            p[XY_TO_ARRAY(i, j)] = (1 - w) * p[XY_TO_ARRAY(i, j)] + w * p_new;
        END_FOR
        set_bounds(0, p);
    }

    // Subtract gradient of pressure from velocity to make it mass-conserving
    FOR_EACH_CELL
        u[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
        v[XY_TO_ARRAY(i, j)] -= 0.5f * N * (p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
    END_FOR
    set_bounds(1, u);
    set_bounds(2, v);
}
