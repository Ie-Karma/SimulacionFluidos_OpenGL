#ifndef _SOLVER_H_
#define SOLVER_H

#define XY_TO_ARRAY(i,j) ((i)+(N+2)*(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}
#define SWAP(x0,x) {float * tmp=x0;(x0)=x;(x)=tmp;}
#define CORRECT_DATA (u != nullptr || v != nullptr || dens != nullptr || u_prev != nullptr || v_prev != nullptr || dens_prev != nullptr)

class solver
{
	float dt_ = 0, diff_ = 0, visc_ = 0, grav_ = 0;
	unsigned N;
	float * u_prev = nullptr, * v_prev = nullptr, * dens_prev = nullptr;

public:

	typedef enum { JACOBI, GAUSS_SEIDEL, SOR } SolverMethod;

	SolverMethod method;
	float * u, * v, * dens;
	void init(unsigned n, float dt, float diff, float visc, SolverMethod method);
	void free_data(void) const;
	void clear_data(void) const;
	bool allocate_data(void);
	void clear_prev_data(void) const;
	void add_density(unsigned i, unsigned j, float source) const;
	void add_velocity(unsigned i, unsigned j, float force_x, float force_y) const;
	void solve(void);
private:
	void dens_step(void);
	void vel_step(void);
	void add_gravity(void) const;

	void add_source(float * x, const float * s) const;
	void set_bounds(int b, float * x) const;
	void lin_solve( int b, float * x, const float * x0, float a, float c) const;
	void diffuse(int b, float * x, const float * x0) const;
	void advect(int b, float * d, const float * d0, const float * u, const float * v) const;
	void jacobi_project(float* p, float* v, float* p1, float* div) const;
	void gauss_seidel_project(float* p, float* v, float* p1, float* div) const;
	void sor_project(float* p, float* v, float* p1, float* div) const;
	void project(float * u, float * v, float * p, float * div) const;
};

#endif
