#ifndef hogks1d_H
#define hogks1d_H

extern int K;
extern double r;
extern double Mu;
extern double Nu;
extern double c1_euler;
extern double c2_euler;
extern double T_inf; //the real temperature
extern double R_gas; //the gas constant for real ideal air
extern double S;
enum TAU_TYPE { Euler, NS, Sutherland };
extern TAU_TYPE tau_type;
extern bool Smooth;
enum Solver_type{gks,hllc,exactRS};
extern Solver_type solver_choose;
enum GKS1d_type{nothing,kfvs1st,kfvs2nd,gks1st,gks2nd};
extern GKS1d_type gks1dsolver;
enum BC_type{ freebc, wall, periodic };
extern BC_type leftbc, rightbc;
enum G0_construct_type { collisionn, collisionnless };
extern G0_construct_type g0type;
enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
//basic data structure
// the block instore the global geometry information 
typedef class Block1d
{
public:
	bool uniform;
	int ghost;
	int nodex;
	int nx;
	double dx;
	double left;
	double right;
	int stages;
	double timecoefficient[5][5][3];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
};

// to remeber the cell avg values
typedef class Fluid1d
{
public:
	double primvar[3];
	double convar[3];
	double convar_old[3];
	double cx; //center coordinate in x direction
	double dx; //the mesh size dx
};

// to remember the fluid information in a fixed point, 
// such as reconstructioned value, or middle point value
typedef class Point1d
{
public:
	double convar[3];
	double convar_old[3];
	double der1[3];
	double der2[3];
	double x; // coordinate of a point 
};

// remember the flux in a fixed interface, 
// for RK method, we only need f
// for 2nd der method, we need derf
// for 3rd der method, we need der2f
typedef class Flux1d
{
public:
	double F[3]; //total flux in dt time
	double f[3]; // the f0 in t=0
	double derf[4]; // the f_t in t=0
	double der2f[4]; // the f_tt in t=0 only active when 3rd order flux
};

// every interfaces have left center and right value.
// and center value for GKS
// and a group of flux for RK method or multi-stage GKS method
typedef class Interface1d
{
public:
	Point1d left;
	Point1d center;
	Point1d right;
	Flux1d *flux;
	double x; // coordinate of the interface, equal to point1d.x
};



// some basic macro variable function
double U(double density, double densityu);
double Temperature(double density, double pressure);
double entropy(double density, double pressure);
double Soundspeed(double density, double pressure);

void Copy_variable_1D(double *newprim, double oldprim[3]);
void Convar_to_primvar(Fluid1d *fluids, Block1d block);
void Convar_to_primvar_1D(double *primvar, double convar[3]);
double Pressure(double density, double densityu, double densityE);

void Primvar_to_convar_1D(double *convar, double primvar[3]);
double DensityU(double density, double u);
double DensityE(double density, double u, double pressure);
double Lambda(double density, double u, double densityE);
void Convar_to_ULambda_1d(double *primvar, double convar[3]);
void Convar_to_char1D(double *character, double primvar[3], double convar[3]);
void Char_to_convar1D(double *convar, double primvar[3], double charvar[3]);
//get the explicit time step through Dtx function
double Dtx(double dtx, double dx, double CFL, double convar[3]);

double Get_MAX(double a1, double a2);
double Get_MIN(double a1, double a2);
//basic gks function
double Alpha(double lambda, double u);
double Beta(double lambda, double u);
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;

public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in);
	void calcualte_MMDF1d();
};

// to calculate the microsolpe moment
void G(int no_u, int no_xi, double *psi, double a[3], MMDF1d m);
void GL(int no_u, int no_xi, double *psi, double a[3], MMDF1d m);
void GR(int no_u, int no_xi, double *psi, double a[3], MMDF1d m);
void Microslope(double *a, double der[3], double prim[3]);

void sod1d();
void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d *interfaces, Flux1d **fluxes);
void ICfor1dRM(Fluid1d *fluids, Fluid1d zone1, Fluid1d zone2, Block1d block);
void output1d(Fluid1d *fluids, Block1d block);
void output1d_checking(Fluid1d *fluids,Interface1d *interfaces,Flux1d **fluxes, Block1d block);
void CopyFluid_new_to_old(Fluid1d *fluids, Block1d block);
//compute the cfl number
double Get_CFL(Block1d &block, Fluid1d *fluids, double tstop);

typedef void(*BoundaryCondition) (Fluid1d *fluids, Block1d block, Fluid1d bcvalue);

void free_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);


void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void(*Reconstruction_within_Cell)(Point1d &left, Point1d &right, Fluid1d *fluids);
extern Reconstruction_within_Cell cellreconstruction;
void Vanleer(Point1d &left, Point1d &right, Fluid1d *fluids);
void WENO5Z(Point1d &left, Point1d &right, Fluid1d *fluids);
void WENO5Z_left(double &var, double &der1, double der2, double wn2, double wn1, double w, double wp1, double wp2, double h);
void WENO5Z_right(double &var, double &der1, double der2, double wn2, double wn1, double w, double wp1, double wp2, double h);
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_3rd(Interface1d &interfaces, Fluid1d *fluids);
void Center_5th(Interface1d &interfaces, Fluid1d *fluids);
void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d &block, int stage);
typedef void(*Flux_function)(Flux1d &flux, Interface1d& interface, double dt);
extern Flux_function flux_function;
double Get_Tau_NS(double density0, double lambda0);
double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt);
double TauNS_Sutherland(double density0, double lambda0);
void GKS(Flux1d &flux, Interface1d& interface, double dt);

void ustarforHLLC(double d1, double u1, double p1, double s1, double star1, double *ustar);
void get_Euler_flux(double p[3], double *flux);
void HLLC(Flux1d &flux, Interface1d& interface, double dt);

typedef void(*TimeMarchingCoefficient)(Block1d &block);
extern TimeMarchingCoefficient timecoe_list;
void S1O1(Block1d &block);
void S1O2(Block1d &block);
void S1O3(Block1d &block);
void S2O4(Block1d &block);
void RK4(Block1d &block);
void Initial_stages(Block1d &block);
void Update(Fluid1d *fluids, Flux1d **fluxes, Block1d block);
void Update(Fluid1d *fluids, Flux1d **fluxes, Block1d block, int stage);

int main();
#endif