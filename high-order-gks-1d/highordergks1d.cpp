//this is a light version of high order GKS
//under fintie volume framework
//with high order sptial discrtization, WENO5Z
//and high order time marching method, like RK or Mulistage Method
//The GKS flux solver contains from kfvs1st to gks2nd so far
//and was packaged in one function, 
//by the author's partial understanding of GKS
//also a omp parallel is inside the program,
//It could be Compiled under WINDOWS or Liunx SYSTEM, 
//the author thinks it would also work on MAC, 
//if you want to start learning coding for high order GKS
//it is not a bad idea to following the comments step by step, 

//first it comes some Standrad lib for compile GKS
#include"highordergks1d.h"//of course the header itself
#include"omp.h" //omp parallel lib
#include<cmath> //bais math operation
#include <cstdlib>//for using exit()
//the followings are for Output ISSUE
#include<iostream>
#include <fstream>
#include <sstream>
#include<iomanip>
//Output end
using namespace std;

//some global parameters are difined as follows, 
//maybe you can just jump to void main() for beginning
#define pi 3.14159265358979323846 //since C++ No Pi

int K = 4;  //interal degree of freedom, default diatomic gas
double r = (K+3)/(K+1); //specific heat ratio
//viscous coefficient
double Mu = -1.0;
double Nu = -1.0;
//artifical viscosity for euler problem
double c1_euler = -1.0;
double c2_euler = -1.0;
double T_inf = 285; //the real temperature, default 20 Celcius degree 
double R_gas = 8.31441 / (28.959e-3); //the gas constant for real ideal air
double S = 110.4; //reference entropy value in Sutherland law
//tau_type, euler==0,ns==1,sutherland==2
TAU_TYPE tau_type;
//to set the flow is smooth or not
bool Smooth = false;

//C type Function Pointer, which is convenient for you
//to switch different 
//solvers
Solver_type solver_choose=gks;
GKS1d_type gks1dsolver=nothing;
Flux_function flux_function = GKS;
//and reconstructions
Reconstruction_within_Cell cellreconstruction=Vanleer;
Reconstruction_forG0 g0reconstruction=Center_3rd;
G0_construct_type g0type=collisionn;
Reconstruction_variable reconstruction_variable = conservative;
//and timemarching scheme
TimeMarchingCoefficient timecoe_list=S1O1;
//and boundary condtions
BC_type leftbc = freebc, rightbc = freebc;

//global parameters end


// the main program begin here
int main()
{
	//a example for sodtest case is shown here,
	//one can create another case, like blastwave(), to run other test case
	sod1d();
	return 0;
}

void sod1d() // a sod test case example
{
	//for a specific problem, we need specify some
	//start parameters
	//the folloing block1d contains most of important simulation parameters in 1D computation
	Block1d block;
	//let's first do mesh
	block.uniform = true; //is it uniform mesh
	block.nodex = 100; //how many mesh points (cell numbers)
	block.ghost = 3; //how ghost cell needed? low order suggest 2, weno5 suggest 3
	block.CFL = 0.5; //CFL number

	leftbc = freebc;  //left boundary type
	rightbc = freebc; //right boundary type
	//in this exampple free means just do simple extrapolation
	// for fixed (or Dirichlet) boundary, you might need 
	// specify the boundary values
	Fluid1d *bcvalue = new Fluid1d[2];  //left index 0, right index 1
	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	if (leftbc == freebc)
	{
		leftboundary = free_boundary_left;
	}
	if (rightbc == freebc)
	{
		rightboundary = free_boundary_right;
	}

	// which kind of gas you want to simulate
	K = 4;
	r = 1.4;
	// gas kind specify end

	//prepare the reconstruction function
	cellreconstruction = WENO5Z;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_5th;


	//prepare the flux function
	flux_function = GKS;
	//if you choose GKS, then, the following parameters
	//will take effect, ohterwise, if you use RM solver
	//this parameter won't influence your program
	gks1dsolver = gks2nd;
	// for a sod example, a numerical tau for euler problem is choose
	tau_type = Euler;
	Smooth = false; //g0 is constructed by a local collision of gl and gr
	c1_euler = 0.05;
	c2_euler = 1;

	//prepare time marching stratedgy
	timecoe_list = S2O4; 
	Initial_stages(block);
	
	//prepare the omp parallel environment
	int nthreads = 2;
	omp_set_num_threads(nthreads);

	//the follwings are Memory-Allocating and Geometry Relation Bulid 
	//it is auto-executed after specifying the above parameters

	// allocate memory for 1-D fluid field
    // in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	Fluid1d *fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d *interfaces = new Interface1d[block.nx+1];
	// then the flux, which the number is identical to interfaces
	Flux1d **fluxes = new Flux1d *[block.nx + 1];
	for (int i = 0; i <= block.nx;i++)
	{
	// for m th step time marching schemes, m subflux needed 	
		fluxes[i] = new Flux1d[block.stages]; 
		for (int j = 0; j < block.stages; j++)
		{
			for (int k = 0; k < 3;k++)
			{
				fluxes[i][j].f[k] = 0.0;
				fluxes[i][j].derf[k] = 0.0;
				fluxes[i][j].der2f[k] = 0.0;
			}
		}
	}
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching
	if (block.uniform == true)
	{
		block.left = 0.0;
		block.right = 1.0;
		block.dx = (block.right - block.left) / block.nodex;
        //set the uniform geometry information
		SetUniformMesh(block, fluids, interfaces, fluxes);
	}
	else
	{
		//set non-uniform geomertry information, should be added later
	}
	//ended mesh part
	//Memory-Allocating and Geometry Relation Bulid complete


	//Then it is General Process for solving a given IC/BC problem,
	//with a governing Euler or NS equation
	//under Fintie Volume Framework

	// first is about initializing, lets first initialize a sod test case
	//you can initialize whatever kind of 1d test case as you like
	Fluid1d zone1; zone1.primvar[0] = 1.0; zone1.primvar[1] = 0.0; zone1.primvar[2] = 1.0;
	Fluid1d zone2; zone2.primvar[0] = 0.125; zone2.primvar[1] = 0.0; zone2.primvar[2] = 0.1;
	ICfor1dRM(fluids, zone1, zone2, block);	
	double tstop = 0.2; //the expected simulation time
	//initializing part end

	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
	                  //initialize inputstep=1, to avoid a 0 result
	while (block.t < tstop)
	{
		// assume you are using command window,
		// you can specify a running step conveniently
		if (block.step%inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit "<<endl;
			cin >> inputstep;
			if (inputstep== 0)
			{
				output1d(fluids, block);
				break;
			}
		}
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt=Get_CFL(block, fluids, tstop);

		for (int istage = 0; istage < block.stages; istage++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);
			Reconstruction_forg0(interfaces, fluids, block);

			//then is solver part
			Calculate_flux(fluxes, interfaces, block,istage);

			//then is update flux part
			Update(fluids, fluxes, block,istage);
			//output1d_checking(fluids, interfaces, fluxes, block);
		}

		block.step++; //update step
		block.t = block.t + block.dt; //update simulation time
		if ((block.t - tstop) > 0)
		{
			output1d(fluids, block); //when it is the required ouput time, Output
		}
	}

	//you can free the memory, but i'm lazy

}

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d *interfaces, Flux1d **fluxes)
{
	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		fluids[i].dx = block.dx; //cell size
		fluids[i].cx = block.left + (i + 0.5-block.ghost)*block.dx; //cell center location
	}
	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		interfaces[i].x = block.left + (i-block.ghost)*block.dx;
		interfaces[i].left.x = interfaces[i].x;
		interfaces[i].right.x = interfaces[i].x;
		interfaces[i].center.x = interfaces[i].x;
		interfaces[i].flux = fluxes[i];
	}
}

void ICfor1dRM(Fluid1d *fluids, Fluid1d zone1, Fluid1d zone2, Block1d block)
{
	for (int i = 0; i < block.nx; i++)
	{
		if (i < 0.5*block.nx)
		{
			Copy_variable_1D(fluids[i].primvar, zone1.primvar);
		}
		else
		{
			Copy_variable_1D(fluids[i].primvar, zone2.primvar);
		}
	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
	}
}

void output1d(Fluid1d *fluids, Block1d block)
{
	//for create a output file, pls include the follow standard lib
	//#include <fstream>
    //#include <sstream>
    //#include<iomanip>
	stringstream name;
	name << "Result1D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "variables = x,density,u,pressure,temperature,entropy,Ma" << endl;
	out << "zone i = " << block.nodex << ", F=POINT" << endl;
	
	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data
	for (int i = block.ghost; i < block.ghost+block.nodex; i++)
		{
			out << fluids[i].cx << " "
				<< setprecision(10)
				<< fluids[i].primvar[0] << " "
				<< fluids[i].primvar[1] << " "
				<< fluids[i].primvar[2] << " "
				<< Temperature(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
				<< entropy(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
				<< sqrt(pow(fluids[i].primvar[1], 2)) / Soundspeed(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
				<< endl;
		}

	//close the file
	out.close();

}
void output1d_checking(Fluid1d *fluids, Interface1d *interfaces, Flux1d **fluxes, Block1d block)
{
	
	
	stringstream name;
	name << "Checking-cell-1D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "variables = x,density,u,pressure,temperature,entropy,Ma" << endl;
	out << "zone i = " << block.nodex << ", F=POINT" << endl;

	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		out << fluids[i].cx << " "
			<< setprecision(10)
			<< fluids[i].primvar[0] << " "
			<< fluids[i].primvar[1] << " "
			<< fluids[i].primvar[2] << " "
			<< Temperature(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< entropy(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< sqrt(pow(fluids[i].primvar[1], 2)) / Soundspeed(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< endl;
	}

	//close the file
	out.close();

	stringstream name1;
	name << "Checking-reconstruction-1D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s1;
	name >> s1;
	ofstream out1(s1);
	out1 << "variables = x,rholeft,rhouleft,rhoeleft"
		           <<   ",rhoright,rhouright,rhoeright"
	               <<   ",rhocenter,rhoucenter,rhoecenter" << endl;
	out1 << "zone i = " << block.nodex+1 << ", F=POINT" << endl;
	//output the data
	for (int i = block.ghost; i < block.ghost + block.nodex+1; i++)
	{
		out1 << interfaces[i].x << " "
			<< setprecision(10)
			<< interfaces[i].left.convar[0] << " "
			<< interfaces[i].left.convar[1] << " "
			<< interfaces[i].left.convar[2] << " "
			<< interfaces[i].right.convar[0] << " "
			<< interfaces[i].right.convar[1] << " "
			<< interfaces[i].right.convar[2] << " "
			<< interfaces[i].center.convar[0] << " "
			<< interfaces[i].center.convar[1] << " "
			<< interfaces[i].center.convar[2] << " "
			<< endl;
	}

	//close the file
	out1.close();



}

double U(double density, double densityu)
{
	return densityu / density;
}

double Temperature(double density, double pressure)
{
	//return pressure / density / R*nitrogen;
	return pressure / density;
}
double entropy(double density, double pressure)
{
	return pressure / pow(density, r);
}
double Soundspeed(double density, double pressure)
{
	return sqrt(r*pressure / density);
}

void Copy_variable_1D(double *newvar, double oldvar[3])
{
	for (int i = 0; i < 3; i++)
	{
		newvar[i] = oldvar[i];
	}
}

void Convar_to_primvar(Fluid1d *fluids, Block1d block)
{
#pragma omp parallel  for 
	for (int i = 0; i < block.nodex; i++)
	{
			Convar_to_primvar_1D(fluids[i].primvar, fluids[i].convar);
	}
}

void Convar_to_primvar_1D(double *primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Pressure(convar[0], convar[1], convar[2]);

}
double Pressure(double density, double densityu, double densityE)
{
	return (r - 1)*(densityE - 0.5*densityu*densityu / density);
}

void Primvar_to_convar_1D(double *convar, double primvar[3])
{
	convar[0] = primvar[0];
	convar[1] = DensityU(primvar[0], primvar[1]);
	convar[2] = DensityE(primvar[0], primvar[1], primvar[2]);
}
double DensityU(double density, double u)
{
	return density*u;
}
double DensityE(double density, double u, double pressure)
{
	return density*(pressure / density / (r - 1) + 0.5*(u*u));	
}

double Lambda(double density, double u, double densityE)
{
	return (K + 1.0) / 4.0*(density / (densityE - 0.5*density*(u*u)));
}
void Convar_to_ULambda_1d(double *primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Lambda(convar[0], primvar[1], convar[2]);

}

void Convar_to_char1D(double *character, double primvar[3], double convar[3])
{
	double c = sqrt(r*primvar[2] / primvar[0]);
	double	alfa = (r - 1.0) / (2.0*c*c);
	double u = primvar[1];
	double s[3][3];
	s[0][0] = alfa*(0.5*u*u +u*c/(r-1.0));
	s[0][1] = alfa*(-u-c/(r-1));
	s[0][2] = alfa;
	s[1][0] = alfa*(-u*u + 2.0*c*c / (r - 1.0));
	s[1][1] = alfa*2.0*u;
	s[1][2] = -2.0*alfa;
	s[2][0] = alfa*(0.5*u*u - u*c / (r - 1.0));
	s[2][1] = alfa*(-u+c/(r-1));
	s[2][2] = alfa;

	for (int i = 0; i < 3; i++)
	{
		character[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			character[i] = character[i] + s[i][j] * convar[j];
		}
	}


}

void Char_to_convar1D(double *convar, double primvar[3], double charvar[3])
{
	double c = sqrt(r*primvar[2] / primvar[0]);
	
	double u = primvar[1];
	double	h = 0.5*u*u + c*c / (r - 1.0);
	double s[3][3];
	s[0][0] = 1.0;
	s[0][1] = 1.0;
	s[0][2] = 1.0;
	s[1][0] = u-c;
	s[1][1] = u;
	s[1][2] = u + c;
	s[2][0] = h-u*c;
	s[2][1] = u*u/2.0;
	s[2][2] = h + u*c;

	for (int i = 0; i < 3; i++)
	{
		convar[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			convar[i] = convar[i] + s[i][j] * charvar[j];
		}
	}

}

void CopyFluid_new_to_old(Fluid1d *fluids,Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost+block.nodex; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fluids[i].convar_old[j] = fluids[i].convar[j];
		}
	}
}

double Get_MAX(double a1, double a2)
{
	return (a1>a2) ? a1 : a2;
}
double Get_MIN(double a1, double a2)
{
	return (a1<a2) ? a1 : a2;
}

double Get_CFL(Block1d &block, Fluid1d *fluids, double tstop)
{
double dt = block.dx;
for (int i = 0; i < block.nodex; i++)
{
	dt = Dtx(dt, block.dx, block.CFL, fluids[i + block.ghost].convar);
}
if (block.t + dt>tstop)
{
	dt = tstop - block.t + 1e-15;
}
//print time step information
cout << "step= " << block.step
<< "time size is " << dt
<< " time= " << block.t << endl;
return dt;
}

double Dtx(double dtx, double dx, double CFL, double convar[3])
{
	double tmp;
	double prim[3];
	Convar_to_primvar_1D(prim, convar);
	tmp = abs(prim[1]) + sqrt(r*prim[2] / prim[0]);
	if (tmp>CFL*dx / dtx)
	{
		dtx = CFL*dx / tmp;
	}
	if (dtx > 0.25*CFL*dx*dx / Mu&&tau_type == NS&&Mu>0)
	{
		dtx = 0.25*CFL*dx*dx / Mu;
	}
	return dtx;
}


void free_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.ghost-1; i >= 0; i--)
	{
			fluids[i].convar[0] = fluids[i + 1].convar[0];
			fluids[i].convar[1] = fluids[i + 1].convar[1];
			fluids[i].convar[2] = fluids[i + 1].convar[2];
	}
}


void free_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
			fluids[i].convar[0] = fluids[i - 1].convar[0];
			fluids[i].convar[1] = fluids[i - 1].convar[1];
			fluids[i].convar[2] = fluids[i - 1].convar[2];

	}
}


void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost-1; i < block.nx-block.ghost+1; i++)
	{
		cellreconstruction(interfaces[i].right, interfaces[i + 1].left, &fluids[i]);
	}
}

void FirsrOrder(Point1d &left, Point1d &right, Fluid1d *fluids)
{
	for (int i = 0; i<3; i++)
	{
		left.convar[i] = fluids[0].convar[i];
		right.convar[i] = fluids[0].convar[i];
			left.der1[i] = 0.0;
			right.der1[i] = 0.0;
			left.der2[i] = 0.0;
			right.der2[i] = 0.0;

	}
}

void Vanleer(Point1d &left, Point1d &right, Fluid1d *fluids)
{

		Fluid1d wn1 = fluids[-1];
		Fluid1d w = fluids[0];
		Fluid1d wp1 = fluids[1];
		double splus[3], sminus[3];

		for (int i = 0; i<3; i++)
		{
			splus[i] = (wp1.convar[i] - w.convar[i]) / ((wp1.dx+w.dx)/2.0);
			sminus[i] = (w.convar[i] - wn1.convar[i]) / ((wn1.dx + w.dx) / 2.0);

			if ((splus[i] * sminus[i]) > 0)
			{
				left.der1[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
				right.der1[i] = left.der1[i];
			}
			else
			{
				left.der1[i] = 0.0;
				right.der1[i] = 0.0;
			}
			left.convar[i] = w.convar[i] - 0.5*w.dx*left.der1[i];
			right.convar[i] = w.convar[i] + 0.5*w.dx*right.der1[i];
		}

		//you shall consider the reduce order part

}

void WENO5Z(Point1d &left, Point1d &right, Fluid1d *fluids)
{
	
	
	Fluid1d wn2 = fluids[-2];
	Fluid1d wn1 = fluids[-1];
	Fluid1d w=fluids[0];
	Fluid1d wp1 = fluids[1];
	Fluid1d wp2 = fluids[2];
	double h = w.dx;
    
	double ren2[3], ren1[3], re0[3], rep1[3], rep2[3];
	double var[3], der1[3], der2[3];

	double base_left[3];
	double base_right[3];
	Convar_to_primvar_1D(wn1.primvar, wn1.convar);
	Convar_to_primvar_1D(w.primvar, w.convar);
	Convar_to_primvar_1D(wp1.primvar, wp1.convar);

	for (int i = 0; i < 3; i++)
	{
		base_left[i] = 0.5*(wn1.primvar[i] + w.primvar[i]);
		base_right[i] = 0.5*(wp1.primvar[i] + w.primvar[i]);
	}

	// cell left
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			ren2[i] = wn2.convar[i];
			ren1[i] = wn1.convar[i];
			re0[i] = w.convar[i];
			rep1[i] = wp1.convar[i];
			rep2[i] = wp2.convar[i];
		}
	}
	else
	{
		Convar_to_char1D(ren2, base_left, wn2.convar);
		Convar_to_char1D(ren1, base_left, wn1.convar);
		Convar_to_char1D(re0, base_left, w.convar);
		Convar_to_char1D(rep1, base_left, wp1.convar);
		Convar_to_char1D(rep2, base_left, wp2.convar);
	}


	for (int i = 0; i < 3; i++)
	{
		WENO5Z_left(var[i], der1[i], der2[i], ren2[i],ren1[i],re0[i],rep1[i],rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			left.convar[i] = var[i];
			left.der1[i] = der1[i];
			left.der2[i] = der2[i];
		}
	}
	else
	{	
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der2, base_left, der2);
		Char_to_convar1D(left.der1, base_left, der1);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			ren2[i] = wn2.convar[i];
			ren1[i] = wn1.convar[i];
			re0[i] = w.convar[i];
			rep1[i] = wp1.convar[i];
			rep2[i] = wp2.convar[i];
		}
	}
	else
	{
		Convar_to_char1D(ren2, base_right, wn2.convar);
		Convar_to_char1D(ren1, base_right, wn1.convar);
		Convar_to_char1D(re0, base_right, w.convar);
		Convar_to_char1D(rep1, base_right, wp1.convar);
		Convar_to_char1D(rep2, base_right, wp2.convar);
	}

	for (int i = 0; i < 3; i++)
	{
		WENO5Z_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			right.convar[i] = var[i];
			right.der1[i] = der1[i];
			right.der2[i] = der2[i];
		}
	}
	else
	{
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der2, base_right, der2);
		Char_to_convar1D(right.der1, base_right, der1);
	}

}
void WENO5Z_left(double &var, double &der1, double der2, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qleft[3];
	double qright[3];
	double dleft[3];
	double dright[3];

	qleft[0] = 11.0 / 6.0 * w - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	qleft[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w - 1.0 / 6.0 * wp1;
	qleft[2] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w;


	dleft[0] = 0.1;
	dleft[1] = 0.6;
	dleft[2] = 0.3;

	double 	beta[3];

	beta[0] = 13.0 / 12.0*pow((w - 2 * wp1 + wp2), 2) + 0.25*pow((3 * w - 4 * wp1 + wp2), 2);
	beta[1] = 13.0 / 12.0*pow((wn1 - 2 * w + wp1), 2) + 0.25*pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0*pow((wn2 - 2 * wn1 + w), 2) + 0.25*pow((wn2 - 4 * wn1 + 3 * w), 2);

	double epsilon = 1e-6;

	double tau5 = abs(beta[0] - beta[2]);

	double alphaleft[3];
	for (int i = 0; i < 3; i++)
	{
		alphaleft[i] = dleft[i] * (1 + tau5 / (beta[i] + epsilon));
	}

	double alphal = alphaleft[0] + alphaleft[1] + alphaleft[2];

	double omegaleft[3];
	omegaleft[0] = alphaleft[0] / alphal;
	omegaleft[1] = alphaleft[1] / alphal;
	omegaleft[2] = alphaleft[2] / alphal;

	double left = omegaleft[0] * qleft[0] + omegaleft[1] * qleft[1] + omegaleft[2] * qleft[2];

	qright[0] = 1.0 / 3.0 * w + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	qright[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w + 1.0 / 3.0 * wp1;
	qright[2] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w;

	dright[0] = 0.3;
	dright[1] = 0.6;
	dright[2] = 0.1;

	double alpharight[3];

	for (int i = 0; i < 3; i++)
	{
		alpharight[i] = dright[i] * (1 + tau5 / (beta[i] + epsilon));
	}

	double alphar = alpharight[0] + alpharight[1] + alpharight[2];

	double omegaright[3];

	omegaright[0] = alpharight[0] / alphar;
	omegaright[1] = alpharight[1] / alphar;
	omegaright[2] = alpharight[2] / alphar;

	double right = omegaright[0] * qright[0] + omegaright[1] * qright[1] + omegaright[2] * qright[2];
	var = left;
	der2 = 6 * (left + right- 2 * w) / (h*h);
	der1 = (right - left) / h - 0.5*h*der2;

}
void WENO5Z_right(double &var, double &der1, double der2, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qleft[3];
	double qright[3];
	double dleft[3];
	double dright[3];

	qleft[0] = 11.0 / 6.0 * w - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	qleft[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w - 1.0 / 6.0 * wp1;
	qleft[2] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w;


	dleft[0] = 0.1;
	dleft[1] = 0.6;
	dleft[2] = 0.3;

	double 	beta[3];

	beta[0] = 13.0 / 12.0*pow((w - 2 * wp1 + wp2), 2) + 0.25*pow((3 * w - 4 * wp1 + wp2), 2);
	beta[1] = 13.0 / 12.0*pow((wn1 - 2 * w + wp1), 2) + 0.25*pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0*pow((wn2 - 2 * wn1 + w), 2) + 0.25*pow((wn2 - 4 * wn1 + 3 * w), 2);

	double epsilon = 1e-6;

	double tau5 = abs(beta[0] - beta[2]);

	double alphaleft[3];
	for (int i = 0; i < 3; i++)
	{
		alphaleft[i] = dleft[i] * (1 + tau5 / (beta[i] + epsilon));
	}

	double alphal = alphaleft[0] + alphaleft[1] + alphaleft[2];

	double omegaleft[3];
	omegaleft[0] = alphaleft[0] / alphal;
	omegaleft[1] = alphaleft[1] / alphal;
	omegaleft[2] = alphaleft[2] / alphal;

	double left = omegaleft[0] * qleft[0] + omegaleft[1] * qleft[1] + omegaleft[2] * qleft[2];

	qright[0] = 1.0 / 3.0 * w + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	qright[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w + 1.0 / 3.0 * wp1;
	qright[2] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w;

	dright[0] = 0.3;
	dright[1] = 0.6;
	dright[2] = 0.1;

	double alpharight[3];

	for (int i = 0; i < 3; i++)
	{
		alpharight[i] = dright[i] * (1 + tau5 / (beta[i] + epsilon));
	}

	double alphar = alpharight[0] + alpharight[1] + alpharight[2];

	double omegaright[3];

	omegaright[0] = alpharight[0] / alphar;
	omegaright[1] = alpharight[1] / alphar;
	omegaright[2] = alpharight[2] / alphar;

	double right = omegaright[0] * qright[0] + omegaright[1] * qright[1] + omegaright[2] * qright[2];
	var = right;
	der2 = 6 * (left + right - 2 * w) / (h*h);
	der1 = (right - left) / h + 0.5*h*der2;

}


void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		g0reconstruction(interfaces[i], &fluids[i-1]);
	}
}
void Center_3rd(Interface1d &interfaces, Fluid1d *fluids)
{
	double w[3], wp1[3];

	//assume this is uniform mesh,
	//when mesh is no uniform, accuracy decrease
	double h = (fluids[0].dx + fluids[1].dx)/2.0;
	for (int i = 0; i < 3; i++)
	{
		w[i] = fluids[0].convar[i];
		wp1[i] = fluids[1].convar[i];
	}
	if (g0type == collisionn)
	{
		//have collision
		// to separate the reconstruction part from the flux part,
		// here we do a collision to get g0
		double convar_left[3], convar_right[3];
		for (int i = 0; i < 3; i++)
		{
			convar_left[i] = interfaces.left.convar[i];
			convar_right[i] = interfaces.right.convar[i];
		}
		
		double prim_left[3], prim_right[3];
		Convar_to_ULambda_1d(prim_left, convar_left);
		Convar_to_ULambda_1d(prim_right, convar_right);
		
		MMDF1d ml(prim_left[1], prim_left[2]);
		MMDF1d mr(prim_right[1], prim_right[2]);

		double unit[3]{ 1.0, 0.0, 0.0 };
		
		double gl[3], gr[3];
		GL(0, 0, gl, unit, ml);
		GR(0, 0, gr, unit, mr);

		for (int i = 0; i < 3; i++)
		{
			interfaces.center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
		}
		for (int i = 0; i < 3; i++)
		{
			interfaces.center.der1[i] = (wp1[i] - w[i]) / h;
			interfaces.center.der2[i] = 3 * ((wp1[i] + w[i]) - 2 * interfaces.center.convar[i]) / h / h;
		}
	}
	else
	{
		//collisionless
		//at this time, only second order accuracy obtained
		for (int i = 0; i < 4; i++)
		{
			interfaces.center.convar[i] = (wp1[i] + w[i]) / 2.0;
			interfaces.center.der1[i] = (wp1[i] - w[i]) / h;
			interfaces.center.der2[i] = 0.0;
		}
	}
}

void Center_5th(Interface1d &interfaces, Fluid1d *fluids)
{
	double wn1[3],w[3], wp1[3],wp2[3];

	//assume this is uniform mesh,
	//when mesh is no uniform, accuracy decrease
	double h = (fluids[0].dx + fluids[1].dx) / 2.0;
	for (int i = 0; i < 3; i++)
	{
		wn1[i] = fluids[-1].convar[i];
		w[i] = fluids[0].convar[i];
		wp1[i] = fluids[1].convar[i];
		wp2[i] = fluids[2].convar[i];
	}
	
	if (g0type == collisionn)
	{
		//have collision
		// to separate the reconstruction part from the flux part,
		// here we do a collision to get g0
		double convar_left[3], convar_right[3];
		for (int i = 0; i < 3; i++)
		{
			convar_left[i] = interfaces.left.convar[i];
			convar_right[i] = interfaces.right.convar[i];
		}

		double prim_left[3], prim_right[3];
		Convar_to_ULambda_1d(prim_left, convar_left);
		Convar_to_ULambda_1d(prim_right, convar_right);

		MMDF1d ml(prim_left[1], prim_left[2]);
		MMDF1d mr(prim_right[1], prim_right[2]);

		double unit[3]{ 1.0, 0.0, 0.0 };

		double gl[3], gr[3];
		GL(0, 0, gl, unit, ml);
		GR(0, 0, gr, unit, mr);

		for (int i = 0; i < 3; i++)
		{
			interfaces.center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
		}

		for (int i = 0; i < 4; i++)
		{
			interfaces.center.der1[i] = (-1.0 / 12.0*(wp2[i] - wn1[i]) + 5.0 / 4.0*(wp1[i] - w[i])) / h;
			interfaces.center.der2[i] = (-1.0 / 8.0*(wp2[i] + wn1[i]) + 31.0 / 8.0*(wp1[i] + w[i]) - 15.0 / 2.0*interfaces.center.convar[i]) / (h*h);

		}
	}
	else
	{
		//collisionless
		//at this time, only 4th order accuracy obtained
		for (int i = 0; i < 4; i++)
		{
			interfaces.center.convar[i] = (-1.0 / 12.0*(wp2[i] + wn1[i]) + 7.0 / 12.0*(wp1[i] + w[i]));
			interfaces.center.der1[i] = (-1.0 / 12.0*(wp2[i] - wn1[i]) + 5.0 / 4.0*(wp1[i] - w[i])) / h;
			interfaces.center.der2[i] = (1.0 / 2.0*(wp2[i] + wn1[i]) - 1.0 / 2.0*(wp1[i] + w[i])) / (h*h);
		}
	}
}

// to prepare the basic element for moment calculation
MMDF1d::MMDF1d(){ u = 1.0; lambda = 1.0; };
MMDF1d::MMDF1d(double u_in, double lambda_in)
{
	u = u_in;
	lambda = lambda_in;
	calcualte_MMDF1d();
}
void MMDF1d::calcualte_MMDF1d()
{
	uwhole[0] = 1;
	uwhole[1] = u;
	uplus[0] = 0.5*Alpha(lambda, -u);
	uminus[0] = 0.5*Alpha(lambda, u);
	uplus[1] = u*uplus[0] + 0.5*Beta(lambda, u);
	uminus[1] = u*uminus[0] - 0.5*Beta(lambda, u);
	for (int i = 2; i <= 9; i++)
	{
		uwhole[i] = u*uwhole[i - 1] + 0.5*(i - 1) / lambda*uwhole[i - 2];
	}
	for (int i = 2; i <= 9; i++)
	{
		uplus[i] = u*uplus[i - 1] + 0.5*(i - 1) / lambda*uplus[i - 2];
		uminus[i] = u*uminus[i - 1] + 0.5*(i - 1) / lambda*uminus[i - 2];
	}
	xi2 = 0.5*K / lambda;
	xi4 = 0.25*(K*K + 2 * K) / (lambda * lambda);
	//xi6?? how to calculate
	xi6 = 0.5*(K + 4) / lambda*xi4;

	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					upxi[i][k] = uplus[i];
					unxi[i][k] = uminus[i];
				}
				if (k == 1)
				{
					upxi[i][k] = uplus[i] * xi2;
					unxi[i][k] = uminus[i] * xi2;
				}
				if (k == 2)
				{
					upxi[i][k] = uplus[i] * xi4;
					unxi[i][k] = uminus[i] * xi4;
				}
				if (k == 3)
				{
					upxi[i][k] = uplus[i] * xi6;
					unxi[i][k] = uminus[i] * xi6;
				}
			}
		}
	}
	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					uxi[i][k] = uwhole[i];
				}
				if (k == 1)
				{
					uxi[i][k] = uwhole[i] * xi2;
				}
				if (k == 2)
				{
					uxi[i][k] = uwhole[i] * xi4;
				}
				if (k == 3)
				{
					uxi[i][k] = uwhole[i] * xi6;
				}
			}
		}
	}
}

double Alpha(double lambda, double u)
{
	return erfc(sqrt(lambda)*u);
}
double Beta(double lambda, double u)
{
	return exp(-lambda*u*u) / sqrt(pi*lambda);
}

//a general G function
void G(int no_u, int no_xi, double *psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.uxi[no_u][no_xi] + a[1] * m.uxi[no_u + 1][no_xi] + a[2] * 0.5*(m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.uxi[no_u + 1][no_xi] + a[1] * m.uxi[no_u + 2][no_xi] + a[2] * 0.5*(m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5*(a[0] * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]) +
		a[1] * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5*(m.uxi[no_u + 4][no_xi] + m.uxi[no_u][no_xi + 2] + 2 * m.uxi[no_u + 2][no_xi + 1]));

}
void GL(int no_u, int no_xi, double *psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.upxi[no_u][no_xi] + a[1] * m.upxi[no_u + 1][no_xi] + a[2] * 0.5*(m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.upxi[no_u + 1][no_xi] + a[1] * m.upxi[no_u + 2][no_xi] + a[2] * 0.5*(m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5*(a[0] * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]) +
		a[1] * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5*(m.upxi[no_u + 4][no_xi] + m.upxi[no_u][no_xi + 2] + 2 * m.upxi[no_u + 2][no_xi + 1]));

}
void GR(int no_u, int no_xi, double *psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.unxi[no_u][no_xi] + a[1] * m.unxi[no_u + 1][no_xi] + a[2] * 0.5*(m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.unxi[no_u + 1][no_xi] + a[1] * m.unxi[no_u + 2][no_xi] + a[2] * 0.5*(m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5*(a[0] * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]) +
		a[1] * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5*(m.unxi[no_u + 4][no_xi] + m.unxi[no_u][no_xi + 2] + 2 * m.unxi[no_u + 2][no_xi + 1]));

}

void Microslope(double *a, double der[3], double prim[3])
{
	double R4, R2;

	R4 = der[2] / prim[0] - 0.5*(prim[1] * prim[1] + 0.5*(K + 1) / prim[2])*der[0] / prim[0];
	R2 = (der[1] - prim[1] * der[0]) / prim[0];
	a[2] = 4 * prim[2] * prim[2] / (K + 1)*(2 * R4 - 2 * prim[1] * R2);
	a[1] = 2 * prim[2] * R2 - prim[1] * a[2];
	a[0] = der[0] / prim[0] - prim[1] * a[1] - 0.5*a[2] * (prim[1] * prim[1] + 0.5*(K + 1) / prim[2]);

}


void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d &block, int stage)
{
#pragma omp parallel  for 
	for (int i = block.ghost; i < block.nodex+block.ghost+1; i++)
	{
		flux_function (fluxes[i][stage], interfaces[i], block.dt);
		
	}
}

double Get_Tau_NS(double density0, double lambda0)
{
	if (tau_type == Euler)
	{
		return 0.0;
	}
	else if (tau_type == NS)
	{
		if (Mu > 0.0)
		{
			return 2.0*Mu*lambda0 / density0;
		}
		else if (Nu>0.0)
		{
			return 2 * Nu *lambda0;
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		return TauNS_Sutherland(density0, lambda0);
	}

}
double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt)
{
	if (tau_type == Euler)
	{
		if (c1_euler <= 0 && c2_euler <= 0)
		{
			return 0.0;
		}
		else
		{
			double C = c2_euler*abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
			if (C < 1)
			{
				return c1_euler*dt + dt*C;
			}
			else
			{
				return c1_euler*dt + dt;
			}
		}
	}
	else if (tau_type == NS)
	{
		if (Smooth == false)
		{
			double tau_n = 1.0*abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right)*dt;
			if (tau_n != tau_n)
			{
				tau_n = 0.0;
			}
			if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
			{
				return 0.0;
			}
			else
			{
				if (Mu > 0.0)
				{
					return tau_n + 2.0*Mu*lambda0 / density0;
				}
				else if (Nu>0.0)
				{
					return tau_n + 2 * Nu *lambda0;
				}
				else
				{
					return 0.0;
				}
			}

		}
		else
		{
			if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
			{
				return 0.0;
			}
			else
			{
				if (Mu > 0.0)
				{
					return 2.0*Mu*lambda0 / density0;
				}
				else if (Nu>0.0)
				{
					return 2 * Nu *lambda0;
				}
				else
				{
					return 0.0;
				}
			}
		}

	}
	else if (tau_type == Sutherland)
	{
		if (Smooth == false)
		{

			double tau_n = 1.0*abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right)*dt;
			if (tau_n != tau_n)
			{
				tau_n = 0.0;
			}
			if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
			{

				return 0.0;
			}
			else
			{
				if (Mu > 0.0)
				{
					return tau_n + TauNS_Sutherland(density0, lambda0);

				}
				else if (Nu>0.0)
				{
					return tau_n + 2 * Nu *lambda0;
				}
				else
				{
					return 0.0;
				}
			}

		}
		else
		{
			if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
			{
				return 0.0;
			}
			else
			{
				if (Mu > 0.0)
				{
					return TauNS_Sutherland(density0, lambda0);
				}
				else if (Nu>0.0)
				{
					return 2 * Nu *lambda0;
				}
				else
				{
					return 0.0;
				}
			}
		}
	}
	else
	{
		return 0.0;
	}
}
double TauNS_Sutherland(double density0, double lambda0)
{
	double T0 = 1.0 / 2.0 / (lambda0*R_gas);
	double mu = Mu*pow(T0 / T_inf, 1.5)*(T_inf + S) / (T0 + S);
	return 2.0*mu*lambda0 / density0;
}

void GKS(Flux1d &flux, Interface1d& interface, double dt)
{
	if (gks1dsolver == nothing)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	double Flux[2][3];
	//change conservative variables to rho u lambda
	double convar_left[3], convar_right[3], convar0[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interface.left.convar[i];
		convar_right[i] = interface.right.convar[i];
		convar0[i] = interface.center.convar[i];
		//cout << convar_left[i] << " " << convar_right[i] << " " << convar0[i] << " ";
	}
	//cout << endl;
	double prim_left[3], prim_right[3], prim0[3];
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);
	Convar_to_ULambda_1d(prim0, convar0);

	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[2]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	//t[0] = dt;
	//t[1] = dt;
	// non equ part time coefficient for gks_2nd algorithm
	t[0] = tau_num*(1 - eta); // this refers glu, gru part
	t[1] = tau_num*(eta*(dt + tau_num) - tau_num) + tau*tau_num*(eta - 1); //this refers aluu, aruu part
	t[2] = tau*tau_num*(eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd
	t[3] = tau_num*eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num*(tau_num - eta*(dt + tau_num) - tau*(eta - 1)) - dt*tau; //this refers a0uu part
	t[5] = 0.5*dt*dt - tau*tau_num*(eta - 1) - tau*dt; //this refers A0u part

	if (gks1dsolver == kfvs1st)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks1dsolver == kfvs2nd)
	{
		t[0] = dt;
		t[1] = -dt*dt / 2.0;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}
	//cout << "non-equ t " << t[0] <<" equ t "<<t[3]<< endl;
	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);

	double unit[3] = { 1, 0.0, 0.0 };

	double glu[3], gru[3];
	GL(1, 0, glu, unit, ml);
	GR(1, 0, gru, unit, mr);

	//only one part, the kfvs1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}

	if (gks1dsolver == kfvs1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i]=Flux[0][i]/dt;
		}
		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF1d m0(prim0[1], prim0[2]);

	double g0u[3];
	G(1, 0, g0u, unit, m0);

	//the equ g0u part, the gks1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
	}


	//cout << flux[0] << " " << flux[1] << " " << flux[3] << endl;

	if (gks1dsolver == gks1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i] / dt;
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1left[3], der1right[3];
	for (int i = 0; i < 3; i++)
	{
		der1left[i] = interface.left.der1[i];
		der1right[i] = interface.right.der1[i];
	}

	double alx[3];
	Microslope(alx, der1left, prim_left);

	double alxuul[3];
	GL(2, 0, alxuul, alx, ml);

	double arx[3];
	Microslope(arx, der1right, prim_right);
	double arxuur[3];
	GR(2, 0, arxuur, arx, mr);
	//cout << alx[0] << " " << alx[1] << " " << alx[2] << endl;
	for (int i = 0; i < 3; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
	}
	if (gks1dsolver == kfvs2nd)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i] / dt;
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd

	//for t[2] Aru,Alu part

	double alxu[4];
	double arxu[4];

	//take <u> moment for al, ar
	G(1, 0, alxu, alx, ml);
	G(1, 0, arxu, arx, mr);

	double Al[3], Ar[3];
	double der_AL[3], der_AR[3];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 3; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i]);
		der_AR[i] = -prim_right[0] * (arxu[i]);
	}
	// solve the coefficient martix b=ma
	Microslope(Al, der_AL, prim_left);
	Microslope(Ar, der_AR, prim_right);

	//to obtain the Alu and Aru
	double Alul[3];
	double Arur[3];
	GL(1, 0, Alul, Al, ml);
	GR(1, 0, Arur, Ar, mr);

	for (int i = 0; i < 3; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
	}

	// for t[4] a0xuu part

	double a0x[3];
	double der1[3];

	for (int i = 0; i < 3; i++)
	{
		der1[i] = interface.center.der1[i];
	}

	//solve the microslope
	Microslope(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[3];
	G(1, 0, a0xu, a0x, m0);
	//a0x <u^2> moment
	double a0xuu[3];
	G(2, 0, a0xuu, a0x, m0);

	for (int i = 0; i < 3; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i]);
	}
	

	// for t[5] A0u part
	double derA0[3];

	for (int i = 0; i < 3; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i]);
	}
	double A0[3];
	Microslope(A0, derA0, prim0);
	double A0u[3];
	G(1, 0, A0u, A0, m0);
	for (int i = 0; i < 3; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (gks1dsolver == gks2nd&&timecoe_list==S1O1)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i] / dt;
		}
		return;
	}
	if (gks1dsolver == gks2nd&&timecoe_list == S2O4)
	{
		double dt2 = 0.5*dt;
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num*(1 - eta); // this refers glu, gru part
		t[1] = tau_num*(eta*(dt2 + tau_num) - tau_num) + tau*tau_num*(eta - 1); //this refers aluu, aruu part
		t[2] = tau*tau_num*(eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num*eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num*(tau_num - eta*(dt2 + tau_num) - tau*(eta - 1)) - dt2*tau; //this refers a0uu part
		t[5] = 0.5*dt2*dt2 - tau*tau_num*(eta - 1) - tau*dt2; //this refers A0u part
		
		for (int i = 0; i < 3; i++)
		{
			// t0 part
			Flux[1][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
			// t1 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
			// t2 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}
		
		for (int i = 0; i < 3; i++)
		{		
			flux.f[i] = (4.0*Flux[1][i] - Flux[0][i]) / dt;
			flux.derf[i] = 4.0*(Flux[0][i] - 2.0*Flux[1][i]) / dt / dt;
		}
		//cout << "here" << endl;
		return;
	}
	else
	{
		cout << "no valid solver specify" << endl;
		exit(0);
	}
}

void ustarforHLLC(double d1, double u1, double p1, double s1, double star1, double *ustar)
{
	double tmp1, tmp2, tmp3;
	tmp1 = d1*(s1 - u1) / (s1 - star1);
	tmp2 = 0.5*(u1*u1) + p1 / ((r - 1.0)*d1);
	tmp3 = star1 + p1 / (d1*(s1 - u1));

	ustar[0] = tmp1;
	ustar[1] = tmp1*star1;
	ustar[2] = tmp1*(tmp2 + (star1 - u1)*tmp3);
}

void get_Euler_flux(double p[3], double *flux)
{
	flux[0] = p[0] * p[1];
	flux[1] = p[0] * p[1] * p[1] + p[2];
	double ENERGS = 0.5*(p[1] * p[1])*p[0] + p[2] / (r - 1);
	flux[2] = p[1] * (ENERGS + p[2]);
}



void HLLC(Flux1d &flux, Interface1d& interface, double dt)
{
	double pl[3], pr[3];
	Convar_to_primvar_1D(pl, interface.left.convar);
	Convar_to_primvar_1D(pr, interface.right.convar);

	double al, ar, pvars, pstar, tmp1, tmp2, qk, sl, sr, star;
	al = sqrt(r*pl[2] / pl[0]); //sound speed
	ar = sqrt(r*pr[2] / pr[0]);
	tmp1 = 0.5*(al + ar);         //avg of sound and density
	tmp2 = 0.5*(pl[0] + pr[0]);

	pvars = 0.5*(pl[2] + pr[2]) - 0.5*(pr[1] - pl[1])*tmp1*tmp2;
	pstar = Get_MAX(0.0, pvars);

	double flxtmp[3], qstar[3];

	if (pstar <= pl[2])
	{
		qk = 1.0;
	}
	else
	{
		tmp1 = (r + 1.0) / (2.0*r);
		tmp2 = (pstar / pl[2] - 1.0);
		qk = sqrt(1.0 + tmp1*tmp2);
	}
	sl = pl[1] - al*qk;

	if (pstar <= pr[2])
	{
		qk = 1.0;
	}
	else
	{
		tmp1 = (r + 1.0) / (2.0*r);
		tmp2 = (pstar / pr[2] - 1.0);
		qk = sqrt(1.0 + tmp1*tmp2);
	}
	sr = pr[1] + ar*qk;
	tmp1 = pr[2] - pl[2] + pl[0] * pl[1] * (sl - pl[1]) - pr[0] * pr[1] * (sr - pr[1]);
	tmp2 = pl[0] * (sl - pl[1]) - pr[0] * (sr - pr[1]);
	star = tmp1 / tmp2;

	if (sl >= 0.0) {
		get_Euler_flux(pl, flux.f);
	}
	else if (sr <= 0.0)
	{
		get_Euler_flux(pr, flux.f);
	}
	else if ((star >= 0.0) && (sl <= 0.0))
	{
		get_Euler_flux(pl, flxtmp);
		ustarforHLLC(pl[0], pl[1], pl[2], sl, star, qstar);

		for (int m = 0; m< 3; m++)
		{
			flux.f[m] = flxtmp[m] + sl*(qstar[m] - interface.left.convar[m]);
		}
	}
	else if ((star <= 0.0) && (sr >= 0.0))
	{
		get_Euler_flux(pr, flxtmp);
		ustarforHLLC(pr[0], pr[1], pr[2], sr, star, qstar);
		for (int m = 0; m < 4; m++)
		{
			flux.f[m] = flxtmp[m] + sr*(qstar[m] - interface.right.convar[m]);
		}
	}

}

void Initial_stages(Block1d &block)
{
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage 
		{
			for (int k = 0; k < 3; k++) //refers f, derf, der2f
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	timecoe_list(block);
}

void S1O1(Block1d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
}

void S1O2(Block1d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][0] = 0.5;
}

void S1O3(Block1d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][0] = 1.0/6.0;
}

void S2O4(Block1d &block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0/8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0/6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;
}

void RK4(Block1d &block)
{
	block.stages = 4;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;
	block.timecoefficient[2][2][0] = 1.0;
	block.timecoefficient[3][0][0] = 1.0 / 6.0;
	block.timecoefficient[3][1][0] = 1.0 / 3.0;
	block.timecoefficient[3][2][0] = 1.0 / 3.0;
	block.timecoefficient[3][3][0] = 1.0 / 6.0;
}

void Update(Fluid1d *fluids, Flux1d **fluxes, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i <block.nodex+block.ghost; i++)
	{
		for (int j = 0; j <3; j++)
		{
			fluids[i].convar[j] = fluids[i].convar_old[j] + 1.0 / fluids[i].dx*(fluxes[i][0].F[j] - fluxes[i + 1][0].F[j]);
		}
	}
}

void Update(Fluid1d *fluids, Flux1d **fluxes, Block1d block, int stage)
{
	if (stage>block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost+1; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double Flux=0.0;
			for (int k = 0; k < stage+1; k++)
			{
				Flux = Flux
					+ block.timecoefficient[stage][k][0] * dt*fluxes[i][k].f[j]
					+ block.timecoefficient[stage][k][1] * dt*dt*fluxes[i][k].derf[j]
					+ block.timecoefficient[stage][k][2] * dt*dt*dt*fluxes[i][k].der2f[j];
				//cout << block.timecoefficient[stage][k][0] << block.timecoefficient[stage][k][1] << block.timecoefficient[stage][k][2] << endl;

			}
			fluxes[i][stage].F[j] = Flux;
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i <block.nodex + block.ghost; i++)
	{
		for (int j = 0; j <3; j++)
		{
			fluids[i].convar[j] = fluids[i].convar_old[j] + 1.0 / fluids[i].dx*(fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
		}
	}
}