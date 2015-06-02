#include "H_eq_solver.h"
#include "R_eq_solver.h"

int main(int argc, char** argv ){

	///////////////////////////////////////
	//set params
	///////////////////////////////////////

	//distance of nuclei
	double r_ab=2.0;

	//check
	if( argc != 5){
		std::cerr<<"4 command line arguments, check sorce,  sample notebook"<<std::endl;
		exit(1);
	}

	//quantum numbers
	int m=atoi(argv[1]);
	int u_i=atoi(argv[2]);
	//int E_i=atoi(argv[3]);

	//important param	
	std::string symmstr(argv[4]);
	bool symm;
	if(symmstr=="symm"){
		symm=true;
	}
	else if ( symmstr=="antisymm" ){
		symm=false;
	}
	else{
		exit(1);
	}

	///////////////////////////////////////
	// calibrate lam2-u relationship
	///////////////////////////////////////

	//initialize H equation solver
	double lam2_0=0;
	double u_0=0;
	int N=100; //number of points in integration
	H_eq_solver my_heq_solver(symm,lam2_0,u_0,m,N);

	//get lam u curve
	double u_min=-30;
	double u_max=150;
	double u_grid_size=200;
	bool verbose=true;
	bool write_res=true;
	my_heq_solver.find_u_range( u_min,u_max,u_grid_size, verbose,write_res,"u_err.csv");	





	//calibrate lam2-u relationship
	double lam_min=-20;
	double lam_max=50;
	double lam_grid_size=40;
	u_min=-30;
	u_max=150;
	u_grid_size=40;
	verbose=true;
	write_res=true;
	std::vector<std::vector<double> > calib;
	calib = my_heq_solver.lam_u_i(lam_min,lam_max,lam_grid_size,
					u_min,u_max,u_grid_size,
					u_i,write_res,"calib.csv",verbose);	


	///////////////////////////////////////
	// Find Energy
	///////////////////////////////////////

	//initialize R-eq solver
	double E_0=0;
	N=1100;
	R_eq_solver my_req_solver(lam2_0,u_0,m,r_ab,E_0,calib[0],calib[1],N);

	//find energy
	double E_min=-3.0;
	double E_max=4.0;
	int E_grid_size=1000;
	my_req_solver.find_E_range(E_min,E_max,E_grid_size,verbose,write_res,"E_err.csv");


}

