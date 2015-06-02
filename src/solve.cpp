#include "H_eq_solver.h"
#include "R_eq_solver.h"

int main(int argc, char** argv ){

	///////////////////////////////////////
	//set params
	///////////////////////////////////////

	//check
        if( argc != 5){
                std::cerr<<"4 command line arguments, check sorce,  sample notebook"<<std::endl;
                exit(1);
        }


	//distance of nuclei
	double r_ab=2.0;

	//quantum numbers
	int m=atoi(argv[1]);
	int u_i=atoi(argv[2]);
	int E_i=atoi(argv[3]);

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

	//calibrate lam2-u relationship
	double lam_min=-30;
	double lam_max=30;
	double lam_grid_size=101;
	double u_min=-50;
	double u_max=150;
	double u_grid_size=100;
	bool verbose=true;
	bool write_result=true;
	std::vector<std::vector<double> > calib;
	calib = my_heq_solver.lam_u_i(lam_min,lam_max,lam_grid_size,
					u_min,u_max,u_grid_size,
					u_i,write_result,"calib.csv",verbose);	


	///////////////////////////////////////
	// Find Energy
	///////////////////////////////////////

	//initialize R-eq solver
	double E_0=0;
	N=1100; //number of points in integration
	R_eq_solver my_req_solver(lam2_0,u_0,m,r_ab,E_0,calib[0],calib[1],N);

	//find energy
	double E_min=-3.0;
	double E_max=10.0;
	int E_grid_size=1000;
	double E= my_req_solver.find_E_i(E_min,E_max,E_grid_size,E_i,verbose);


	///////////////////////////////////////
	//Solve equations
	///////////////////////////////////////

	//solve, and save R equation
	my_req_solver.set_E(E);	
	my_req_solver.solve(verbose,write_result,"r_res.csv");	

	//set appropriate params to H-eq solver
	my_heq_solver.set_lam2(my_req_solver.get_lam2());
	my_heq_solver.set_u(my_req_solver.get_u());

	//solve, and save H equation
	my_heq_solver.solve(verbose,write_result,"h_res.csv");	
}

