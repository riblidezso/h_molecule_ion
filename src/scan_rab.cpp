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
	double lam_min=-100;
	double lam_max=100;
	double lam_grid_size=301;
	double u_min=-100;
	double u_max=100;
	double u_grid_size=201;
	bool verbose=false;
	bool write_result=true;
	std::vector<std::vector<double> > calib;
	calib = my_heq_solver.lam_u_i(lam_min,lam_max,lam_grid_size,
					u_min,u_max,u_grid_size,
					u_i,write_result,"calib.csv",verbose);	


	///////////////////////////////////////
	// Find Energy-r_ab relationship
	///////////////////////////////////////

	//initialize R-eq solver
	double r_ab_0=2.0;
	double E_0=0;
	N=1100;
	R_eq_solver my_req_solver(lam2_0,u_0,m,r_ab_0,E_0,calib[0],calib[1],N);

	//find energy
	double E_min=-3.0;
	double E_max=0.0;
	int E_grid_size=100;
	int r_grid=50;
	double r_ab_min=0.8;
	double r_ab_max=10.0;

	std::vector<double> r_abs,Es,E_fulls;
	double r_ab,E;

	for (int i=0;i<r_grid;i++){
		r_ab=r_ab_min+i*(r_ab_max-r_ab_min)/r_grid;
		my_req_solver.set_r_ab(r_ab);
		E = my_req_solver.find_E_i(E_min,E_max,E_grid_size,E_i,verbose);
		r_abs.push_back(r_ab);
		Es.push_back(E);
		E_fulls.push_back(E+2.0/r_ab);
	}
	
	std::cerr<<"writng r_ab-E to: "<<"r_ab_E.csv"<<"\n"<<std::endl;
	std::ofstream myfile;
	myfile.open("r_ab_E.csv");
	myfile << "r_ab\tE_e\tE_full\n";
	for (int i=0;i<r_grid;i++){
	        myfile << r_abs[i]<<"\t"<<Es[i]<<"\t"<<E_fulls[i]<<"\n";
	}
	myfile.close();

}

