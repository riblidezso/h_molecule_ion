/*

	Solves the "R" equation of the hidrogen molecule ion

*/

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef R_EQ_SOLVER 
#define R_EQ_SOLVER

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <random>
#include <chrono>
#include <thread>
#include <string.h>

#include <initializer_list>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "H_eq_solver.h"

class R_eq_solver{

//////////////////////////////////////////////////////////
//Data
//////////////////////////////////////////////////////////
public:
	int N; //number of steps
	double step; //step size

	//equation eigenvalues
	int m; 
	double lam2,u;

	//not eigenvalue paramaters
	double r_ab,E;

	//lam-mu calibration vector
	std::vector<double> lam_vec,u_vec;


//////////////////////////////////////////////////////////
//some helper functions
//////////////////////////////////////////////////////////

	//find root with linear interpolation
	std::vector<double> root_lin_interpol(std::vector<double> x, std::vector<double> y){
		std::vector<double> roots;
		for (size_t i=1;i<x.size();i++){
        		if((y[i-1]>0 && y[i]< 0) || (y[i-1]<0 && y[i]> 0)){
            			roots.push_back(x[i-1]+(x[i]-x[i-1])*(-1*y[i-1])/(y[i]-y[i-1]));
			}
		}
    		return roots;
	}


	//find value with linear interpolation
	double lin_interpol(double x0, std::vector<double> x, std::vector<double> y){
		for (size_t i=0;i<x.size();i++){
			if(x[i-1] == x0){
				return y[i-1];
			}
        		if( (x[i-1]< x0 && x[i]> x0)){
            			return y[i-1]+(x0-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1]);
			}
		}
		std::cerr<<"Error cant extrapolate "<<x0<<std::endl;
		exit(1);
	}



//////////////////////////////////////////////////////////
//API
//////////////////////////////////////////////////////////
public:

	//////////////////////////////////////////////////////////
	// constructor
	//////////////////////////////////////////////////////////
	R_eq_solver(double lam2, double u, int m, double r_ab, double E, std::vector<double> lam_vec, std::vector<double> u_vec, int N){
		this->lam2=lam2;
		this->u=u;
		this->m=m;
		this->r_ab=r_ab;
		this->E=E;
		this->lam_vec=lam_vec;
		this->u_vec=u_vec;
		this->N=N;

		this->step=11.0/ (double) N;

		return;
	}


	//////////////////////////////////////////////////////////
	// setters 
	//////////////////////////////////////////////////////////
	void set_r_ab(double r_ab){
		this->r_ab=r_ab;
		this->lam2=(this->E * this->r_ab*this->r_ab) /4.0;
		this->u=lin_interpol(this->lam2,lam_vec,u_vec);
		return;
	}
	void set_E(double E){
		this->E=E;
		this->lam2=(this->E * this->r_ab*this->r_ab) /4.0;
		this->u=lin_interpol(this->lam2,lam_vec,u_vec);
		return;
	}





        //////////////////////////////////////////////////////////
        // getters
        //////////////////////////////////////////////////////////
        double get_lam2(){
                return this->lam2;
        }
        double get_u(){
                return this->u;
        }



	
	//////////////////////////////////////////////////////////
	// solves the equation
	//////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////
	// solves the equation
	//////////////////////////////////////////////////////////
	double solve(bool verbose, bool write_res, std::string filename){
                if(verbose){
                        std::cerr<<"solving R equation..."<<std::endl;
                }

		std::vector<double> x(this->N);
		for (int i=0;i<this->N;i++){
			x[i]=1+i*this->step+this->step;
		}
		std::vector<double> y0(this->N);
		std::vector<double> y1(this->N);
		for (int i=0;i<this->N;i++){
			y0[i]=0;
			y1[i]=0;
		}
		
		std::vector<double> temp;
		temp = this->init_val();
		y0[0]=temp[0];
		y1[0]=temp[1];

		for (int i=1;i<this->N;i++){
			temp = this->runge_kutta_4_fix(x[i-1],temp);
			y0[i]=temp[0];
			y1[i]=temp[1];
		}

		std::vector<double> R(this->N);
		for (int i=0;i<this->N;i++){
			R[i]= pow((x[i]*x[i]-1),this->m/2.0) * y0[i];
		}

		double err = R[this->N-1];
		if(verbose){
			std::cerr<<"Solution accurate to: "<<err<<std::endl;
		}

		if(write_res){
			std::cerr<<"writng solution to: "<<filename<<"\n"<<std::endl;
			std::ofstream myfile;
			myfile.open(filename);
  			myfile << "x\tR\n";
			for (int i=0;i<this->N;i++){
  				myfile << x[i]<<"\t"<<R[i]<<"\n";
			}
  			myfile.close();	
		}
		return err;
	}
		

	//////////////////////////////////////////////////////////
	// find E in E range
	//////////////////////////////////////////////////////////
	std::vector<double> find_E_range(double E_min, double E_max, int E_grid_size,
					bool verbose,bool write_res, std::string filename){
		if(verbose){
                        std::cerr<<"finding E in a range...\n"<<std::endl;
                }

		std::vector<double> Es(E_grid_size);
		std::vector<double> errs(E_grid_size);
		for (size_t i=0;i<(size_t)E_grid_size;i++){
			Es[i]=E_min+i*(E_max-E_min)/E_grid_size;
			this->set_E(Es[i]);
			errs[i]=this->solve(false,false,"");
			//if(verbose){
			//	std::cout<<Es[i]<<" "<<errs[i]<<std::endl;
			//}
		}

		std::vector<double> roots=root_lin_interpol(Es,errs);
		if( roots.size() < 1){
			std::cerr<<"Solution not found for E in this E range!"<<std::endl;
			exit(1);
		}

		if(write_res){
                        std::cerr<<"writng errors to: "<<filename<<"\n"<<std::endl;
                        std::ofstream myfile;
                        myfile.open(filename);
                        myfile << "E\terror\n";
                        for (size_t i=0;i<Es.size();i++){
                                myfile << Es[i]<<"\t"<<errs[i]<<"\n";
                        }
                        myfile.close();
                }


		return roots;
	}

	//////////////////////////////////////////////////////////
	// find E around original guess 
	//////////////////////////////////////////////////////////
	double find_E(double guess0, double eps_err,double d){
		double err1,err2;
		err1=err2=1.0;
		double guess1=guess0;
		double guess2=guess1+d;
		int i=0;
		while(pow(err1*err1,0.5)>eps_err && i<100){
			i++;
			this->set_E(guess1);
			err1=this->solve(false,false,"");
			this->set_E(guess2);
			err2=this->solve(false,false,"");
			guess1=guess1 - d*err1/(err2-err1);
			guess2=guess1+d;
			//std::cout<<guess1<<" "<<err1<<"\n";
			//std::cout<<guess2<<" "<<err2<<"\n\n";
		}
		if ( i==100 ){
			std::cout<<"last iteration in E calculation!\n";
			//exit(1);
		}

		//std::cout<<"done\n";
		
		return guess1;
	}



	//////////////////////////////////////////////////////////
	// find E_i
	//////////////////////////////////////////////////////////
	double find_E_i(double E_min, double E_max, int grid_size,
			int E_i, bool verbose){
		//get i_i approx root with linear interpol
		double root_guess = find_E_range( E_min, E_max,  grid_size, verbose,false,"")[E_i];
		//calculate more precise root with Newton-Rapshon
		double root_accurate = find_E(root_guess,1e-4,1e-5);
		if( verbose){
			std::cout<<"E found: "<<root_accurate<<std::endl;
		}
		return root_accurate;
		
	}



//////////////////////////////////////////////////////////
// Inside job
//////////////////////////////////////////////////////////
private:

	//derivative function
	std::vector<double> fun(double x, std::vector<double> y){
		std::vector<double> res;
		res.push_back(y[1]);
		res.push_back( ( 1.0/ (x*x - 1.0) ) * (
				 - 2.0*x*(this->m+1)*y[1] - 
				(2.0*this->r_ab*x + this->m*(this->m+1) - this->u +this->lam2*x*x)*y[0] ) );
		return res;
	}

	//initial slope
	double slope(){
		return -1*( 2*this->r_ab + this->m*(this->m+1) + this->lam2 - this->u ) / (2 * (this->m+1 ) );
	}
	
	
	//initial value
	std::vector<double>  init_val(){
		std::vector<double> res;
		res.push_back( 1 + this->step * this->slope());
		res.push_back( this->slope());
		return res;
	}
	
	//4th order runge kutta, fix stepsize	
	std::vector<double> runge_kutta_4_fix(double x, std::vector<double> y){
		std::vector<double> k1=this->fun(x,y);

		std::vector<double> y2=std::vector<double>(y);
		y2[0]+= k1[0] * (this->step/2.0);
		y2[1]+= k1[1] * (this->step/2.0);
		std::vector<double> k2=this->fun(x+this->step/2.0,y2);

		std::vector<double> y3=std::vector<double>(y);
		y3[0]+= k2[0] * (this->step/2.0);
		y3[1]+= k2[1] * (this->step/2.0);
		std::vector<double> k3=this->fun(x+this->step/2.0,y3);

		std::vector<double> y4=std::vector<double>(y);
		y4[0]+= k3[0] * (this->step);
		y4[1]+= k3[1] * (this->step);
		std::vector<double> k4=this->fun(x+this->step,y4);

		std::vector<double> res=std::vector<double>(y);;
		res[0]+= (this->step/6) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0] );
		res[1]+= (this->step/6) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1] );

		return res;
	}


};



#endif
