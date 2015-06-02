/*

	Solves the "H" equation of the hidrogen molecule ion

*/

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef H_EQ_SOLVER 
#define H_EQ_SOLVER

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

class H_eq_solver{

//////////////////////////////////////////////////////////
//Data
//////////////////////////////////////////////////////////
public:
	int N; //number of steps
	double step; //step size

	bool symm; //symmetry of solution

	//equation eigenvalues
	int m; 
	double lam2,u;


//////////////////////////////////////////////////////////
//some helper functions
//////////////////////////////////////////////////////////

	//find root with linear interpolation
	std::vector<double> root_lin_interpol(std::vector<double> x, std::vector<double> y){
		std::vector<double> roots(0);
		for (size_t i=1;i<x.size();i++){
        		if((y[i-1]>0 && y[i]< 0) || (y[i-1]<0 && y[i]> 0)){
            			roots.push_back(x[i-1]+(x[i]-x[i-1])*(-1*y[i-1])/(y[i]-y[i-1]));
			}
		}
    		return roots;
	}

//////////////////////////////////////////////////////////
//API
//////////////////////////////////////////////////////////
public:

	//////////////////////////////////////////////////////////
	// constructor
	//////////////////////////////////////////////////////////
	H_eq_solver(bool symm, double lam2, double u, int m, int N){
		this->symm=symm;
		this->lam2=lam2;
		this->u=u;
		this->m=m;
		this->N=N;

		this->step=1.0/ (double) N;

		return;
	}


	//////////////////////////////////////////////////////////
	// setters
	//////////////////////////////////////////////////////////
	void set_lam2(double lam2){
		this->lam2=lam2;
		return;
	}
	void set_u(double u){
		this->u=u;
		return;
	}
	
	//////////////////////////////////////////////////////////
	// solves the equation
	//////////////////////////////////////////////////////////
	double solve(bool verbose, bool write_res, std::string filename){
		if(verbose){
			std::cerr<<"solving H equation..."<<std::endl;
		}

		std::vector<double> x(this->N);
		for (int i=0;i<this->N;i++){
			x[i]=-1+i*this->step+this->step;
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

		std::vector<double> H(this->N);
		for (int i=0;i<this->N;i++){
			H[i]= pow((1-x[i]*x[i]),this->m/2.0) * y0[i];
		}

		double err;
		if(this->symm){
			err=(H[this->N-1]-H[this->N-2])/this->step;
		}
		else{
			err=H[this->N-1];
		}
		if(verbose){
			std::cerr<<"Solution accurate to: "<<err<<std::endl;
		}

		if(write_res){
			std::cerr<<"writng solution to: "<<filename<<"\n"<<std::endl;
			std::ofstream myfile;
			myfile.open(filename);
  			myfile << "x\tH\n";
			for (int i=0;i<this->N;i++){
  				myfile << x[i]<<"\t"<<H[i]<<"\n";
			}
  			myfile.close();	
		}

		return err;
	}
		
	//////////////////////////////////////////////////////////
	// find all u in range given lam2
	//////////////////////////////////////////////////////////
	std::vector<double> find_u_range(double u_min, double u_max, int grid_size,
					 bool verbose, bool write_res, std::string filename){
                if(verbose){
                        std::cerr<<"finding u in a range...\n"<<std::endl;
                }
		std::vector<double> us(grid_size);
		std::vector<double> errs(grid_size);
		for (size_t i=0;i<(size_t)grid_size;i++){
			us[i]=u_min+i*(u_max-u_min)/grid_size;
			this->u=us[i];
			errs[i]=this->solve(false,false,"");
		}

		std::vector<double> roots=root_lin_interpol(us,errs);
		if( roots.size() < 1){
			std::cerr<<"Solution not found for u in this u range!"<<std::endl;
			exit(1);
		}

		if(write_res){
                        std::cerr<<"writng errors to: "<<filename<<"\n"<<std::endl;
                        std::ofstream myfile;
                        myfile.open(filename);
                        myfile << "u\terror\n";
                        for (size_t i=0;i<us.size();i++){
                                myfile << us[i]<<"\t"<<errs[i]<<"\n";
                        }
                        myfile.close();
                }

		return roots;
	}

	//////////////////////////////////////////////////////////
	// find u around original guess given lam2
	//////////////////////////////////////////////////////////
	double find_u(double guess0, double eps_err,double d){
		double err1,err2;
		err1=err2=1.0;
		double guess1=guess0;
		double guess2=guess1+d;
		int i=0;
		while(pow(err1*err1,0.5)>eps_err && i<100){
			i++;
			this->u=guess1;
			err1=this->solve(false,false,"");
			this->u=guess2;
			err2=this->solve(false,false,"");
			guess1=guess1 - d*err1/(err2-err1);
			guess2=guess1+d;
			//std::cout<<guess1<<" "<<err1<<"\n";
			//std::cout<<guess2<<" "<<err2<<"\n\n";
		}
		//std::cout<<"done\n";
		
		return guess1;
	}


	//////////////////////////////////////////////////////////
	// find u given lam2
	//////////////////////////////////////////////////////////
	double find_u_i(double u_min, double u_max, int grid_size,
			int u_i, bool verbose){
		//get i_i approx root with linear interpol
		double root_guess = find_u_range( u_min, u_max,  grid_size, verbose,false,"")[u_i];
		//calculate more precise root with Newton-Rapshon
		double root_accurate = find_u(root_guess,1e-4,1e-5);
		return root_accurate;
	}


	//////////////////////////////////////////////////////////
	// find lam2-u relationship
	//////////////////////////////////////////////////////////
	std::vector<std::vector<double> > lam_u_i(double lam_min, double lam_max, int lam_grid_size, 
			double u_min, double u_max, int u_grid_size,
			int u_i, bool write_res, std::string filename, bool verbose){

                if(verbose){
                        std::cerr<<"calibrating lam2-u relationship...\n"<<std::endl;
                }

		std::vector<double> lams(lam_grid_size);
		std::vector<double> us(lam_grid_size);

		for (size_t i=0;i<(size_t)lam_grid_size;i++){
			lams[i]=lam_min+i*(lam_max-lam_min)/lam_grid_size;
			this->lam2=lams[i];
			us[i]=this->find_u_i(u_min,u_max,u_grid_size,u_i,false);
			//std::cout<<lams[i]<<" "<<us[i]<<std::endl;
		}

		if(write_res){
			std::cerr<<"writng calibration to: "<<filename<<"\n"<<std::endl;
			std::ofstream myfile;
			myfile.open(filename);
  			myfile << "lam2\tu\n";
			for (int i=0;i<lam_grid_size;i++){
  				myfile << lams[i]<<"\t"<<us[i]<<"\n";
			}
  			myfile.close();	
		}
		return {lams,us};
	}
		

//////////////////////////////////////////////////////////
// Inside job
//////////////////////////////////////////////////////////
private:

	//derivative function
	std::vector<double> fun(double x, std::vector<double> y){
		std::vector<double> res;
		res.push_back(y[1]);
		res.push_back( ( 1.0/ (1.0-x*x) ) * (
				(2.0*x*(this->m+1)*y[1]) +
				(this->m*(this->m+1) - this->u +this->lam2*x*x)*y[0] ) );
		return res;
	}

	//initial slope
	double slope(){
		if (this->symm){
			return ( this->m*(this->m+1) + this->lam2 - this->u ) / (2 * (this->m+1 ) );
		}
		else{
			return -1*( this->m*(this->m+1) + this->lam2 - this->u ) / (2 * (this->m+1 ) );
		}
	}
	
	
	//initial value
	std::vector<double>  init_val(){
		std::vector<double> res;
		if (this->symm){
			res.push_back( 1 + this->step * this->slope());
			res.push_back( this->slope());
		}
		else{
			res.push_back( - 1 + this->step * this->slope());
			res.push_back( this->slope());
		}
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
