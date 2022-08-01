//Author: Animesh Rastogi
//Email: animesh.rastogi@utexas.edu
//M.S. student - Civil Engineering [Computational Mechanics]
//The University of Texas at Austin

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>
using namespace std;
using namespace Eigen;

//Define the FEA class
class oned_finite_element{
    double h;   //Element Size
    double J;   //Jacobian
    double kwinkler;     //Winkler constant
    int nnodes;         //Number of nodes
    MatrixXd K_global;  //Global stiffness matrix
    MatrixXd F_global;  //Global Force Vector
    MatrixXd u_all;     //Solution Vector
public:
    double Force;   //Applied Force
    int nel;        //Number of elements
    double L;       //Lenght of the bar
    double E;       //Young's modulus
    double A;       //Area of cross-section
    oned_finite_element();
    void initialize();
    double shape_fun(double, int);
    double shape_fun(double, int, int);
    double der_shape_fun(double, int);
    double der_shape_fun(double, int, int);
    double gauss_point(int);
    void assemble();
    void solve();
    void write_output();
};

//Constructor
oned_finite_element::oned_finite_element(){
    Force = 1.0;
    nel = 2;
    L = 1.0;
    E = 1.0;
    A = 1.0;
    kwinkler = 1.0;
}

//Function to initialize the dependent parameters and the matrices/vectors
void oned_finite_element::initialize(){
    
    cout<<"Initializing FEA system\n";
    h = L/nel;
    J = h/2;
    nnodes = 2*nel+1;
    K_global.resize(nnodes,nnodes);
    F_global.resize(nnodes,1);
    K_global.setZero();
    F_global.setZero();
    cout<<"Number of elements: "<<nel<<"\t"<<"Degrees of Freedom: "<< nnodes<<"\n";
}

//Quadratic Shape Function
double oned_finite_element::shape_fun(double x, int j){
    if(j == 0){
        return -0.5*(x*(1-x));
    }
    else if(j == 1){
        return 1-(x*x);
    }
    else{
        return 0.5*(x*(1+x));
    }
}

double oned_finite_element::shape_fun(double x, int j, int k){
    
    Matrix<double, 3, 1> shape_function;
    shape_function(0) = shape_fun(x, 0);
    shape_function(1) = shape_fun(x, 1);
    shape_function(2) = shape_fun(x, 2);
    Matrix<double, 3, 3> shape_function_matrix = shape_function*shape_function.transpose();
    return shape_function_matrix(j,k);
    
}

//Quadratic Shape Function Derivatives
double oned_finite_element::der_shape_fun(double x, int j){
    if(j == 0){
        return x-0.5;
    }
    else if(j == 1){
        return -2*x;
    }
    else{
        return x+0.5;
    }
}

double oned_finite_element::der_shape_fun(double x, int j, int k){
    Matrix<double, 3, 1> der_shape_function;
    der_shape_function(0) = der_shape_fun(x, 0);
    der_shape_function(1) = der_shape_fun(x, 1);
    der_shape_function(2) = der_shape_fun(x, 2);
    Matrix<double, 3, 3> der_shape_function_matrix = der_shape_function*der_shape_function.transpose();
    return der_shape_function_matrix(j,k);
}


//Gauss-Point coordinate [2 point quadrature rule]
double oned_finite_element::gauss_point(int index){
    double gp;
    if(index == 0)
        gp = -1.0/sqrt(3);
    else
        gp = 1/sqrt(3);
    return gp;
}

//Assemble the global stiffness matrix
void oned_finite_element::assemble(){
    cout<<"Assembling the global stiffness and force vector...\n";
    
    double gp;  //Gauss-point coordinate
    double fval;
    
    //Element loop
    for(int i = 0; i<nel; i++){
        //Gauss-point loop
        for(int l = 0; l<2; l++){
            gp = gauss_point(l);    //Get the gauss point coordinate
        for(int j=0;j<3;j++){
            for (int k=0;k<3;k++){
                //Assemble the global stiffness
                    K_global(2*i+j,2*i+k) = K_global(2*i+j,2*i+k) + (E*A*der_shape_fun(gp,j,k))/J + kwinkler*J*shape_fun(gp, j, k);
            }
            //Assemble the global force vector
            fval = h/2*(shape_fun(gp,0)*(2*i) + shape_fun(gp,1)*(2*i+1)+
                        shape_fun(gp,2)*(2*i+2));
            F_global(2*i+j) = F_global(2*i+j) + (Force/L)*J*fval*shape_fun(gp,j);
            }
        }
    }
}

void oned_finite_element::solve(){
    cout<<"Solving the system...\n";
    //Solve for unconstraint degree of freedom
    MatrixXd K_constraint;
    MatrixXd F_constraint;
    MatrixXd u_constraint;
    K_constraint = K_global.block(1,1,nnodes-2,nnodes-2);
    F_constraint = F_global.block(1,0,nnodes-2,1);
    u_constraint = K_constraint.colPivHouseholderQr().solve(F_constraint);
    cout<<"u_constraint\n"<<u_constraint<<"\n";
}

int main(int argc, const char * argv[]) {


    time_t tstart, tend;
    tstart = time(0);
    oned_finite_element fea;
    fea.nel = 2;
    fea.initialize();
    fea.assemble();
    fea.solve();
    tend = time(0);
    cout << "Total wall time "<< difftime(tend, tstart) <<" second(s)."<< endl;
    return 0;
}
