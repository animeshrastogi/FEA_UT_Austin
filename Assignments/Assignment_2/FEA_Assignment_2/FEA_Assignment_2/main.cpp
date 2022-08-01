#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//Define the FEA class
class oned_finite_element{
    double h;   //Element Size
    double J;   //Jacobian
    double kspring;     //Spring constant
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
    void initialize();
    double shape_fun(double, int);
    double gauss_point(int);
    void assemble();
    void solve();
    void write_output();
};

//Function to initialize the dependent parameters and the matrices/vectors
void oned_finite_element::initialize(){
    h = L/nel;
    J = h/2;
    kspring = 4*E*A*L;
    nnodes = nel+1;
    K_global.resize(nnodes,nnodes);
    F_global.resize(nnodes,1);
    K_global.setZero();
    F_global.setZero();
}

//Linear Shape Function
double oned_finite_element::shape_fun(double x, int index){
    if(index == 0){
        return 0.5 - 0.5*x;
    }
    else{
        return 0.5 + 0.5*x;
    }
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
    
    //Shape function derivatives [for global stiffness]
    Matrix<double, 2, 1> der_shape_function;
    der_shape_function(0) = -0.5;
    der_shape_function(1) = 0.5;
    
    double gp;  //Gauss-point coordinate
    
    Matrix<double, 2, 2> der_shape_function_matrix;
    der_shape_function_matrix = der_shape_function*der_shape_function.transpose();
    
    //Element loop
    for(int i = 0; i<nel; i++){
        //Gauss-point loop
        for(int l = 0; l<2; l++){
            gp = gauss_point(l);    //Get the gauss point coordinate
        for(int j=0;j<2;j++){
            for (int k=0;k<2;k++){
                //Assemble the global stiffness
                    K_global(i+j,i+k) = K_global(i+j,i+k) + (E*A*der_shape_function_matrix(j,k))/J;
            }
            //Assemble the global force vector
            F_global(i+j) = F_global(i+j) + Force*J*shape_fun(gp,j);
            }
        }
    }
}
void oned_finite_element::solve(){
    cout<<"Solving the system...\n";
    //Apply boundary condition
    K_global(nnodes-1,nnodes-1) = K_global(nnodes-1,nnodes-1) +kspring;
    //Solve the system using direct solver [QR decomposition]
    u_all = K_global.colPivHouseholderQr().solve(F_global);
    
    //Display the solution
    cout<<K_global<<"\n";
    cout<<F_global<<"\n";
    cout<<u_all<<"\n";
}

int main(int argc, const char * argv[]) {
    oned_finite_element fea;
    fea.nel = 2;
    fea.Force = 1;
    fea.A = 1.0;
    fea.E = 1.0;
    fea.L = 1.0;
    fea.initialize();
    fea.assemble();
    fea.solve();
    return 0;
}
