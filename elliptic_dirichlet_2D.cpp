#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

// defining pi
double pi = 2*asin(1);

// function for the Poisson equation
double f(double x, double y) {
    return 0;
}
    
// functions for the boundary conditions 
double gy_xRight(double y) {
    return sin(y*pi);
}

double gy_xLeft(double x) {
    return sin(x*pi);
}

double gx_yBottom(double x) {
    return sin(x*pi);
}

double gx_yTop(double x) {
    return sin(x*pi);
}

// discretization data
const int m = 3; // number of nodes in the x-axis
const int n = 3; // number of nodes in the y-axis
const int size = m*n; 


// LU solver function definition
// this is needed for solving the later linear system of equations
// this function implements the LU method to do that

double* LU_solver(double * M, double b[]){
    double L[m*n][m*n], U[m*n][m*n], z[m*n];
    static double x[m*n];
    double sum0,sumL,sumU,sumz,sumx;
    
    for(int i=0; i<=m*n-1;i++)for(int j=0; j<=m*n-1;j++){
        if (i<=j){   // Conditional for constructing U
            if (i<j) // Conditional for construcing the upper triangle of L
                {L[i][j]=0;}
            
            if (i==j) // Conditional for contrcuting the main diagonal of L
                {L[i][j]=1;}
            sumU=0;
            for(int k=0; k<=i-1; k++){ 
                sumU=sumU+L[i][k]*U[k][j];}
            U[i][j] = *(M + (j+i*m*n))-sumU;
            }
            
        if (i>j){ // Conditional for constructing U
            U[i][j]=0;
            sumL=0;
            for(int k=0; k<=j-1; k++){
                sumL=sumL+L[i][k]*U[k][j];
                }
            L[i][j]=(*(M + (j+i*m*n))-sumL)/U[j][j];
            } 
    }

    for(int i=0; i<=m*n-1;i++){
        sumz=0;
        for(int j=0; j<=i-1;j++){ 
            sumz=sumz+L[i][j]*z[j];
            }
        z[i]=b[i]-sumz;
    }
    for(int i=m*n-1; i>=0; i--){
        double cx=0, sumx=0;
        for (int j=i+1; j<=m*n-1; j++){ //Loop utilizado para hallar el sumatoriot
            sumx=sumx+U[i][j]*x[j];
        }
        x[i]=(z[i]-sumx)/U[i][i];
        }

    return x;
    }


int main()

{   
    // geometrical limits of the problem
    double x_left = 0;
    double x_right = 1;

    double y_bottom = 0;
    double y_top = 1;

    double h = (x_right - x_left)/(m-1); // space for the x-axis discretization
    double k = (y_top - y_bottom)/(n-1); // space for the y-axis discretization
    double t = -2*(1/(h*h) + 1/(k*k));

    // vectors of X and Y discretized positions
    double X[m], Y[n];
    for (int i=0; i < m; i++){
        X[i] = x_left + i*h;
    } 

    for (int i=0; i < n; i++){
        Y[i] = y_bottom + i*h;
    } 

    // main matrix calculation
    double A[m*n][m*n];

    // initializing it as a diagonal matrix
    for (int i=0;i<m*n; i++){
        for (int j=0;j<m*n; j++){
            if (i==j){
                A[i][j] = 1;
            }
            else {
                A[i][j] = 0;
            }
        }
    }

    // filling the non-boundary elements
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            int p = i + (j)*m;
            A[p][p] = t;
            A[p][i + 1 + (j)*m] = 1/(h*h);
            A[p][i - 1 + (j)*m] = 1/(h*h);
            A[p][i + (j+1)*m] = 1/(k*k);
            A[p][i + (j-1)*m] = 1/(k*k);            
        }
    } 

    double b[m*n];
    // creo que el [1] sobra
    // first and last 
    for (int i=0;i<m;i++){
        b[i] = gx_yBottom(X[i]);
        b[(n-1)*m +i] = gx_yTop(X[i]);
    }
    // intermediate sub-arrays
    for (int j=1;j<n-1;j++){
        for(int i=0;i<m;i++){
            b[j*m + i] = f(X[i],Y[j]);
        }
    }

    for (int j=1;j<n-1;j++){
        b[j*m] = gy_xLeft(Y[j]);
        b[(j+1)*m -1] = gy_xRight(Y[j]);
    }

    // calling LU solver
    double* sol_pointer = LU_solver(*A,b);
    
    // saving the data for a later plot with python
    ofstream out_file ("solution_data.txt");
    
    out_file << "X values" << endl;
    for (int i=0;i<m;i++){
        out_file << X[i] << endl;
    }

    out_file << endl;
    out_file << "Y values" << endl;
    for (int i=0;i<n;i++){
        out_file << Y[i] << endl;
    }

    out_file << endl;
    out_file << "U values" << endl;
    for (int i=0;i<m*n;i++){
        out_file << *(sol_pointer + i) << endl;
    }
    out_file << endl;

}