#include <Eigen/Sparse>
#include <iostream>
#include "writer.hpp"
#include <cmath>
#include <Eigen/SparseCholesky>
#include <stdexcept>

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;

//! Vector type
typedef Eigen::VectorXd Vector;

//! Our function pointer, typedef'd to make it easier to use
typedef double(*FunctionPointer)(double, double);

//----------------poissonBegin----------------

//! Create the Poisson matrix for 2D finite difference.
//! @param[out] A will be the Poisson matrix (as in the exercise)
//! @param[in] N number of elements in the x-direction
//! @param[in] dx the cell width

void createPoissonMatrix2D(SparseMatrix& A, int N, double dx)
{
     // Fill the matrix A using setFromTriplets - method

// (write your solution here)

    /// Start of my solution ///

    // determine the size of the matrix A
    A.resize(N*N, N*N);

    // vector of triplet to fill the matrix A
    std::vector<Triplet> triplets;

    // number of triplets needed: B-matrices:  N * (N + 2 * N - 2)
    //                            I-matrices:  (2 * N*N - 2) * N
    triplets.reserve(N * (N + 2 * N - 2) + (2 * N - 2) * N);


    for (int i = 0; i < N*N; i++)
    {
        // prepare the triplets to fill the matrix A => Triplet(row, column, value)

        // 4 on the diagonal
        triplets.push_back(Triplet(i, i, 4 + dx*dx));

        // -1 of B-matrices
        if (i%N != 0   &&   (i-N+1)%N != 0)
        {
            triplets.push_back(Triplet(i, i-1, -1));
            triplets.push_back(Triplet(i, i+1, -1));
        }
        else if (i%N == 0)
        {
            triplets.push_back(Triplet(i, i+1, -1));
        }
        else if ((i-N+1)%N == 0)
        {
            triplets.push_back(Triplet(i, i-1, -1));
        }

        // -1 of I-matrices
        if (i > (N-1))
        {
            triplets.push_back(Triplet(i-N, i, -1));
            triplets.push_back(Triplet(i, i-N, -1));
        }
    }

    // fill the matrix with the triplets
    A.setFromTriplets(triplets.begin(), triplets.end());

    /// End of my solution ///

}

//----------------poissonEnd----------------


//----------------RHSBegin----------------

//! Create the Right hand side for the 2D finite difference
//! @param[out] rhs will at the end contain the right hand side
//! @param[in] f the right hand side function f
//! @param[in] N the number of points in the x direction
//! @param[in] dx the cell width
//! @param[in] g the boundary condition function g

void createRHS(Vector& rhs, FunctionPointer f, int N, double dx, FunctionPointer g)
{
    rhs.resize(N * N);

    // fill up RHS
    // remember that the index (i,j) corresponds to j*N+i

// (write your solution here)

    /// Start of my solution ///

    for (int j = 0; j < N; j++)
    {
        // y_j  with  j = 1, ..., N
        double y_j = (j + 1) * dx;

        for (int i = 0; i < N; i++)
        {
            // x_i  with  i = 1, ..., N
            double x_i = (i + 1) * dx;

            // f_(i,j) = dx^2 * f(x_i, y_i)
            rhs[j*N + i] = dx * dx * f(x_i, y_j);

            // add the function g on the boundaries
            if (i == 0) // if i = 1
            {
                rhs[j*N + i] += g(0, y_j);
            }
            if (j == 0) // if j = 1
            {
                rhs[j*N + i] += g(x_i, 0);
            }
            if (i == N-1) // if i = N
            {
                rhs[j*N + i] += g(N+1, y_j);
            }
            if (j == N-1) // if j = N
            {
                rhs[j*N + i] += g(x_i, N+1);
            }
        }
    }

    /// End of my solution ///

}

//----------------RHSEnd----------------


//----------------solveBegin----------------

//! Solve the Poisson equation in 2D
//! @param[out] u will contain the solution u
//! @param[in] f the function pointer to f
//! @param[in] N the number of points to use (in x direction)

void poissonSolve(Vector& u, FunctionPointer f, int N, FunctionPointer g)
{
    // Solve Poisson 2D here
// (write your solution here)

    /// Start of my solution ///

    // step size
    double dx = 1.0 / (N + 1);

    // matrix A
    SparseMatrix A;
    createPoissonMatrix2D(A, N, dx);

    // vector F
    Vector rhs;
    createRHS(rhs, f, N, dx, g);

    // matrix to solve AU=F
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(A);

    // print an error message if it doesn't work
    if ( solver.info() !=  Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }

    // resize the vector u to include the points at the boundary
    u.resize((N+2) * (N+2));
    u.setZero();

    Vector innerU = solver.solve(rhs);

    // copy vector to inner u
    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            u[i * (N + 2) + j] = innerU[(i - 1) * N + j - 1];
        }
    }

    /// End of my solution ///

}

//----------------solveEnd----------------


double F(double x, double y) {
     return (1 + 8*M_PI*M_PI)*sin(2*M_PI*x)*cos(2*M_PI*y);
}

double G(double x, double y) {
     return sin(2*M_PI*x);
}



int main(int, char**) {
    Vector u;
    poissonSolve(u, F, 100, G);
    writeToFile("u_fd.txt", u);
}
