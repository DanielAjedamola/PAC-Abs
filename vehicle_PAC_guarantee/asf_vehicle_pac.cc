/*
* Computing an ASF for the vehicle model abstraction with PAC guarantee
* 
* created on: Feb 2025
*     author: Daniel A.
*/

#include <array>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "gurobi_c++.h"
#include "scots.hh"
#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicSets.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "RungeKutta4 _distb.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

// set up state-space initial parameters
/* lower bounds of the hyper rectangle */
double lb[sDIM]={0,0,-M_PI-0.4};  
/* upper bounds of the hyper rectangle */
double ub[sDIM]={10,10,M_PI+0.4}; 
// state grid parameter
double eta_x[sDIM]={0.1,0.1,0.2};

/* sampling time */
const double tau = 0.25;

/* vehicle dynamics as black-box simulator */
state_type vehicleDyn(state_type x, input_type u){

    state_type dx;
    state_type xn;
    double alpha=std::atan(std::tan(u[1])/2.0);
    dx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
    dx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
    dx[2] = u[0]*std::tan(u[1]);
    for (int j = 0; j < sDIM; ++j) {
        xn[j] = x[j] + tau * dx[j]; // forward Euler discretization
    }
    return xn;
}

// Define a struct to hold the quantized simulator for abstract state,input
struct AbstractVec {
    std::vector<state_type> X_hat;
    std::vector<input_type> U_hat;
};

state_type quant(const state_type& x, std::vector<state_type> X_hatt) {
    // Ensure x is within bounds
    for (std::size_t i = 0; i < sDIM; ++i) {
        assert(x[i] >= lb[i] && x[i] <= ub[i]);
    }

    // Search for a close enough quantized value
    for (const auto& x_hat : X_hatt) {
        bool within_threshold = true;
        for (std::size_t i = 0; i < sDIM; ++i) {
            if (std::round(std::abs(x_hat[i] - x[i]) * 10000.0) / 10000.0 > eta_x[i]) {
                within_threshold = false;
                break;
            }
        }
        if (within_threshold) {
            return x_hat;
        }
    }

    throw std::runtime_error("No quantized state found"); 
}

/* forward declaration of the functions to setup the state space and input space of the vehicle */
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr);

// converting the grid encoded as BDDs to a vectors of points 
std::vector<double> bdd_to_grid_point(Cudd& manager, scots::SymbolicSet& symbolicSet) {
    BDD bdd = symbolicSet.getSymbolicSet();

    size_t dim = symbolicSet.getDimension();
    const double* eta = symbolicSet.getEta();
    const double* firstGridPoint = symbolicSet.getFirstGridPoint();

    // Find the number of grid points in the BDD
    size_t no_gp = symbolicSet.getSize();

    // Initialize the vector of grid points to be returned
    std::vector<double> gp(no_gp * dim);

    // Set up iteration to iterate over BDD cubes
    DdManager* dd = manager.getManager();
    int* cube;
    CUDD_VALUE_TYPE value;
    DdGen* gen;
    size_t counter = 0;

    // Iterate over BDD cubes
    Cudd_ForeachCube(dd, bdd.getNode(), gen, cube, value) {
        size_t offset = 1;
        for (size_t i = 0; i < dim; i++) {
            size_t no_vars = symbolicSet.getNofBddVars()[i];
            size_t* indBddVars = symbolicSet.getIndBddVars()[i];
            for (size_t j = 0; j < no_vars; j++) {
                size_t id = indBddVars[j];
                if (cube[id] == 1) {
                    for (size_t k = 0; k < offset; k++) {
                        gp[(counter + k) * dim + i] = firstGridPoint[i] + (size_t{1} << (no_vars - 1 - j)) * eta[i];
                    }
                }
                /* take care of don't care */
                if (cube[id] == 2) {
                    for (size_t k = 0; k < offset; k++) {
                        for (size_t l = 0; l <= i; l++) {
                            gp[(counter + k + offset) * dim + l] = gp[(counter + k) * dim + l];
                        }
                    }
                    for (size_t k = 0; k < offset; k++) {
                        gp[(counter + k + offset) * dim + i] = firstGridPoint[i] + (size_t{1} << (no_vars - 1 - j)) * eta[i];
                    }
                    offset = (offset << 1);
                }
            }
        }
        counter += offset;
    }
    // The vector 'gp' should contain the grid points, check if it's empty
    if (gp.empty()) {
        std::cout << "Error: Grid points vector is empty!" << std::endl;
    }
    return gp;
}

int main() {
    /* to measure time */
    TicToc tt;
    /* there is one unique manager to organize the bdd variables */
    Cudd mgr;

    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    tt.tic();
    scots::SymbolicSet ss = vehicleCreateStateSpace(mgr);
    std::vector<double> ss_vector = bdd_to_grid_point(mgr, ss);
    
    /****************************************************************************/
    /* construct SymbolicSet for the input space */
    /****************************************************************************/
    scots::SymbolicSet is = vehicleCreateInputSpace(mgr);
    std::vector<double> is_vector = bdd_to_grid_point(mgr, is);

    std::vector<state_type> X_hat;
    std::vector<input_type> U_hat;

    
    // Nested loop to iterate over each pair (s, i) in ss_vector x is_vector
    for (size_t s = 0; s < ss_vector.size(); s += 3) {
        for (size_t i = 0; i < is_vector.size(); i += 2) {
            // Extract the current (x, u) pair from ss_vector and is_vector
            // std::vector<double> x(ss_vector.begin() + s, ss_vector.begin() + s + 3);
            // std::vector<double> u(is_vector.begin() + i, is_vector.begin() + i + 2);
            // state_type xx;
            // input_type uu;

            // std::copy(ss_vector.begin() + s, ss_vector.begin() + s + 3, xx.begin());
            // std::copy(is_vector.begin() + i, is_vector.begin() + i + 2, uu.begin());
            X_hat.push_back({ss_vector[s], ss_vector[s + 1], ss_vector[s + 2]});
            U_hat.push_back({is_vector[i], is_vector[i + 1]});
            
            // // to print the values of pair (r,u)
            // std::cout << "r: ";
            // for (const auto& elem : r) {
            //     std::cout << elem << " ";
            // }
            // std::cout << std::endl;

            // // Print the elements of u
            // std::cout << "u: ";
            // for (const auto& elem : u) {
            //     std::cout << elem << " ";
            // }
            // std::cout << std::endl;
        }
    }

    size_t Nx = X_hat.size();
    size_t M = U_hat.size();
    std::vector<std::vector<state_type>> X_U_hat(Nx, std::vector<state_type>(M));
    for (size_t j = 0; j < Nx; j++) {
        for (size_t i = 0; i < M; i++) {
            // Calculate successor of (x, u) pair
            state_type xn = vehicleDyn(X_hat[j], U_hat[i]);
            state_type xn_hat = quant(xn, X_hat);
            // Store the calculated successor in the vectors
            X_U_hat[j][i] = xn_hat;
        }
    }

    // // Print M1_vector
    // std::cout << "M1_vector:" << std::endl;
    // for (const auto& matrix : M1_vector) {
    //     for (const auto& row : matrix) {
    //         for (double value : row) {
    //             std::cout << value << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // // Print M2_vector
    // std::cout << "M2_vector:" << std::endl;
    // for (const auto& row : M2_vector) {
    //     for (double value : row) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::cout << "Abstract x-u-x' concluded" << std::endl;
    std::cout << "Now optimising" << std::endl;
    tt.toc();
    
    // try {
    //     // Set up Gurobi environment
    //     GRBEnv env = GRBEnv(true);
    //     env.set(GRB_IntParam_OutputFlag, 0);
    //     env.set(GRB_IntParam_LogToConsole, 0);
    //     // Start the environment
    //     env.start();

    //     // optimising the grid parameter for the least abstraction size.
    //     size_t J = M1_vector.size(); // Set the size of M1 and M2 here
    //     GRBModel model = GRBModel(env);
    //     // Decision variables bound
    //     double x_bd = 20.5;//0.50; 
    //     double n = -1;
    //     double y_bd = std::exp(x_bd);
    //     double z_bd = std::exp(2.0*x_bd);
    //     double s_bd = std::exp(3.0*x_bd);

    //     GRBVar x1 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x1");
    //     GRBVar x2 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x2");
    //     GRBVar x3 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x3");
    //     //auxilliary variable y=exp(x)
    //     GRBVar y1 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y1");
    //     GRBVar y2 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y2");
    //     GRBVar y3 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y3");

    //     std::vector<GRBVar> z(6);
    //     for (int i = 0; i < 6; i++) {
    //         z[i] = model.addVar(0, z_bd, 0, GRB_CONTINUOUS, "z" + std::to_string(i + 1));
    //     }

    //     std::vector<GRBVar> s(10);
    //     for (int i = 0; i < 10; i++) {
    //         s[i] = model.addVar(0, s_bd, 0, GRB_CONTINUOUS, "s" + std::to_string(i + 1));
    //     }
    //     // // Create an auxiliary variable for constants * x_i
    //     GRBVar x11 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x11");
    //     GRBVar x12 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x12");
    //     GRBVar x13 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x13");
    //     GRBVar x22 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x22");
    //     GRBVar x23 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x23");
    //     GRBVar x33 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x33");
    //     GRBVar x111 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x111");
    //     GRBVar x112 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x112");
    //     GRBVar x113 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x113");
    //     GRBVar x122 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x122");
    //     GRBVar x123 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x123");
    //     GRBVar x222 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x222");
    //     GRBVar x133 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x133");
    //     GRBVar x233 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x233");
    //     GRBVar x223 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x223");
    //     GRBVar x333 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x333");

    //     // Update the model
    //     model.update();

    //     // Objective function
    //     GRBQuadExpr objective = 0.0;
    //     // Calculate the exponential term
    //     double exp_term = std::exp(-1 * gamma);
    //     for (size_t j = 0; j < J; j++) {
    //         double m11 = M1_vector[j][0][0];
    //         double m12 = M1_vector[j][0][1];
    //         double m13 = M1_vector[j][0][2];
    //         double m21 = M1_vector[j][1][0];
    //         double m22 = M1_vector[j][1][1];
    //         double m23 = M1_vector[j][1][2];
    //         double m31 = M1_vector[j][2][0];
    //         double m32 = M1_vector[j][2][1];
    //         double m33 = M1_vector[j][2][2];

    //         double p1 = M2_vector[j][0];
    //         double p2 = M2_vector[j][1];
    //         double p3 = M2_vector[j][2];

    //         // Add the quadratic terms to the objective expression
    //         objective += p1 * p2 * p3 +
    //                     p1 * p2 * (m31 * y1) + p1 * p2 * (m32 * y2) + p1 * p2 * (m33 * y3) +
    //                     p1 * (m21 * y1) * p3 + p1 * m21 * (m31 * z[0]) + p1 * m21 * (m32 * z[1]) + p1 * m21 * (m33 * z[2]) +
    //                     p1 * (m22 * y2) * p3 + p1 * (m22 * z[1]) * m31 + p1 * m22 * (m32 * z[3]) + p1 * m22 * (m33 * z[4]) +
    //                     p1 * (m23 * y3) * p3 + p1 * (m23 * z[2]) * m31 + p1 * m23 * (m32 * z[4]) + p1 * m23 * (m33 * z[5]) +
    //                     (m11 * y1) * p2 * p3 + m11 * p2 * (m31 * z[0]) + m11 * p2 * (m32 * z[1]) + m11 * p2 * (m33 * z[2]) +
    //                     m11 * (m21 * z[0]) * p3 + m11 * m21 * (m31 * s[0]) + m11 * m21 * (m32 * s[1]) + m11 * m21 * (m33 * s[2]) +
    //                     m11 * (m22 * z[1]) * p3 + m11 * (m22 * s[1]) * m31 + m11 * m22 * (m32 * s[3]) + m11 * m22 * (m33 * s[4]) +
    //                     m11 * (m23 * z[2]) * p3 + m11 * (m23 * s[2]) * m31 + m11 * (m23 * s[4]) * m32 + (m11 * s[5]) * m23 * m33 +
    //                     (m12 * y2) * p2 * p3 + (m12 * z[1]) * p2 * m31 + m12 * p2 * (m32 * z[3]) + m12 * p2 * (m33 * z[4]) +
    //                     (m12 * z[1]) * m21 * p3 + (m12 * s[1]) * m21 * m31 + m12 * (m21 * s[3]) * m32 + m12 * (m21 * s[4]) * m33 +
    //                     m12 * (m22 * z[3]) * p3 + m12 * (m22 * s[3]) * m31 + m12 * (m22 * s[6]) * m32 + m12 * m22 * (m33 * s[7]) +
    //                     m12 * (m23 * z[4]) * p3 + (m12 * s[4]) * m23 * m31 + m12 * m23 * (m32 * s[7]) + m12 * (m23 * s[8]) * m33 +
    //                     (m13 * y3) * p2 * p3 + (m13 * z[2]) * p2 * m31 + (m13 * z[4]) * p2 * m32 + m13 * p2 * (m33 * z[5]) +
    //                     (m13 * z[2]) * m21 * p3 + (m13 * s[2]) * m21 * m31 + (m13 * s[4]) * m21 * m32 + m13 * (m21 * s[5]) * m33 +
    //                     m13 * (m22 * z[4]) * p3 + m13 * (m22 * s[4]) * m31 + (m13 * s[7]) * m22 * m32 + (m13 * s[8]) * m22 * m33 +
    //                     m13 * (m23 * z[5]) * p3 + m13 * (m23 * s[5]) * m31 + m13 * (m23 * s[8]) * m32 + m13 * m23 * (m33 * s[9]);

    //     }
    //     objective *= exp_term;
    //     model.setObjective(objective, GRB_MINIMIZE);

    //     // Constraint: x1 + x2 + x3 = gamma
    //     GRBLinExpr constraint123 = x1 + x2 + x3;
    //     model.addConstr(constraint123 == gamma);
    //     // constraints due to the auxilliary var
    //     model.addConstr(x11 == 2.0 * x1); // Add the constraint for x1_double = 2 * x1
    //     model.addConstr(x22 == 2.0 * x2);
    //     model.addConstr(x33 == 2.0 * x3);
    //     model.addConstr(x111 == 3.0 * x1);
    //     model.addConstr(x12 == x1 + x2);
    //     model.addConstr(x13 == x1 + x3);
    //     model.addConstr(x23 == x2 + x3);
    //     model.addConstr(x112 == x1 + x1 + x2);
    //     model.addConstr(x113 == x1 + x1 + x3);
    //     model.addConstr(x122 == x1 + x2 + x2);
    //     model.addConstr(x123 == x1 + x2 + x3);
    //     model.addConstr(x222 == x2 + x2 + x2);
    //     model.addConstr(x133 == x1 + x3 + x3);
    //     model.addConstr(x223 == x2 + x2 + x3);
    //     model.addConstr(x233 == x2 + x3 + x3);
    //     model.addConstr(x333 == x3 + x3 + x3);
    //     model.addGenConstrExp(x1, y1);
    //     model.addGenConstrExp(x2, y2);
    //     model.addGenConstrExp(x3, y3);
    //     model.addGenConstrExp(x11, z[0]);
    //     model.addGenConstrExp(x12, z[1]);
    //     model.addGenConstrExp(x13, z[2]);
    //     model.addGenConstrExp(x22, z[3]);
    //     model.addGenConstrExp(x23, z[4]);
    //     model.addGenConstrExp(x33, z[5]);
    //     model.addGenConstrExp(x111, s[0]);
    //     model.addGenConstrExp(x112, s[1]);
    //     model.addGenConstrExp(x113, s[2]);
    //     model.addGenConstrExp(x122, s[3]);
    //     model.addGenConstrExp(x123, s[4]);
    //     model.addGenConstrExp(x133, s[5]);
    //     model.addGenConstrExp(x222, s[6]);
    //     model.addGenConstrExp(x223, s[7]);
    //     model.addGenConstrExp(x233, s[8]);
    //     model.addGenConstrExp(x333, s[9]);
        
    //     tt.tic();
    //     // Optimize the model
    //     model.optimize();

    //     tt.toc();
    //     // Check the optimization status and retrieve the optimized values
    //     if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
    //         double opt_x1 = x1.get(GRB_DoubleAttr_X);
    //         double opt_x2 = x2.get(GRB_DoubleAttr_X);
    //         double opt_x3 = x3.get(GRB_DoubleAttr_X);
    //         double eta_x_1 = y1.get(GRB_DoubleAttr_X);
    //         double eta_x_2 = y2.get(GRB_DoubleAttr_X);
    //         double eta_x_3 = y3.get(GRB_DoubleAttr_X);

    //         std::cout << "Optimal values: x1 = " << opt_x1 << ", x2 = " << opt_x2 << ", x3 = " << opt_x3 << std::endl;
    //         std::cout << "Optimal values: eta_x_1 = " << eta_x_1 << ", eta_x_2 = " << eta_x_2 << ", eta_x_3 = " << eta_x_3 << std::endl;
    //         std::cout << "for eta_x" << eta_x[0] << "," << eta_x[1] << "," << eta_x[2] << std::endl;
    //     } else {
    //         std::cout << "Optimization failed." << std::endl;
    //     }
    // } catch (GRBException e) {
    //     std::cout << "Error code = " << e.getErrorCode() << std::endl;
    //     std::cout << e.getMessage() << std::endl;
    // } catch (...) {
    //     std::cout << "Exception during optimization" << std::endl;
    // }

    return 1;
}


scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr) {

    /* setup the workspace of the uniform grid */   
    // using pre-defined parameters
    scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta_x);

    /* add the grid points to the SymbolicSet ss */
    ss.addGridPoints();

 return ss;
}

scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[iDIM]={-1,-1};  
  /* upper bounds of the hyper rectangle */
  double ub[iDIM]={1,1}; 
  /* grid node distance diameter */
  double eta[iDIM]={.3,.3};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}