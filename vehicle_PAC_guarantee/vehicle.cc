/*
 * vehicle.cc
 *
 *  created on: 10.02.2025
 *      author: Daniel A.
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound_for_PAC.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif


/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* state space interval*/
const double hn1 = 0.18;
const double hn2 = 0.18;
/* sampling time */
const double tau = 0.25;
/* number of samples */
double const N = 700;
/* lower bounds of the hyper rectangle */
double lbn[sDIM]={0,0,-M_PI-hn2*2};  
/* upper bounds of the hyper rectangle */
double ubn[sDIM]={10,10,M_PI+hn2*2}; 

/* we discretize the continuous time vehicle ode by forward Euler (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, input_type &u) -> void {

  /* the difference inclusion describing the vehicle */
  auto rhs =[](state_type& dx,  const state_type &x, input_type &u) {
      double alpha=std::atan(std::tan(u[1])/2.0);
      dx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      dx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      dx[2] = u[0]*std::tan(u[1]);
  };
    state_type dx;
    rhs(dx, x, u);
    for (int j = 0; j < sDIM; ++j) {
        x[j] += tau * dx[j]; 
    }
};

/* dynamics for computing N iid successors */
state_type dyn(state_type x, input_type u){

    state_type dx;
    state_type xn;
    double alpha=std::atan(std::tan(u[1])/2.0);
    dx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
    dx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
    dx[2] = u[0]*std::tan(u[1]);
    for (int j = 0; j < sDIM; ++j) {
        xn[j] = x[j] + tau * dx[j];
    }
    return xn;
}

/* for over-approximating one-time step reachable set (the result is stored in r)  */
/*in the case of using subgrid as it is in the paper*/
auto radius_post = [](state_type &r, state_type &x, input_type &u) {
    double eta_x[sDIM] = {hn1*0.5,hn1*0.5,hn2*0.5};

    std::vector<state_type> successors;

    // Determine grid points per dimension
    int grid_pts_per_dim = static_cast<int>(std::round(std::pow(N, 1.0 / sDIM)));
    // int total_samples = std::pow(grid_pts_per_dim, sDIM);

    // Build subcell-centered grid per dimension
    std::vector<std::vector<double>> grid(sDIM);
    for (int i = 0; i < sDIM; ++i) {
        double step = (2 * eta_x[i]) / grid_pts_per_dim;
        for (int j = 0; j < grid_pts_per_dim; ++j) {
            double center = x[i] - eta_x[i] + (j + 0.5) * step;
            grid[i].push_back(center);
        }
    }

    // Cartesian product to generate center points of subcells
    std::function<void(int, state_type&)> recurse;
    recurse = [&](int dim, state_type& current_sample) {
        if (dim == sDIM) {
            successors.push_back(dyn(current_sample, u));
            return;
        }
        for (double val : grid[dim]) {
            current_sample[dim] = val;
            recurse(dim + 1, current_sample);
        }
    };

    state_type x_sample;
    recurse(0, x_sample);

    // Truncate to N if needed
    if (successors.size() > N) {
        successors.resize(N);
    }

    // Compute mean
    state_type mean;
    mean.fill(0.0);
    for (const auto &s : successors) {
        for (int i = 0; i < sDIM; ++i) {
            mean[i] += s[i];
        }
    }
    for (int i = 0; i < sDIM; ++i) {
        mean[i] /= successors.size();
    }

    // Compute max deviation
    state_type max_dev;
    max_dev.fill(0.0);
    for (const auto &s : successors) {
        for (int i = 0; i < sDIM; ++i) {
            max_dev[i] = std::max(max_dev[i], std::abs(s[i] - mean[i]));
        }
    }

    for (int i = 0; i < sDIM; ++i) {
        r[i] = max_dev[i];
    }
    
    // // check out the radius 
    // for(int i=0; i<sDIM; i++){
    //   std::cout << r[i] << "\n";
    // } 
};

/* in the case of using N iid samples within a cell */
// auto radius_post = [](state_type &r, state_type &x, input_type &u) {

//     double eta_x[sDIM] = {hn1*0.5,hn1*0.5,hn2*0.5};
//     std::vector<state_type> successors(N);
//     std::random_device rd;
//     std::mt19937 gen(rd());

//     // Generate N uniform samples around x
//     for (int j = 0; j < N; ++j) {
//         state_type x_sample;
        
//         // Generate a full state sample
//         for (int i = 0; i < sDIM; ++i) {
//             std::uniform_real_distribution<double> dist(x[i] - eta_x[i], x[i] + eta_x[i]);
//             x_sample[i] = dist(gen);
//         }

//         // Apply dynamics and store the successor state
//         successors[j] = dyn(x_sample, u);
//     }

//     // Compute mean successor
//     state_type mean;
//     mean.fill(0.0);
//     for (const auto &s : successors) {
//         for (int i = 0; i < sDIM; ++i) {
//             mean[i] += s[i];
//         }
//     }
//     for (int i = 0; i < sDIM; ++i) {
//         mean[i] /= N;
//     }

//     // Compute max deviation
//     state_type max_dev;
//     max_dev.fill(0.0);
//     for (const auto &s : successors) {
//         for (int i = 0; i < sDIM; ++i) {
//             max_dev[i] = std::max(max_dev[i], std::abs(s[i] - mean[i]));
//         }
//     }

//     for (int i=0; i<sDIM; i++) {
//         r[i] = max_dev[i];
//     }

//     // // check out the radius 
//     // for(int i=0; i<3; i++){
//     //   std::cout << r[i] << "\n";
//     // } 
// };



/* forward declaration of the functions to setup the state space
 * input space and obstacles of the vehicle example */
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr);

void vehicleCreateObstacles(scots::SymbolicSet &obs);


int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  scots::SymbolicSet ss=vehicleCreateStateSpace(mgr);
  ss.writeToFile("vehicle_ss.bdd");

  /****************************************************************************/
  /* construct SymbolicSet for the obstacles */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet obs(ss);
  vehicleCreateObstacles(obs);
  obs.writeToFile("vehicle_obst.bdd");

  /****************************************************************************/
  /* we define the target set */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet ts(ss);
  /* define the target set as a symbolic set */
  double H[4*sDIM]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
  /* compute inner approximation of P={ x | H x<= h1 }  */
  double h[4] = {-9,9.51,-0, .51};
  ts.addPolytope(4,H,h, scots::INNER);
  ts.writeToFile("vehicle_target.bdd");

  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  scots::SymbolicSet is=vehicleCreateInputSpace(mgr);

  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  /* first create SymbolicSet of post variables 
   * by copying the SymbolicSet of the state space and assigning new BDD IDs */
  scots::SymbolicSet sspost(ss,1);
  /* instantiate the SymbolicModel */
  scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(vehicle_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

  /****************************************************************************/
  /* we continue with the controller synthesis */
  /****************************************************************************/
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T = ts.getSymbolicSet();
  BDD O = obs.getSymbolicSet();
  tt.tic();
  /* compute controller */
  BDD C=fp.reachAvoid(T,O,1);
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller(ss,is);
  controller.setSymbolicSet(C);
  std::cout << "Controller size: " << controller.getSize() << std::endl;
  controller.writeToFile("vehicle_controller.bdd");

  return 1;
}

scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr) {

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  double eta[sDIM]={hn1,hn1,hn2};   


  scots::SymbolicSet ss(mgr,sDIM,lbn,ubn,eta);

  /* add the grid points to the SymbolicSet ss */
  ss.addGridPoints();

 return ss;
}

void vehicleCreateObstacles(scots::SymbolicSet &obs) {

  /* add the obstacles to the symbolic set */
  /* the obstacles are defined as polytopes */
  /* define H* x <= h */
  double H[4*sDIM]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
  /* add outer approximation of P={ x | H x<= h1 } form state space */
  double h1[4] = {-1,1.2,-0, 9};
  obs.addPolytope(4,H,h1, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h2 } form state space */
  double h2[4] = {-2.2,2.4,-0,5};
  obs.addPolytope(4,H,h2, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h3 } form state space */
  double h3[4] = {-2.2,2.4,-6,10};
  obs.addPolytope(4,H,h3, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h4 } form state space */
  double h4[4] = {-3.4,3.6,-0,9};
  obs.addPolytope(4,H,h4, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h5 } form state space */
  double h5[4] = {-4.6 ,4.8,-1,10};
  obs.addPolytope(4,H,h5, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h6 } form state space */
  double h6[4] = {-5.8,6,-0,6};
  obs.addPolytope(4,H,h6, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h7 } form state space */
  double h7[4] = {-5.8,6,-7,10};
  obs.addPolytope(4,H,h7, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h8 } form state space */
  double h8[4] = {-7,7.2,-1,10};
  obs.addPolytope(4,H,h8, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h9 } form state space */
  double h9[4] = {-8.2,8.4,-0,8.5};
  obs.addPolytope(4,H,h9, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h10 } form state space */
  double h10[4] = {-8.4,9.3,-8.3,8.5};
  obs.addPolytope(4,H,h10, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h11 } form state space */
  double h11[4] = {-9.3,10,-7.1,7.3};
  obs.addPolytope(4,H,h11, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h12 } form state space */
  double h12[4] = {-8.4,9.3,-5.9,6.1};
  obs.addPolytope(4,H,h12, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h13 } form state space */
  double h13[4] = {-9.3,10,-4.7,4.9};
  obs.addPolytope(4,H,h13, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h14 } form state space */
  double h14[4] = {-8.4,9.3,-3.5,3.7};
  obs.addPolytope(4,H,h14, scots::OUTER);
  /* add outer approximation of P={ x | H x<= h15 } form state space */
  double h15[4] = {-9.3,10 ,-2.3,2.5};
  obs.addPolytope(4,H,h15, scots::OUTER);

}

scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={-1,-1};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={1,1}; 
  /* grid node distance diameter */
  double eta[sDIM]={.3,.3};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}


