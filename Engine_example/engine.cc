/*
 * engine.cc
 *
 *  created on: Nov.05.2025
 *      author: Daniel A.
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <iostream>
#include <cmath>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound_for_PAC.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

/* state space dim */
#define sDIM 2
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,2> state_type;
typedef std::array<double,2> input_type;

/* sampling time */
const double tau = 0.1;
/* state space discretization parameter*/
const double hn1 = 16e-6;
/* input space discretization parameter*/
const double hu1 = 1.0/30;
const double hu2 = 0.1;
/* number of samples */
const double N = 100; 
/* dynamics parameter */
const double a = 1.0/3.5;
const double H = 0.18;
const double l = 8;
const double B = 2;
const double W = 0.25;

/* we discretize the continuous time engine ode by forward Euler (the result is stored in x)  */
auto  engine_post = [](state_type &x, input_type &u) -> void {

  /* the difference equation describing the engine */
  auto rhs =[](state_type& dx,  const state_type &x, input_type &u) {
    double psi = a + H*(1 + 1.5*(x[0]/W-1) - 0.5*pow(x[0]/W-1, 3));
    dx[0] = 1/l*(psi - x[1]) + u[0];
    dx[1] = 1/(4*l*B*B)*(x[0] - u[1]*sqrt(x[1]));
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
    double psi = a + H*(1 + 1.5*(x[0]/W-1) - 0.5*pow(x[0]/W-1, 3));
    dx[0] = 1/l*(psi - x[1]) + u[0];
    dx[1] = 1/(4*l*B*B)*(x[0] - u[1]*sqrt(x[1]));

    for (int j = 0; j < sDIM; ++j) {
        xn[j] = x[j] + tau * dx[j];
    }
    return xn;
}

/* for over-approximating one-time step reachable set (the result is stored in r)  */
/*in the case of using subgrid as it is in the paper*/
auto radius_post = [](state_type &r, state_type &x, input_type &u) {
    double eta_x[sDIM] = {hn1*0.5,hn1*0.5};

    std::vector<state_type> successors;

    // Determine grid points per dimension
    int grid_pts_per_dim = static_cast<int>(std::round(std::pow(N, 1.0 / sDIM)));

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


int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0.4479,0.6473};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={0.4559,0.6553}; 
  /* grid node distance diameter */
  double eta[sDIM]={hn1,hn1};   
  scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);
  ss.addGridPoints();

  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  double ilb[iDIM]={-0.05, 0.5};  
  double iub[iDIM]={0.05, 0.8}; 
  double ieta[iDIM]={hu1, hu2};   
  scots::SymbolicSet is(mgr,iDIM,ilb,iub,ieta);
  is.addGridPoints();

  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  /* first create SymbolicSet of post variables 
   * by copying the SymbolicSet of the state space and assigning new BDD IDs */
  scots::SymbolicSet sss(ss,1);
  /* instantiate the SymbolicModel */
  scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sss);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(engine_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

  /****************************************************************************/
  /* we continue with the specification */
  /****************************************************************************/
  /* construct SymbolicSet for the safe set */
  double eps = hn1; /*approximate value of the ASF, obtained from the .ipynb file*/
  double H[4*sDIM]={-1, 0,
                     1, 0,
                     0,-1,
                     0, 1};
  /* add outer approximation of P={ x | H x<= h } form state space */
  double h[4] = {-0.4479+eps, 0.4559-eps, -0.6473+eps, 0.6553-eps};
  /* initialize the safe set with the ss 
   * in order to obtain all the necessary information */
  scots::SymbolicSet safe(ss);
  safe.addPolytope(4,H,h, scots::INNER);
  std::cout << "Safe set details:" << std::endl;
  safe.writeToFile("engine.bdd");

  /****************************************************************************/
  /* we continue with the controller synthesis */
  /****************************************************************************/
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD Z = safe.getSymbolicSet();
  /* contains the state space grid points and associated inputs */
  BDD C; 
  tt.tic();
  C = fp.safe(Z, 1);
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet dom(ss);
  dom.setSymbolicSet(C);
  std::cout << "Domain size: " << dom.getSize() << std::endl;

  scots::SymbolicSet controller(ss,is);
  controller.setSymbolicSet(C);
  std::cout << "Controller size: " << controller.getSize() << std::endl;
  controller.writeToFile("engine_controller.bdd");

  /*newly added*/
  ss.writeToFile("engine_state_grid.bdd");
  is.writeToFile("engine_input_grid.bdd");


  return 1;
}