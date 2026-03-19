/*
 * dcdc.cc
 *
 *  created on: 09.29.2025
 *      author: Daniel A.
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <iostream>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound_for_PAC.hh"
// #include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

/* state space dim */
#define sDIM 2
#define iDIM 1

/* data types for the ode solver */
typedef std::array<double,2> state_type;
typedef std::array<double,1> input_type;

/* sampling time */
const double tau = 0.5;
/* state space interval*/
const double hn1 = 0.0005;
/* number of samples */
double const N = 300;
/* dynamics parameters */
const double r0=1.0 ; 
const double vs = 1.0 ;
const double rl = 0.05 ;
const double rc = rl / 10 ;
const double xl = 3.0 ;
const double xc = 70.0 ;


/* we discretize the continuous time dcdc ode by forward Euler (the result is stored in x)  */
auto  dcdc_post = [](state_type &x, input_type &u) -> void {

  /* the difference equation describing the vehicle */
  auto rhs =[](state_type& dx,  const state_type &x, input_type &u) {
    const double b[2]={vs/xl, 0};

    double a[2][2];
    if(u[0]==1) {
      a[0][0] = -rl / xl;
      a[0][1] = 0;
      a[1][0] = 0;
      a[1][1] = (-1 / xc) * (1 / (r0 + rc));
    } else {
      a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
      a[0][1] =  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
      a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
      a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
    }
    dx[0] = a[0][0]*x[0]+a[0][1]*x[1] + b[0];
    dx[1] = a[1][0]*x[0]+a[1][1]*x[1] + b[1];
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
    const double b[2]={vs/xl, 0};

    double a[2][2];
    if(u[0]==1) {
      a[0][0] = -rl / xl;
      a[0][1] = 0;
      a[1][0] = 0;
      a[1][1] = (-1 / xc) * (1 / (r0 + rc));
    } else {
      a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
      a[0][1] =  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
      a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
      a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
    }
    dx[0] = a[0][0]*x[0]+a[0][1]*x[1] + b[0];
    dx[1] = a[1][0]*x[0]+a[1][1]*x[1] + b[1];

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
    
    // check out the radius 
    for(int i=0; i<sDIM; i++){
      std::cout << r[i] << "\n";
    } 
};

// /* computation of the growth bound (the result is stored in r)  */
// auto radius_post = [](state_type &r, input_type &u) -> void {

//   /* the ode to determine the radius of the cell which over-approximates the
//    * attainable set see: http://arxiv.org/abs/1503.03715v1 */
//   auto growth_bound_ode = [](state_type &drdt,  const state_type &r, const input_type &u) {
//     /* for the dcdc boost converter the growth bound is simply given by the metzler matrix of the system matrices */ 
//     const double r0=1.0 ; 
//     const double rl = 0.05 ;
//     const double rc = rl / 10 ;
//     const double xl = 3.0 ;
//     const double xc = 70.0 ;

//     double a[2][2];
//     if(u[0]==1) {
//       a[0][0] = -rl / xl;
//       a[0][1] = 0;
//       a[1][0] = 0;
//       a[1][1] = (-1 / xc) * (1 / (r0 + rc));
//     } else {
//       a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
//       a[0][1] =  ((1 / xl) * (r0 / (r0 + rc))) / 5 ;
//       a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
//       a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
//     }
//     drdt[0] = a[0][0]*r[0]+a[0][1]*r[1];
//     drdt[1] = a[1][0]*r[0]+a[1][1]*r[1];
//   };
//   ode_solver(growth_bound_ode,r,u);
// };

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
  double lb[sDIM]={1.15,5.45};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={1.55,5.85}; 
  /* grid node distance diameter */
  double eta[sDIM]={hn1,hn1};   
  scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);
  ss.addGridPoints();

  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  double ilb[iDIM]={1};  
  double iub[iDIM]={2}; 
  double ieta[iDIM]={1};   
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
  abstraction.computeTransitionRelation(dcdc_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

  /****************************************************************************/
  /* we continue with the specification */
  /****************************************************************************/
  /* construct SymbolicSet for the safe set */
  double eps = 2*0.0005; /*approximate value of the ASF, obtained from the .ipynb file*/
  double H[4*sDIM]={-1, 0,
                     1, 0,
                     0,-1,
                     0, 1};
  /* add outer approximation of P={ x | H x<= h } form state space */
  double h[4] = {-1.15+eps,1.55-eps,-5.45+eps, 5.85-eps};
  /* initialize the safe set with the ss 
   * in order to obtain all the necessary information */
  scots::SymbolicSet safe(ss);
  safe.addPolytope(4,H,h, scots::INNER);
  std::cout << "Safe set details:" << std::endl;
  safe.writeToFile("dcdc.bdd");

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
  controller.writeToFile("dcdc_controller.bdd");

  /*newly added*/
  ss.writeToFile("dcdc_state_grid.bdd");
  is.writeToFile("dcdc_input_grid.bdd");


  return 1;
}

// int main() {
//   /* to measure time */
//   TicToc tt;
//   /* there is one unique manager to organize the bdd variables */
//   Cudd mgr;

//   /****************************************************************************/
//   /* construct SymbolicSet for the state space */
//   /****************************************************************************/
//   /* setup the workspace of the synthesis problem and the uniform grid */
//   /* lower bounds of the hyper rectangle */
//   double lb[sDIM]={0.65,4.95};  
//   /* upper bounds of the hyper rectangle */
//   double ub[sDIM]={1.65,5.95};  
//   /* grid node distance diameter */
//   double eta[sDIM]={hn1,hn1};  //0.005
//   scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);
//   ss.addGridPoints();
//   /****************************************************************************/
//   /* construct SymbolicSet for the input space */
//   /****************************************************************************/
//   double ilb[iDIM]={1};  
//   double iub[iDIM]={2}; 
//   double ieta[iDIM]={1};   
//   scots::SymbolicSet is(mgr,iDIM,ilb,iub,ieta);
//   is.addGridPoints();

//   /****************************************************************************/
//   /* setup class for symbolic model computation */
//   /****************************************************************************/
//   scots::SymbolicSet sspost(ss,1); /* create state space for post variables */ 
//   scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
//   /* compute the transition relation */
//   tt.tic();
//   abstraction.computeTransitionRelation(dcdc_post, radius_post);
//   std::cout << std::endl;
//   tt.toc();
//   /* get the number of elements in the transition relation */
//   std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

//   /****************************************************************************/
//   /* we continue with the controller synthesis for FG (target) */
//   /****************************************************************************/
//   /* construct SymbolicSet for target (it is a subset of the state space)  */
//   scots::SymbolicSet target(ss);
//   /* add inner approximation of P={ x | H x<= h } form state space */
//   double H[4*sDIM]={-1, 0,
//                      1, 0,
//                      0,-1,
//                      0, 1};
//   double h[4] = {-1.1,1.6,-5.4, 5.9};
//   target.addPolytope(4,H,h, scots::INNER);
//   std::cout << "Target set details:" << std::endl;
//   target.writeToFile("dcdc_target.bdd");

//   /* we setup a fixed point object to compute reach and stay controller */
//   scots::FixedPoint fp(&abstraction);
//   /* the fixed point algorithm operates on the BDD directly */
//   BDD T = target.getSymbolicSet();

//   /* we implement the nested fixed point algorithm
//    *
//    * mu X. nu Y. ( pre(Y) & T ) | pre(X)
//    *
//    */
//   tt.tic();
//   size_t i,j;
//   /* outer fp*/
//   BDD X=mgr.bddOne();
//   BDD XX=mgr.bddZero();
//   /* inner fp*/
//   BDD Y=mgr.bddZero();
//   BDD YY=mgr.bddOne();
//   /* the controller */
//   BDD C=mgr.bddZero();
//   BDD U=is.getCube();
//   /* as long as not converged */
//   for(i=1; XX != X; i++) {
//     X=XX;
//     BDD preX=fp.pre(X);
//     /* init inner fp */
//     YY = mgr.bddOne();
//     for(j=1; YY != Y; j++) {
//       Y=YY;
//       YY= ( fp.pre(Y) & T ) | preX;
//     }
//     XX=YY;
//     std::cout << "Iterations inner: " << j << std::endl;
//     /* remove all (state/input) pairs that have been added
//      * to the controller already in the previous iteration * */
//     BDD N = XX & (!(C.ExistAbstract(U)));
//     /* add the remaining pairs to the controller */
//     C=C | N;
//     //std::cout << C.CountMinterm(17) << std::endl;
//   }
//   std::cout << "Iterations outer: " << i << std::endl;
//   tt.toc();

//   /****************************************************************************/
//   /* last we store the controller as a SymbolicSet 
//    * the underlying uniform grid is given by the Cartesian product of 
//    * the uniform gird of the space and uniform gird of the input space */
//   /****************************************************************************/
//   scots::SymbolicSet controller(ss,is);
//   controller.setSymbolicSet(C);
//   std::cout << "Domain size: " << controller.getSize() << std::endl;
//   controller.writeToFile("dcdc_controller.bdd");

//   ss.writeToFile("dcdc_state_grid.bdd");
//   is.writeToFile("dcdc_input_grid.bdd");


//   return 1;
// }

