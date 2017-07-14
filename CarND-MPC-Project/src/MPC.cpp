#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the timestep length and duration
size_t N = 15;
double dt = 0.15;

// Start-indices for the various values
size_t x_start     = 0;                    // N values
size_t y_start     = x_start + N;	   // N values
size_t psi_start   = y_start + N;	   // N values
size_t v_start     = psi_start + N;	   // N values
size_t cte_start   = v_start + N;	   // N values
size_t epsi_start  = cte_start + N;	   // N values
size_t delta_start = epsi_start + N;       // N-1 values
size_t a_start     = delta_start + N - 1;  // N-1 values

// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set between 40 - 100 mph.
double ref_cte  = 0;
double ref_epsi = 0;
double ref_v    = 80;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // MPC Implementation (mainly repurposed from Quiz solution)
    // `fg` a vector of the cost constraints,
    // `vars` is a vector of variable values (state & actuators)

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.


    // Weights for the cost function
    const double w_cte    = 1000;
    const double w_epsi   = 1000;
    const double w_dv     = 1;
    const double w_delta  = 100;
    const double w_a      = 10;
    const double w_ddelta = 10;
    const double w_da     = 10;
    
    // Initialize
    fg[0] = 0;

    for (size_t t = 0; t < N; t++) {
      // Minimize cross-track error
      fg[0] += w_cte  * CppAD::pow(vars[cte_start + t], 2);

      // Minimize error in direction
      fg[0] += w_epsi * CppAD::pow(vars[epsi_start + t], 2);

      // Minimize deviation from reference velocity
      fg[0] += w_dv * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    for (unsigned int t = 0; t < N - 1; t++) {
      // Minimize use of steering
      fg[0] += w_delta * CppAD::pow(vars[delta_start + t], 2);

      // Minimize use of throttle
      fg[0] += w_a * CppAD::pow(vars[a_start + t], 2);
    }

    for (unsigned int t = 0; t < N - 2; t++) {
      // Minimize sudden turns
      fg[0] += w_ddelta * CppAD::pow(vars[delta_start + t + 1] -
				     vars[delta_start + t], 2);

      // Minimize sudden accelerations or braking
      fg[0] += w_da * CppAD::pow(vars[a_start + t + 1] -
				 vars[a_start + t], 2);
    }

    // The state at time 0
    fg[x_start    + 1] = vars[x_start];
    fg[y_start    + 1] = vars[y_start];
    fg[psi_start  + 1] = vars[psi_start];
    fg[v_start    + 1] = vars[v_start];
    fg[cte_start  + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    // The rest of the constraints
    for (unsigned int t = 0; t < N-2; t++) {
      // The state at time t.
      AD<double> x0    = vars[x_start    + t];
      AD<double> y0    = vars[y_start    + t];
      AD<double> psi0  = vars[psi_start  + t];
      AD<double> v0    = vars[v_start    + t];
      AD<double> cte0  = vars[cte_start  + t];
      AD<double> epsi0 = vars[epsi_start + t];


      // The state at time t+1 .
      AD<double> x1    = vars[x_start    + t + 1];
      AD<double> y1    = vars[y_start    + t + 1];
      AD<double> psi1  = vars[psi_start  + t + 1];
      AD<double> v1    = vars[v_start    + t + 1];
      AD<double> cte1  = vars[cte_start  + t + 1];
      AD<double> epsi1 = vars[epsi_start + t + 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t];
      AD<double> a0     = vars[a_start     + t];

      // Apply polynomial equation for CTE
      AD<double> f0 =             
	coeffs[0] +		  
	coeffs[1] * x0 +	  
	coeffs[2] * x0 * x0 +     
	coeffs[3] * x0 * x0 * x0;
      
      // Tangent for psi
      AD<double> psides0 = CppAD::atan(                    
				       1*coeffs[1] +       
				       2*coeffs[2]*x0 +
				       3*coeffs[3]*x0*x0);

      // Cost variables by application of predictive model from time t0 to t1
      fg[x_start    + t + 2] = x1    - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[y_start    + t + 2] = y1    - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[psi_start  + t + 2] = psi1  - (psi0 + v0 * delta0 / Lf * dt);
      fg[v_start    + t + 2] = v1    - (v0 + a0 * dt);
      fg[cte_start  + t + 2] = cte1  - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[epsi_start + t + 2] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];


  size_t n_vars = N*6 + (N-1)*2;
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_start   ] = x;
  vars[y_start   ] = y;
  vars[psi_start ] = psi;
  vars[v_start   ] = v;
  vars[cte_start ] = cte;
  vars[epsi_start] = epsi;
  
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] =  1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] =  0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] =  1.0;
  }
 
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // The bounds are fixed at starting point
  constraints_lowerbound[x_start   ] = x;
  constraints_lowerbound[y_start   ] = y;
  constraints_lowerbound[psi_start ] = psi;
  constraints_lowerbound[v_start   ] = v;
  constraints_lowerbound[cte_start ] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start   ] = x;
  constraints_upperbound[y_start   ] = y;
  constraints_upperbound[psi_start ] = psi;
  constraints_upperbound[v_start   ] = v;
  constraints_upperbound[cte_start ] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);
  for (size_t i = 0; i<N-1; ++i) {
    result.push_back(solution.x[x_start + i]);
    result.push_back(solution.x[y_start + i]);  
  }
  return result;
}
