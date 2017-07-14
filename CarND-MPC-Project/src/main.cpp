#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta = j[1]["steering_angle"];
	  double a = j[1]["throttle"];

	  delta *= -1; // Adjust for negative steering angle

	  // Way-points from the car's perspective
	  Eigen::VectorXd xvals(ptsx.size());
	  Eigen::VectorXd yvals(ptsy.size());

	  // Let's fill in the waypoints information by translating to origin
	  // and rotating the coordinate system to make the xvals horizontal
	  for (size_t i = 0; i < ptsx.size(); ++i) {
	    // Translate
	    double x_trans = ptsx[i]-px;
	    double y_trans = ptsy[i]-py;

	    // Rotate
	    xvals[i] = x_trans*cos(-psi) - y_trans*sin(-psi);
	    yvals[i] = x_trans*sin(-psi) + y_trans*cos(-psi);
	  }

	  auto coeffs = polyfit(xvals, yvals, 3);

	  double cte = polyeval(coeffs, 0);
	  double epsi = atan(coeffs[1]);

	  // Latency adjustment
	  double dt_lat = 0.1;
	  double x0 = 0;
	  double y0 = 0;
	  double psi0 = 0;
	  double v0 = v;
	  double Lf = 2.67;

	  if (dt_lat>0) {
	    // Updated using the initial steering angle
	    x0   += v*cos(delta)*dt_lat;
	    y0   += v*sin(delta)*dt_lat;
	    psi0 += v*delta/Lf*dt_lat;
	    v0   += a*dt_lat;

	    // Using the kinematic update for the cte and epsi
	    cte += v*sin(delta)*dt_lat;
	    epsi += v*delta/Lf*dt_lat;

	    // While we can use the kinematic update equations to update CTE and EPSE (as
	    // suggested by the first reviewer), re-avaluating using the polynomial using
	    // the updated x0 (from latency) is better as it provides a nonlinear and
	    // more accurate update.
	    /*
	    cte = polyeval(coeffs, x0);
	    epsi = atan(1*coeffs[1] +
			2*coeffs[2]*x0 +
			3*coeffs[3]*x0*x0);
	    */
	      }

	  // Fill the state and solve for vars
	  Eigen::VectorXd state(6);
	  state << x0, y0, psi0, v0, cte, epsi;

	  auto vars = mpc.Solve(state, coeffs);

          double steer_value    =  -vars[0]; // Steering angle is negative in rotated coordinates
          double throttle_value =  vars[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the waypoints/reference line (Yellow line)
          vector<double> next_x_vals;
          vector<double> next_y_vals;
	  int npoints = 10;
	  double dn = 5.0;
	  for (int i = 1; i < npoints+2; ++i) {
	    next_x_vals.push_back(dn*i);
	    next_y_vals.push_back(polyeval(coeffs, dn*i));
	  }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

	  //Display the MPC predicted trajectory (Green line)
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
	  for (size_t i = 2; i < vars.size(); ++i) {
	    if (i%2 == 0) {
	      mpc_x_vals.push_back(vars[i]);
	    } else {
	      mpc_y_vals.push_back(vars[i]);
	    }
	  }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(int(dt_lat*1000)));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
