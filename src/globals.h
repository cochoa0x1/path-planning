#ifndef GLOBALS_H
#define GLOBALS_H

#include <vector>


using namespace std;

typedef struct SimulationVars {
    double car_x;
    double car_y;
    double car_s;
    double car_d;
    double car_yaw;
    double car_speed;

    vector<double> previous_path_x;
    vector<double> previous_path_y;

    double end_path_s;
    double end_path_d;
}SimulationVars;


typedef struct MapVars {
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

}MapVars;

typedef struct OtherCar {

  int _id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;
  int lane;
  double v;
  double next_s;
  vector<double> path_x;
  vector<double> path_y;
}OtherCar;

#endif
