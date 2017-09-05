#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "globals.h"

using namespace std;

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
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

/* coordinate transform: shift then rotate */
vector<double> transform(double x, double y, double theta, double x0, double y0){
  vector<double> new_point;

  //shift
  double dx = x - x0;
  double dy = y - y0;

  //rotate
  double new_x = dx*cos(0-theta)-dy*sin(0-theta);
  double new_y = dx*sin(0-theta)+dy*cos(0-theta);

  new_point.push_back(new_x);
  new_point.push_back(new_y);

  return new_point;
}
/* rotate opposite direction then shift */
vector<double> inv_transform(double x, double y, double theta, double x0, double y0){

  //rotate
  vector<double> rotated = transform(x, y, -1.0*theta, 0.0, 0.0);

  //shift back
  rotated[0] = rotated[0]+x0;
  rotated[1] = rotated[1]+y0;
  return rotated;
}

void transform_vector(vector<double> &xs, vector<double> &ys, double theta, double x0, double y0){
  /* transforms in place */

  for(int i=0; i < xs.size(); i++){
    vector<double> new_v = transform(xs[i],ys[i], theta,x0,y0);
    xs[i] = new_v[0];
    ys[i] = new_v[1];
  }
}

void inv_transform_vector(vector<double> &xs, vector<double> &ys, double theta, double x0, double y0){
  /* transforms in place backwards direction*/

  for(int i=0; i < xs.size(); i++){
    vector<double> new_v = inv_transform(xs[i],ys[i], theta,x0,y0);
    xs[i] = new_v[0];
    ys[i] = new_v[1];
  }
}

double min_distance_between_paths(vector<double> x1, vector<double> y1,vector<double> x2, vector<double> y2){
  //kinda cheat here and just check the min distance for the vertices of the path
  double min_distance = 1e7;

  for( int i=0; i< x1.size(); i++){
    for( int j =0; j < x2.size(); j++){
      double dx = x1[i]-x2[j];
      double dy = y1[i]-y2[j];
      double d = sqrt(dx*dx+dy*dy);
      if( d < min_distance){
        min_distance = d;
      }
    }
  }

  return min_distance;
}

/* construct a path ending in lane that drives the car at target_speed */
void generate_path(vector<double> &xs, vector<double> &ys, double aggro, double target_speed, int lane, SimulationVars smv, MapVars mapv){


  vector<double> anchor_x;
  vector<double> anchor_y;

  int N = smv.previous_path_x.size();

  if( N < 2){
    double previous_x = smv.car_x - cos(smv.car_yaw);
    double previous_y = smv.car_y - sin(smv.car_yaw);

    anchor_x.push_back(previous_x);
    anchor_y.push_back(previous_y);

    anchor_x.push_back(smv.car_x);
    anchor_y.push_back(smv.car_y);

    smv.end_path_s = smv.car_s;

  }
  else{
    anchor_x.push_back(smv.previous_path_x[N-2]);
    anchor_y.push_back(smv.previous_path_y[N-2]);

    anchor_x.push_back(smv.previous_path_x[N-1]);
    anchor_y.push_back(smv.previous_path_y[N-1]);
  }

  /* add some target points for the spline */
  for(int i=0; i<3;i++){
    vector<double> xy  = getXY(smv.end_path_s+(i+1)*aggro,4*lane+2, mapv.map_waypoints_s, mapv.map_waypoints_x, mapv.map_waypoints_y);
    anchor_x.push_back(xy[0]);
    anchor_y.push_back(xy[1]);
  }

  /* transform the anchor points to "car" coordinates so that the first point is at 0,0
  and the path faces toward the x direction. In otherwords do not use car_yaw use ref_angle
  */
  double ref_angle = atan2(anchor_y[1]-anchor_y[0],anchor_x[1]-anchor_x[0]);

  /* transform to future car-at-path-end coordinates */
  double x0 = anchor_x[1];
  double y0 = anchor_y[1];
  transform_vector(anchor_x,anchor_y,ref_angle,x0,y0);

  /* build the spline */
  tk::spline spline_path;
  spline_path.set_points(anchor_x,anchor_y);

  vector<double> next_x_vals = smv.previous_path_x;
  vector<double> next_y_vals = smv.previous_path_y;

  transform_vector(next_x_vals,next_y_vals,ref_angle,x0,y0);

  double last_x = anchor_x[1];

  //add the new path elements
  //double dist_inc = .2;
  for( int i = 0; i < 50-N ; i++){
    double _x = last_x + (i+1)*.02*target_speed*0.44704; //0.44704 mph to meter per second
    double _y = spline_path(_x);

    next_x_vals.push_back(_x);
    next_y_vals.push_back(_y);
  }


  inv_transform_vector(next_x_vals,next_y_vals,ref_angle,x0,y0);

  xs = next_x_vals;
  ys = next_y_vals;
}

void path_dt(vector<double> x, vector<double> y, vector<double> &dx, vector<double> &dy){
  //returns the discrete time dx, dy

  vector<double> _dx;
  vector<double> _dy;

  const double dt = .02;

  if(x.size() <=1){
      return;
  }

  for(int i = 1; i < x.size(); i++){
    _dx.push_back( (x[i] - x[i-1])/dt );
    _dy.push_back( (y[i] - y[i-1])/dt );
  }
  dx = _dx;
  dy = _dy;
}

vector<double> path_norm_stats(vector<double> x, vector<double> y){
  //returns average norm, max_norm, min_norm
  double _avg =0;
  double _max = -1e7;
  double _min = 1e7;

  for( int i =0; i< x.size(); i++){
    double norm = sqrt(x[i]*x[i] + y[i]*y[i]);
    _avg += norm;
    if( norm > _max){
      _max = norm;
    }
    if( norm < _min){
      _min = norm;
    }

  }
  if(x.size() !=0){
    _avg /= x.size();
  }
  return {_avg, _min, _max};
}
double path_score(vector<double> xs
                ,vector<double> ys
                ,double target_speed
                ,int target_lane
                ,double current_speed
                ,int current_lane
                ,vector<OtherCar> cars
                ,SimulationVars smv){
  //lower score is better

  vector<double> vx, vy;
  path_dt(xs,ys, vx, vy);

  vector<double> v_stats = path_norm_stats(vx,vy);
  cout << "v stats: " << v_stats[0] << ", " << v_stats[1] << ", " << v_stats[2] << endl;


  vector<double> ax, ay;
  path_dt(vx,vy, ax, ay);

  vector<double> a_stats = path_norm_stats(ax,ay);
  cout << "a stats: " << a_stats[0] << ", " << a_stats[1] << ", " << a_stats[2] << endl;

  vector<double> jx, jy;
  path_dt(ax,ay, jx, jy);

  vector<double> j_stats = path_norm_stats(jx,jy);
  cout << "j stats: " << j_stats[0] << ", " << j_stats[1] << ", " << j_stats[2] << endl;

  //in general don't like changing lanes
  double lane_change_penalty = 10*abs(target_lane-current_lane);

  //slow speed penalty
  double speed_penalty = 10.0*(49.5 - target_speed);
  if(target_speed > 49.5 || target_speed < 0){
    speed_penalty = 2000;
  }

  //acceleration penalty, in general don't want to accelerate
  double accel_penalty = 1.0*abs(target_speed-current_speed);



  //if there is going to be a head on collision, have a large penalty

  //if there is a car ahead of us change lanes

  double collision_penalty =0;
  for(auto car : cars){

      //see if paths intersect. get the min distance between the paths_x
      double min_distance = min_distance_between_paths(car.path_x, car.path_y, xs, ys);
      cout << "car id: " << car._id << " lane: " << car.lane << " distance: " << min_distance << endl;;

      //if the car is in the same lane, and it is currently in front
      if(car.lane == target_lane && min_distance < 12){
          collision_penalty += 1e6/(min_distance*min_distance+1);
      }

      //if the car is between the lanes we are changing, add a penalty if we would be close
      if(car.lane != target_lane && car.lane != current_lane && abs(current_lane-target_lane) ==2 && min_distance< 10){
        collision_penalty += 1e6/(min_distance*min_distance+1);
      }
      //
      // //if the car is in front of us either in this lane or when we change
      // double d = fabs(car.next_s - smv.car_s);
      //
      // if(car.lane == target_lane && d < 60 && car.next_s >= smv.car_s){
      //   collision_penalty+=1e7/(d*d+1);
      // }
      //
      // //if we are changing lanes, we don't want to cut another car off
      // if(car.lane == target_lane && d < 45 && (car.next_s < smv.car_s) && current_lane!=target_lane){
      //   collision_penalty+=1000;
      // }
      //
      // //guard against the double lange change
      // if(car.lane != target_lane && car.lane != current_lane && d < 45 && abs(current_lane-target_lane) ==2){
      //   collision_penalty+=1000;
      // }
  }

  //don't drive slowly in the fast lane. or ie if you are gonna go fast
  //drive in the fast lane
  double fast_lane = 0;
  if(target_lane != 0 && current_speed > 30){
    fast_lane = 100;
  }

  return lane_change_penalty+speed_penalty+collision_penalty+accel_penalty;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  //double target_speed = 5.0;

  // Load up map values
  MapVars mapv; //this will make it easier to pass these to functions

  mapv.map_waypoints_x=map_waypoints_x;
  mapv.map_waypoints_y=map_waypoints_y;
  mapv.map_waypoints_s=map_waypoints_s;
  mapv.map_waypoints_dx=map_waypoints_dx;
  mapv.map_waypoints_dy=map_waypoints_dy;


  h.onMessage([&mapv,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            int lane = int(floor(car_d/4));
            /* unpack the variables */
            SimulationVars smv;
            smv.car_x = car_x;
            smv.car_y = car_y;
            smv.car_s = car_s;
            smv.car_d = car_d;
            smv.car_yaw = car_yaw;
            smv.car_speed = car_speed;
            smv.end_path_s = end_path_s;
            smv.end_path_d = end_path_d;

            //convert whatever the previous_path_* vars are into std::vector<double>
            vector<double> prev_path_x;
            vector<double> prev_path_y;
            for(int i= 0; i < previous_path_x.size(); i++){
              prev_path_x.push_back( previous_path_x[i]);
              prev_path_y.push_back( previous_path_y[i]);
            }

            smv.previous_path_x = prev_path_x;
            smv.previous_path_y = prev_path_y;

            /* end */

            /* process the sensor fusion var */
            vector<OtherCar> cars;
            for(int i=0; i < sensor_fusion.size(); i++){
              OtherCar car;
              car._id = sensor_fusion[i][0];
              car.x = sensor_fusion[i][1];
              car.x = sensor_fusion[i][2];
              car.vx  = sensor_fusion[i][3];
              car.vy = sensor_fusion[i][4];
              car.s = sensor_fusion[i][5];
              car.d = sensor_fusion[i][6];
              //0-4:2 4-8-1 8-12 2
              car.lane = int(floor(car.d/4));
              car.v = sqrt(car.vx*car.vx+car.vy*car.vy);
              car.next_s = car.s + .02*car.v;

              //calculate the car's path for the next time interval
              vector<double> car_path_x;
              vector<double> car_path_y;
              for(int i=0; i < 20; i++){
                double next_s = car.s + (i)*.02*car.v;
                double next_d = car.d;
                vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                car_path_x.push_back(xy[0]);
                car_path_y.push_back(xy[1]);
              }

              car.path_x = car_path_x;
              car.path_y = car_path_y;

              cars.push_back(car);
            }

          	json msgJson;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            vector<int> lane_options;
            lane_options.push_back(0);
            lane_options.push_back(1);
            lane_options.push_back(2);

            vector<double> dv_options;
            dv_options.push_back(-2);
            dv_options.push_back(-1);
            dv_options.push_back(0);
            dv_options.push_back(1);
            dv_options.push_back(2.0);
            dv_options.push_back(3.0);

            vector<vector<double>> paths_x;
            vector<vector<double>> paths_y;
            vector<double> path_scores;

            cout << "current car speed: " << car_speed << endl;
            for(auto new_lane : lane_options){
              for(auto dv : dv_options){
                double target_speed = car_speed+dv;

                vector<double> next_x_vals;
                vector<double> next_y_vals;
                generate_path(next_x_vals, next_y_vals,30.0, target_speed, new_lane, smv, mapv);
                double score = path_score(next_x_vals, next_y_vals, target_speed, new_lane, car_speed, lane,cars,smv);
                cout << " lane: "<< new_lane << " , dv: " << dv << " score: " << score << endl;

                paths_x.push_back(next_x_vals);
                paths_y.push_back(next_y_vals);
                path_scores.push_back(score);

              }
            }

            //find the lowest score path_scores
            double lowest =1e8;
            int lowest_index =0;
            for(int i =0; i< path_scores.size(); i++){
              if(path_scores[i] < lowest){
                lowest = path_scores[i];
                lowest_index=i;
              }
            }

            cout <<"best score: " << lowest << endl;

          	msgJson["next_x"] = paths_x[lowest_index];
          	msgJson["next_y"] = paths_y[lowest_index];


          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
