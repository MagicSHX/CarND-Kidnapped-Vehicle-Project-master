/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (is_initialized) {
	return;
	}
	// Initializing the number of particles
	num_particles = 150;

	// Set standard deviations for x, y, and theta
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; i++) {

	Particle particle;
	particle.id = i;
	particle.x = dist_x(gen);
	particle.y = dist_y(gen);
	particle.theta = dist_theta(gen);
	particle.weight = 1.0f;
	particles.push_back(particle);
	}
	is_initialized = true;
	cout<<"initialized already"<<endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	normal_distribution<double> dist_x(0.0, std_x);
	normal_distribution<double> dist_y(0.0, std_y);
	normal_distribution<double> dist_theta(0.0, std_theta);

	for (int i = 0; i < num_particles; i++) {
		double theta = particles[i].theta;
	if ( fabs(yaw_rate) < 0.00001 ) { 
	  particles[i].x += velocity * delta_t * cos( theta );
	  particles[i].y += velocity * delta_t * sin( theta );
	} else {
	  particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) )+ dist_x(gen);
	  particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) )+ dist_y(gen);
	  particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
	}
	}
	cout<<"x, y : "<<particles[0].x<<'\t'<<particles[0].y<<endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	unsigned int Observations_num = observations.size();
	unsigned int Predictions_num = predicted.size();
	for (unsigned int i = 0; i < Observations_num; i++) { 
	double Distance_min = numeric_limits<double>::max();
	int indicator = -1;
	for (unsigned j = 0; j < Predictions_num; j++ ) {
	  double distance = (observations[i].x - predicted[j].x)*(observations[i].x - predicted[j].x)+(observations[i].y - predicted[j].y)*(observations[i].y - predicted[j].y);
	  if ( distance < Distance_min ) {
		Distance_min = distance;
		indicator = predicted[j].id;
	  }
	}
	observations[i].id = indicator;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  for(Particle& p : particles){    

    //check range, get close landmarks
    vector<LandmarkObs> predictions;    
    for(Map::single_landmark_s& landmark : map_landmarks.landmark_list){
      double d = dist(landmark.x_f, landmark.y_f, p.x, p.y);      
      if(d < sensor_range){	
	LandmarkObs new_obs = {};
	new_obs.id = landmark.id_i;
	new_obs.x  = landmark.x_f;
	new_obs.y  = landmark.y_f;
	predictions.push_back(new_obs);
      }
    }
        
		
		
		
    //transform observations to map's coordinate
    //helper variables
    const double x = p.x;
    const double y = p.y;
    const double theta = p.theta;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    //temp vector for transformed observations
    vector<LandmarkObs> particle_observations;    
    for(LandmarkObs& o : observations){
      LandmarkObs no = {};
      no.id = o.id;
      //Udacity Transformation (not working)
      //no.x +=  x*cos_theta + y*sin_theta;
      //no.y += -x*sin_theta + y*cos_theta;
      //http://planning.cs.uiuc.edu/node99.html Transformation
      no.x = o.x*cos_theta - o.y*sin_theta + x;
      no.y = o.x*sin_theta + o.y*cos_theta + y;     
      particle_observations.push_back(no);
    }
    
    //data association between landmarks and observations
    dataAssociation(predictions, particle_observations);   
    
    //record associations
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    //compute weights
    double weight=1;
    double std_2_pi = 2.0*M_PI*std_landmark[0]*std_landmark[1];
    double std_x_2  = 2.0*std_landmark[0]*std_landmark[0];
    double std_y_2  = 2.0*std_landmark[1]*std_landmark[1];
    for(LandmarkObs& o : particle_observations){
      //recover associated landmark
      Map::single_landmark_s m =  map_landmarks.landmark_list.at(o.id -1);//map_id->second];
      //compute Multivariate-Gaussian probability
      double e1 = pow(o.x - m.x_f, 2);
      double e2 = pow(o.y - m.y_f, 2);
      double e = (e1/std_x_2 + e2/std_y_2);
      double ee = exp(-e);
      double w = ee/std_2_pi;
      //prod of all weights
      weight *= w;
      //record association
      associations.push_back(o.id);
      sense_x.push_back(o.x);
      sense_y.push_back(o.y);      
    }
    //update particle with final weight
    p.weight = weight;
    //insert into weight vector
    weights.push_back(weight);
    //update particle's associations
    SetAssociations(p, associations, sense_x, sense_y);    
  }
  cout << "UPDATE WEIGHTS OK\n";
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Get weights and max weight.
	default_random_engine gen;

	std::discrete_distribution<> dist_weighted(weights.begin(),weights.end());
	std::vector<Particle> particles_sampled;
	for (int i = 0; i < particles.size() ; ++i) {
		int sample_index = dist_weighted(gen);
		particles_sampled.push_back(particles[sample_index]);
	}
	particles = particles_sampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();
	
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
