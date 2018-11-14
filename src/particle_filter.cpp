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

void ParticleFilter::init(double x, double y, double theta, double sigma[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (!is_initialized){
		// cout << "Initializing";
		// Number of particles
		num_particles = 200;

		double std_x = sigma[0];
		double std_y = sigma[1];	
		double std_theta = sigma[2];

		// This line creates a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(x, std_x);
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta, std_theta);

		for (int i = 0; i < num_particles; ++i) {
			Particle particle;
			particle.id = i;
			particle.x = dist_x(gen);
			particle.y = dist_y(gen);
			particle.theta = dist_theta(gen); 
			// Also init weight
			particle.weight = 1.0;

			particles.push_back(particle);
		}
		is_initialized = true;
	}


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	double std_x = std_pos[0];
	double std_y = std_pos[1];	
	double std_theta = std_pos[2];

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	// 
	// Calculate new state.
	for (int i = 0; i < num_particles; i++) {

		double theta = particles[i].theta;
		// is Too low, close to zero
		if (fabs(yaw_rate) < EPS) { 
		  particles[i].x += velocity * delta_t * cos(theta);
		  particles[i].y += velocity * delta_t * sin(theta);
		} else { 
		 // Yaw is 
		  particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t ) - sin(theta));
		  particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
		  particles[i].theta += yaw_rate * delta_t;
		}

		// Adding noise.
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	
	// Observations
	unsigned int n_observations = observations.size();
	// Close landmarks
  	unsigned int n_predictions = predicted.size();
	for (int i = 0; i < n_observations; i++){
		    // Initialize min distance as a really big number.
    	double min_distance = numeric_limits<double>::max();

    	// Initialize the found map in something not possible.
    	int map_id = -1;
		for (unsigned j = 0; j < n_predictions; j++) { // For each predition.

			double x_distance = observations[i].x - predicted[j].x;
			double y_distance = observations[i].y - predicted[j].y;
			// Eucladian distance
			double distance = sqrt(pow((x_distance), 2) + pow((y_distance) ,2));

			// If the "distance" is less than min, stored the id and update min.
			if ( distance < min_distance ) {
				min_distance = distance;
				map_id = predicted[j].id;
			}
	    }
		// Update the observation identifier.
    	observations[i].id = map_id;
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
	double sig_landmark_x = std_landmark[0];
	double sig_landmark_y = std_landmark[1];
	// For each particle 
	for (int  i = 0; i < num_particles; i++){
		// Particle measurement vehicle space 
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		vector<LandmarkObs> close_landmarks;
		// Find landmark close to the particle
		for (int j=0; j < map_landmarks.landmark_list.size(); j++){
			// Get a landmark
			float landmark_x = map_landmarks.landmark_list[j].x_f;
			float landmark_y = map_landmarks.landmark_list[j].y_f;
			// 
			int id = map_landmarks.landmark_list[j].id_i;
			// Diff 
			double d_x = x - landmark_x;
			double d_y = y - landmark_y;
			// if the sum of the difference squared is within sensor range 
			// then take this landmark 
			
			if (pow(d_x, 2) + pow(d_y, 2) <= pow(sensor_range,2)) {
				close_landmarks.push_back(LandmarkObs{id, landmark_x, landmark_y});
      		}
		}
		/*##############################################################################
									        Transform
		################################################################################*/
		vector<LandmarkObs> mapped_observations;
	    for(unsigned int j = 0; j < observations.size(); j++) {
			double x_map = cos(theta) * observations[j].x - sin(theta) * observations[j].y + x;
			double y_map = sin(theta) * observations[j].x + cos(theta) * observations[j].y + y;
			mapped_observations.push_back(LandmarkObs{observations[j].id, x_map, y_map});
    	}



		/*##############################################################################
									        Association
		################################################################################*/

    	
		dataAssociation(close_landmarks, mapped_observations);


		/*##############################################################################
									        Weight
		################################################################################*/

		// Reinitializing weight
		particles[i].weight = 1.0;
		
		for(unsigned int j = 0; j < mapped_observations.size(); j++) {
			double x_map = mapped_observations[j].x;
      		double y_map = mapped_observations[j].y;
      		int id = mapped_observations[j].id;
      		
      		double landmark_x, landmark_y;
		    unsigned int k = 0;
		    unsigned int n_landmarks = close_landmarks.size();
		    
	      	bool found = false;

	      	while(!found && k < n_landmarks) {
	      		// cout << "Found" << endl;
	        	if (close_landmarks[k].id == id) {
	          		found = true;
	          		landmark_x = close_landmarks[k].x;
	          		landmark_y = close_landmarks[k].y;
	        	}
	        	k++;
	      	}
	      	// Difference
	      	double x_diff_squared = pow(x_map - landmark_x, 2);
	      	double y_diff_squared = pow(y_map - landmark_y, 2);
	      	// //calculate exponent
	      	// double exponent = exp(-(x_diff_squared / (2 * pow(sig_landmark_x, 2)) + (y_diff_squared /( 2 * pow(sig_landmark_y, 2)))));
	      	// // calculate normalization term
	      	// double gauss_norm = (1 / (2 * M_PI * sig_landmark_x * sig_landmark_y));
	      	// // calculate weight using normalization terms and exponent
	      	// double weight = gauss_norm * exponent;
	      	double weight = ( 1/(2*M_PI*sig_landmark_x*sig_landmark_y)) * exp( -(x_diff_squared/(2*sig_landmark_x*sig_landmark_x) + (y_diff_squared/(2*sig_landmark_y*sig_landmark_y)) ) );

	      	
	      	// ZERO check
	      	if (weight == 0){
	      		particles[i].weight *= EPS;
		     } else {
		     	particles[i].weight *= weight;
		     }
		}

		// observations = mapped_observations;
		// map_landmarks = 

		// //calculate normalization term
		// double gauss_norm = (1/(2 * M_PI * sig_x * sig_y))

		// //calculate exponent
		// double exponent = ((x - mu_x)**2)/(2 * sig_x**2) + ((y - mu_y)**2)/(2 * sig_y**2)

		// //calculate weight using normalization terms and exponent
		// double weight = gauss_norm * math.exp(-exponent)

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<double> weights;
	//initialize with a small value 
  	double max_weight = numeric_limits<double>::min();

  	// Get max weight 
  	for(int i = 0; i < num_particles; i++) {
    	weights.push_back(particles[i].weight);
    	if (particles[i].weight > max_weight) {
      		max_weight = particles[i].weight;
    	}
  	}

	uniform_real_distribution<double> dist_cont(0.0, 2 * max_weight);
	uniform_int_distribution<int> dist_index(0, num_particles - 1);
	
	// Wheel or pie chart
	// Pick an index from the distribution
	int index = dist_index(gen);

	double beta = 0.;
	vector<Particle> resampled_particles;


	for (int i = 0; i < num_particles; i++){
		beta += dist_cont(gen);
		while (beta > weights[index]){
			beta -= weights[index];
			// To make sure we dont go out of bounds
			index = (index + 1) % num_particles;
			// index++;
		}
		resampled_particles.push_back(particles[index]);
	}

	// Set the particles to the newly weighted resampled particles 
	particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
