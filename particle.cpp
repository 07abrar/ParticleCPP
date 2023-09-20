#include<iostream>
#include <cstring>
#include<cmath>
#include<vector>
#include<algorithm>
#include<limits>
#include<fstream>
#include<ctime>

#define THREAD 4
#include<omp.h>

//define some variables
const int radius = 5;
const int fr = 1;
const double speed = 1.0;
const int fps = 20;
const double dt = 1.0/fps;
const int final_time = 200;
const int steps = final_time/dt;

// Define a struct to represent a particle
struct Particle {
    double x;
    double y;
    double vx;
    double vy;
};

//initialize particle location
std::vector<Particle> initialize_particles(int history){
    std::vector<Particle> particles(history);
    
    //#pragma omp parallel for num_threads(THREAD)
    for(int i=0; i<history; i++){
        double teta = ((double)rand() / RAND_MAX) * 2 * M_PI;
        double position = ((double)rand() / RAND_MAX) * radius;
        double teta2 = ((double)rand() / RAND_MAX) * 2 * M_PI;
        particles[i] = {position*cos(teta), position*sin(teta), speed*cos(teta2), speed*sin(teta2)};
    }
    return particles;
}

//initialize fuel
double fuel[1][2];
void initialize_fuel(double fuel[1][2]){
    fuel[0][0] = 0;
    fuel[0][1] = 0;
}

void update_particles(std::vector<Particle>& particles, int i) {
    std::vector<Particle> new_particles;
    std::vector<Particle> new_particles_to_add;

    //#pragma omp parallel for num_threads(THREAD)
    for(auto& particle : particles) {
        if (std::isnan(particle.x)){
            continue;
        }

        else if (particle.x >= fuel[0][0] - fr && particle.x <= fuel[0][0] + fr &&
                 particle.y >= fuel[0][1] - fr && particle.y <= fuel[0][1] + fr &&
                 sqrt(particle.x * particle.x + particle.y * particle.y) <= fr) {
            // Delete the particle as NaN
            particle.x = std::numeric_limits<double>::quiet_NaN();
            particle.y = std::numeric_limits<double>::quiet_NaN();
            particle.vx = std::numeric_limits<double>::quiet_NaN();
            particle.vy = std::numeric_limits<double>::quiet_NaN();

            new_particles = initialize_particles(2);
            for(auto& new_particle : new_particles) {
                new_particles_to_add.push_back(new_particle);
            }
        }
        else if (sqrt(particle.x * particle.x + particle.y * particle.y) < radius) {
            // Particle is inside the domain, update position
            particle.x += particle.vx * dt;
            particle.y += particle.vy * dt;
        }
        else {
            // Particle hits boundary, update velocity and position
            double distance_to_origin = sqrt(particle.x * particle.x + particle.y * particle.y);
            double normal_vector_x = particle.x / distance_to_origin;
            double normal_vector_y = particle.y / distance_to_origin;
            double dot_product = particle.vx * normal_vector_x + particle.vy * normal_vector_y;
            particle.vx -= 2 * dot_product * normal_vector_x;
            particle.vy -= 2 * dot_product * normal_vector_y;
            particle.x += particle.vx * dt;
            particle.y += particle.vy * dt;
        }
    }
    particles.insert(particles.end(), new_particles_to_add.begin(), new_particles_to_add.end());
    
    /*
    std::ofstream outputFile("particles_time_step_" + std::to_string(i) + ".csv", std::ios::app);
    if (outputFile.is_open()) {
        outputFile << "x" << "," << "y" << "\n";
        for (const auto& particle : particles) {
            if (!std::isnan(particle.x)) {
                outputFile << particle.x << "," << particle.y << "\n";
            }
        }
        outputFile.close();
    }
    else {
        std::cerr << "Error: Unable to open file for writing.\n";
        }
    */
}

//main loop of the program
int main(){
    std::clock_t c_start = std::clock();
    srand(time(0));

    int history = 5000;
    std::vector<Particle> particles = initialize_particles(history);
    int c=0;
    for(auto i:particles){
        c++;
    }

    std::cout<<"The length of array is "<<c<<std::endl;
    std::cout<<"Position 1"<<std::endl;
    std::cout<<particles[299].x<<std::endl;
    std::cout<<particles[299].y<<std::endl;
    
    initialize_fuel(fuel);
    for(int i=0; i<steps; i++){
        update_particles(particles, i);
    }


    int e=0;
    for(auto i:particles){
        e++;
    }
    std::cout<<"The length of array is "<<e<<std::endl;
    std::cout<<"Position 3"<<std::endl;
    std::cout<<particles[299].x<<std::endl;
    std::cout<<particles[299].y<<std::endl;

    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms / 1000.0 << " s\n";
    return 0;
}
