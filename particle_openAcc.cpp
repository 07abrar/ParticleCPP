#include<iostream>
#include <cstring>
#include<cmath>
#include<vector>
#include<algorithm>
#include<limits>
#include<fstream>
#include<ctime>
#include<iomanip>
#include <sstream>
#ifdef _OPENACC
#include <openacc.h>
#endif



//define some variables
const int radius = 5;
const int fr = 1;
const double speed = 1.0;
const int fps = 20;
const double dt = 1.0/fps;
const int final_time = 10;
const int steps = final_time/dt;

// Define a struct to represent a particle
struct Particle {
    double x;
    double y;
    double vx;
    double vy;
};

//initialize particle
std::vector<Particle> initialize_particles(int history){
    std::vector<Particle> particles(history);

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

void update_particles(std::vector<Particle>& particles, int history, std::vector<int>& collision, int j) {

    #pragma acc kernels
    for(int i=0; i<history; i++) {
        if (std::isnan(particles[i].x)){
            continue;
        }

        else if (particles[i].x >= fuel[0][0] - fr && particles[i].x <= fuel[0][0] + fr &&
                 particles[i].y >= fuel[0][1] - fr && particles[i].y <= fuel[0][1] + fr &&
                 sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y) <= fr) {
            // Delete the particle as NaN
            particles[i].x = std::numeric_limits<double>::quiet_NaN();
            particles[i].y = std::numeric_limits<double>::quiet_NaN();
            particles[i].vx = std::numeric_limits<double>::quiet_NaN();
            particles[i].vy = std::numeric_limits<double>::quiet_NaN();
            collision[i] = 1;
        }
        else if (sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y) < radius) {
            // Particle is inside the domain, update position
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
        }
        else {
            // Particle hits boundary, update velocity and position
            double distance_to_origin = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
            double normal_vector_x = particles[i].x / distance_to_origin;
            double normal_vector_y = particles[i].y / distance_to_origin;
            double dot_product = particles[i].vx * normal_vector_x + particles[i].vy * normal_vector_y;
            particles[i].vx -= 2 * dot_product * normal_vector_x;
            particles[i].vy -= 2 * dot_product * normal_vector_y;
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
        }
    }
    std::ostringstream filename;
    filename << "particles_time_step_" << std::setfill('0') << std::setw(3) << j << ".csv";
    std::ofstream outputFile(filename.str(), std::ios::app);
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
}

//main loop of the program
int main(){
    std::clock_t c_start = std::clock();
    srand(time(0));

    int history = 500;
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
    for(int j=0; j<steps; j++){
        std::vector<int> collision(history, 0);
        #pragma acc data copy(particles, collision)
        update_particles(particles, history, collision, j);
        std::vector<Particle> new_particles;
        std::vector<Particle> new_particles_to_add;
        for(int k=0; k<history; k++){
            if(collision[k]==1){
                new_particles = initialize_particles(2);
                for(auto& new_particle : new_particles) {
                    new_particles_to_add.push_back(new_particle);
                    }
                history += 2;
            }
        }
        particles.insert(particles.end(), new_particles_to_add.begin(), new_particles_to_add.end());
        new_particles.clear();
        new_particles_to_add.clear();
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
