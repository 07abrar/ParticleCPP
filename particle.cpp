#include<iostream>
#include<cmath>
#include<vector>

//define some variables
const int history = 5000;
const int radius = 5;
const int fr = 1;
const double speed = 1.0;
const int fps = 20;
const double dt = 1.0/fps;
const int final_time = 200;
const int steps = final_time/dt;

//initialize particle location
double r[history][2];
void initialize_particles(double r[][2], int history){
    for(int i=0; i<history; i++){
        double teta = ((double)rand() / RAND_MAX) * 2 * M_PI;
        double position = ((double)rand() / RAND_MAX) * radius;
        r[i][0] = position*cos(teta);
        r[i][1] = position*sin(teta);
    }
}

//initialize fuel
double fuel[1][2];
void initialize_fuel(double fuel[1][2]){
    double fuel[0][0] = 0;
    double fuel[0][1] = 0;
}

//initialize particle speed
double v[history][2];
void initialize_speed(double v[][2], int history){
    for(int i=0; i<history; i++){
        double teta =((double)rand() / RAND_MAX) *2 * M_PI;
        V[i][0] = speed * cos(teta)
        v[i][1] = speed * sin(teta)
    }
}

//main loop of the program
int main(){
    srand(time(0));
    
    double record[steps][history][2];
    memcpy(record[0], r, sizeof(r));
    
    initialize_particles(r,history);
    initialize_fuel();
    initialize_speed(v, history);

    for(int i=1; i<steps; i++){
        for(int j=0; j<history; j++){
            if(r[j][0] == 11000){
                continue;
            }
            else if (r[j][0] >= fuel[0][0]-fr && r[j][0] <= fuel[0][0]+fr &&
                     r[j][0] >= fuel[0][0]-fr && r[j][0] <= fuel[0][0]+fr &&
                     sqrt(r[j][0]*r[j][0] + r[j][0]*r[j][0]) <= fr){
                r[j][0] = 11000;
                v[j][0] = 0;
                double r2[2][2];
                initialize_particles(r2, 2);
                double v2[2][2];
                initialize_particles(v2, 2);
                for (int k = 0; i<2, i++){
                    r[history+k][0] = r2[k][0]
                    r[history+k][1] = r2[k][1]
                    v[history+k][0] = v2[k][0]
                    v[history+k][1] = v2[k][1]
                }
                history += 2;
                double record2[steps][history][2];
                memcpy(record2, record, sizeof(record));
                memcpy(record2[i], r, sizeof(r));
                memcpy(record, record2, sizeof(record2));
                continue;
            }
            else if (sqrt(r[j][0]*r[j][0] + r[j][1]*r[j][1]) < radius) {
                r[j][0] += v[j][0]*dt;
                r[j][1] += v[j][1]*dt;
                continue;
            }
            else {
                double distance_to_origin = sqrt(r[j][0]*r[j][0] + r[j][1]*r[j][1]);
                double normal_vector[2] = {r[j][0]/distance_to_origin, r[j][1]/distance_to_origin};
                double dot_product = v[j][0]*normal_vector[0] + v[j][1]*normal_vector[1];
                v[j][0] -= 2 * dot_product * normal_vector[0];
                v[j][1] -= 2 * dot_product * normal_vector[1];
                r[j][0] += v[j][0] * dt;
                r[j][1] += v[j][1] * dt;
            }
        }
        memcpy(record[i], r, sizeof(r));
    }
    return 0;
}