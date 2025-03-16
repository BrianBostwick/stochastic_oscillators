#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 9.81     // Gravitational acceleration (m/s^2)
#define L 1.0      // Length of pendulum (m)
#define T_MAX 30.0 // Maximum simulation time
#define DT 0.01    // Time step

// Structure to store state variables
typedef struct {
    double theta;
    double omega;
} State;

// Global variables for stochastic kicks
double *kick_times;
int *kick_directions;
int num_kicks;
double amplitude;

// Function to initialize random kicks
void initialize_kicks() {
    kick_times = (double *)malloc(num_kicks * sizeof(double));
    kick_directions = (int *)malloc(num_kicks * sizeof(int));
    srand(time(NULL));

    for (int i = 0; i < num_kicks; i++) {
        kick_times[i] = ((double)rand() / RAND_MAX) * T_MAX;
        kick_directions[i] = (rand() % 2) ? 1 : -1;
    }
}

// Function defining the nonlinear pendulum ODE with kicks
void pendulum_equation(double t, State s, double gamma, double *dtheta_dt, double *domega_dt) {
    *dtheta_dt = s.omega;
    *domega_dt = - (G / L) * sin(s.theta) - gamma * s.omega;

    for (int i = 0; i < num_kicks; i++) {
        if (fabs(t - kick_times[i]) < DT) {
            *domega_dt += kick_directions[i] * amplitude;
        }
    }
}

// Runge-Kutta 4th order integration
State rk4_step(State s, double t, double dt, double gamma) {
    State s_new;
    double k1_theta, k1_omega, k2_theta, k2_omega, k3_theta, k3_omega, k4_theta, k4_omega;

    pendulum_equation(t, s, gamma, &k1_theta, &k1_omega);
    State s2 = {s.theta + 0.5 * dt * k1_theta, s.omega + 0.5 * dt * k1_omega};

    pendulum_equation(t + 0.5 * dt, s2, gamma, &k2_theta, &k2_omega);
    State s3 = {s.theta + 0.5 * dt * k2_theta, s.omega + 0.5 * dt * k2_omega};

    pendulum_equation(t + 0.5 * dt, s3, gamma, &k3_theta, &k3_omega);
    State s4 = {s.theta + dt * k3_theta, s.omega + dt * k3_omega};

    pendulum_equation(t + dt, s4, gamma, &k4_theta, &k4_omega);

    s_new.theta = s.theta + (dt / 6.0) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
    s_new.omega = s.omega + (dt / 6.0) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);

    return s_new;
}

// Function to simulate the pendulum and return data to Python
void simulate_pendulum(double gamma, double amplitude_input, int num_kicks_input, double *t_out, double *theta_out, double *omega_out, int num_steps) {
    num_kicks = num_kicks_input;
    amplitude = amplitude_input;

    initialize_kicks();

    State s = { M_PI/1.0, 0.0 }; // Initial angle slightly off 180 degrees, initial angular velocity = 0
    double t = 0.0;

    for (int i = 0; i < num_steps; i++) {
        t_out[i] = t;
        theta_out[i] = s.theta;
        omega_out[i] = s.omega;

        s = rk4_step(s, t, DT, gamma);
        t += DT;
    }

    free(kick_times);
    free(kick_directions);
}
