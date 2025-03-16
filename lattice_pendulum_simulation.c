#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 9.81     // Gravitational acceleration (m/s^2)
#define L 1.0      // Length of pendulum (m)
#define T_MAX 30.0 // Maximum simulation time
#define DT 0.01    // Time step
#define N 10       // Number of pendulums in the lattice
#define K 0.1      // Coupling strength

// Structure to store state variables
typedef struct {
    double theta;
    double omega;
} State;

// Global variables for stochastic kicks
double **kick_times;
int **kick_directions;
int num_kicks;
double amplitude;

// Function to initialize random kicks for each pendulum
void initialize_kicks() {
    kick_times = (double **)malloc(N * sizeof(double *));
    kick_directions = (int **)malloc(N * sizeof(int *));
    srand(time(NULL));

    for (int i = 0; i < N; i++) {
        kick_times[i] = (double *)malloc(num_kicks * sizeof(double));
        kick_directions[i] = (int *)malloc(num_kicks * sizeof(int));

        for (int j = 0; j < num_kicks; j++) {
            kick_times[i][j] = ((double)rand() / RAND_MAX) * T_MAX;
            kick_directions[i][j] = (rand() % 2) ? 1 : -1;
        }
    }
}

// Function defining the coupled nonlinear pendulum ODE with kicks
void pendulum_equation(double t, State lattice[], double gamma, double *dtheta_dt, double *domega_dt) {
    for (int i = 0; i < N; i++) {
        dtheta_dt[i] = lattice[i].omega;
        double coupling_term = 0.0;

        // Nearest-neighbor coupling
        if (i > 0) {
            coupling_term += lattice[i - 1].theta;
        }
        if (i < N - 1) {
            coupling_term += lattice[i + 1].theta;
        }
        coupling_term -= 2 * lattice[i].theta;

        domega_dt[i] = - (G / L) * sin(lattice[i].theta) - gamma * lattice[i].omega + K * coupling_term;

        // Apply stochastic kicks
        for (int j = 0; j < num_kicks; j++) {
            if (fabs(t - kick_times[i][j]) < DT) {
                domega_dt[i] += kick_directions[i][j] * amplitude;
            }
        }
    }
}

// Runge-Kutta 4th order integration
void rk4_step(State lattice[], double t, double dt, double gamma) {
    State temp_lattice[N];
    double k1_theta[N], k1_omega[N], k2_theta[N], k2_omega[N], k3_theta[N], k3_omega[N], k4_theta[N], k4_omega[N];

    pendulum_equation(t, lattice, gamma, k1_theta, k1_omega);

    for (int i = 0; i < N; i++) {
        temp_lattice[i].theta = lattice[i].theta + 0.5 * dt * k1_theta[i];
        temp_lattice[i].omega = lattice[i].omega + 0.5 * dt * k1_omega[i];
    }
    pendulum_equation(t + 0.5 * dt, temp_lattice, gamma, k2_theta, k2_omega);

    for (int i = 0; i < N; i++) {
        temp_lattice[i].theta = lattice[i].theta + 0.5 * dt * k2_theta[i];
        temp_lattice[i].omega = lattice[i].omega + 0.5 * dt * k2_omega[i];
    }
    pendulum_equation(t + 0.5 * dt, temp_lattice, gamma, k3_theta, k3_omega);

    for (int i = 0; i < N; i++) {
        temp_lattice[i].theta = lattice[i].theta + dt * k3_theta[i];
        temp_lattice[i].omega = lattice[i].omega + dt * k3_omega[i];
    }
    pendulum_equation(t + dt, temp_lattice, gamma, k4_theta, k4_omega);

    for (int i = 0; i < N; i++) {
        lattice[i].theta += (dt / 6.0) * (k1_theta[i] + 2 * k2_theta[i] + 2 * k3_theta[i] + k4_theta[i]);
        lattice[i].omega += (dt / 6.0) * (k1_omega[i] + 2 * k2_omega[i] + 2 * k3_omega[i] + k4_omega[i]);
    }
}

// Function to simulate the lattice of pendulums and return data to Python
void simulate_pendulums(double gamma, double amplitude_input, int num_kicks_input, double *t_out, double *theta_out, double *omega_out, int num_steps) {
    num_kicks = num_kicks_input;
    amplitude = amplitude_input;

    initialize_kicks();

    State lattice[N];
    for (int i = 0; i < N; i++) {
        lattice[i].theta = M_PI / 1.0;  // Initial angle
        lattice[i].omega = 0.0;         // Initial angular velocity
    }

    double t = 0.0;

    for (int i = 0; i < num_steps; i++) {
        t_out[i] = t;

        for (int j = 0; j < N; j++) {
            theta_out[i * N + j] = lattice[j].theta;
            omega_out[i * N + j] = lattice[j].omega;
        }

        rk4_step(lattice, t, DT, gamma);
        t += DT;
    }

    for (int i = 0; i < N; i++) {
        free(kick_times[i]);
        free(kick_directions[i]);
    }
    free(kick_times);
    free(kick_directions);
}
