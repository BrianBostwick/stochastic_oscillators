# stochastic_oscillators

- Code to simulate non-linear coupled oscillators with couplings to noise via ramdon kicks.
- To compile the C code use: gcc -shared -o libpendulum.so -fPIC -arch x86_64 pendulum_simulation.c -lm
- Same but to compile lattice code: gcc -shared -o liblatticependulum.so -fPIC -arch x86_64 lattice_pendulum_simulation.c -lm
