#include "ising.h"
using namespace std;







int main()
{   

    Ising prueba;
    bool  spins[SIZE]; //Stores the spins
    
	int i,j; //Counters

	int N; //Number of averages done into the system

	double h[5]; //Values of the exp(-J/kT)

	double energy; //Value of the energy of the system

    double m[DATA]; //This array contains several moments of the magnetization and energy

	int neighs[SIZE][4]; //To store nearest neighbours

    mt19937 gen(958431198); //Mersenne Twister RNG
	uniform_int_distribution<int> brandom(0, 1); //Get any random integer
	uniform_int_distribution<int> ran_pos(0, SIZE-1); //Get any random integer
	uniform_real_distribution<double> ran_u(0.0, 1.0); //Our uniform variable generator

	double tstar; //Control parameter
	double deltat, deltat_crit; //Change in control parameter by iteration
	double tmax, tmin, tcrit_up, tcrit_down; //Max and min temperature, and interval where we apply Wolff

    tmax = 5.0;
    tmin = 0.1;
    tcrit_up = 2.4;
    tcrit_down = 2.2;

    deltat = 0.1;
    deltat_crit = 0.01;

    ofstream output; //Output of the stream

    prueba.initialize(spins, gen, brandom); //Init randomly
    prueba.get_neighbors(neighs); //Get neighbour table
	energy = prueba.get_energy(spins, neighs); //Compute initial energy

    N = 1e5;

	output.open("test_alt1281.txt");
	for (tstar = tmax; tstar > tcrit_up; tstar -= deltat)
	{
	    prueba.do_step(spins, neighs, tstar, N, h, energy, gen, ran_u, ran_pos, m);
	    prueba.w_output(output, tstar, N, m);
	    cout << tstar << endl;
	}
    for (tstar = tcrit_up - deltat_crit; tstar > tcrit_down; tstar -= deltat_crit)
	{
	    prueba.do_step_wolff(spins, neighs, tstar, N, energy, gen, ran_u, ran_pos, m);
	    prueba.w_output(output, tstar, N, m);
	    cout << tstar << endl;
        
        }
    for (tstar = tcrit_down - deltat; tstar >= tmin; tstar -= deltat)
	{
	    prueba.do_step(spins, neighs, tstar, N, h, energy, gen, ran_u, ran_pos, m);
	    prueba.w_output(output, tstar, N, m);
	    cout << tstar << endl;
        
	}
	output.close();

	return 0;

}


