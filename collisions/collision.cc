#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

#define N 16 // number of particles
#define X 100 // number of bins for histogram

struct particle // defines a particle and the properties it should have
{
  double x ; // position
  double v ; // velocity
  double im ; // inverse mass
  double a ; // radius of particle
  double T ; // kinetic energy
  double Ttime; // kinetic energy x time
  double p ; // momentum
};

struct result // used as an output
{
  int arecollisionshappening ; // 0 if collisions aren't happening, 1 if they are
  int collisionnumber ; // give the pair of particles that will collide
  double time ; // gives the time before the next collision
};

void showState (particle particles[], int num) // outputs all properties of a specified particle, used in debugging
{
  cout << "particle " << num << "\t";
  cout << "x " << particles[num].x << "\t";
  cout << "v " << particles[num].v << "\t";
  cout << "im " << particles[num].im << "\t";
  cout << "a " << particles[num].a << "\t";
  cout << "T " << particles[num].T << "\t";
  cout << "p " << particles[num].p << endl;
}

void outputtotalEnergy(particle particles[], double time) // outputs the total kinetic energy of a system
{
  double totalenergy = 0;

  for (int i = 0; i < N; i++)
  {
    totalenergy += particles[i].T;
  }

  cout << time << "\t" << totalenergy << endl;
}

void showSystemPositions (particle particles[], double time) // outputs the time and positions of all particles, used to plot particles
{
  cout << time << "\t";
  for (int i = 0; i < N; i++)
  {
    cout << particles[i].x << "\t";
  }
  cout << endl;
}

void updateenergyandmomentum (particle particles[], int num) // updates the energy and momentum of a specified particle pair
{
    if (particles[num].im == 0) // if inverse mass is zero, sets T and p to zero to avoid breaking everything
    {
      particles[num].T = 0;
      particles[num].p = 0;
    }
    if (particles[num].im != 0) // computes values if non zero
    {
      particles[num].T = 0.5 * (1/ particles[num].im) * pow(particles[num].v, 2);
      particles[num].p = (1/particles[num].im) * particles[num].v;
    }
}

void generateParticles (particle particles[]) // used in debugging to specify individual particles
{

  // It's a bit boring to specify them individually it turns out.

  particles[0].x = -10;
  particles[0].v = 1;
  particles[0].im = 0.7;
  particles[0].a = 1;

  particles[1].x = 10;
  particles[1].v = -1;
  particles[1].im = 1;
  particles[1].a = 2;


  particles[2].x = 15;
  particles[2].v = 0;
  particles[2].im = 0;
  particles[2].a = 3;

/*
  particles[3].x = -0;
  particles[3].v = 2;
  particles[3].im = 1;
  particles[3].a = 0;

  particles[4].x = 5;
  particles[4].v = 3;
  particles[4].im = 3;
  particles[4].a = 0;

  particles[5].x = 8;
  particles[5].v = 0.5;
  particles[5].im = 0.1;
  particles[5].a = 0;

  particles[6].x = 10;
  particles[6].v = 0;
  particles[6].im = 0;
  particles[6].a = 0;
*/

  for (int i = 0; i < N; i++) {
    updateenergyandmomentum(particles, i);
  }

}

void generateRandomParticles (particle particles[], double maximum_speed, double maximum_mass) // generates particles with random masses and velocities
{
  srand(time(0)); // generates a random seed
  double randomdoubles[N], randominverses[N]; // initiates arrays which are used later

/*
Positions particles uniformly distributed between 0 and 1, and gives them a radius of zero
it is important that they are positioned in order, otherwise the rest of the program will break.
However, it turns out that positioning them in order isn't too difficult.
*/
  for (int i = 0; i < N; i++)
  {
    particles[i].x = (1.0 * i) / (1.0 * (N-1));
    particles[i].a = 0;
    particles[i].Ttime = 0;
  }


  for (int i = 0; i < N; i++) // gives particles a random speed, uniformly distributed between defined values
  {
    randomdoubles[i] = 2 * maximum_speed * (((1.0 * rand()) / RAND_MAX) - 0.5);
    particles[i].v = randomdoubles[i];
  }

  for (int i = 0; i < N; i++) // gives particles random masses, again uniformly distributed
  {
    randomdoubles[i] = 0.5 * maximum_mass * (((1.0 * rand()) / RAND_MAX) + 0.1);
    randominverses[i] = (1 / randomdoubles[i]);
    particles[i].im = randominverses[i];

  }

  for (int i = 0; i < N; i++) // for all the particles that have been produced, given them a Ekin and momentum
  {
    updateenergyandmomentum(particles, i);
  }

}

void generateStructure (particle particles[], int structure, double maximum_speed, double maximum_mass) // generates specific features, for example walls and a piston
{
  if (structure == 1 || structure == 2) // generate walls of infinite mass, they can easily be made to move
  {
    particles[0].v = 0;
    particles[0].im = 0;
    particles[N-1].v = 0;
    particles[N-1].im = 0;
  }

  if (structure == 2) // generate piston halfway between walls, again, can be give starting values etc
  {
    int halfway = int(N/2);
    double mass = 50 * maximum_mass;
    particles[halfway].im = 1 / mass;
    particles[halfway].v = maximum_speed;
    particles[halfway].a = 0.05;
  }
}

bool dotheycollide (particle particles[], int n1, int n2) // determines whether two particles collide
{
  double deltav = particles[n1].v - particles[n2].v ; // they collide if the left hand particle has a more positive velocity
  if (deltav > 0) {
    return true;
  }
  else {
    return false;
  }
}

double timetoparticlecollision (particle particles[], int n1, int n2) // determines the time between two particle collisions
{
  double deltax = particles[n2].x - particles[n1].x - particles[n1].a - particles[n2].a;
  double deltav = particles[n1].v - particles[n2].v ;
  double deltat = deltax / deltav;
  return deltat;
}

result particlecollisions (particle particles[]) // outputs which, if any, particles collide and the time that takes
{
  result results;

  // initialises variables
  double deltat, deltatmin = 1000;
  int numberofcollisions = 0, particlenumber;

  for (int k = 0; k < N - 1; k++) // loops over particles
  {
    if (dotheycollide(particles, k, k + 1) == true) // if the particles collide
    {
      deltat = timetoparticlecollision(particles, k, k+1); // determines the time that collision takes
      numberofcollisions += 1;
      if (deltat <= deltatmin) // and if this is the smallest time so far, says that they will be the next pair to collide
      {
        deltatmin = deltat;
        particlenumber = k;
      }
    }
  }
  if (numberofcollisions == 0) // if there are no more collisions
  {
    results.arecollisionshappening = 0;
    results.collisionnumber = 0;
    results.time = 0;
  }
  if (numberofcollisions != 0) // if there are still collisions we can have some fun
  {
    results.arecollisionshappening = 1;
    results.collisionnumber = particlenumber;
    results.time = deltatmin;
  }
  return results;
}

void collide (particle particles[], result results) // collides a pair of particles, and updates their momenta and kinetic energies
{
  int n1 = results.collisionnumber, n2 = results.collisionnumber + 1; // to keep later code concise
  double v1init = particles[n1].v ;
  double v2init = particles[n2].v ;

  // collision equations, determined by hand
  particles[n1].v = ((particles[n2].im - particles[n1].im)/(particles[n1].im + particles[n2].im)) * v1init + ((2 * particles[n1].im)/(particles[n1].im + particles[n2].im)) * v2init;
  particles[n2].v = ((2 * particles[n2].im)/(particles[n1].im + particles[n2].im)) * v1init + ((particles[n1].im - particles[n2].im)/(particles[n1].im + particles[n2].im)) * v2init;

  // since these are the only particles to have done anything in this step, it is only necessary to update their energies and momenta
  updateenergyandmomentum(particles, n1);
  updateenergyandmomentum(particles, n2);
}

void moveparticles (particle particles[], result results) // moves particles to position of next collision
{

  for (int i = 0; i < N; i++) // moves all the particles by their speed x time before next collision
  {
    particles[i].x += particles[i].v * results.time;
    particles[i].Ttime += particles[i].T * results.time;
  }
}

void printmovingparticles (particle particles[], result results, double timestep) // an attempt to dynamically print particles positions so that simulation could occur, given up upon
{
  particle movingparticles[N];
  for (int i = 0; i < N; i++)
  {
    movingparticles[i] = particles[i]; // copies the array of particles to the new array
  }

  double time = 0;
  while (time <= results.time)
  {
    for (int i = 0; i < N; i++)
    {
      cout << movingparticles[i].x << "\t" << 0 << "\t" << movingparticles[i].a << "\t";
      movingparticles[i].x += movingparticles[i].v * timestep;
    }
    cout << endl;
    time += timestep;
  }

}

void histogram(particle particles[]) // used to plot a histogram of particle energies
{
  double Tav = 0;

  for (int i = 0; i < N; i++)
  {
    Tav += particles[i].T;
  }

  Tav = Tav / (1.0 * N);

  double binborders[X];
  for (int i = 0; i < X; i++)
  {
    binborders[i] = Tav * (1.0 * i) / (1.0 * X);
  }
  double binwidth = (1.0 / (1.0 * X));

  double binvalues[X];
  for (int i = 0; i < X; i++)
  {
    binvalues[i] = 0;
    for (int k = 0; k < N; k++)
    {
      if (particles[k].T > binborders[i] & particles[k].T < binborders[i] + binwidth)
      {
        binvalues[i] += 1;
      }
    }
  }

  for (int i = 0; i < X; i++)
  {
    cout << binborders[i] << "\t" << binvalues[i] << endl;
  }

}

void outputenergies(particle particles[], int output, double time) // outputs either average Ekin, instantaneous Ekin, or a histogram of instantaneous Ekins
{
  if (output == 1)
  {
    double Taverage[N];
    double particlemass[N];
    for (int i = 0; i < N; i++)
    {
      Taverage[i] = particles[i].Ttime / time ;

      if (particles[i].im == 0)
      {
        particlemass[i] = 0;
      }
      if (particles[i].im != 0)
      {
        particlemass[i] = 1 / particles[i].im ;
      }
      cout << particlemass[i] << "\t" << Taverage[i] << endl;
    }
  }

  if (output == 2)
  {
    double particlemass[N];
    for (int i = 0; i < N; i++)
    {
      if (particles[i].im == 0)
      {
        particlemass[i] = 0;
      }
      if (particles[i].im != 0)
      {
        particlemass[i] = 1 / particles[i].im ;
      }
      cout << particlemass[i] << "\t" << particles[i].T << endl;
    }
  }

  if (output ==3)
  {
    histogram(particles);
  }

}

int main()
{
  double time = 0; // initiates time
  double maximum_mass = 1; // defines the maximum mass (not inverse mass) of a particle
  double maximum_speed = 10; // defines the maximum speed of a particle

  int structure = 1; // 0 -> nothing, 1 -> walls, 2 -> walls and piston
  int output = 0; // 0 -> visualise collisions, 1 -> Taverage against particle mass, 2 -> T at an instant against particle mass, 3 -> T at an instant histogram, 4 -> Kinetic energy against time

  int numberofcollisions = 0;
  int maxnumberofcollisions = 1000; // governs after what number of collisions the simulation ends

  // initiates structures
  result results;
  particle particles[N];

  generateRandomParticles(particles, maximum_speed, maximum_mass); // generates an array of random particles
  generateStructure(particles, structure, maximum_speed, maximum_mass); // and, if required, adds walls or a piston
  results = particlecollisions(particles); // necessary to determine whether any particles collide initially


  while (numberofcollisions <= maxnumberofcollisions) // main loop
  {
    if (output == 0) // which you choose initially
    {
      showSystemPositions(particles, time); // outputting time and x coordinates
    }
    if (output == 4)
    {
      outputtotalEnergy(particles, time); // outputs Total Ekin vs time, used for moving walls
    }

    results = particlecollisions(particles);

    if (results.arecollisionshappening == 0) // breaks the loop if no more collisions happen
    {
      break;
    }
    if (results.arecollisionshappening == 1) // performs the collision
    {
      moveparticles(particles, results);
      collide(particles, results);
      time += results.time;
      numberofcollisions += 1;
    }
  }

  outputenergies(particles, output, time); // used if other outputs are desired, contains logic in function

}
