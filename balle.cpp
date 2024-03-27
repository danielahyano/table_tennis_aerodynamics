// balle - Program to compute the trajectory of a baseball
//         using the Euler method.
#include "NumMeth.h"
using namespace std;

int main() {

  //* Set initial position and velocity of the baseball
  double y1, speed, theta;
  double r1[2+1], v1[2+1], r[2+1], v[2+1], accel[2+1]; 
  cout << "Enter initial height (meters): "; cin >> y1;
  r1[1] = 0;  r1[2] = y1;     // Initial vector position
  cout << "Enter initial speed (m/s): "; cin >> speed; 
  cout << "Enter initial angle (degrees): "; cin >> theta;
  const double pi = 3.141592654; 
  v1[1] = speed*cos(theta*pi/180);   // Initial velocity (x)
  v1[2] = speed*sin(theta*pi/180);   // Initial velocity (y)
  r[1] = r1[1];  r[2] = r1[2];  // Set initial position and velocity
  v[1] = v1[1];  v[2] = v1[2];             

  //* Set physical parameters (mass, Cd, etc.)
  double Cd = 0.35;      // Drag coefficient (dimensionless)
  double area = 4.3e-3;  // Cross-sectional area of projectile (m^2)
  double grav = 9.81;    // Gravitational acceleration (m/s^2)
  double mass = 0.145;   // Mass of projectile (kg)
  double airFlag, rho;
  airFlag = 1;
  if( airFlag == 0 )
    rho = 0;      // No air resistance
  else
    rho = 1.2;    // Density of air (kg/m^3)
  double air_const = -0.5*Cd*rho*area/mass;  // Air resistance constant

  //* Loop until ball hits table or max steps completed
  double tau = 0.001;
  int iStep, maxStep = 1000;   // Maximum number of steps
  double *xplot, *yplot, *xNoAir, *yNoAir;
  xplot = new double [maxStep+1];  yplot = new double [maxStep+1];
  xNoAir = new double [maxStep+1]; yNoAir = new double [maxStep+1];
  for( iStep=1; iStep<=maxStep; iStep++ ) {

    //* Record position (computed and theoretical) for plotting
    xplot[iStep] = r[1];   // Record trajectory for plot
    yplot[iStep] = r[2];
    double t = (iStep-1)*tau;     // Current time
    xNoAir[iStep] = r1[1] + v1[1]*t;
    yNoAir[iStep] = r1[2] + v1[2]*t - 0.5*grav*t*t;
  
    //* Calculate the acceleration of the ball 
    double normV = sqrt( v[1]*v[1] + v[2]*v[2] );
    accel[1] = air_const*normV*v[1];   // Air resistance
    accel[2] = air_const*normV*v[2];   // Air resistance
    accel[2] -= grav;                  // Gravity
  
    //* Calculate the new position and velocity using Euler method
    r[1] += tau*v[1];                 // Euler step
    r[2] += tau*v[2];                 
    v[1] += tau*accel[1];     
    v[2] += tau*accel[2];     
  
    //* If ball reaches table (y<0.75), use conservation of energy to change direction
    if( r[2] < 0.75 )  {
      v[2] = -v[2];             
    } 

    //* If ball reaches end of the table (x>2.7), break
    if( r[1] > 2.7 )  {
      xplot[iStep+1] = r[1];  // Record last values computed
	    yplot[iStep+1] = r[2];
      break;                  // Break out of the for loop
    } 
  }

  //* Print maximum range and time of flight
  cout << "Maximum range is " << r[1] << " meters" << endl;
  cout << "Time of flight is " << iStep*tau << " seconds" << endl;

  //* Print out the plotting variables: 
  //    xplot, yplot, xNoAir, yNoAir
  ofstream xplotOut("xplot.txt"), yplotOut("yplot.txt"), 
	       xNoAirOut("xNoAir.txt"), yNoAirOut("yNoAir.txt");
  int i;
  for( i=1; i<=iStep+1; i++ ) {
    xplotOut << xplot[i] << endl;
    yplotOut << yplot[i] << endl;
  }
  for( i=1; i<=iStep; i++ ) {
    xNoAirOut << xNoAir[i] << endl;
    yNoAirOut << yNoAir[i] << endl;
  }

  delete []  xplot, yplot, xNoAir, yNoAir; // Release memory

  return 0;

}