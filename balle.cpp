// balle - Program to compute the trajectory of a baseball
//         using the Euler method.
#include "NumMeth.h"
#include <cmath>
using namespace std;

// Need to calculate cross product: https://www.tutorialspoint.com/cplusplus-program-to-compute-cross-product-of-two-vectors
void crossProduct(double v_A[], double v_B[], double c_P[]) {
   c_P[1] = v_A[2] * v_B[3] - v_A[3] * v_B[2];
   c_P[2] = -(v_A[1] * v_B[3] - v_A[3] * v_B[1]);
   c_P[3] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
}

void magnitude(double v[], double mag){
  mag = sqrt((v[1]*v[1]) + (v[2]*v[2] + v[3]*v[3]));
}

int main() {

  //* Set initial position and velocity of the baseball
  double z1, x1, y1, speed, theta;
  double r1[3+1], v1[3+1], r[3+1], v[3+1], accel[3+1], w[3+1]; 
  // cout << "Enter initial height (meters): "; cin >> y1;
  x1 = 0;
  y1 = 0.5;
  z1 = 0.8;
  r1[1] = x1;  r1[2] = y1;     r1[3]=z1;// Initial vector position
  // cout << "Enter initial speed (m/s): "; cin >> speed; 
  speed = 5;
  // cout << "Enter initial angle (degrees): "; cin >> theta;
  theta = -20;
  const double pi = 3.141592654; 
  v1[1] = speed*cos(theta*pi/180);   // Initial velocity (x)
  v1[2] = 0;   // Initial velocity (y)
  v1[3] = speed*sin(theta*pi/180);
  r[1] = r1[1];  r[2] = r1[2];  r[3]=r1[3]; // Set initial position and velocity
  v[1] = v1[1];  v[2] = v1[2];  v[3]=v1[3];          
  // set initial spin 
  w[1] = 20000; w[2] = -2000; w[3]=0;

  //* Set physical parameters (mass, Cd, etc.)
  double Cd = 0.155;      // Drag coefficient (dimensionless)
  double radius = 0.02;  // Ball radius in meters
  double area = pi*radius*radius;  // Cross-sectional area of projectile (m^2)
  double grav = 9.81;    // Gravitational acceleration (m/s^2)
  double mass = 0.145;   // Mass of projectile (kg)
  double airFlag, rho;
  double s0;//spin coefficient
  double alpha = 2/3;

  airFlag = 1;
  if( airFlag == 0 )
    rho = 0;      // No air resistance
  else
    rho = 1.2;    // Density of air (kg/m^3)
  double air_const = -0.5*Cd*rho*area/mass;  // Air resistance constant

  s0 = -0.7 * pi * pow(radius, 3)*rho; 
  double magn_v; 
  magnitude(v, magn_v);
  Cd = .4127 * radius/magn_v;
  cout << s0 << " ";

  //* Loop until ball hits table or max steps completed
  double tau = 0.001;
  int iStep, maxStep = 1000;   // Maximum number of steps
  double *xplot, *yplot, *xNoAir, *yNoAir, *zplot, *zNoAir;
  xplot = new double [maxStep+1];  yplot = new double [maxStep+1]; zplot = new double [maxStep+1];
  xNoAir = new double [maxStep+1]; yNoAir = new double [maxStep+1]; zNoAir = new double [maxStep+1];
  for( iStep=1; iStep<=maxStep; iStep++ ) {

    //* Record position (computed and theoretical) for plotting
    xplot[iStep] = r[1];   // Record trajectory for plot
    yplot[iStep] = r[2];
    zplot[iStep] = r[3];
    double t = (iStep-1)*tau;     // Current time
    xNoAir[iStep] = r1[1] + v1[1]*t;
    yNoAir[iStep] = r1[2] + v1[2]*t;
    zNoAir[iStep] = r1[3] + v1[3]*t - 0.5*grav*t*t;
  
    //* Calculate the acceleration of the ball 
    double normV = sqrt( v[1]*v[1] + v[2]*v[2] );
    double c_P[3+1]; 
    crossProduct(w, v, c_P);
    // cout << c_P[1] << c_P[2] << c_P[3] <<" ";
    

    accel[1] = air_const*normV*v[1];   // Air resistance
    accel[1] += s0/mass * c_P[1]; // Magnus Force
    accel[2] = air_const*normV*v[2];   // Air resistance
    accel[2] += s0/mass * c_P[2]; // Magnus Force
    accel[3] = air_const*normV*v[3];   // Air resistance
    accel[3] += s0/mass * c_P[3]; // Magnus Force
    accel[3] -= grav;                  // Gravity
  
    //* Calculate the new position and velocity using Euler method
    r[1] += tau*v[1];                 // Euler step
    r[2] += tau*v[2]; 
    r[3] += tau*v[3];            
    v[1] += tau*accel[1];     
    v[2] += tau*accel[2];     
    v[3] += tau*accel[3];
  
    //* If ball reaches table (y<0.75), use conservation of energy to change direction
    if( r[3] < 0.75 )  {
      double v2[3+1], w2[3+1]; 
      v2[1] = v[1]*(1-alpha)/(alpha + 1) - alpha*(1+1)/(alpha + 1)*radius*w[2];
      v2[2] = v[2]*(1-alpha)/(alpha + 1) + alpha/(alpha + 1)*radius*w[1];
      v2[3] = -v[3];
      w2[1] = v[3]*(1+1)/(radius*(alpha + 1)) + (alpha -1 )/(alpha + 1)*radius*w[1];
      w2[2] = -v[1]*(1+1)/(radius*(alpha + 1)) + (alpha-1)/(alpha + 1)*radius*w[2];
      w2[3] = w[3];
      v[1] = v2[1]; v[2] = v2[2]; v[3]=v2[3];
      w[1] = w2[1]; w[2] = w2[2]; w[3]=w2[3];
      cout<<"v1"<<v[3]<< endl;
      cout<< "v2"<<v2[3]<< endl;

    } 

    //* If ball reaches end of the table (x>2.7), break
    if( r[1] > 2.7 )  {
      xplot[iStep+1] = r[1];  // Record last values computed
	    yplot[iStep+1] = r[2];
      zplot[iStep+1] = r[3];
      break;                  // Break out of the for loop
    } 
  }

  //* Print maximum range and time of flight
  cout << "Time of flight is " << iStep*tau << " seconds" << endl;

  //* Print out the plotting variables: 
  //    xplot, yplot, xNoAir, yNoAir
  ofstream xplotOut("xplot.txt"), yplotOut("yplot.txt"), zplotOut("zplot.txt"),
	       xNoAirOut("xNoAir.txt"), yNoAirOut("yNoAir.txt"), zNoAirOut("zNoAir.txt") ;
  int i;
  for( i=1; i<=iStep+1; i++ ) {
    xplotOut << xplot[i] << endl;
    yplotOut << yplot[i] << endl;
    zplotOut << zplot[i] << endl;
  }
  for( i=1; i<=iStep; i++ ) {
    xNoAirOut << xNoAir[i] << endl;
    yNoAirOut << yNoAir[i] << endl;
    zNoAirOut << zNoAir[i] << endl;
  }

  delete [] xplot; // Release memory
  delete [] yplot;
  delete [] zplot;
  delete [] xNoAir;
  delete [] yNoAir;
  delete [] zNoAir;
  return 0;

}