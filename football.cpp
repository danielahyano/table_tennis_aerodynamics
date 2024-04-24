// balle - Program to compute the trajectory of a baseball
//         using the Euler method.
#include "NumMeth.h"
#include <cmath>
using namespace std;

// Fotball ball, to compare results'

// Need to calculate cross product: https://www.tutorialspoint.com/cplusplus-program-to-compute-cross-product-of-two-vectors
void crossProduct(double v_A[], double v_B[], double c_P[]) {
   c_P[1] = v_A[2] * v_B[3] - v_A[3] * v_B[2];
   c_P[2] = -(v_A[1] * v_B[3] - v_A[3] * v_B[1]);
   c_P[3] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
}

double magnitude(double v[]){
  double mag = sqrt((v[1]*v[1]) + (v[2]*v[2]) + (v[3]*v[3]));
  return mag;
}

int main() {

  //* Set initial position and velocity of the ball
  double z1, x1, y1, speed, theta;
  double r1[3+1], v1[3+1], r[3+1], v[3+1], accel[3+1], w[3+1]; 
  // cout << "Enter initial height (meters): "; cin >> y1;
  x1 = 2.5;
  y1 = 0;
  z1 = 1.5;
  r1[1] = x1;  r1[2] = y1;     r1[3]=z1;// Initial vector position

  v1[1] = -5.5;   // Initial velocity (x)
  v1[2] = 0;   // Initial velocity (y)
  v1[3] = -5;  // Initial velocity (z)
  r[1] = r1[1];  r[2] = r1[2];  r[3]=r1[3]; // Set initial position and velocity
  v[1] = v1[1];  v[2] = v1[2];  v[3]=v1[3];          
  // set initial spin 
  w[1] = 0; w[2] = 0; w[3]= 0;

  //* Set physical parameters (mass, Cd, etc.)
  const double pi = 3.141592654; 
  double Cd = 0.15;      // Drag coefficient (dimensionless)
  double Cl = 0.9;      // Spin coeffiecient
  double radius = 0.11;  // Ball radius in meters
  double area = pi*radius*radius;  // Cross-sectional area of projectile (m^2)
  double grav = 9.81;    // Gravitational acceleration (m/s^2)
  double mass = 0.430;   // Mass of projectile (kg)
  double airFlag, rho;
  double s0;//spin coefficient
  double alpha = 2/3;
  double vc=12.19; 
  double vs=1.309;

  airFlag = 0;
  if( airFlag == 0 )
    rho = 0;      // No air resistance
  else
    rho = 1.2;    // Density of air (kg/m^3)
 
  // update Cd value
  double magn_w; 
  double magn_v1; 
  double magn_v; 
  magn_v1 = magnitude(v1);
  magn_w = magnitude(w);
  double S = magn_w/(radius*magn_v1);
  cout << "magn_v1 " <<  magn_v1 << endl;
  cout << "w[2] " <<  w[2] << endl;
  cout << "magn_w " << magn_w << endl;
  cout << "S " << S << endl;
  if (S > 0.05){
    Cd = .4127 * pow(S,0.3056);
  }
  else{
    Cd = 0.155 + 0.36/(exp((magn_v1-vc)/vs));
  }

  bool bounce = true; // whether the ball kicks or not.
  int bounceCount = 0;
  int bounceTime = 0;


  //* Loop until ball hits floor or max steps completed
  double tau = 0.001;
  int iStep, maxStep = 10000;   // Maximum number of steps
  double *xplot, *yplot, *xNoAir, *yNoAir, *zplot, *zNoAir;
  double *vxplot, *vyplot,*vzplot; 
  double *wxplot, *wyplot,*wzplot; 
  xplot = new double [maxStep+1];  yplot = new double [maxStep+1]; zplot = new double [maxStep+1];
  vxplot = new double [maxStep+1];  vyplot = new double [maxStep+1]; vzplot = new double [maxStep+1];
  wxplot = new double [maxStep+1];  wyplot = new double [maxStep+1]; wzplot = new double [maxStep+1];
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
    vxplot[iStep] = v[1];   // Record velocity for plot
    vyplot[iStep] = v[2];
    vzplot[iStep] = v[3];
    wxplot[iStep] = w[1];   // Record angular velocity for plot
    wyplot[iStep] = w[2];
    wzplot[iStep] = w[3];
  
    //* Calculate the acceleration of the ball 
    double normV;
    normV = magnitude(v);
    double c_P[3+1]; 
    crossProduct(w, v, c_P);
    S = magn_w * radius / normV;
    
    if (S > 0.05 && normV>vc){
      Cd = .4127 * pow(S,0.3056);
    }  
    else{
      Cd = 0.155 + 0.36/(exp((magn_v1-vc)/vs));
    }


    double coeff_drag = -(Cd*area*rho)/(2*mass);
    double coeff_mag = Cl * area * radius * rho/mass; 
    accel[1] = coeff_drag*normV*v[1];   // Air resistance
    accel[1] += coeff_mag * (w[3]*v[2] - w[2]*v[3]); // Magnus Force
    // cout << "accel x" << accel[1] << endl;
    accel[2] = coeff_drag*normV*v[2];   // Air resistance
    accel[2] += coeff_mag * (w[1]*v[3] - w[3]*v[1]); // Magnus Force
    // cout << "accel y" << accel[2] << endl;
    accel[3] = coeff_drag*normV*v[3];   // Air resistance
    accel[3] += coeff_mag * (w[2]*v[1] - w[1]*v[2]); // Magnus Force
    accel[3] -= grav;                  // Gravity
    // cout << "accel z" << accel[3] << endl;
  
    //* Calculate the new position and velocity using Euler midpoint method
    r[1] += (tau/2)*(2*v[1] + tau*accel[1]);                 // Euler midpoint step
    r[2] += (tau/2)*(2*v[2] + tau*accel[2]); 
    r[3] += (tau/2)*(2*v[3] + tau*accel[3]);             
    v[1] += tau*accel[1];     
    v[2] += tau*accel[2];     
    v[3] += tau*accel[3];
  
    //* If ball reaches floor, use conservation of energy to change direction
    if( r[3] < 0 && bounce)  {
      
      double diff = iStep - bounceTime;
      cout << "diff " << diff <<  endl;
      if(diff  > 10){
        double v2[3+1], w2[3+1]; 
        bounceTime = iStep;
        bounceCount = bounceCount + 1;
        // here I am assuming total conservation of energy.
        double ex = -0.75; //-1 corresponds to frictionless gliding
        double ey = -0.75; //-1 corresponds to frictionless gliding
        double ez = 0.57;
        v2[1] = v[1]*(1-alpha*ex)/(alpha + 1) - alpha*(1+ex)/(alpha + 1)*radius*w[2];
        v2[2] = v[2]*(1-alpha*ey)/(alpha + 1) + alpha*(1+ey)/(alpha + 1)*radius*w[1];
        v2[3] = -ez*v[3];
        w2[1] = v[2]*(1 + ey)/(radius*(alpha + 1)) + (alpha - ey )/(alpha + 1)*w[1];
        w2[2] = -v[1]*(1 + ex)/(radius*(alpha + 1)) + (alpha - ex)/(alpha + 1)*w[2];
        w2[3] = w[3];
        v[1] = v2[1]; v[2] = v2[2]; v[3]=v2[3];
        w[1] = w2[1]; w[2] = w2[2]; w[3]=w2[3];
        cout << "vx: "<<v2[1] << endl;
        cout << "vz: "<<v2[3] << endl;
        cout << "z: "<<r[3] << endl;
      }
    } 

    //* If ball reaches end of the table (x>3), break
    if( bounceCount > 5)  {
      cout << "END "<< bounceCount << endl;
      xplot[iStep+1] = r[1];  // Record last values computed
	    yplot[iStep+1] = r[2];
      zplot[iStep+1] = r[3];
      vxplot[iStep+1] = v[1]; 
	    vyplot[iStep+1] = v[2];
      vzplot[iStep+1] = v[3];
      break;                  // Break out of the for loop
    } 
  }

  //* Print maximum range and time of flight
  cout << "Time of flight is " << iStep*tau << " seconds" << endl;

  //* Print out the plotting variables: 
  //    xplot, yplot, xNoAir, yNoAir
  ofstream xplotOut("xplot.txt"), yplotOut("yplot.txt"), zplotOut("zplot.txt"),
	       xNoAirOut("xNoAir.txt"), yNoAirOut("yNoAir.txt"), zNoAirOut("zNoAir.txt"),
         vxplotOut("vxplot.txt"), vyplotOut("vyplot.txt"), vzplotOut("vzplot.txt"),
         wxplotOut("wxplot.txt"), wyplotOut("wyplot.txt"), wzplotOut("wzplot.txt");
  int i;
  for( i=1; i<=iStep+1; i++ ) {
    xplotOut << xplot[i] << endl;
    yplotOut << yplot[i] << endl;
    zplotOut << zplot[i] << endl;
    vxplotOut << vxplot[i] << endl;
    vyplotOut << vyplot[i] << endl;
    vzplotOut << vzplot[i] << endl;
    wxplotOut << wxplot[i] << endl;
    wyplotOut << wyplot[i] << endl;
    wzplotOut << wzplot[i] << endl;
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
  delete [] vxplot; 
  delete [] vyplot;
  delete [] vzplot;
  delete [] wxplot; 
  delete [] wyplot;
  delete [] wzplot;
  return 0;

}