#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include <time.h>
//Constants
 double rshaft = 0.15; // radius of the cylindrical shaft
 int revrate = 377; // angular velocity of the cylindrical shaft
 double nu = 0.0000148; // kinematic viscosity of fluid
 double rho = 1.225; // density of fluid
 double pi = 3.141592653589793;
 #define DERIV(u, ul, ur) (nu * ((ur - 2.0 * (u) + ul) / (dy * dy))) // Defining the derivative function f = nu * (d2u/dy2) for Runge-Kutta later.
 double timespace[1000][16] = {0}; // Timespace Matrix Initialization 
 double torque_data_storage[1000] = {0}; // Matrix Initialization for the storage (and later output) of torque on the bar across time
 double torque_analytical_storage[1000] = {0}; // Matrix Initialization for the storage (and later output) of analytical torque values on the bar across time

int main()
{
 double torque = 0; // Variable for tracking the torque applied to the shaft
 double output = 0; // function value
 double u_surface = (rshaft * revrate); // function value at shaft surface
 double dy = 0.0107; // spacial step (optimal = 0.0107)
 double dt = 0.001; // time step (optimal = 0.130)
 char method [3]; // storage for the method choice of the operator
 int t = 0; // initial value of time
 FILE *fp = NULL;

 
 printf("Please Select a Method of Discretization,\nAB for Adams-Bashforth or RK for Runge-Kutta: "); // User imput for Method selection
 scanf("%s", method); // recieves a selection from the User
 int selection = 0;

 if (strcmp(method, "AB") == 0) // Selects the AB solution
 {
   printf("\nAdams-Bashforth Selected, Please Wait"); 
   void waitFor (unsigned int secs); { // Adds a 5 second time delay for the User to read the above message
    unsigned int retTime = time(0) + 2;   
    while (time(0) < retTime);              
   }
   while (t <= 999 && timespace[1000][15] < 0.01*u_surface) // Far Field 1% shaft velocity condition
      {
         int y = 0; 
         selection = 0;
         while (y*0.01 <= rshaft) // Boundary Layer = shaft radius condition (max y of 15, 15 = 0.15*0.02)
         {
            if (t == 0) // Initial Conditions
         {
            timespace[t][y] = 0;  
            selection = 1;
         }
            else if (t > 0 && y == 0) // Boundary Conditions (Shaft)
         {
            u_surface = (rshaft * revrate);
            timespace[t][y] = u_surface;
            selection = 2;
         }
            else if (t > 0 && y == 1) //Explicit Euler Function
         {
            timespace[t][y] = timespace[t-1][y] + nu*((timespace[t-1][y-1] - 2*timespace[t-1][y] + timespace[t-1][y+1])/(dy*dy))*dt;
            selection = 4;
         }
            else if (t >= 2 && y > 0 && y < 15) // Adams Bashforth Function
         {
            double rhs_now = nu * ((timespace[t-1][y+1] - 2 * timespace[t-1][y] + timespace[t-1][y-1]) / (dy*dy));
            double rhs_past = nu * ((timespace[t-2][y+1] - 2 * timespace[t-2][y] + timespace[t-2][y-1]) / (dy*dy));
            timespace[t][y] = timespace[t-1][y] + (dt/2.0) * (3.0*rhs_now - rhs_past);
            selection = 5;
         }
         // Block calculates the analytical velocity and stores it
         double u_analytical = 0;
         if (t > 0) {
            double eta = (y * dy) / (2.0 * sqrt(nu * (t * dt))); // Calculate the similarity variable "eta"
            u_analytical = u_surface * erfc(eta);
         } else {
            u_analytical = 0; // At t=0, velocity is 0 everywhere
         }
         output = timespace[t][y];
         printf("\nt = {%d} || y = {%d} || u = {%.6e} || selection = {%i}", t, y, output, selection);
         y = y + 1;

         } 
         // Block calculates the numerical torque and stores it
         if (t > 0) {
          double mu = nu * rho; //declaring the dynamic viscosity of the fluid
          double du_dy_wall = (timespace[t][1] - timespace[t][0]) / dy; // Using first-order forward differencing scheme at the shaft surface
          double tau_w = mu * du_dy_wall;
          torque = 2.0 * pi * (rshaft * rshaft) * tau_w;
          torque_data_storage[t] = torque;
          // Block calculates the analytical torque and stores it separately
          double current_time = t*dt; // may cause issues
          double tau_w_analytical = (nu * rho) * (u_surface / sqrt(pi * nu * current_time));
          double torque_analytical = 2.0 * pi * (rshaft * rshaft) * tau_w_analytical;

          torque_analytical_storage[t] = -torque_analytical;

          printf("\nTime Step %d | Torque: %.5e N*m/m", t, torque);
         }
         t = t + 1;
         printf("\n ");
      }
 }

 else if (strcmp(method, "RK") == 0) // Selects the RK solution
 {
   printf("\nRunge-Kutta Selected, Please Wait");
   void waitFor (unsigned int secs); { // Adds a 5 second time delay to read the above message
    unsigned int retTime = time(0) + 2;   
    while (time(0) < retTime);              
   }
       while (t <= 999 && timespace[1000][15] < 0.01*u_surface) // Boundary Layer 2% shaft radius condition
      {
         int y = 0; // initial value of distance
         selection = 0;
         while (y*0.01 <= rshaft)
         {
            if (t == 0) // Initial Conditions
         {
            timespace[t][y] = 0;  
            selection = 1;
         }
            else if (t > 0 && y == 0) // Boundary Conditions (Shaft)
         {
            u_surface = (rshaft * revrate);
            timespace[t][y] = u_surface;
            selection = 2;
         }
            else if (t >= 2 && y > 0 && y < 15) // 4th Order Runge-Kutta Function
         {
            double u_center = timespace[t-1][y];
            double u_left = timespace[t-1][y-1];
            double u_right = timespace[t-1][y+1];
            double k1 = DERIV(u_center, u_left, u_right);
            double k2 = DERIV(u_center + (dt/2.0)*k1, u_left, u_right);
            double k3 = DERIV(u_center + (dt/2.0)*k2, u_left, u_right);
            double k4 = DERIV(u_center + dt*k3, u_left, u_right);

            timespace[t][y] = u_center + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4);
            selection = 5;
         }
         output = timespace[t][y];
         printf("\nt = {%d} || y = {%d} || u = {%.6e} || selection = {%i}", t, y, output, selection);
         y = y + 1;
         } 
         // Block calculates the numerical torque and stores it
         if (t > 0) {
          double mu = nu * rho; //declaring the dynamic viscosity of the fluid
          double du_dy_wall = (timespace[t][1] - timespace[t][0]) / dy; // Using first-order forward differencing scheme at the shaft surface
          double tau_w = mu * du_dy_wall;
          torque = 2.0 * pi * (rshaft * rshaft) * tau_w;
          torque_data_storage[t] = torque;
          // Block calculates the analytical torque and stores it separately
          double current_time = t*dt; // may cause issues
          double tau_w_analytical = (nu * rho) * (u_surface / sqrt(pi * nu * current_time));
          double torque_analytical = 2.0 * pi * (rshaft * rshaft) * tau_w_analytical;

          torque_analytical_storage[t] = -torque_analytical;

          printf("\nTime Step %d | Torque: %.5e N*m/m", t, torque);
         }
         t = t + 1;
         printf("\n ");
      }
}
    // Section below produces and opens a CSV file, then writes the values of our timespace array
    fp = fopen("C:\\Users\\J-for\\Downloads\\4304_CBA_1_output.csv", "w");
    if (fp == NULL) {
        printf("Error: Could not create file!\n");
        return 1;
    }

    
    fprintf(fp, "Time_Step"); // Loop which writes the header within the file, first the time index, then the y positions
    for (int j = 0; j < 16; j++) {
        fprintf(fp, ",y%d", j);
    }
    fprintf(fp, ",Torque_Nm_m"); // Writes a single torque column header
    fprintf(fp, ",Ana Torque_Nm_m\n"); // Writes a single analytical torque column header
    
    for (int i = 0; i < t; i++) { //Loop which writes the values created in the simulation
        fprintf(fp, "%d", i); 
        for (int j = 0; j < 16; j++) {
            fprintf(fp, ",%.6e", timespace[i][j]);
        }
        fprintf(fp, ",%.6e", torque_data_storage[i]);
        fprintf(fp, ",%.6e\n", torque_analytical_storage[i]);
    }
    fclose(fp);
    printf("\nThe code has run in its entirety. Results saved to 4304_CBA_1_output.csv\n");
    return 0;
}
