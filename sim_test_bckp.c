#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// constants
//#define L 500 //boxsize
#define N 100  //number of atoms
#define timesteps 500000 //runntime
#define dt 0.01  //time step size
#define Dim 3  //dimensions
#define pi 3.14159 
#define nhis 250  //total no. of bins for rdf
//#define gamma 0.01 //frictional coefficient
#define Kbend 1  //for bend potential

//Lennard-Jones potential parameter
#define epsilon 1.0  // potential well depth 
#define sigma 1.0   // constant distance between 2 atom when potential is 0 -- van der Waal's radius
double K=50;  //for bond potential
int bondR=1;
double t;
double pos[N][Dim];
double pos_m[N][Dim], vel_m[N][Dim];
double vel[N][Dim];  
double sumv_x, sumv_y, sumv_z, sumv2;
double r = 0 ;
double g[nhis];
int count;//=ceil(cbrt(N));
double L; //=count*sigma;
double rho;// = N / (L * L * L) ;  //number density of the particles in the system
double delg; // bin size for rdf
double temp=1.0; // instantaneous temperature
double etot; // total energy per particle 
int s;
int ngr;int ig;
int sw, it0;
double nsamp = 0.5; //frequency at which the dif() is called
int ntel, delt, tt0=0, t0 = 0;
double dtime=0; 
double ntime[timesteps], vacf[timesteps], r2t[timesteps];
double vx0[N][timesteps], x0[N][timesteps], time0[timesteps];
double tim = 0;  // time for MSD
int t0max = 500;  // ???
double x[N], vx[N];
int it0;
int Nbond = N-1;
int Nbend = N-2;
int a2b[N][N];  //angle to bond
double bx[N], by[N], bz[N];  //pseudo bond vectors
double b[N];  //pseudo bond length
double angle[N-2];
double a1[N], a2[N], a3[N];  //beads in angle


//Function to compute positions and velocities

void init()
{
    sumv_x = 0, sumv_y = 0, sumv_z = 0, sumv2 = 0;
    
    //position assigning in lattice
    int n = 0;
    double scale=1.2;
    count=ceil(cbrt(N));
    L=count*sigma*scale; 
    //printf("L %f", L);
    rho = N / (L * L * L) ;  //number density of the particles in the system
    for (int x = 0; x < count; x++)
    {
        for (int y = 0; y < count; y++)
        {
            for (int z = 0; z < count; z++)
            {
                if (n < N)
                {
                    pos[n][0] = (double)x ;
                    pos[n][1] = (double)y ;
                    pos[n][2] = (double)z ;
                    n++;
                }
            }
        }
    }
                    pos[0][0] = 0.0;//(double)x ;
                    pos[0][1] = 0.0;//(double)y ;
                    pos[0][2] = 0.0;//(double)z ;

    FILE* read_pos=fopen("last_snap","r");
    for (int i=0; i < N; i++) 
    {
        for (int j=0; j < 3; j++)
        {
     //       fscanf(read_pos,"%lf",&pos[i][j]);
            if(j==2 && i>0)pos[i][2]=pos[i-1][2]+1.20;else pos[i][j]=0;
            //pos[i][j] = pos[i][j]*sigma*scale;//L * (rand() / (double)RAND_MAX ) ;           // random position
            vel[i][j] = 0.1 * (rand() / (double)RAND_MAX - 0.5 ) ;    // random velocity
            
            if (j==0) {sumv_x += vel[i][j] ;}
            if (j==1) {sumv_y += vel[i][j] ;}
            if (j==2) {sumv_z += vel[i][j] ;}

            sumv2 +=  (vel[i][j] * vel[i][j]) ;  
            //printf("%lf\t",pos[i][j]);       
        }
     //   printf("%lf\t%lf\t%lf\n",pos[i][0],pos[i][1],pos[i][2]);
       // printf("\n");
    }
    fclose(read_pos);

    //exit(-1);
    sumv_x = sumv_x / N ;
    sumv_y = sumv_y / N ;
    sumv_z = sumv_z / N ;
        //printf("sumv_x %f ", sumv_x);
        //exit(-1);
    sumv2 = sumv2 / N ;
    double fs = sqrt(3 * temp / sumv2) ;

    //printf("%f %f %f %f fs: %f", sumv_x,sumv_y,sumv_z,sumv2,fs);
//printf("%d\t%lf\t%lf\n",count,L,rho);exit(-1);
    //exit(-1);
    for (int i=0; i < N; i++)
    {
       for (int j=0; j < 3; j++)
       {
          if (j==0) {vel[i][j] = (vel[i][j] - sumv_x) * fs ;} 
          if (j==1) {vel[i][j] = (vel[i][j] - sumv_y) * fs ;}
          if (j==2) {vel[i][j] = (vel[i][j] - sumv_z) * fs ;}
          pos_m[i][j] = pos[i][j] - vel[i][j] * dt ;
       }
    }


    //Initialization of particles in bending interactions
    for (int i = 0; i < Nbend; i++) 
    {
        a1[i] = i;
        a2[i] = i + 1;
        a3[i] = i + 2;
    }

    //initialization of rdf
    if (s == 0)
    {
        delg=L/nhis;
        for (int i=0; i < nhis; i++)
        {
            g[i] = 0;
        }

    }


   //initialization of dif 
   if (sw == 0)
   {
     ntel = 0; // time counter
     dtime =dt *nsamp;  // time b/w 2 samples
     for (int i=0; i<timesteps; i++)
     {
       ntime[i] = 0;  //no. of samples for time i
       vacf[i] = 0;
       r2t[i] = 0;
     }
   }
     
}



//function to compute forces

double compute_forces(double force[N][Dim]) {
    for (int i = 0; i < N; i++) {
        for (int d = 0; d < Dim; d++) {
            force[i][d] = 0.0;
        }
    }
    
    double f = 0 ;
    double ftot = 0;
    
    double rc = L/2;
    double rc2 = rc * rc ;
    double rc6 = rc2 * rc2 * rc2 ;
    double rc12 = rc6 * rc6 ;
    double sigma2=sigma*sigma;
    double rcut2=sigma2*2.5*2.5;//pow(2.0,1.0/3.0);
    double ecut = 4 * ((1/rc12) - (1/rc6)) ;
    double en=0.0;
    for (int i = 0; i < N-1-bondR ; i++)    // for atom i
    {
        for (int j = i + 1+bondR; j< N ; j++)   // for force b/w atom i with atom i+1 
        {
        
            double dx = pos[i][0] - pos[j][0];
            double dy = pos[i][1] - pos[j][1];
            double dz = pos[i][2] - pos[j][2];
            // Periodic boundary conditions
            //printf("before :  %lf\t%lf\t%lf\n",dx,dy,dz);
            //while (dx<-L/2)dx += L ;while(dx>L/2)dx-=L;//* ((int) (dx/L)) ;
            //while (dy<-L/2)dy += L ;while(dy>L/2)dy-=L;//* ((int) (dx/L)) ;
            //while (dz<-L/2)dz += L ;while(dz>L/2)dz-=L;//* ((int) (dx/L)) ;
            //dy -= L * ((int) (dy/L)) ;
            //dz -= L * ((int) (dz/L)) ;
            //printf("after :  %lf\t%lf\t%lf\n",dx,dy,dz);

            double r2 = dx*dx + dy*dy + dz*dz ;
             //printf("%lf\t%lf\t%lf\t%lf\n",sqrt(r2),dx,dy,dz);
            if (r2 < rcut2) 
            {
               
              double r2i = sigma2/r2 ;
              double r6i = r2i * r2i * r2i ;
             // printf("%lf\n",r2i);
              f = 48 * r6i * (r6i - 0.5)/r2 ;  // Lennard - Jones potential
              
              
              
              //fbond = sqrt(pow(fbond_x,2)+pow(fbond_y,2)+pow(fbond_z,2));

              ftot = f ;                          // Total force 
              
              //printf("LJ force, Total force: %lf %lf\n", f, ftot);

              force[i][0] += (f * dx) ;   // update force
              force[i][1] += (f * dy) ;
              force[i][2] += (f * dz) ;
              force[j][0] -= (f * dx) ;
              force[j][1] -= (f * dy) ;
              force[j][2] -= (f * dz) ;
              en += 4 * r6i * (r6i - 1) - ecut ;
             // printf("en: %lf\n", en);
            }
              //else {printf("not working\t"); exit(1); }
        }
    }
    return en/N;
    //FILE *energy = fopen("energy.txt", "w");
    //fprintf(energy, "%lf %lf \n", en, t);
    //fclose(energy);
}



//Function to calculate bond potential

double bond(double force[N][Dim])
{

    double fbond = 0;
    double rb = 1.2;
    double en_bond = 0;

    for(int i=0; i < Nbond; i++)
    {
    
     double dx = pos[i][0] - pos[i+1][0];
     double dy = pos[i][1] - pos[i+1][1];
     double dz = pos[i][2] - pos[i+1][2];
           // while (dx<-L/2)dx += L ;while(dx>L/2)dx-=L;//* ((int) (dx/L)) ;
           // while (dy<-L/2)dy += L ;while(dy>L/2)dy-=L;//* ((int) (dx/L)) ;
           // while (dz<-L/2)dz += L ;while(dz>L/2)dz-=L;//* ((int) (dx/L)) ;
    
     double r = sqrt((dx * dx)+(dy * dy)+(dz * dz)) ;

     //if(r<rb+0.5 && r>rb-0.5)K=10;else K=1;
     fbond = -K *(r - rb)/r ;  // bond potential   
     en_bond+=0.5*K*(r-rb)*(r-rb);//fbond = -K *(r - rb)/r ;  // bond potential   
    //printf("%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",fbond,r,dx,dy,dz,i,i+1);
     force[i][0] += fbond * dx ;  // force update
     force[i][1] += fbond * dy ;
     force[i][2] += fbond * dz ;
     force[i+1][0] -= fbond * dx ;
     force[i+1][1] -= fbond * dy ;
     force[i+1][2] -= fbond * dz ;
     //en_bond += fbond;

     }
     return en_bond/N;
}

void benddef()
{
  int i,j,k,p1, p2, p3;
  double b1x, b2x, b1y, b2y, b1z, b2z;
  double tmp1;

  // Initialize a2b[][] to -1
  for (int i = 0; i < N; i++) 
  {
    for (int j = 0; j < N; j++) 
    {
        a2b[i][j] = -1;  // Initialize to -1 to indicate no mapping
    }
  }  


  /* bond vectors */
  for (i=0; i<Nbond; i++) 
  {
    j= i; //Central atom 
    k= i+1; //First atom

    a2b[k][j] = i ; // Map pair (j, k) to index i
    
    bx[i]=pos[k][0]-pos[j][0];
    by[i]=pos[k][1]-pos[j][1];
    bz[i]=pos[k][2]-pos[j][2];
    
    //while (bx[i]<-L/2)bx[i] += L ;while(bx[i]>L/2)bx[i]-=L;//* ((int) (dx/L)) ;
    //while (by[i]<-L/2)by[i] += L ;while(by[i]>L/2)by[i]-=L;//* ((int) (dx/L)) ;
    //while (bz[i]<-L/2)bz[i] += L ;while(bz[i]>L/2)bz[i]-=L;//* ((int) (dx/L)) ;
    
    b[i]=sqrt(bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i]);
    
    bx[i]=bx[i]/b[i]; 
    by[i]=by[i]/b[i]; 
    bz[i]=bz[i]/b[i];

    
   //printf("Bond vector for pair (%d, %d): bx = %lf, by = %lf, bz = %lf, b = %lf\n", j, k, bx[i], by[i], bz[i], b[i]);
   }

  /* bond angles */
  for (i=0; i<Nbend; i++) 
  {
    p1=a1[i]; //First atom 
    p2=a2[i]; //Central atom
    p3=a3[i]; //Second atom
   
    b1x=bx[a2b[p2][p1]];
    b1y=by[a2b[p2][p1]];
    b1z=bz[a2b[p2][p1]];

    b2x=bx[a2b[p3][p2]];
    b2y=by[a2b[p3][p2]];
    b2z=bz[a2b[p3][p2]];

    tmp1=b1x*b2x+b1y*b2y+b1z*b2z;
    if (tmp1<-1.0) {tmp1=-1.0; }
    if (tmp1>1.0) {tmp1=1.0; }

    angle[i]=pi-acos(tmp1);
   
    //printf("Angle for triplet (%d, %d, %d): %lf\n", p1, p2, p3, angle[i]);
   }

}



//Function to calculate bend potential

double bend(double force[N][Dim])
{
  int i, j, k, l;
  double fbend;
  double e=0;
  double theta, dtheta=0;
  double cth, sth, a, b1;
  double dix, diy, diz, dkx, dky, dkz;
  double en_bend;
  //double tmp;


  for(i=0; i<Nbend; i++)
  {
   j=a1[i]; 
   k=a2[i]; 
   l=a3[i]; 
   
   theta = angle[i];
   cth=cos(theta); 
   sth=sin(theta);
   
   dtheta=theta-pi;

   e=dtheta*dtheta;

   a=-2*Kbend*dtheta/(sth*b[a2b[k][j]]); 
   b1=-2*Kbend*dtheta/(sth*b[a2b[l][k]]);

   dix=a*(-bx[a2b[k][j]]*cth - bx[a2b[l][k]]); 
   diy=a*(-by[a2b[k][j]]*cth - by[a2b[l][k]]); 
   diz=a*(-bz[a2b[k][j]]*cth - bz[a2b[l][k]]);
  
   //printf("Force components for atom %d: dix = %lf, diy = %lf, diz = %lf\n", j, dix, diy, diz);

   dkx=b1*(bx[a2b[l][k]]*cth + bx[a2b[k][j]]); 
   dky=b1*(by[a2b[l][k]]*cth + by[a2b[k][j]]); 
   dkz=b1*(bz[a2b[l][k]]*cth + bz[a2b[k][j]]);

   force[j][0] +=dix;  
   force[j][1] +=diy;  
   force[j][2] +=diz;
   
   force[k][0] +=(-dix-dkx);  
   force[k][1] +=(-diy-dky);  
   force[k][2] +=(-diz-dkz);
   
   force[l][0] +=dkx;  
   force[l][1] +=dky;  
   force[l][2] +=dkz;

   en_bend+=a;
   }

   return (en_bend/N) ;

}

//Function to integrate equations of motion

void integrate(double force[N][Dim], double en)
{
 
 
  sumv_x = 0, sumv_y = 0, sumv_z = 0;
  sumv2 = 0;
  double vi = 0;
  double xx = 0;

    for(int i=0; i<N; i++)
        {   
            for(int d=0; d<Dim; d++)
            {
                //double xx = pos[i][d] + (vel[i][d] * dt) + (0.5 * (dt*dt) * (force[i][d] + eta - (gamma * vi))) ;
                double xx = 2 * pos[i][d] - pos_m[i][d] + (dt * dt) * force[i][d]; //Verlet Algorithm
                vi = (xx - pos_m[i][d]) / (2 * dt) ;
                //vi = vel[i][d] + (dt * (force[i][d] + eta - (gamma * vi))) ;
                //sumv_x += vi ;    //velocity center of mass
                //sumv_y += vi ;
                //sumv_z += vi ;
                sumv2 += vi * vi ;   //total kinetic energy
                pos_m[i][d] = pos[i][d] ;  //update positions previous time
                pos[i][d] = xx ;    //update positions current time

                
               // vel[i][d] = vi ;   //velocity update
               
               // printf("Pos before: %f\n", pos_m[i][d]);
               // printf("Pos After: %f\n", pos[i][d]);

            } 
        }
	//printf("%lf\n", sumv2);
        temp = sumv2/ (3*N) ;  //instantaneous temperature
        sumv2 *= 0.5/N ;
        //printf("%f\t%f\t%f\n ", pos[i][0], pos[i][1], pos[i][2]) ;
        //printf("\n");
        etot = (en + sumv2);// / N ;  //total energy per particle
        //printf("%f", etot) ;   
       // }
 //printf("en, sumv2: %lf, %lf\n", en, sumv2);
       
    
}



//Function to calculate the Radial Distribution Function

void gr()
{ 
    //int ngr; 
    //int ig;
    //FILE *rdf = fopen("rdf.txt", "w");
     

        for (int i=0; i < N-1; i++)
        {
            for (int j=i+1; j < N; j++)
            {
                double dx = pos[i][0] - pos[j][0] ;
                double dy = pos[i][1] - pos[j][1] ;
                double dz = pos[i][2] - pos[j][2] ;
                
                //printf("%lf %lf %lf\n", dx, dy, dz);

                //periodic boundary conditions
     //           while (dx<-L/2)dx += L ;while(dx>L/2)dx-=L; //* (L * (int) (dx/L)) ;
     //           while (dy<-L/2)dy += L ;while(dy>L/2)dy-=L; //* (L * (int) (dx/L)) ;
     //           while (dz<-L/2)dz += L ;while(dz>L/2)dz-=L; //* (L * (int) (dx/L)) ;

                double r = sqrt((dx*dx) + (dy*dy) + (dz*dz)) ;

                //printf("r: %lf \n", r);
                //exit(1);

                ig = (int) (r/delg);
                if (ig < nhis)
                {
                    //printf("ig %d\n", ig);
                    //exit(1);
                    g[ig] = g[ig] + 2.0 ;   // contribution of particles i and j
     //               printf("ig, g[ig] %d  %lf\n", ig, g[ig]);
                }
            }
        }

}


//Function to calculate diffusion (MSD)

 void diffusion()
{
    static int ntel;
    if(sw==0){
    ntel=0;}
   if (sw == 1)
   {
     ntel++;
     //int it0 = 1;  // ???  (is it "t"?)
     //it0 /= 10;
     //printf("it0, ntel %d %d\n", it0, ntel);
     if (ntel%it0 == 0)
     {
       t0 += 1;
       tt0 = fmod((t0 - 1), t0max) + 1;
       time0[tt0] = ntel;
       for (int i=0; i<N; i++)
       {
         x[i] = sqrt(pow(pos[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2));
         vx[i] = sqrt(pow(vel[i][0],2) + pow(vel[i][1],2) + pow(vel[i][2],2));
         x0[i][tt0] = x[i];  //store position for given t=0
         vx0[i][tt0] = vx[i];  //store velocity for given t=0
       }
     }
     int m = ((t0<t0max)? t0 : t0max) ;
     for (int t=0; t < m; t++)  //update vacf and r2, for t=0
     {
       delt = ntel - time0[t] + 1; //actual time minus t=0
       //printf("delt: %d\n", delt);
       if (delt < timesteps)
       {
         //printf("delt: %d\n", delt);
         ntime[delt] += 1;
         //printf("ntime[delt]: %lf\n", ntime[delt]);
         for (int i=0; i<N; i++)
         {
           vacf[delt] += vx[i] * vx0[i][t];  //update velocity autocorrelation
           //printf("%lf %d\n", vacf[delt], delt);
           r2t[delt] += pow((x[i] - x0[i][t]), 2);  //update mean-squared displacement
         }
       }
     }
   }



} 





// Main function

int main() {
    double force[N][Dim], en, en_bond, en_bend;  
    double netpos, sumnetpos, stddev;
    
    s=0;  // rdf initialization
    sw=0;  // diffusion initialization

    init();   // initialization
    

    FILE *energy = fopen("energy.txt", "w");
    FILE *position = fopen("position.pdb", "w");
    FILE *bonds = fopen("bonds.pml", "w");
    FILE *sd = fopen("SD.txt", "w"); 
    FILE *temperature = fopen("temperature.txt", "w");
    FILE *rdf = fopen("rdf.txt", "w");
    FILE *dif = fopen("dif.txt", "w");
    FILE *vac = fopen("vac.txt", "w");


    for (double t = 0; t < timesteps; t++) 
    {
        if((int)t%500==0){
        printf("%lf\n", t);} 
        en=compute_forces(force);   // determine the forces
        if(bondR==1)en_bond=bond(force);               // determine bond forces
        //benddef();  //define bond vectors and bond angles
        //en_bend=bend(force);     //determine bend forces

        integrate(force, en);   // integrate equations of motion
        
        
   //     eta = generate_gaussian(0.0, 1.0) ;
   //     fprintf(gaussian, "%lf %lf\n", eta, t);

        t += dt ;
      /*  for (int i=0; i<N; i++)
        {
                fprintf (forces, "%lf %lf \n", force[i][0]+force[i][1]+force[i][2], t) ;
        }  */


    // Output positions
        
        char *a = "CA"; 
        char *b = "GLY";
        if((int)t%500==0){

        double e_potential=en+en_bond+en_bend;
        double e_total=e_potential+sumv2;
        // output energy
        fprintf(energy, " %d %lf %lf %lf %lf %lf %lf %lf\n", (int)t, en, en_bond, en_bend, e_potential,sumv2, e_total,temp); // LJ, bond, bend, KE, total
        // output temperature
        fprintf(temperature, "%lf %lf\n", temp, t); 
        
        for (int i = 0; i < N; i++) 
        {
            fprintf(position, "ATOM   %4u  %3s %3s  %4u    %8.3f%8.3f%8.3f  1.00\n", i+1, a, b, i+1, pos[i][0], pos[i][1], pos[i][2]);
           // fprintf(position, "ATOM\t%d\tCA\tGLY\t%d\t%lf\t%lf\t%lf\t1.00\n", i+1, i+1, pos[i][0], pos[i][1], pos[i][2]);
        }
        fprintf(position, "ENDMDL\n");  
        }
    
    //Output standard deviation

        sumnetpos = 0; stddev = 0; 
       
        if((int)t%500==0){
        for (int i = 0; i < N; i++)
        {
            sumnetpos += sqrt(pow(pos[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2));
        }

        double mean = sumnetpos / N ;

        for (int i=0 ; i < N; i++)
        {
            netpos = 0;
            netpos = sqrt(pow(pos[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2));
            stddev += pow(netpos - mean, 2) ;
        }
	      stddev/=N;
	      stddev=sqrt(stddev);
        fprintf(sd, "%lf %lf\n", stddev, t); 
        }
    
    
    it0 = ((int) t +1) / 10 ; // no. of time a subroutine is being called
    //printf("it0= %d, it0 mod 10 = %lf\n", it0, fmod(it0, 10));

    // rdf sample 
    
    if((int)t>=100000){
    if ((int)t%500==0)
    {
      ngr++;
      gr();
    }}

    // dif 

    //sw=0;
    //diffusion();

     int i = it0 ; // / 10;
     
     //if ((int)t>=10000){
     if (it0%10==0.0)
      {  
     

        sw=1;
        diffusion(); 
       // printf("%lf", ntime[10]); 
    
       sw=2;
       

      // printf("i: %d\t", i);  
       //printf("before: %lf\n", ntime[4]);
       tim = dtime * (i+0.5);  //time

       if(it0>=1000 && ntel%it0==0)  // ntel = no. of times sw=1 called |  it0 = no. of times sampling done
       {
        vacf[i] = vacf[i] / (N * ntime[i] * 100); //volume velocity autocorrelation
        r2t[i] = r2t[i] / (N * ntime[i] * 100);  //mean-squared displacement 

        //printf("vacf[%d]: %lf\n",i, vacf[i]);
       
        fprintf(dif, "%lf %lf\n", r2t[i], tim);
        fprintf(vac, "%lf %lf\n", vacf[i], tim);

       }

      }//}

  }



    
    // rdf results initialisation    NOTE:  (g(r) = no. density at a distance r)
    
  
        printf("ngr %d\n", ngr);
        for (int i=0; i < nhis; i++)
        {
            r = delg * (i + 0.5);    // distance r
            double vb = (pow(i+1,3) - pow(i,3)) * pow(delg, 3);  // volume b/w bin i+1 and i
            //printf("vb: %lf \n", vb);
            double nid = (4.0/3.0) * pi * vb * rho ;
            //printf("%lf %lf\t",r, g[i]); 
            g[i] = g[i] / (ngr * N * nid);   // normalize g(r)
            //printf("%lf %d %d %lf\n",g[i], ngr,N,nid); 

            fprintf(rdf, "%lf %lf\n", g[i], r); 
        }

    



/*   sw=2;

   if (sw == 2)
   {
     for (int i=1; i<timesteps; i++)
     {
       time = dtime * (i+0.5);  //time

       if(ntime[i]!=0.0 && time < 2)
       {
        vacf[i] = vacf[i] / (N * ntime[i]);  //volume velocity autocorrect
        r2t[i] = r2t[i] / (N * ntime[i]);  //mean-squared displacement

        fprintf(dif, "%lf %lf\n", r2t[i], time);
       }
     }
   }  */


        //output diffusion
//        int sw;

//        sw=0;
//        diffusion(sw);


//        sw=2;
//        diffusion(sw);



        for (int i = 1; i < N; i++)
        {
          fprintf(bonds, "bond resi %d, resi %d\n", i, i+1);
        }


    fclose(energy);
    fclose(position);
    fclose(bonds);
    fclose(sd);
    fclose(temperature);
    fclose(rdf);
    fclose(dif);
    fclose(vac);

    return 0; 

}
