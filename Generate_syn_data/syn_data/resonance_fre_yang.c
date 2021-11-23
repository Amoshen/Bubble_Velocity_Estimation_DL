#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<complex.h>


#define V0 1535.0
#define K 3.89*1000000000
#define G 2.52*1000000
#define G_prime 123*1000



#define pi 3.1415926
//#define r 0.01
#define row 1612.0 /* density of sediment saturated with water */
#define gamma 1.31
#define rhog 0.717

#define sp 2.19
#define Cg 0.0311
#define ng 0.001

#define time_loop 1  /* pressure loop */
#define freq_loop 2500 /* number of frequency sampling*/
#define dtime 1000  /* time sampling for pressure change (s) */
#define dfreq 1 /* step interval of frequency (Hz) */
#define waterdepth 20
#define porosity 0.63
#define rhos 2851.2    /* density of sediment grain */
#define rhow 1000.0   /* density of water */

#define dts 0.0002

#define Tnumb 5000  /* number of time sampling */
#define fnumb Tnumb/2 /* number of frequency sampling*/
#define fstep (1.0/dts)/Tnumb /* step interval of frequency (Hz) *//*1hz*/

#define N 2 /* dimension of matrix */
#define layer 40 /* number of layers*/
#define seafloor 20 /* number of seafloor*/
#define layer_interval 1.0 /* layer interval (m)*/

#define amplifier 10.0 /* small velue for stablizing the calculation */

Mul(double _Complex [][N], double _Complex [][N], double _Complex [][N]);
Copy(double _Complex [][N], double _Complex [][N]);

void resonance_fre_generator(char strA[], char strB[], char strC[], char strD[], char strA1[], char strB1[], char strD1[], float r){
    long i,j;
FILE *A,*B,*C,*D,*A1,*B1,*C1,*D1;

static float as,Pa,Ph,Ps,P0,tidalheight,times[time_loop+1],X_value,frequency[freq_loop+1],X_ast,Y_ast,d_ast,f_ast,A_value, B_value,Xm,Ym,ans,final_ans,k,f0,dr,dt,d,df,alpha,density,Sg,P0ave;
static int depth = 10;
    A=fopen(strcat(strA,"resonance_velocity_tokyobaymm.txt"),"w");
    B=fopen(strcat(strB,"resonance_attenuation_tokyobaymm.txt"),"w");
    C=fopen(strcat(strC,"resonance_density_tokyobaymm.txt"),"w");
    D=fopen(strcat(strD,"resonance_qval_tokyobaymm.txt"),"w");

    A1=fopen(strcat(strA1,"resonance_velocity_tokyobay_freq.txt"),"w");
    B1=fopen(strcat(strB1,"resonance_attenuation_tokyobay_freq.txt"),"w");
    D1=fopen(strcat(strD1,"resonance_qval_tokyobay_freq.txt"),"w");

  for (i=1; i<=time_loop; i++)
    {
        times[i]=(i-1)*dtime;
        tidalheight = waterdepth - 0.75*cos(2*pi*times[i]/12/3600);
        Pa = 101325*exp((-9.8*tidalheight*0.02896968)/(15.01*8.314462618));
        Ph = tidalheight*1000*9.8;
        Ps = row*depth*9.8;
        P0 = Ph + Pa + Ps;
        P0ave=101325*exp((-9.8*waterdepth*0.02896968)/(15.01*8.314462618))+waterdepth*1000*9.8+row*depth*9.8;
        
        for (j=1; j<=freq_loop; j++)
        {
            frequency[j]=j*dfreq;
    
    
 
   X_value= r*sqrt(2*2*pi*frequency[j]*rhog*sp/Cg);
            
    B_value=3*(gamma-1)*(X_value*(sinh(X_value)+sin(X_value))-2*(cosh(X_value)-cos(X_value)))/(pow(X_value,2.0)*(cosh(X_value)-cos(X_value))+3*(gamma-1)*(sinh(X_value)-sin(X_value)));
    
    as=K/(gamma*P0+4/3*G);
            A_value=(1+pow(B_value,2.0))*(1+3*(gamma-1)/X_value*(sinh(X_value)-sin(X_value))/(cosh(X_value)-cos(X_value)));
    f0=(1.0/(2*pi*r))*sqrt((3*gamma*P0/(A_value*row))+(4*G/row));
    printf("Resonance_fre = %f",f0);
    k=2*pi*f0/V0;
    dt=B_value;
    dr=k*r;
    df=4*G_prime/row/(2*pi*f0)/(2*pi*f0)/r/r;
    d=dt+dr+df;
    
    
    f_ast=frequency[j]/f0;
    d_ast=d*pow(f_ast,2.0);
    X_ast=ng*(1.0-pow(f_ast,2.0))/((1.0-pow(f_ast,2.0))*(1.0-pow(f_ast,2.0))+d_ast*d_ast);
    Y_ast=ng*d_ast/((1.0-pow(f_ast,2.0))*(1.0-pow(f_ast,2.0))+d_ast*d_ast);
    Xm=X_ast;
    Ym=Y_ast;
            
            
          ans=(1.0+as*Xm)/2.0*(1+sqrt(1+as*as*Ym*Ym/(1+as*Xm)/(1+as*Xm)));
   
            if(ans<0)
            {
                ans=(1.0+as*Xm)/2.0*(1-sqrt(1+as*as*Ym*Ym/(1+as*Xm)/(1+as*Xm)));
            }
            

         
        printf("value=%f %f %f\n",ans,f0,frequency[j]);
            
   
    final_ans=sqrt(V0*V0/fabs(ans));
    
    

            fprintf(A," %f",final_ans);
            fprintf(A1," %f %f\n",frequency[j]/1000.0,final_ans);
            
            alpha=pi*frequency[j]/V0*final_ans/V0*as*Ym;
            fprintf(B," %f",alpha);
            fprintf(B1," %f %f\n",frequency[j]/1000.0,alpha);
            
            fprintf(D," %f",pi*frequency[j]/final_ans/alpha);
            fprintf(D1," %f %f\n",frequency[j]/1000.0,pi*frequency[j]/final_ans/alpha);
            
            
        
    }
        
        Sg=(row-(1-porosity)*rhos-porosity*rhow)/(-rhow*porosity+rhog*porosity);
        
        density=rhos*(1-porosity)+rhow*(1-Sg)*porosity+rhog*P0/P0ave*Sg*porosity;
     
                    fprintf(A,"\n");
                    fprintf(B,"\n");
                    fprintf(C," %f\n",density);
                    fprintf(D,"\n");
        fclose(A);
        fclose(B);
        fclose(C);
        fclose(D);
        fclose(A1);
        fclose(B1);
        fclose(D1);
       
    }
};

void yang_model_generator(char strA[], char strB[], char strC[], char strD[], char strv[]){    
    FILE *Af,*Bf,*fp1,*fp2,*fp3,*fp4,*fp5;
static double traveltime,depth,alpha[layer+1][fnumb+1],beta[layer+1][fnumb+1],freq[fnumb+1],hval[layer+1],qval[layer+1],vel[layer+1],qval_f[time_loop+1][fnumb+1],vel_f[time_loop+1][fnumb+1],density[layer+1],density_t[time_loop+1],cval[layer+1],tval[layer+1],Real_DFT_total[layer+1][Tnumb+1],Ima_DFT_total[layer+1][Tnumb+1],spike_total_real[layer+1][Tnumb+1],spike_total_imag[layer+1][Tnumb+1];
static long i,j,l,k,t,f;
static double _Complex w_value[layer+1][fnumb+1],S11[layer+1][fnumb+1],S12[layer+1][fnumb+1],S21[layer+1][fnumb+1],S22[layer+1][fnumb+1],F11[fnumb+1],F12[fnumb+1],F21[fnumb+1],F22[fnumb+1],A[N][N],B[N][N],C[N][N],R[fnumb+1],upgoing[layer+1][fnumb+1],downgoing[layer+1][fnumb+1],VA[N],VC[N],spike_total[layer+1][Tnumb+1];
    
    
    double P = 2.0*pi/Tnumb;

    Af=fopen(strcat(strA,"resonance_velocity_tokyobaymm.txt"),"r"); /* input file for velocity*/
    Bf=fopen(strcat(strB,"resonance_qval_tokyobaymm.txt"),"r");  /* input file for Q-value*/
    fp2=fopen(strcat(strC,"resonance_density_tokyobaymm.txt"),"r");  /* input file for density value*/
    fp1=fopen(strcat(strD,"reflectivity_yang_tokyobaymm.txt"),"w"); /* ouput file for reflectivity*/

    fp3=fopen(strcat(strv,"2-D_velocity_tokyobay.txt"),"w"); /* output file for velocity label about tokyo bay*/


    // reading velocity and Q-value data from file
   

        for (t=1; t<=time_loop; t++)
        {
            for (f=1; f<=fnumb; f++)
        {

        fscanf(Af,"%le",&vel_f[t][f]);
        fscanf(Bf,"%le",&qval_f[t][f]);
       
            qval_f[t][f]=qval_f[t][f]*amplifier;


        }
            fscanf(fp2,"%le",&density_t[t]);
    }
    
    for (f=0; f<=fnumb; f++)
    {
        freq[f]=f*fstep;
    }
    

    ///
    for (t=1; t<=time_loop; t++)
    {
        for (f=1; f<=fnumb; f++)

    {

// layer definition
        for (i=1; i<=layer; i++)
        {

                hval[i]=layer_interval;

            // 1st layer (sea)
            if (i<seafloor)
            {
                vel[i]=1500.0;
                qval[i]=10000000000.0;
                density[i]=1000.0;
            }
            
            
            // 2nd layer (sand layer containing gas bubbles)
            else if (i>=seafloor && i<=layer)
            {
               
                vel[i]=vel_f[t][f];
                qval[i]=qval_f[t][f];
                density[i]=density_t[t];

            }
            fprintf(fp3," %lf",vel[i]);
            

            
           }
           fprintf(fp3,"\n");





   
   

      for (i=1; i<=layer; i++)
      {

              alpha[i][f]=pi*freq[f]/vel[i]/qval[i];
              beta[i][f]=2.0*pi*freq[f]/vel[i];
          
              w_value[i][f]=cexp(I*beta[i][f]*hval[i])*exp(-alpha[i][f]);
          
          
      }
      

      
      for (i=1; i<=layer-1; i++)
      {
         
          cval[i]=(density[i+1]*vel[i+1]-density[i]*vel[i])/(density[i+1]*vel[i+1]+density[i]*vel[i]);
          tval[i]=2.0*density[i+1]*vel[i+1]/(density[i+1]*vel[i+1]+density[i]*vel[i]);
          
      }
      
      for (i=1; i<=layer-1; i++)
      {

              S11[i][f]=1.0/w_value[i][f]/tval[i];
              S12[i][f]=cval[i]*w_value[i][f]/tval[i];
              S21[i][f]=cval[i]/w_value[i][f]/tval[i];
              S22[i][f]=w_value[i][f]/tval[i];

      }
      
      
      
      /* derivation of R */
      

          A[0][0]=S11[1][f];
          A[0][1]=S12[1][f];
          A[1][0]=S21[1][f];
          A[1][1]=S22[1][f];
          
          
      for (i=1; i<=layer-2; i++)
      {
          
          B[0][0]=S11[i+1][f];
          B[0][1]=S12[i+1][f];
          B[1][0]=S21[i+1][f];
          B[1][1]=S22[i+1][f];
          
          Mul(B, A, C);
          
          Copy(C, A);
          
      }
          F11[f]=C[0][0];
          F12[f]=C[0][1];
          F21[f]=C[1][0];
          F22[f]=C[1][1];


          R[f]=F12[f]/(F11[f]-F12[f]);
   
      

      
      /* derivation of Rk */
      

          VA[0]=-R[f];
          VA[1]=1.0+R[f];

          
          for (i=1; i<=layer-1; i++)
          {
              VC[0]=S11[i][f]*VA[0]+S12[i][f]*VA[1];
              VC[1]=S21[i][f]*VA[0]+S22[i][f]*VA[1];
              
              upgoing[i][f]=VC[0];
              downgoing[i][f]=VC[1];
              
              VA[0]=VC[0];
              VA[1]=VC[1];
           
          }
          
      }
      
      for (i=0; i<=layer-1; i++)
      {
          for (k=1; k<=Tnumb; k++)
          {
              
              spike_total[i][k]=0.0;
          }
      }
   
      
      for (i=0; i<=layer-1; i++)
      {
          for (k=1; k<=Tnumb; k++)
          {
              for (j=1; j<=fnumb; j++)
              {
                  
                  spike_total[i][k] += (upgoing[i][j]+downgoing[i][j])*cexp(-I*2.0*pi*freq[j]*k*dts);

              }
          }
      }
      
   
    
    traveltime=0.0;
    depth=0.0;
    
    for (i=1; i<=layer-2; i++)
    {
        if(i>0)
        {
            traveltime+=hval[i]/vel[i];
            depth+=hval[i];
        }
        
    for (k=1; k<=Tnumb; k++)
    {
        fprintf(fp1," %lf",creal(spike_total[i][k]));
    }
        fprintf(fp1,"\n");

        if(traveltime > Tnumb*dts)
        {
            printf("traveltime is larger than Tnumb*dt !! \n");
            printf("traveltime=%f Tnumb*dt=%f \n",traveltime,Tnumb*dts);
        }
    }

    
    }
    fclose(Af);
        fclose(Bf);
        fclose(fp2);
        fclose(fp1);
        fclose(fp3);



   };

int main(int argc, char *argv[])
{
    //srand((unsigned)time(NULL));
    FILE *radius;
    float r = atof(argv[1])*1e-2;
    //float r = (rand()%1000+1)*1e-4;/*0.1mm-100mm*/
    char strA[50],strB[50],strC[50],strD[50],strE[50],strF[50],strG[50],strH[50],strA1[50],strB1[50],strD1[50],strv[50];
    radius=fopen("radius_tokyobay.txt","a"); /* input file for velocity*/
    fprintf(radius,"%.2f\n",r*1e3);
    fclose(radius);
    sprintf(strA,"%.2f",r*1e3);
    sprintf(strB,"%.2f",r*1e3);
    sprintf(strC,"%.2f",r*1e3);
    sprintf(strD,"%.2f",r*1e3);
    sprintf(strE,"%.2f",r*1e3);
    sprintf(strF,"%.2f",r*1e3);
    sprintf(strG,"%.2f",r*1e3);
    sprintf(strH,"%.2f",r*1e3);
    sprintf(strA1,"%.2f",r*1e3);
    sprintf(strB1,"%.2f",r*1e3);
    sprintf(strD1,"%.2f",r*1e3);
    sprintf(strv,"%.2f",r*1e3);
    resonance_fre_generator(strA, strB, strC, strD, strA1, strB1, strD1, r);
    yang_model_generator(strE, strF, strG, strH,strv);
    return 0;
          		
}

Mul(double _Complex x[][N], double _Complex y[][N], double _Complex z[][N]){
    int i,j,k;
    
    for(i = 0; i<N; i++){
        for(j = 0; j<N; j++){
            z[i][j] = 0.0;
            for(k = 0; k<N; k++){
                z[i][j] = z[i][j] + x[i][k] * y[k][j];
            }
        }
    }
}
Copy(double _Complex x[][N], double _Complex y[][N]){
    int i,j;
    
    for(i = 0; i<N; i++){
        for(j = 0; j<N; j++){
            y[i][j] = x[i][j];
        }
    }
}
