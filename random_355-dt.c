   #include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define		pi      3.1415926535897932384626433832795028841971693993751058

#define		N	6 //Número osciladores 
#define		ni	0 //Inicio del ciclo para los osciladores


//Para guardar los datos

FILE    *mag_x,
        *mag_y,
        *mag_z,
        *CI_mag_x,
        *CI_mag_y,
        *CI_mag_z;
        
void dM( int n); 	//La funcion
double dMX[N]; //Las derivadas
double dMY[N];
double dMZ[N];

double k1x[N],k2x[N],k3x[N],k4x[N],k1y[N],k2y[N],k3y[N],k4y[N],k1z[N],k2z[N],k3z[N],k4z[N],mx[N],mx_1[N],mx_rk[N],my[N],my_1[N],my_rk[N],mz[N],mz_1[N],mz_rk[N],d,t,t_rk,dt,tf; //Declarar las cositas para el RK
//Las cosas con "_rk[n]" se llaman así porque tuve que renombrar esas cosas en las funciones dMX[n], dMY[n] y dMZ[n] para poder hacer el Runge Kutta

double rx[N], ry[N], rz[N]; //Posiciòn de cada oscilador

int n, n_otro; //Para hacer el ciclo de los osciladores, el contador de los osciladores

//Cositas para no guardar tantos datos

long int B1_red; //Beta 1 multiplicado por e4 y redondeado a un entero


long int B1_div = 10; //Cada cuantos B1 guardo beta1

long int t_div = 10; //Cada cuantos t guardo tiempo

//Otros parametros

double  B1, B1F, dB1;

double  w=0.08;                  //  Voltage frequency
double  h=0.1, phi=0.3;               // Parameters asociated to External Fields

double  B0x=0.0,B0z=0.05;              // Normalized constants of Anisotropy and Demagnetization
double  alpha=0.005; 

double V0 = 1256.64; //Volumen del material

/***********************************************************************/

int main()
{    
     srand48(355);
    /*===================================*/                  
    /*-------------PARAMETROS----------- */   
    /*===================================*/   
    

//Las cosas del tiempo. 
    dt = 2.0*pi/(w*1000.0);
    tf = 18000.0*2.0*pi/w;   
    
    B1=0.16; // ELl beta1 inicial
    B1F=0.16; // El beta1 final
    dB1=0.0001; // El paso de beta 1
    
    
/***********************************************************************/   
   
    /*===================================*/                  
    /*--------------GUARDAR------------- */   
    /*===================================*/             
    
    
    char	Name_Mx[150];
    char	Name_My[150];
    char	Name_Mz[150];
    
/***********************************************************************/
    
    	/*===================================*/      
    	/*-----------ABRIR CSV CI----------- */   
    	/*===================================*/  
    		
    CI_mag_x=fopen("CI_MX.csv","w");
    CI_mag_y=fopen("CI_MY.csv","w");
    CI_mag_z=fopen("CI_MZ.csv","w");
      		
    //Esto para el titulo del espacio que necesito vacío
    fprintf(CI_mag_x,"Beta1,"); 
    fprintf(CI_mag_y,"Beta1,");
    fprintf(CI_mag_z,"Beta1,");

        
	/*===================================*/            
	/*-----------------CI--------------- */   
	/*===================================*/    
    
     mx[0]= -0.391009;
     my[0]= 0.763533;
     mz[0]= -0.513936;
          
     mx[1]= -0.947049;
     my[1]= -0.298996;
     mz[1]= -0.117043;

     mx[2]= 0.092641;
     my[2]= -0.706958;
     mz[2]= -0.701162; 
     
     mx[3]= -0.578792;
     my[3]= 0.682681;
     mz[3]= -0.446034;
          
     mx[4]= 0.440036;
     my[4]= 0.28852;
     mz[4]= 0.850367;

     mx[5]= 0.066988;
     my[5]= -0.203915;
     mz[5]= -0.976694;
          
    // mx[6]= -0.997497;
    // my[6]= -0.031585;
    // mz[6]= -0.06326;
                                                   
    /*--------------POSICIONES------------- */   
    
    
    for (n = ni; n < N; n = n+1)
    {
        rx[n] = 0.0;
        ry[n] = 0.0;
        rz[n] = 50.0*n;
    }

    /*======================================*/
    
    for (n = ni; n < N; n = n+1)
    {	
    	// Primera fila CI
    	fprintf(CI_mag_x,"%d,",n);
       	fprintf(CI_mag_y,"%d,",n);
       	fprintf(CI_mag_z,"%d,",n);
    }

    /*=====================================================*/      
    /*--------------------CICLO DE BETA1-------------------*/   
    /*=====================================================*/ 
    while(B1<=B1F)
    {    
   
    	B1_red = lround(B1*10000);
    	

	if(B1_red % B1_div == 0)
	{
	
		/*===================================*/      
    		/*-----------ABRIR CSV DATOS-------- */   
    		/*===================================*/ 
    		
		sprintf(Name_Mx,"_B1(e4)_%ld_MX.csv",B1_red);
		sprintf(Name_My,"_B1(e4)_%ld_MY.csv",B1_red);
		sprintf(Name_Mz,"_B1(e4)_%ld_MZ.csv",B1_red);

        	mag_x=fopen(Name_Mx,"w");
        	mag_y=fopen(Name_My,"w");
      		mag_z=fopen(Name_Mz,"w");
      		
      		//Esto para el titulo del espacio que necesito vacío
      		fprintf(mag_x,"Oscilador,");
    		fprintf(mag_y,"Oscilador,");
    		fprintf(mag_z,"Oscilador,");
    		
    		
        	//Primera Columna CI
    		fprintf(CI_mag_x,"\n%ld,",B1_red);
        	fprintf(CI_mag_y,"\n%ld,",B1_red);
        	fprintf(CI_mag_z,"\n%ld,",B1_red);
        		
		for (n = ni; n < N; n = n+1)
    		{     	
    		
    			//Primera fila mag
       			fprintf(mag_x,"%d,",n);
       			fprintf(mag_y,"%d,",n);
       			fprintf(mag_z,"%d,",n);
       			    		
       			//Dentro CI
       			fprintf(CI_mag_x,"%f,",mx[n]);
        		fprintf(CI_mag_y,"%f,",my[n]);
        		fprintf(CI_mag_z,"%f,",mz[n]);
       		}
       			
       			/*===================================*/            
			/*-----------------RK--------------- */   
			/*===================================*/ 
			
    		for (t = 0.0; t <= tf; t = t+dt)
    		{
	
			//K1
        		for (n = ni; n < N; n = n+1)
        		{ 
					
        			mx_rk[n] = mx[n];
        			my_rk[n] = my[n];
       				mz_rk[n] = mz[n];
       				t_rk = t;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       				
			for (n = ni; n < N; n = n+1)
        		{         			
                       		k1x[n] = dt*dMX[n];
        			k1y[n] = dt*dMY[n];
            			k1z[n] = dt*dMZ[n];
        		}
        		            			
            		//k2
            		for (n = ni; n < N; n = n+1)
        		{ 
            		
            		        mx_rk[n] = mx[n] + 0.5*k1x[n];
        			my_rk[n] = my[n] + 0.5*k1y[n];
        			mz_rk[n] = mz[n] + 0.5*k1z[n];
        			t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       			
        		for (n = ni; n < N; n = n+1)
        		{        			
				k2x[n] = dt*dMX[n];
				k2y[n] = dt*dMY[n];
				k2z[n] = dt*dMZ[n];
			}
            			
            		//k3
            		for (n = ni; n < N; n = n+1)
        		{ 
				mx_rk[n] = mx[n]+0.5*k2x[n];
				my_rk[n] = my[n]+0.5*k2y[n];
				mz_rk[n] = mz[n]+0.5*k2z[n];
				t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       			
			for (n = ni; n < N; n = n+1)
 			{        			
				k3x[n] = dt*dMX[n];
				k3y[n] = dt*dMY[n];
				k3z[n] = dt*dMZ[n];
			}
			//k4
			for (n = ni; n < N; n = n+1)
			{ 
				mx_rk[n] = mx[n]+k3x[n];
				my_rk[n] = my[n]+k3y[n];
 				mz_rk[n] = mz[n]+k3z[n];
				t_rk = t+dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       			        			
			for (n = ni; n < N; n = n+1)
			{ 
				k4x[n] = dt*dMX[n];
				k4y[n] = dt*dMY[n];
				k4z[n] = dt*dMZ[n];
			}
        			
			//Final
			for (n = ni; n < N; n = n+1)
 			{ 
				mx[n] = mx[n] + (k1x[n] + 2.0*k2x[n] + 2.0*k3x[n] + k4x[n])/6.0;
				my[n] = my[n] + (k1y[n] + 2.0*k2y[n] + 2.0*k3y[n] + k4y[n])/6.0;
				mz[n] = mz[n] + (k1z[n] + 2.0*k2z[n] + 2.0*k3z[n] + k4z[n])/6.0;
			}            			
     			if(t>600000.0 )
			{
    				fprintf(mag_x,"\n%f,",t+dt);
        			fprintf(mag_y,"\n%f,",t+dt);
        			fprintf(mag_z,"\n%f,",t+dt);
     				for (n = ni; n < N; n = n+1)
     				{
     					/*===========================*/      
    					/*-----------GUARDAR-------- */   
     					/*===========================*/  
     					
		                        fprintf (mag_x, "%f,",mx[n]);
      					fprintf (mag_y, "%f,",my[n]);
      					fprintf (mag_z, "%f,",mz[n]);

				}
			}
			else
			{

			}
			
		}
            	fclose(mag_x);
    		fclose(mag_y);
    		fclose(mag_z);
	}
	
	/*-----------------LO QUE NO SE GUARDA--------------- */   
	
	else
	{
    		for (t = 0.0; t <= tf; t = t+dt)
    		{
    			//K1
        		for (n = ni; n < N; n = n+1)
        		{ 
					
        			mx_rk[n] = mx[n];
        			my_rk[n] = my[n];
        			mz_rk[n] = mz[n];
        			t_rk = t;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       				
			for (n = ni; n < N; n = n+1)
        		{         			
				k1x[n] = dt*dMX[n];
            			k1y[n] = dt*dMY[n];
            			k1z[n] = dt*dMZ[n];
        		}
        			            			
            		//k2
            		for (n = ni; n < N; n = n+1)
        		{ 
            			
                		mx_rk[n] = mx[n] + 0.5*k1x[n];
        			my_rk[n] = my[n] + 0.5*k1y[n];
        			mz_rk[n] = mz[n] + 0.5*k1z[n];
        			t_rk = t + 0.5*dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{        			
            			k2x[n] = dt*dMX[n];
           			k2y[n] = dt*dMY[n]; 
            			k2z[n] = dt*dMZ[n];
            		}
            			
            		//k3
            		for (n = ni; n < N; n = n+1)
        		{ 
            			mx_rk[n] = mx[n]+0.5*k2x[n];
        			my_rk[n] = my[n]+0.5*k2y[n];
        			mz_rk[n] = mz[n]+0.5*k2z[n];
        			t_rk = t + 0.5*dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       			
        		for (n = ni; n < N; n = n+1)
        		{        			
            			k3x[n] = dt*dMX[n];
            			k3y[n] = dt*dMY[n];
            			k3z[n] = dt*dMZ[n];
        		}
        		//k4
        		for (n = ni; n < N; n = n+1)
        		{ 
       				mx_rk[n] = mx[n]+k3x[n];
        			my_rk[n] = my[n]+k3y[n];
        			mz_rk[n] = mz[n]+k3z[n];
        			t_rk = t+dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       			}
       						
            		for (n = ni; n < N; n = n+1)
        		{ 
            			k4x[n] = dt*dMX[n];
            			k4y[n] = dt*dMY[n];
            			k4z[n] = dt*dMZ[n];
        		}
        			
        		//Final
        		for (n = ni; n < N; n = n+1)
        		{ 
            			mx[n] = mx[n] + (k1x[n] + 2.0*k2x[n] + 2.0*k3x[n] + k4x[n])/6.0;
            			my[n] = my[n] + (k1y[n] + 2.0*k2y[n] + 2.0*k3y[n] + k4y[n])/6.0;
            			mz[n] = mz[n] + (k1z[n] + 2.0*k2z[n] + 2.0*k3z[n] + k4z[n])/6.0;
			}            			
		
		}
	}
	printf("Beta1=%f\n",B1);
	B1 = B1 + dB1;	
    }
    fclose(CI_mag_x);
    fclose(CI_mag_y);
    fclose(CI_mag_z);
}

/***********************************************************************/      


/*=====================================*/      
/*-----------LA COSA A INTEGRAR--------*/   
/*=====================================*/ 


//En X
void dM(int n)
{
    	double G,mx2, my2, mz2;
    	double hext_x,hext_y,hext_z,ha_x,ha_y,ha_z, Heff_x, Heff_y, Heff_z;
    	double drx, dry, drz, dr1, dr2, dr5, coef;
    	double rnm, hd_x=0.0 ,hd_y=0.0, hd_z=0.0;// (debe ir fuera)
  
	
    	G=B1*cos(w*t_rk);
    
    	//rx[n_otro]-rx[n] es la diferencia de las posiciones de ambos osciladores
    
    	//Campo externo
   
   	hext_x= h*sin(phi);
   
   	hext_y= 0.0;
   
   	hext_z= h*cos(phi);
   
   	//Anisotropia
   
   	ha_x= 0.0;
   	
	ha_y= 0.0;
   
   	ha_z= (B0z + G)*mz_rk[n];
   
   	//Campo de un m sobre el otro:
   
   
   
   	
   	   // Osciladores después de n
       for(n_otro = ni; n_otro < n; n_otro = n_otro+1)
    	{
            //diferencia de las posiciones de ambos osciladores elevado a 5
            	drx = rx[n_otro]-rx[n];
            	dry = ry[n_otro]-ry[n];
            	drz = rz[n_otro]-rz[n];
                
            	dr2 = drx*drx+dry*dry+drz*drz;
            	dr1 = sqrt(dr2);
                
            	dr5 = dr2*dr2*dr1;
                
              
            	coef = V0/(4.0*pi*dr5);
                
		hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
	
		hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
	
		hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}
    	
    	// Osciladores antes de n	
        for(n_otro = n+1; n_otro < N; n_otro = n_otro+1)
    	{
            //diferencia de las posiciones de ambos osciladores elevado a 5
            	drx = rx[n_otro]-rx[n];
            	dry = ry[n_otro]-ry[n];
            	drz = rz[n_otro]-rz[n];
                    
            	dr2 = drx*drx+dry*dry+drz*drz;
            	dr1 = sqrt(dr2);
            	        
            	dr5 = dr2*dr2*dr1;
            	        
            	coef = V0/(4.0*pi*dr5);

            	hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
        
            	hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
        
            	hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}

   	//Campo efectivo
   
      	Heff_x = hext_x + ha_x + hd_x;
   
    	Heff_y = hext_y + ha_y + hd_y;
   
    	Heff_z = hext_z + ha_z + hd_z;
   
   	//m al cuadrado
   	mx2=mx_rk[n]*mx_rk[n];
    	my2=my_rk[n]*my_rk[n];
    	mz2=mz_rk[n]*mz_rk[n];

    	dMX[n]=Heff_x*(my2 + mz2)*alpha + Heff_y*(mz_rk[n] - mx_rk[n]*my_rk[n]*alpha) - Heff_z*(my_rk[n] + mx_rk[n]*mz_rk[n]*alpha);
    

    	dMY[n]=Heff_y*(mx2 + mz2)*alpha - Heff_x*(mz_rk[n] + mx_rk[n]*my_rk[n]*alpha) +  Heff_z*(mx_rk[n] - my_rk[n]*mz_rk[n]*alpha);

	
    	dMZ[n]= Heff_z*(mx2+my2)*alpha+Heff_x*(my_rk[n]-mx_rk[n]*mz_rk[n]*alpha)-Heff_y*(mx_rk[n]+my_rk[n]*mz_rk[n]*alpha);

}
