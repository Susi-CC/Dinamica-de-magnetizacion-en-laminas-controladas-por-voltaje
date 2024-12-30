#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define		pi      3.1415926535897932384626432132795028841971693993751058

#define		N	6 //Número osciladores 
#define		ni	0 //Inicio del ciclo para los osciladores


//Para guardar los datos

FILE    *LLE,
	*dm_x,
        *dm_y,
        *dm_z;
         
void dDelta(int n); 	//La funcion para Lyapunov
double ddMX[N], ddMY[N], ddMZ[N]; //Las derivadas de las devirvadas
double mx_1[N], my_1[N],mz_1[N];
double mx_i[N], my_i[N],mz_i[N]; //otra magnetización para el calculo de Lyapunov
double delta_mx[N], delta_my[N], delta_mz[N]; //Para Lyapunov
double Norma_deltaM, deltaM2, deltaM2_1; // Para normalizar delta 
double delta_mx_rk[N], delta_my_rk[N], delta_mz_rk[N],k1lx[N],k2lx[N],k3lx[N],k4lx[N],k1ly[N],k2ly[N],k3ly[N],k4ly[N],k1lz[N],k2lz[N],k3lz[N],k4lz[N];

double ExpLyapunov, CoefExp, lle;
int n_e; //Para sacar prom del exponente

        
void dM( int n); 	//La funcion
double dMX[N]; //Las derivadas
double dMY[N];
double dMZ[N];

double k1x[N],k2x[N],k3x[N],k4x[N],k1y[N],k2y[N],k3y[N],k4y[N],k1z[N],k2z[N],k3z[N],k4z[N],mx[N],mx_rk[N],my[N],my_rk[N],mz[N],mz_rk[N],d,t,t_rk,dt,tf; //Declarar las cositas para el RK
//Las cosas con "_rk[n]" se llaman así porque tuve que renombrar esas cosas en las funciones dMX[n], dMY[n] y dMZ[n] para poder hacer el Runge Kutta

double rx[N], ry[N], rz[N]; //Posiciòn de cada oscilador

int n, n_otro; //Para hacer el ciclo de los osciladores, el contador de los osciladores

//Cositas para no guardar tantos datos

long int B1_red; //Beta 1 multiplicado por e4 y redondeado a un entero


long int B1_div = 100; //Cada cuantos B1 guardo beta1

long int t_div = 10; //Cada cuantos t guardo tiempo

//Otros parametros

double  B1, B1F, dB1;

double  w=0.08;                  //  Voltage frequency
double  h=0.1, phi=0.3;               // Parameters asociated to External Fields

double  B0x=0.0,B0z=0.05;              // Normalized constants of Anisotropy and Dedmnetization
double  alpha=0.005; 

double V0 = 1256.64; //Volumen del material

int j = 355;

/***********************************************************************/

int main()
{    
for (j = 300; j < 350; j = j+1)
    {  
    /*===================================*/                  
    /*-------------PARAMETROS----------- */   
    /*===================================*/   
    

//Las cosas del tiempo. 
    dt = 2.0*pi/(w*1000.0);
    tf = 18000.0*2.0*pi/w;   
    
    B1=0.16; // El beta1 inicial
    B1F=0.16; // El beta1 final
    dB1=0.0001; // El paso de beta 1
    
    
/***********************************************************************/   
   
    /*===================================*/                  
    /*--------------GUARDAR------------- */   
    /*===================================*/             
    
    
    char	Name_Mx[150];
    char	Name_My[150];
    char	Name_Mz[150];
    
         /*===================================*/      
    	/*-----------ABRIR CSV CI----------- */   
    	/*===================================*/  
    		
    LLE=fopen("LLE.csv","w");
   		
    //Esto para el titulo del espacio que necesito vacío
    fprintf(LLE,"Beta1,"); 

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
     
  
     srand48(j);
     	/*===================================*/            
	/*-------------CI-deltaM------------ */   
	/*===================================*/    

    /*--------------Mi------------- */   
    for (n = ni; n < N; n = n+1)
    {
        mx_i[n] = -1.0+2.0*drand48();
        my_i[n] = -1.0+2.0*drand48();
        mz_i[n] = -1.0+2.0*drand48();
        
        double Norma=sqrt(mx_i[n]*mx_i[n]+my_i[n]*my_i[n]+mz_i[n]*mz_i[n]);
        
        mx_i[n] = mx_i[n]/Norma;
        my_i[n] = my_i[n]/Norma;
        mz_i[n] = mz_i[n]/Norma;

    }

    /*--------------deltaM------------- */
    
    deltaM2 = 0.0;//El valor inicial
    
    //Le quito la proyección sobre M y calcular norma 
    for (n = ni; n < N; n = n+1)
    {   
    	delta_mx[n] = (mx_i[n]-(mx_i[n]*mx[n]+my_i[n]*my[n]+mz_i[n]*mz[n])*mx[n]);
    	delta_my[n] = (my_i[n]-(mx_i[n]*mx[n]+my_i[n]*my[n]+mz_i[n]*mz[n])*my[n]);
    	delta_mz[n] = (mz_i[n]-(mx_i[n]*mx[n]+my_i[n]*my[n]+mz_i[n]*mz[n])*mz[n]);
    
        deltaM2_1 = delta_mx[n]*delta_mx[n]+delta_my[n]*delta_my[n]+delta_mz[n]*delta_mz[n];
        deltaM2 = deltaM2 + deltaM2_1; //Quiero que se vayan sumando para cada n
    }
    
    Norma_deltaM = sqrt(deltaM2);

    //Normalizar
    
    for(n = ni; n < N; n = n+1)
    {    
    	delta_mx[n] = delta_mx[n]/Norma_deltaM;
    	delta_my[n] = delta_my[n]/Norma_deltaM;
    	delta_mz[n] = delta_mz[n]/Norma_deltaM;
    	
    	printf("\n dm_x = %f",delta_mx[n]);
    }
                                            
    /*--------------POSICIONES------------- */   
    
    
    for (n = ni; n < N; n = n+1)
    {
        rx[n] = 0.0;
        ry[n] = 0.0;
        rz[n] = 50.0*n;
    }

    /*=====================================================*/      
    /*--------------------CICLO DE BETA1-------------------*/   
    /*=====================================================*/ 
    while(B1<=B1F)
    {    
   
    	B1_red = lround(B1*10000);
    	

	if(B1_red % B1_div == 0)
	{
	
		//Primera Columna datos
    		fprintf(LLE,"\n%d,",j);

       		/*===================================*/      
    		/*-----------ABRIR CSV DATOS-------- */   
    		/*===================================*/ 
    		
		sprintf(Name_Mx,"SEED-%d-_B1(e4)_%ld_Delta-MX_.csv",j,B1_red);
		sprintf(Name_My,"SEED-%d-_B1(e4)_%ld_Delta-MY_.csv",j,B1_red);
		sprintf(Name_Mz,"SEED-%d-_B1(e4)_%ld_Delta-MZ_.csv",j,B1_red);

        	dm_x=fopen(Name_Mx,"w");
        	dm_y=fopen(Name_My,"w");
      		dm_z=fopen(Name_Mz,"w");
      		
      		//Esto para el titulo del espacio que necesito vacío
      		fprintf(dm_x,"Oscilador,");
    		fprintf(dm_y,"Oscilador,");
    		fprintf(dm_z,"Oscilador,");
    		
    		
		for (n = ni; n < N; n = n+1)
    		{     	
    			//Primera fila dm
       			fprintf(dm_x,"%d,",n);
       			fprintf(dm_y,"%d,",n);
       			fprintf(dm_z,"%d,",n);
       		}    		
       			/*===================================*/            
			/*-----------------RK--------------- */   
			/*===================================*/ 

		//Iniciar contadores
		CoefExp = 0.0;			 
                ExpLyapunov = 0.0;				
                n_e = 0.0;
		
    		for (t = 0.0; t <= tf; t = t+dt)
    		{
    		
			//K1
        		for (n = ni; n < N; n = n+1)
        		{ 
					
        			mx_rk[n] = mx[n];
        			my_rk[n] = my[n];
       				mz_rk[n] = mz[n];
       				
       				
       				delta_mx_rk[n] = delta_mx[n];			
        			delta_my_rk[n] = delta_my[n];        			
        			delta_mz_rk[n] = delta_mz[n];
        			
       				
       				t_rk = t;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       				
			for (n = ni; n < N; n = n+1)
        		{         			
                       		k1x[n] = dt*dMX[n];
        			k1y[n] = dt*dMY[n];
            			k1z[n] = dt*dMZ[n];
            			
            			k1lx[n] = dt*ddMX[n];
        			k1ly[n] = dt*ddMY[n];
            			k1lz[n] = dt*ddMZ[n];
        			
        		}
        		            			
            		//k2
            		for (n = ni; n < N; n = n+1)
        		{ 
            		
            		        mx_rk[n] = mx[n] + 0.5*k1x[n];
        			my_rk[n] = my[n] + 0.5*k1y[n];
        			mz_rk[n] = mz[n] + 0.5*k1z[n];
        			
        			delta_mx_rk[n] = delta_mx[n]+ 0.5*k1lx[n];			
        			delta_my_rk[n] = delta_my[n]+ 0.5*k1ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ 0.5*k1lz[n];
        			
        			
        			t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				dDelta(n);
        			
       			}
       			
        		for (n = ni; n < N; n = n+1)
        		{        			
				k2x[n] = dt*dMX[n];
				k2y[n] = dt*dMY[n];
				k2z[n] = dt*dMZ[n];
				
				k2lx[n] = dt*ddMX[n];
        			k2ly[n] = dt*ddMY[n];
            			k2lz[n] = dt*ddMZ[n];
        			
				
			}
            			
            		//k3
            		for (n = ni; n < N; n = n+1)
        		{ 
				mx_rk[n] = mx[n]+0.5*k2x[n];
				my_rk[n] = my[n]+0.5*k2y[n];
				mz_rk[n] = mz[n]+0.5*k2z[n];
				
				delta_mx_rk[n] = delta_mx[n]+ 0.5*k2lx[n];			
        			delta_my_rk[n] = delta_my[n]+ 0.5*k2ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ 0.5*k2lz[n];
        			
				
				t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       			
			for (n = ni; n < N; n = n+1)
 			{        			
				k3x[n] = dt*dMX[n];
				k3y[n] = dt*dMY[n];
				k3z[n] = dt*dMZ[n];
				
				k3lx[n] = dt*ddMX[n];
        			k3ly[n] = dt*ddMY[n];
            			k3lz[n] = dt*ddMZ[n];
        			
			}
			//k4
			for (n = ni; n < N; n = n+1)
			{ 
				mx_rk[n] = mx[n]+k3x[n];
				my_rk[n] = my[n]+k3y[n];
 				mz_rk[n] = mz[n]+k3z[n];
 				
 				delta_mx_rk[n] = delta_mx[n]+ k3lx[n];			
        			delta_my_rk[n] = delta_my[n]+ k3ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ k3lz[n];
        			
 				
				t_rk = t+dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       			        			
			for (n = ni; n < N; n = n+1)
			{ 
				k4x[n] = dt*dMX[n];
				k4y[n] = dt*dMY[n];
				k4z[n] = dt*dMZ[n];
				
				k4lx[n] = dt*ddMX[n];
        			k4ly[n] = dt*ddMY[n];
            			k4lz[n] = dt*ddMZ[n];
        			
			}
        			
			//Final
			for (n = ni; n < N; n = n+1)
 			{ 
				mx[n] = mx[n] + (k1x[n] + 2.0*k2x[n] + 2.0*k3x[n] + k4x[n])/6.0;
				my[n] = my[n] + (k1y[n] + 2.0*k2y[n] + 2.0*k3y[n] + k4y[n])/6.0;
				mz[n] = mz[n] + (k1z[n] + 2.0*k2z[n] + 2.0*k3z[n] + k4z[n])/6.0;
				
				delta_mx[n] = delta_mx[n] + (k1lx[n] + 2.0*k2lx[n] + 2.0*k3lx[n] + k4lx[n])/6.0;
				delta_my[n] = delta_my[n] + (k1ly[n] + 2.0*k2ly[n] + 2.0*k3ly[n] + k4ly[n])/6.0;
				delta_mz[n] = delta_mz[n] + (k1lz[n] + 2.0*k2lz[n] + 2.0*k3lz[n] + k4lz[n])/6.0;
				
			}   

			if(t<600000.0)
			{
			        //Le quito la proyección sobre M y calcular norma 
   				deltaM2 = 0.0;//El valor inicial
				for (n = ni; n < N; n = n+1)
    				{   
    					delta_mx[n] = (delta_mx[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mx[n]);
    					delta_my[n] = (delta_my[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*my[n]);
    					delta_mz[n] = (delta_mz[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mz[n]);
    
        				deltaM2_1 = delta_mx[n]*delta_mx[n]+delta_my[n]*delta_my[n]+delta_mz[n]*delta_mz[n];
        				deltaM2 = deltaM2 + deltaM2_1; //Quiero que se vayan sumando para cada n
    				}
    
   				
   				Norma_deltaM = sqrt(deltaM2);

    				//Normalizar
    				for(n = ni; n < N; n = n+1)
    				{    
    					delta_mx[n] = delta_mx[n]/Norma_deltaM;
    					delta_my[n] = delta_my[n]/Norma_deltaM;
    					delta_mz[n] = delta_mz[n]/Norma_deltaM;
    				}

			}	
     			
			else
			{
                		
                		//Le quito la proyección sobre M y calcular norma 
   				deltaM2 = 0.0;//El valor inicial
				for (n = ni; n < N; n = n+1)
    				{   
    					delta_mx[n] = (delta_mx[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mx[n]);
    					delta_my[n] = (delta_my[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*my[n]);
    					delta_mz[n] = (delta_mz[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mz[n]);
    
        				deltaM2_1 = delta_mx[n]*delta_mx[n]+delta_my[n]*delta_my[n]+delta_mz[n]*delta_mz[n];
        				deltaM2 = deltaM2 + deltaM2_1; //Quiero que se vayan sumando para cada n
    				}
    
   				
   				Norma_deltaM = sqrt(deltaM2);
                		
                /*/----------------------------Calculo del exponente-----------------------------/*/   
                /*/ 				       /*/					 /*/
                /*/					 					 /*/
        	/*/		CoefExp = (1/(dt)) * log(Norma_deltaM);			 	 /*/
                /*/		ExpLyapunov = ExpLyapunov + CoefExp;				 /*/
                /*/		n_e = n_e + 1;							 /*/
            	/*/										 /*/
            	/*////////////////////////////////////////////////////////////////////////////////*/   
				
				if(lround(t*1000*w/(2*pi))% t_div == 0)
				{
					fprintf(dm_x,"\n%f,",t+dt);
        				fprintf(dm_y,"\n%f,",t+dt);
        				fprintf(dm_z,"\n%f,",t+dt);
        				
     					for (n = ni; n < N; n = n+1)
     					{
     						/*===========================*/      
    						/*-----------GUARDAR-------- */   
     						/*===========================*/  
     					
		                        	fprintf (dm_x, "%f,",delta_mx[n]);
      						fprintf (dm_y, "%f,",delta_my[n]);
      						fprintf (dm_z, "%f,",delta_mz[n]);

					}
					
				//Normalizar
    				for(n = ni; n < N; n = n+1)
    				{    
    					delta_mx[n] = delta_mx[n]/Norma_deltaM;
    					delta_my[n] = delta_my[n]/Norma_deltaM;
    					delta_mz[n] = delta_mz[n]/Norma_deltaM;
    				}
    				
				}
				else
				{
    				//Normalizar
    				for(n = ni; n < N; n = n+1)
    				{    
    					delta_mx[n] = delta_mx[n]/Norma_deltaM;
    					delta_my[n] = delta_my[n]/Norma_deltaM;
    					delta_mz[n] = delta_mz[n]/Norma_deltaM;
    				}
    				}
		
			}
			
		}
	lle = ExpLyapunov/n_e;
       	fprintf(LLE,"%f,",lle);
       	printf("Lyapunov=%f\n", lle);
       	
       	fclose(dm_x);
    	fclose(dm_y);
    	fclose(dm_z);
	
	}
	
	/*-----------------NO SE CALCULA EXPONENTE--------------- */   
	
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
       				
       				
       				delta_mx_rk[n] = delta_mx[n];			
        			delta_my_rk[n] = delta_my[n];        			
        			delta_mz_rk[n] = delta_mz[n];
        			
       				
       				t_rk = t;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       				
			for (n = ni; n < N; n = n+1)
        		{         			
                       		k1x[n] = dt*dMX[n];
        			k1y[n] = dt*dMY[n];
            			k1z[n] = dt*dMZ[n];
            			
            			k1lx[n] = dt*ddMX[n];
        			k1ly[n] = dt*ddMY[n];
            			k1lz[n] = dt*ddMZ[n];
        			
        		}
        		            			
            		//k2
            		for (n = ni; n < N; n = n+1)
        		{ 
            		
            		        mx_rk[n] = mx[n] + 0.5*k1x[n];
        			my_rk[n] = my[n] + 0.5*k1y[n];
        			mz_rk[n] = mz[n] + 0.5*k1z[n];
        			
        			delta_mx_rk[n] = delta_mx[n]+ 0.5*k1lx[n];			
        			delta_my_rk[n] = delta_my[n]+ 0.5*k1ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ 0.5*k1lz[n];
        			
        			
        			t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				dDelta(n);
        			
       			}
       			
        		for (n = ni; n < N; n = n+1)
        		{        			
				k2x[n] = dt*dMX[n];
				k2y[n] = dt*dMY[n];
				k2z[n] = dt*dMZ[n];
				
				k2lx[n] = dt*ddMX[n];
        			k2ly[n] = dt*ddMY[n];
            			k2lz[n] = dt*ddMZ[n];
        			
				
			}
            			
            		//k3
            		for (n = ni; n < N; n = n+1)
        		{ 
				mx_rk[n] = mx[n]+0.5*k2x[n];
				my_rk[n] = my[n]+0.5*k2y[n];
				mz_rk[n] = mz[n]+0.5*k2z[n];
				
				delta_mx_rk[n] = delta_mx[n]+ 0.5*k2lx[n];			
        			delta_my_rk[n] = delta_my[n]+ 0.5*k2ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ 0.5*k2lz[n];
        			
				
				t_rk = t + 0.5*dt;
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       			
			for (n = ni; n < N; n = n+1)
 			{        			
				k3x[n] = dt*dMX[n];
				k3y[n] = dt*dMY[n];
				k3z[n] = dt*dMZ[n];
				
				k3lx[n] = dt*ddMX[n];
        			k3ly[n] = dt*ddMY[n];
            			k3lz[n] = dt*ddMZ[n];
        			
			}
			//k4
			for (n = ni; n < N; n = n+1)
			{ 
				mx_rk[n] = mx[n]+k3x[n];
				my_rk[n] = my[n]+k3y[n];
 				mz_rk[n] = mz[n]+k3z[n];
 				
 				delta_mx_rk[n] = delta_mx[n]+ k3lx[n];			
        			delta_my_rk[n] = delta_my[n]+ k3ly[n];        			
        			delta_mz_rk[n] = delta_mz[n]+ k3lz[n];
        			
 				
				t_rk = t+dt;
        			
       			}
       			
       			for (n = ni; n < N; n = n+1)
        		{ 

       				dM(n);
       				
       				dDelta(n);
        			
       			}
       			        			
			for (n = ni; n < N; n = n+1)
			{ 
				k4x[n] = dt*dMX[n];
				k4y[n] = dt*dMY[n];
				k4z[n] = dt*dMZ[n];
				
				k4lx[n] = dt*ddMX[n];
        			k4ly[n] = dt*ddMY[n];
            			k4lz[n] = dt*ddMZ[n];
        			
			}
        			
			//Final
			for (n = ni; n < N; n = n+1)
 			{ 
				mx[n] = mx[n] + (k1x[n] + 2.0*k2x[n] + 2.0*k3x[n] + k4x[n])/6.0;
				my[n] = my[n] + (k1y[n] + 2.0*k2y[n] + 2.0*k3y[n] + k4y[n])/6.0;
				mz[n] = mz[n] + (k1z[n] + 2.0*k2z[n] + 2.0*k3z[n] + k4z[n])/6.0;
				
				delta_mx[n] = delta_mx[n] + (k1lx[n] + 2.0*k2lx[n] + 2.0*k3lx[n] + k4lx[n])/6.0;
				delta_my[n] = delta_my[n] + (k1ly[n] + 2.0*k2ly[n] + 2.0*k3ly[n] + k4ly[n])/6.0;
				delta_mz[n] = delta_mz[n] + (k1lz[n] + 2.0*k2lz[n] + 2.0*k3lz[n] + k4lz[n])/6.0;
				
			}   

			//Le quito la proyección sobre M y calcular norma 
   			deltaM2 = 0.0;//El valor inicial
			for (n = ni; n < N; n = n+1)
    			{   
    				delta_mx[n] = (delta_mx[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mx[n]);
    				delta_my[n] = (delta_my[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*my[n]);
    				delta_mz[n] = (delta_mz[n]-(delta_mx[n]*mx[n]+delta_my[n]*my[n]+delta_mz[n]*mz[n])*mz[n]);
    
        			deltaM2_1 = delta_mx[n]*delta_mx[n]+delta_my[n]*delta_my[n]+delta_mz[n]*delta_mz[n];
        			deltaM2 = deltaM2 + deltaM2_1; //Quiero que se vayan sumando para cada n
    			}
    
   				
   			Norma_deltaM = sqrt(deltaM2);

    			//Normalizar
    			for(n = ni; n < N; n = n+1)
    			{    
    				delta_mx[n] = delta_mx[n]/Norma_deltaM;
    				delta_my[n] = delta_my[n]/Norma_deltaM;
    				delta_mz[n] = delta_mz[n]/Norma_deltaM;
    			}  
		}
	}
	
	B1 = B1 + dB1;	
    }
	printf("j=%d\n",j);
}
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
    
    /*--------------CAMPO EXTERNO------------- */   
   
   	hext_x= h*sin(phi);
   
   	hext_y= 0.0;
   
   	hext_z= h*cos(phi);
   
    /*--------------ANISOTROPIA------------- */   
   
   	ha_x= 0.0;
   	
	ha_y= 0.0;
   
   	ha_z= (B0z + G)*mz_rk[n];
   
    /*--------------CAMPO DE UN OSC SOBRE EL OTRO------------- */   
   
	// Osciladores antes de n
       for(n_otro = ni; n_otro < n; n_otro = n_otro+1)
    	{
            
            	drx = rx[n_otro]-rx[n];
            	dry = ry[n_otro]-ry[n];
            	drz = rz[n_otro]-rz[n];
                
            	dr2 = drx*drx+dry*dry+drz*drz;
            	dr1 = sqrt(dr2);
                
            	dr5 = dr2*dr2*dr1;//diferencia de las posiciones de ambos osciladores elevado a 5
                
              
            	coef = V0/(4.0*pi*dr5);
                
		hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
	
		hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
	
		hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}
    	
    	// Osciladores después de n	
        for(n_otro = n+1; n_otro < N; n_otro = n_otro+1)
    	{
            
            	drx = rx[n_otro]-rx[n];
            	dry = ry[n_otro]-ry[n];
            	drz = rz[n_otro]-rz[n];
                    
            	dr2 = drx*drx+dry*dry+drz*drz;
            	dr1 = sqrt(dr2);
            	        
            	dr5 = dr2*dr2*dr1; //diferencia de las posiciones de ambos osciladores elevado a 5
            	        
            	coef = V0/(4.0*pi*dr5);

            	hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
        
            	hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
        
            	hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}

    /*--------------CAMPO EFECTIVO------------- */   
   
      	Heff_x = hext_x + ha_x + hd_x;
   
    	Heff_y = hext_y + ha_y + hd_y;
   
    	Heff_z = hext_z + ha_z + hd_z;
   
    /*--------------M AL CUADRADO------------- */   
   	mx2=mx_rk[n]*mx_rk[n];
    	my2=my_rk[n]*my_rk[n];
    	mz2=mz_rk[n]*mz_rk[n];
    	
    /*--------------DM------------- */   
    	dMX[n]=Heff_x*(my2 + mz2)*alpha + Heff_y*(mz_rk[n] - mx_rk[n]*my_rk[n]*alpha) - Heff_z*(my_rk[n] + mx_rk[n]*mz_rk[n]*alpha);
    

    	dMY[n]=Heff_y*(mx2 + mz2)*alpha - Heff_x*(mz_rk[n] + mx_rk[n]*my_rk[n]*alpha) +  Heff_z*(mx_rk[n] - my_rk[n]*mz_rk[n]*alpha);

	
    	dMZ[n]= Heff_z*(mx2+my2)*alpha+Heff_x*(my_rk[n]-mx_rk[n]*mz_rk[n]*alpha)-Heff_y*(mx_rk[n]+my_rk[n]*mz_rk[n]*alpha);

}


/*=====================================*/      
/*--------------Jacobiano--------------*/   
/*=====================================*/ 

void dDelta(int n)
{
	double G;
	double hext_x,hext_y,hext_z,ha_x,ha_y,ha_z, Heff_x, Heff_y, Heff_z;
	double drx, dry, drz, dr1, dr2, dr5, coef;
	double rnm, hd_x ,hd_y, hd_z;
	double Int11[N], Int12[N], Int13[N], Int21[N], Int22[N], Int23[N], Int31[N], Int32[N], Int33[N];
	int n_dm; //Contador derivada
	int n_m; //Contador dm
	
	
	G=B1*cos(w*t_rk);
    
    	//rx[n_otro]-rx[n] es la diferencia de las posiciones de ambos osciladores
    
    /*--------------CAMPO EXTERNO------------- */   
   
   	hext_x= h*sin(phi);
   
   	hext_y= 0.0;
   
   	hext_z= h*cos(phi);
   
    /*--------------ANISOTROPIA------------- */   
   	
   	ha_x= 0.0;
   	
	ha_y= 0.0;
   
   	ha_z= (B0z + G)*mz_rk[n];
   
    	/*--------------CAMPO DE UN OSC SOBRE EL OTRO------------- */   
	//Iniciar el contador para cada n   		
	hd_x=0.0;
	hd_y=0.0;
	hd_z=0.0;

	// Osciladores después de n
       	for(n_otro = ni; n_otro < n; n_otro = n_otro+1)
    	{
            
            		drx = rx[n_otro]-rx[n];
            		dry = ry[n_otro]-ry[n];
            		drz = rz[n_otro]-rz[n];
                
            		dr2 = drx*drx+dry*dry+drz*drz;
            		dr1 = sqrt(dr2);
                
            		dr5 = dr2*dr2*dr1;//diferencia de las posiciones de ambos osciladores elevado a 5

              
            		coef = V0/(4.0*pi*dr5);
                
			hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
	
			hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
	
			hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}
    	
    		// Osciladores antes de n	
        for(n_otro = n+1; n_otro < N; n_otro = n_otro+1)
    	{
            
            		drx = rx[n_otro]-rx[n];
            		dry = ry[n_otro]-ry[n];
            		drz = rz[n_otro]-rz[n];
                    
            		dr2 = drx*drx+dry*dry+drz*drz;
            		dr1 = sqrt(dr2);
            	        
            		dr5 = dr2*dr2*dr1; //diferencia de las posiciones de ambos osciladores elevado a 5
            	        
            		coef = V0/(4.0*pi*dr5);

            		hd_x+= coef*(3.0*( drx )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mx_rk[n_otro]*dr2);
	
			hd_y+= coef*(3.0*( dry )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- my_rk[n_otro]*dr2);
	
			hd_z+= coef*(3.0*( drz )*( (mx_rk[n_otro]*(drx)) + (my_rk[n_otro]*(dry)) + (mz_rk[n_otro]*(drz)) )- mz_rk[n_otro]*dr2);

    	}

    	/*--------------CAMPO EFECTIVO------------- */   
   
      	Heff_x = hext_x + ha_x + hd_x;
   
    	Heff_y = hext_y + ha_y + hd_y;
   
    	Heff_z = hext_z + ha_z + hd_z;
    	
    	
    	/*---------------JACOBIANO------------- */       
    	//matriz de 18 x18-> En Mathematica
	
	
	//La interacción:
	
	//Primero el valor inicial para que no sume eternamente
		
	Int11[n] = 0.0 ;
	Int12[n] = 0.0 ;
	Int13[n] = 0.0 ;
	Int21[n] = 0.0 ;
	Int22[n] = 0.0 ;
	Int23[n] = 0.0 ;
	Int31[n] = 0.0 ;
	Int32[n] = 0.0 ;
	Int33[n] = 0.0 ;

			// Osciladores antes de n
       	for(n_otro = ni; n_otro < n; n_otro = n_otro+1)
    	{
			drx = rx[n_otro]-rx[n];
            		dry = ry[n_otro]-ry[n];
            		drz = rz[n_otro]-rz[n];
                
            		dr2 = drx*drx+dry*dry+drz*drz;
            		dr1 = sqrt(dr2);
                
            		dr5 = dr2*dr2*dr1;//diferencia de las posiciones de ambos osciladores elevado a 5

              
            		coef = V0/(4.0*pi*dr5);

			//"Int" es la parte de la interacción de los osciladores para el ddm

			Int11[n] += delta_mx_rk[n_otro] * -coef * alpha* (my_rk[n]*my_rk[n]*dr2 + mz_rk[n]*mz_rk[n]*dr2);
			Int12[n] += delta_my_rk[n_otro] * coef * (alpha*mx_rk[n]*my_rk[n]*dr2-mz_rk[n]*dr2);
			Int13[n] += delta_mz_rk[n_otro] * coef * (my_rk[n]*dr2 - 3.0*my_rk[n]*drz*drz - alpha* (-mx_rk[n]*mz_rk[n]*dr2 + 3.0*mx_rk[n] *mz_rk[n]*drz*drz));

			
			Int21[n] += delta_mx_rk[n_otro] * coef * (alpha*mx_rk[n]*my_rk[n]*dr2 + mz_rk[n]*dr2);
			Int22[n] += delta_my_rk[n_otro] * -coef * alpha * (mx_rk[n]*mx_rk[n]*dr2+mz_rk[n]*mz_rk[n]*dr2);
			Int23[n] += delta_mz_rk[n_otro] * coef * (-mx_rk[n]*dr2 + 3.0* mx_rk[n] *drz*drz - alpha*(- my_rk[n]*mz_rk[n]*dr2 + 3.0*my_rk[n]*mz_rk[n]*drz*drz));

			
			Int31[n] += delta_mx_rk[n_otro] * coef * (-my_rk[n]*dr2 + alpha*mx_rk[n]*mz_rk[n]*dr2);
			Int32[n] += delta_my_rk[n_otro] * coef * (mx_rk[n]*dr2 + alpha*my_rk[n]*mz_rk[n]*dr2);
			Int33[n] += delta_mz_rk[n_otro] * -coef * alpha* (mx_rk[n]*mx_rk[n]*dr2 + my_rk[n]*my_rk[n]*dr2 - 3.0*mx_rk[n]*mx_rk[n]*drz*drz - 3.0*my_rk[n]*my_rk[n]*drz*drz);

	}

	//Osciladores despues de n
        for(n_otro = n+1; n_otro < N; n_otro = n_otro+1)
    	{
			drx = rx[n_otro]-rx[n];
            		dry = ry[n_otro]-ry[n];
            		drz = rz[n_otro]-rz[n];
                
            		dr2 = drx*drx+dry*dry+drz*drz;
            		dr1 = sqrt(dr2);
                
            		dr5 = dr2*dr2*dr1;//diferencia de las posiciones de ambos osciladores elevado a 5

              
            		coef = V0/(4.0*pi*dr5);

			//"Int" es la parte de la interacción de los osciladores para el ddm

			Int11[n] += delta_mx_rk[n_otro] * -coef * alpha* (my_rk[n]*my_rk[n]*dr2 + mz_rk[n]*mz_rk[n]*dr2);
			Int12[n] += delta_my_rk[n_otro] * coef * (alpha*mx_rk[n]*my_rk[n]*dr2 - mz_rk[n]*dr2);
			Int13[n] += delta_mz_rk[n_otro] * coef * (my_rk[n]*dr2 - 3.0*my_rk[n]*drz*drz - alpha* (-mx_rk[n]*mz_rk[n]*dr2 + 3.0*mx_rk[n] *mz_rk[n]*drz*drz));

			
			Int21[n] += delta_mx_rk[n_otro] * coef * (alpha*mx_rk[n]*my_rk[n]*dr2 + mz_rk[n]*dr2);
			Int22[n] += delta_my_rk[n_otro] * -coef * alpha * (mx_rk[n]*mx_rk[n]*dr2+mz_rk[n]*mz_rk[n]*dr2);
			Int23[n] += delta_mz_rk[n_otro] * coef * (-mx_rk[n]*dr2 + 3.0* mx_rk[n] *drz*drz - alpha*(- my_rk[n]*mz_rk[n]*dr2 + 3.0*my_rk[n]*mz_rk[n]*drz*drz));

			
			Int31[n] += delta_mx_rk[n_otro] * coef * (-my_rk[n]*dr2 + alpha*mx_rk[n]*mz_rk[n]*dr2);
			Int32[n] += delta_my_rk[n_otro] * coef * (mx_rk[n]*dr2 + alpha*my_rk[n]*mz_rk[n]*dr2);
			Int33[n] += delta_mz_rk[n_otro] * -coef * alpha* (mx_rk[n]*mx_rk[n]*dr2 + my_rk[n]*my_rk[n]*dr2 - 3.0*mx_rk[n]*mx_rk[n]*drz*drz - 3.0*my_rk[n]*my_rk[n]*drz*drz);

	}

       /*--------------Derivada de dM------------- */   
	

    	ddMX[n] = delta_mx_rk[n]*(-alpha*( Heff_y *my_rk[n] + Heff_z*mz_rk[n])) + delta_my_rk[n] * (- Heff_z - alpha*( mx_rk[n] * Heff_y + 2.0 * my_rk[n] * Heff_x )) + delta_mz_rk[n] * (Heff_y - my_rk[n] * (B0z + G) - 
  alpha* ( mx_rk[n] *( hd_z + hext_z ) - 2.0*mz_rk[n]*Heff_x + 2.0*mx_rk[n]*mz_rk[n]*(B0z + G)))+ Int11[n] + Int12[n] + Int13[n];
 
    	ddMY[n] = delta_mx_rk[n]*(Heff_z - 
  alpha * (-2.0*mx_rk[n]*Heff_y + my_rk[n] * Heff_x)) + delta_my_rk[n]* - alpha * (mx_rk[n] * Heff_x + mz_rk[n]*Heff_z) + delta_mz_rk[n]*(-Heff_x + mx_rk[n] *(B0z + G ) - 
 alpha*( my_rk[n] * (hd_z  + hext_z) - 2.0*mz_rk[n]*Heff_y + 2.0* (B0z + G)* my_rk[n]*mz_rk[n])) +  Int21[n] + Int22[n] + Int23[n];

    	ddMZ[n] = delta_mx_rk[n]*(-Heff_y - alpha*(-2.0 * mx_rk[n] * Heff_z + mz_rk[n] * Heff_x)) + delta_my_rk[n]*(Heff_x - alpha*(-2.0*my_rk[n]* Heff_z + mz_rk[n]*Heff_y)) + delta_mz_rk[n]*(-alpha*(mx_rk[n] * Heff_x - mx_rk[n]*mx_rk[n]*(B0z + G) + my_rk[n]*Heff_y - my_rk[n]*my_rk[n]*(B0z + G ))) +  Int31[n] + Int32[n] + Int33[n];

    	
}

