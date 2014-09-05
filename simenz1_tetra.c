/* 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	  

   simen1.c
   Una red en las que hay varios grupos de enzimas trabajando en paralelo.
   * Las enzimas son recorridas aleatoriamente para su actualizacion; de esta 
   manera no se favorece a ninguna reaccion. 
   * El recorrido de las enzimas es realizado leyendo secuencialmente las 
     posiciones (reaccion-enzima) contenidas en un vector previamente generado ad hoc. 
   * Cada substancia tiene un TMAX distinto (afinidad) dependiendo de la direccion 
     de la reaccion en la que esta involucrada. 
   * Cada reaccion se define una vez, el codigo automaticamente toma en cuenta la reaccion inversa.
   * La secuencia de las reacciones fue comprobada
   * El valor x = 0.5 fue excluido de los mapas debido a problemas a bajas probabilidades
     (ver histograma del mapa tienda)
   * Los parametros (Vmax, K) de las enzimas se leen desde un archivo
   * La estructura del programa es modular: La actualizacion de las enzimas y el proceso de toma de sustrato
     y liberacion de producto se realizan desde una funcion (8/12/2005)
   * La cantidad de enzimas varia de una enzima a otra (9/12/2005)
   * Se sigue idea de A Parravano para el arranque de una enzima en direccion directa o reversa
   * Se modifico rutina "mapa_s"
   * Se incluyen oligomeros
   * Se incluye re-definicion de las afinidades para oligomeros
   * Una unidad de substancia = 1 microMolar
   * Se incluye re-definicion de las velocidades para oligomeros (13/02/06)
   * Se cambio la estructura de datos de las afinidades (22/02/06)
   * La cooperatividad de los oligomeros se define individualmente para cada enzima
   * Las velocidades son multiplicadas por DT antes de realizarse cualquier cálculo (09/05/2006)
   * Se corrigio problema con avance de fase en oligomeros (17/06/2006)

   * Probabilidad de salida en mapa I fija
   * Calibracion de las probabilidades, curva teorica hallada con Runge-Kutta, i. competitiva
     pin=7.0*I*vmax[pos1][0]/(1.0*ki[inh[pos1][0]]);	
     pout=.02; 

   * Calibracion de las probabilidades, curva teorica hallada con Runge-Kutta, i. no competitiva
     INHIBICION COMPETITIVA   
     pin=0.75*I*vmax[pos1][0]/ki[inh[pos1][2]];	
     pout=.02; 
     INHIBICION NO COMPETITIVA
     pin=0.75*I*vmax[pos1][0]/ki[inh[pos1][2]];	
     pout=.02; 

   * Se incluye tiempo de recuperacion, error corregido  
   * Version operativa, probabilidades para inhibidores dependen de n[inhib[x][x]]
   * Se incluyen tiempos de recuperacion

   * Explorando sincronizacion!

   * NO SE CONSIDERA INHIBICION EN OLIGOMEROS!
   * Corregido problema con signo en funcion VERIFINV  
   * Corregido problema en SUSTPROD para inhibicion con reaccion inversa
   * Incluye animacion
   * SE CORRIGIO PROBLEMA CON REACCIONES UNIDIRECCIONALES
   * DOS SUBSTANCIAS, DOS REACCIONES!
   * Se suprime salida gnuplot
   * Los parametros son escalados para enzimas individuales
   * Se inyecta sustrato de forma exponencial
   * Se lleva contabilidad de enzimas en espera, procesando y recuperandose
   
   * ESTADO = 0,1: Enzima en espera
   * ESTADO = 2: Enzima procesando producto
   * ESTADO = 3: Enzima en recuperacion
     
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<time.h>
//#include "cpgplot.h"

#define NSUB 17
#define NREAC 10
#define TOP 50		// tiempo de simulacion, en multiplos de tau
#define NENZ 50		// Numero de enzimas
#define DT 1.e-6	// Incremento de tiempo correspondiente a una iteracion (segundos)
#define RES 100		// Resolucion del histograma
#define NA 6.022e23	// Numero de Avogadro
#define VSIM 10.3e-18 	// Volumen de simulacion (litros)
#define NINJ 200.00 	// Concentracion final de sustrato (microMolar)

// Variables globales
FILE *dptr;
double p1[2][NREAC]={0.0},vmax[NREAC][2]={0.0},km[2*NREAC][NSUB]={0.0},bdir[NREAC]={0.0},binv[NREAC]={0.0};
double phicdir[NREAC]={0.0},phicinv[NREAC]={0.0},a,b,p2,ki[NSUB]={0.0},coop[NENZ]={0.0},trec[NREAC][2]={0.0};
float fas[RES]={0.0},hh[RES]={0.0};
int eps[NREAC][NSUB]={0},n[NSUB]={0},tenz,zz;  
int nenz[NREAC]={0},tipenz[NREAC]={0},inh[NREAC][5]={-1};
int *cont;
int ***estado;
int **act;
double ***fase;

void rand_init();
void baraja();
void probab();
void sustprod();
void hist();
void whist();
double mapa_s();
double mapa_i();
double enzima47();
double avandir();
double avaninv();
double avaninh();
int verifdir();
int verifinv();
int drnd();

main()
{
   // Declaraciones 
   FILE *dptr;
   FILE *fhist;
   //FILE *dptt;
   //FILE *dptq;
   int i,j=0,k=0,z,reac,enz,prod=0,m,cdir,cinv,index1,index2,gof,gob,cuenta=0,flag,ntop,A_acc,delta,flag2=0;
   double fa1,fa2,fa3,fb1,fb2,fb3,f0,f1,f2,f3,f4,km0;
   double prods,B,prodp,tiempo=0.0,deltat,dado,tau,sigma,Aext_prev,Aext;
   double mm[NREAC] = {0.0};
   
   rand_init();

   // Abre archivos de salida
   if((dptr=fopen("discrete_tetra.dat","w"))==NULL)
      printf("\nEl archivo de salida no puede ser abierto\n");

   if((fhist=fopen("histo.dat","w"))==NULL)
      printf("\nEl archivo de salida no puede ser abierto\n");

   // Inicia pgplot
   //cpgopen("/XSERVE");	
   //cpgask(0);

   // ****************************** DEFINICION DE PARAMETROS *******************************************
 
   a = 0.1;
   p2 = 0.99;

   // Velocidades y afinidades
   // IMPORTANTE: Las afinidades de los sustratos deben ser ingresadas en el mismo orden en el que los
   // sustratos aparecen en la lista de las sustancias!
   
   // PESOS MOLARES: en Daltons
   mm[0] = 5000.;
   mm[1] = 5000.;
   
   // VELOCIDADES: micromol(mg-min)
   vmax[0][0] = 500.0;		// Vmax reaccion directa
   vmax[0][1] = -1.;		// Vmax reaccion inversa
   vmax[1][0] = 500.0;		// Vmax reaccion directa
   vmax[1][1] = -1.;		// Vmax reaccion inversa

   for(j=0; j<NREAC; j++)
   {
      vmax[j][0] *= 1.e-3*mm[j]/60.0;		
      vmax[j][1] *= 1.e-3*mm[j]/60.0;		
      printf("\nLa velocidad de la reaccion %d directa es %lf [sustratos/segundo-enzima]", j, vmax[j][0]);
      printf("\nLa velocidad de la reaccion %d inversa es %lf [sustratos/segundo-enzima]", j, vmax[j][1]);
   }
   
   tau = 1./vmax[0][0];

   printf("\nEl tau es %lf segundos", tau);
   
   printf("\n");

   // AFINIDADES: microMolar	
   km[0][0] = 16.25;		
   km[2][1] = 16.25;
   km0=km[0][0];

   for(j=0; j<2*NREAC; j++)
   {
      for(k=0; k<NSUB; k++)	   
      {	      
         km[j][k] *= NA*1.e-6*VSIM;	      
	 if(km[j][k] != 0.0)
	    printf("\nks[%d][%d] es %lf [sustratos]", j, k, km[j][k]);
      }	 
   }

   printf("\n");

   // TIEMPOS DE RECUPERACION
   trec[0][0] = 0.01;		// Reaccion directa
   trec[0][1] = 1.;		// Reaccion inversa
   trec[1][0] = 0.01;		// Reaccion directa
   trec[1][1] = 1.;		// Reaccion inversa

   // Numero de enzimas que catalizan cada reaccion
   nenz[0] = NENZ;
   nenz[1] = NENZ;

   // Tipo de enzima que cataliza cada reaccion
   tipenz[0] = 4;
   tipenz[1] = 4;

   // Cooperatividad de los oligomeros
   coop[0] = 4.5;
   coop[1] = 4.5;

   printf("\nalpha = %lf",coop[0]);
   printf("\nn_1/2 = %lf",km[0][0]*pow(coop[0],1.5));
   printf("\nS_1/2 = %lf",km0*pow(coop[0],1.5));
   printf("\n");

   // Definicion de inhibicion
   inh[0][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[0][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[0][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[0][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[0][4] = -1;		// Total

   inh[1][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[1][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[1][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[1][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[1][4] = -1;		// Total

   // La reaccion es inhibida?
   for(j=0; j<NREAC; j++)
   {
      inh[j][4] = 0;   
      for(k=0; k<4; k++)
         inh[j][4] += inh[j][k];
   }

   // Afinidades para inhibidores
   ki[0] = 100.0;
   ki[1] = 100.0;
   
   for(j=0; j<NSUB; j++)
      ki[j] *= NA*1.e-6*VSIM;	      

   // Halla numero total de enzimas
   tenz = 0; 
   for(j=0; j<NREAC; j++)
      tenz += nenz[j];

   // Cantidades iniciales de sustrato 
   n[0] = 0;		// S0	// Cantidades minimas para completar la reaccion
   n[1] = 0;		// S1

   // Se definen reacciones
   eps[0][0] = -1;  
   eps[0][1] = 1;  
   eps[1][0] = 1;  
   eps[1][1] = -1;  

   // ************************* RESERVA DE MEMORIA *********************************************

   // Reserva memoria para vector "cont"
   // Tamaño = #tot enz
   if((cont=(int*)calloc(tenz,sizeof(int)))==NULL)
   {
      printf("\nERROR: No hay memoria para el vector 'cont'!");
      exit(1);
   }

   // Reserva memoria para la matriz "estado"
   // Tamaño = #reac x #enz x #subunid (sin inhibicion)
   // Tamaño = #reac x #enz x 2 x #subunid (con inhibicion)
   if((estado=(int***)calloc(NREAC,sizeof(int**)))==NULL)
   {
     printf("\nNo hay memoria para la matriz 'estado'");
     exit(1);
   }
   for(j=0;j<NREAC;j++)
   {
     if((estado[j]=(int**)calloc(nenz[j],sizeof(int*))) == NULL)
     {
       printf("\nNo hay memoria para la matriz 'estado'");
       exit(1);
     }
   }
   for(j=0;j<NREAC;j++)
   {
      if(inh[j][4]==-4)			// No hay inhibicion
      { 
         for(k=0;k<nenz[j];k++)
         {
            if((estado[j][k]=(int*)calloc(tipenz[j],sizeof(int))) == NULL)		
            {
               printf("\nNo hay memoria para la matriz 'estado'");
               exit(1);
            }
         } 
      } 
      if(inh[j][4]!=-4)			// Hay inhibicion
      { 
         for(k=0;k<nenz[j];k++)
         {
            if((estado[j][k]=(int*)calloc(2*tipenz[j],sizeof(int))) == NULL)		
            {
               printf("\nNo hay memoria para la matriz 'estado'");
               exit(1);
            }
         } 
      } 
   }

   // Reserva memoria para la matriz "fase"
   // Tamaño = #reac x #enz x #subunid (sin inhibicion)
   // Tamaño = #reac x #enz x 2 x #subunid (con inhibicion)
   if((fase=(double***)calloc(NREAC,sizeof(double**)))==NULL)
   {
     printf("\nNo hay memoria para la matriz 'fase'");
     exit(1);
   }
   for(j=0;j<NREAC;j++)
   {
     if((fase[j]=(double**)calloc(nenz[j],sizeof(double*))) == NULL)
     {
       printf("\nNo hay memoria para la matriz 'fase'");
       exit(1);
     }
   }
   for(j=0;j<NREAC;j++)
   {
      if(inh[j][4]==-4)			// No hay inhibicion
      { 
         for(k=0;k<nenz[j];k++)
         {
            if((fase[j][k]=(double*)calloc(tipenz[j],sizeof(double))) == NULL)	
            {
               printf("\nNo hay memoria para la matriz 'estado'");
               exit(1);
            }
         }
      }
      if(inh[j][4]!=-4)			// Hay inhibicion
      { 
         for(k=0;k<nenz[j];k++)
         {
            if((fase[j][k]=(double*)calloc(2*tipenz[j],sizeof(double))) == NULL)	// duplicado!
            {
               printf("\nNo hay memoria para la matriz 'estado'");
               exit(1);
            }
         }
      }
   }

   // Reserva memoria para la matriz "act"
   // Tamaño = #reac x #enz
   if((act=(int**)calloc(NREAC,sizeof(int*)))==NULL)
   {
     printf("\nNo hay memoria para la matriz 'act'");
     exit(1);
   }
   for(j=0;j<NREAC;j++)
   {
     if((act[j]=(int*)calloc(nenz[j],sizeof(int))) == NULL)
     {
       printf("\nNo hay memoria para la matriz 'act'");
       exit(1);
     }
   }

   //******************************* CALCULO DE PARAMETROS ADICIONALES **********************************

   // Modificacion de las velocidades de los oligomeros
   /*for(j=0;j<NREAC;j++)
   {
      if(tipenz[j]==2)			// Dimero
      { 
         vmax[j][0] /= 2.0; 
         vmax[j][1] /= 2.0; 
      }  

      if(tipenz[j]==3)			// Trimero 
      { 
         vmax[j][0] /= 3.0;
         vmax[j][1] /= 3.0;
      }

      if(tipenz[j]==4) 			// Tetramero
      {
         vmax[j][0] /= 4.0;		
         vmax[j][1] /= 4.0;		
      }
   }*/

   // Modificacion de las afinidades de los oligomeros
   /*for(j=0; j<2*NREAC; j++)
   {
      if(j%2 == 0) 		// Numero de la reaccion
         cdir = j/2; 

      // Modificacion de km 
      for(k=0; k<NSUB; k++)
      {
         if(tipenz[cdir] == 2)			// Dimero
         {		 
            km[j][k] *= sqrt(coop[cdir]); 
	    printf("\nDimero: ks[%d][%d] es %lf [sustratos]", j, k, km[j][k]);
	 }   

         if(tipenz[cdir] == 3)			// Trimero 
	 {	 
            km[j][k] *= 1.05*coop[cdir];
	    printf("\nTrimero: ks[%d][%d] es %lf [sustratos]", j, k, km[j][k]);
	 }   

         if(tipenz[cdir]==4) 
	 {	 
            km[j][k] *= pow(coop[cdir],1.53);	// Tetramero
	    printf("\nTetramero: ks[%d][%d] es %lf [sustratos]", j, k, km[j][k]);
         } 
      }
   }
   printf("\n");*/

   // Modificacion de las afinidades de los inhibidores
   for(j=0; j<NREAC; j++)
   {
      // Modificacion de ki, inhibicion competitiva
      if(inh[j][0] != -1)				// La reaccion directa es inhibida
      {
         if(tipenz[cdir] == 2)			// Dimero
            ki[inh[j][0]] *= sqrt(coop[cdir]); 

         if(tipenz[cdir] == 3)			// Trimero 
            ki[inh[j][0]] *= 1.05 * coop[cdir];

         if(tipenz[cdir] == 4) 			// Tetramero
            ki[inh[j][0]] *= pow(coop[cdir],1.53);
      }
      if(inh[j][1] != -1)				// La reaccion inversa es inhibida
      {
         if(tipenz[cdir] == 2)			// Dimero
            ki[inh[j][1]] *= sqrt(coop[cdir]); 

         if(tipenz[cdir] == 3)			// Trimero 
            ki[inh[j][1]] *= 1.05 * coop[cdir];

         if(tipenz[cdir] == 4) 			// Tetramero
            ki[inh[j][1]] *= pow(coop[cdir],1.53);
      }

      // Modificacion de ki, inhibicion no competitiva
      if(inh[j][2] != -1)				// Reaccion directa es inhibida
      {
         if(tipenz[cdir] == 2)			// Dimero
            ki[inh[j][2]] *= sqrt(coop[cdir]); 

         if(tipenz[cdir] == 3)			// Trimero 
            ki[inh[j][2]] *= 1.05 * coop[cdir];

         if(tipenz[cdir] == 4) 			// Tetramero
            ki[inh[j][2]] *= pow(coop[cdir],1.53);
      }
      if(inh[j][3] != -1)				// Reaccion inversa es inhibida
      {
         if(tipenz[cdir] == 2)			// Dimero
            ki[inh[j][3]] *= sqrt(coop[cdir]); 

         if(tipenz[cdir] == 3)			// Trimero 
            ki[inh[j][3]] *= 1.05 * coop[cdir];

         if(tipenz[cdir] == 4) 			// Tetramero
            ki[inh[j][3]] *= pow(coop[cdir],1.53);
       } 
   }  

   // Calcula los parametros "b" y "phic" para cada reaccion, ambas direcciones
   for(j=0;j<NREAC;j++)
   {
      bdir[j]=a/(vmax[j][0]*DT)-3.0*a;
      printf("\nEl numero de iteraciones para la r. dir. es %f, el No. de ciclos/iteracion es %lf",(bdir[j]/a)+3,(double)TOP/((bdir[j]/a)+3)); 
      phicdir[j]=1.0+bdir[j]*trec[j][0];
       
      if(vmax[j][1]<0.0)			// Reaccion inversa que no se da!
         binv[j]=0.0;
      else if(vmax[j][1]>0.0)   		// Si la reaccion inversa se da...	
         binv[j]=a/(vmax[j][1]*DT)-3.0*a;
      phicinv[j]=1.0+binv[j]*trec[j][1];
      if(vmax[j][1]>0.0)   		// Si la reaccion inversa se da...	
         printf("\nEl numero de iteraciones para la r. inv. es %f, el No. de ciclos/iteracion es %lf",(binv[j]/a)+3,(double)TOP/((binv[j]/a)+3)); 
   }
   printf("\n");

   //******************************************* INICIALIZACIONES *************************************	

   // Fija fases y estados iniciales de las enzimas
   for(i=0;i<NREAC;i++)
   {
      for(j=0;j<nenz[i];j++)
      {
         if(inh[i][4]==-4)			// No hay inhibicion 
         {
            for(k=0;k<tipenz[i];k++)		 
            {
               fase[i][j][k]=drand48();
               //fase[i][j][k]=0.3186;
               estado[i][j][k]=0; 
            }
         } 
         if(inh[i][4]!=-4)			// Hay inhibicion 
         {
            for(k=0;k<2*tipenz[i];k++)		
            {
               if(k%2==0)
                  fase[i][j][k]=drand48();
                  //fase[i][j][k]=0.3186;
               else
                  fase[i][j][k]=0.005;
               estado[i][j][k]=0;
            }
         } 
      }
   }
   
   // Calcula histograma 
   //for(j=0;j<RES;j++)
      //hh[j]=0.0;
   //hist(bdir[0]+1.);
   //hist(1.);
   
   // Inicializa registro de control de actualizacion de enzimas
   for(j=0;j<tenz;j++)
      cont[j]=j;

   // Inicializa matriz 'act'. Si se trata de un monomero 'act' siempre vale cero!
   for(j=0;j<NREAC;j++)
   {
      for(k=0;k<nenz[j];k++)
         act[j][k]=0;
   }

   A_acc=0.0;
   Aext=0.0;
   // Registra concentraciones de sustancias
   fprintf(dptr,"%lf \t %d \t %d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo/tau, n[0], n[1], n[0]+n[1], n[0]*1.e6/(NA*VSIM), n[1]*1.e6/(NA*VSIM), 1., 0., 0., 1., 0., 0., 0.); 
   tiempo+=DT;
  
   //***************************************** SIMULACION ************************************
   
   // Se agrega sustrato gradualmente
   sigma=13.0*pow(tau,2.0);
   Aext_prev = 0.0;
   A_acc = 0;
   flag = 0;
   flag2=0;
   n[0] = 0;
   ntop = TOP*tau/DT;

   //printf("\nEL VALOR DE NTOP ES %d\n",ntop);

   // Ciclo de iteracion de las enzimas
   for(zz=0;zz<ntop;zz++)
   {
      if(flag==0)
      {   
         Aext = 0.5*NINJ*NA*1.e-6*VSIM*(1. + erf(((double)zz*DT-2.0*tau)/sigma));
	 delta = (int)Aext - (int)Aext_prev;
     
         if(delta >= 1)
         {		 
	    n[0] += delta;
            A_acc += delta;
	    Aext_prev = (int)Aext;

	    if(A_acc >= (int)(NINJ*NA*1.e-6*VSIM))
               flag=1;		  
	 }   
      }	      
      
      if(flag==1)
      {
         printf("\nLa concentracion deseada de %lf Molar se alcanzo en el tiempo %lf tau", NINJ, zz*DT/tau);       
	 printf("\nSe agregaron %d moleculas de sustrato\n",A_acc);
         flag=2;
      }

      /*if((tiempo/tau > 4.899) && (tiempo/tau < 4.901))
      {
         //printf("\nCalculando histograma :)");     
         // Calcula histograma 
         for(j=0;j<RES;j++)
            hh[j]=0.0;
         hist(bdir[1]+1);
         //flag2 = 1;
      }*/ 	 

      // Baraja vector de recorrido de las enzimas
      baraja(zz); 
  
      // Barrido de las enzimas
      for(i=0;i<tenz;i++)
      {
         // Determina reac = reaccion que sera actualizada 
         index1=nenz[0];
         index2=0;
         for(j=1;j<NREAC+1;j++)
         {
            if(cont[i]<index1)
            {  
               reac=j-1; 
               break; 
            } 
            else 
               index1+=nenz[j];
            index2+=nenz[j-1];
         }
        
         // Determina enz = enzima que sera actualizada
         enz=cont[i]-index2;	        		

         // Calcula probabilidad de acuerdo al sustrato disponible
         probab(); 

         // Monomero
         if(tipenz[reac]==1)
         {
            if(inh[reac][4]==-4)					// Sin inhibicion
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)) 		// La reaccion se da
                  fase[reac][enz][0]=mapa_s(reac,enz,0,fase[reac][enz][0]);    	
            if(inh[reac][4]!=-4)					// Con inhibicion
            {
               fase[reac][enz][0]=mapa_s(reac,enz,0,fase[reac][enz][0]);    	// Mapa S
               fase[reac][enz][1]=mapa_i(reac,enz,1,fase[reac][enz][1]);    	// Mapa I
            } 
         } 
          
         // Dimero 
         if(tipenz[reac]==2)
         {
            // Actualiza fases para unidades en la enzima
            for(j=0;j<2;j++)
            {
               // Encuentra numero de unidades ocupadas en la enzima
               act[reac][enz]=0;
               for(k=0;k<2;k++)
               {
                  if(abs(estado[reac][enz][k])>0)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)) 		// La reaccion se da
                  fase[reac][enz][j]=mapa_s(reac,enz,j,fase[reac][enz][j]);    	
            }
         }  

         // Trimero
         if(tipenz[reac]==3)
         {
            // Actualiza fases para unidades en la enzima
            for(j=0;j<3;j++)
            {
               // Encuentra numero de unidades ocupadas en la enzima
               act[reac][enz]=0;
               for(k=0;k<3;k++)
               {
                  if(abs(estado[reac][enz][k])>0)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)) 		// La reaccion se da
                  fase[reac][enz][j]=mapa_s(reac,enz,j,fase[reac][enz][j]);    	
            } 
         }  

         // Tetramero 
         if(tipenz[reac]==4)
         {
            // Actualiza fases para unidades en la enzima
            for(j=0;j<4;j++)
            {
               // Encuentra numero de unidades ocupadas en la enzima
               act[reac][enz]=0;
               for(k=0;k<4;k++)
               {
                  if(abs(estado[reac][enz][k])>0)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)||(estado[reac][enz][j]>=2)) 		// La reaccion se da
                  fase[reac][enz][j]=mapa_s(reac,enz,j,fase[reac][enz][j]);    	
            }
         }  
      }

      // Registra proporciones de enzimas en reposo, procesando y en recuperacion
      fa1 = 0.;
      fa2 = 0.;
      fa3 = 0.;
      fb1 = 0.;
      fb2 = 0.;
      fb3 = 0.;
      for(j=0;j<NENZ;j++)
      {
         if(estado[0][j][0] == 0 || abs(estado[0][j][0]) == 1)	      
            fa1+=1.0;
	 else if(abs(estado[0][j][0]) == 2)
            fa2+=1.0;
	 else if(abs(estado[0][j][0]) == 3)
	    fa3+=1.0;	 
	 
         if(estado[1][j][0] == 0 || abs(estado[1][j][0]) == 1)	      
            fb1+=1.0;
	 else if(abs(estado[1][j][0]) == 2)
            fb2+=1.0;
	 else if(abs(estado[1][j][0]) == 3)
	    fb3+=1.0;	 

         if(estado[0][j][1] == 0 || abs(estado[0][j][1]) == 1)	      
            fa1+=1.0;
	 else if(abs(estado[0][j][1]) == 2)
            fa2+=1.0;
	 else if(abs(estado[0][j][1]) == 3)
	    fa3+=1.0;	 
	 
         if(estado[1][j][1] == 0 || abs(estado[1][j][1]) == 1)	      
            fb1+=1.0;
	 else if(abs(estado[1][j][1]) == 2)
            fb2+=1.0;
	 else if(abs(estado[1][j][1]) == 3)
	    fb3+=1.0;	 

         if(estado[0][j][2] == 0 || abs(estado[0][j][2]) == 1)	      
            fa1+=1.0;
	 else if(abs(estado[0][j][2]) == 2)
            fa2+=1.0;
	 else if(abs(estado[0][j][2]) == 3)
	    fa3+=1.0;	 
	 
         if(estado[1][j][2] == 0 || abs(estado[1][j][2]) == 1)	      
            fb1+=1.0;
	 else if(abs(estado[1][j][2]) == 2)
            fb2+=1.0;
	 else if(abs(estado[1][j][2]) == 3)
	    fb3+=1.0;	 

         if(estado[0][j][3] == 0 || abs(estado[0][j][3]) == 1)	      
            fa1+=1.0;
	 else if(abs(estado[0][j][3]) == 2)
            fa2+=1.0;
	 else if(abs(estado[0][j][3]) == 3)
	    fa3+=1.0;	 
	 
         if(estado[1][j][3] == 0 || abs(estado[1][j][3]) == 1)	      
            fb1+=1.0;
	 else if(abs(estado[1][j][3]) == 2)
            fb2+=1.0;
	 else if(abs(estado[1][j][3]) == 3)
	    fb3+=1.0;	 
      }
      
      fa1 /= (double)NENZ*4.;
      fa2 /= (double)NENZ*4.;
      fa3 /= (double)NENZ*4.;
      fb1 /= (double)NENZ*4.;
      fb2 /= (double)NENZ*4.;
      fb3 /= (double)NENZ*4.;

      // Registra fraccion de subunidades ocupadas
      f0=0.0;
      f1=0.0;
      f2=0.0;
      f3=0.0;
      f4=0.0;
      for(j=0;j<NENZ;j++)
      {
         if(act[0][j]==0)
            f0+=1.0;
         if(act[0][j]==1)
            f1+=1.0;
         if(act[0][j]==2)
	    f2+=1.0;
         if(act[0][j]==3)
	    f3+=1.0;
         if(act[0][j]==4)
	    f4+=1.0;	 
      }

      /*f0/=(double)NENZ;      
      f1/=(double)NENZ;      
      f2/=(double)NENZ;      
      f3/=(double)NENZ;      
      f4/=(double)NENZ;*/      

      // Calcula histograma
      /*if(tiempo/tau==4.9)
      {	      
         hist(bdir[0]+1.);
	 printf("\nHistograma calculado!");
      }*/ 	 

      // Registra concentraciones de sustancias
      fprintf(dptr,"%lf \t %d \t %d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo/tau, n[0], n[1], n[0]+n[1], n[0]*1.e6/(NA*VSIM), n[1]*1.e6/(NA*VSIM), 1., 0., 0., 1., 0., 0., 0.); 
      /*if(zz%30==0){
         fprintf(dptr,"%lf \t %d \t %d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
			 tiempo/tau,n[0],n[1],n[0]+n[1],fa1,fa2,fa3,fb1,fb2,fb3,f0,f1,f2,f3,f4); 
	 //printf("\n fb2 = %lf",fb2);
      }*/
      //fprintf(dptr,"%lf \t %d \t %d \t %d\n", tiempo/tau, n[0], n[1], n[0]+n[1]); 
      tiempo+=DT;

      //if(zz%(10*(int)((bdir[0]/a)+3))==0)
      //{
         // Calcula histograma 
         //for(j=0;j<RES;j++)
         //   hh[j]=0.0;
         //hist(bdir[0]+1);

         // Presenta histograma
         /*cpgscr(1,1.0,0.7,0.5);  
         cpgenv(-0.01,1.01,0.0,30.0,0,0);
         cpglab("Fase","Frecuencia","Evolucion Temporal Histograma de Fase");
         cpgbbuf();
         cpgbin(RES,fas,hh,1);	
         cpgebuf();*/
      //}
   }
   fclose(dptr);

   // Conserva histograma final
   for(j=0;j<RES;j++)
      fprintf(fhist,"%f \t %f \n",fas[j],hh[j]);
   fclose(fhist);

   // Termina pgplot
   //cpgend();
   //return(0);
}

//******************************************************************************************************
//***************************************** SUBRUTINAS *************************************************
//******************************************************************************************************

// Barajado del orden en que las enzimas son recorridas
void baraja(z)
int z;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern int *cont;
   int i,j,k,pos;

   if(z==0)					// Si es la primera vez, baraja 20 veces
   {
      for(k=0;k<20;k++)
      {   
         for(j=0;j<tenz;j++)
         {
            pos=drnd(tenz); 
            i=cont[pos];
            cont[pos]=cont[j];
            cont[j]=i;
         }
      }        
   }
   else if(z%20==0)		// Baraja cada 20 iteraciones
   {
      for(j=0;j<tenz;j++)
      {
         pos=drnd(tenz); 
         i=cont[pos];
         cont[pos]=cont[j];
         cont[j]=i;
      }
   }
}

//***********************************************************************************************************

// Calculo de la probabilidad de arranque de las enzimas
void probab()
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double p1[2][NREAC],km[2*NREAC][NSUB],vmax[NREAC][2];
   // extern int eps[NREAC][NSUB],n[NSUB];
   int j,k,cdir,cinv, yo;
   double pimp;

   for(j=0;j<NREAC;j++)
   {
      p1[0][j]=1.0;		// Reaccion directa
      p1[1][j]=1.0;		// Reaccion inversa	
   }

   for(j=0;j<NREAC;j++)
   {
      for(k=0;k<NSUB;k++)
      { 
         if(eps[j][k]<0) 			// Reaccion directa
         {
            p1[0][j]*=DT*pow((double)n[k]*(vmax[j][0]/km[2*j][k]),(double)abs(eps[j][k]));
	    yo = n[k];
	    if(n[k] == 1)
	       yo = 0;	    
	    pimp = p1[0][j];
	    pimp = km[2*j][k];
            if(p1[0][j]>=1)
               p1[0][j]=0.99999999;   
         }

         if(eps[j][k]>0) 			// Reaccion inversa
         {
            if(vmax[j][1]<0.0)			// La reaccion no se da
               p1[1][j]*=0.0;   
            else				// La reaccion se da
            { 
               if(km[2*j+1][k]==0.0)		// Substancia no participa en la reaccion
                  p1[1][j]*=1.0;   
               if(km[2*j+1][k]!=0.0)
                  p1[1][j]*=DT*pow((double)n[k]*(vmax[j][1]/km[2*j+1][k]),(double)abs(eps[j][k]));   
            } 
            if(p1[1][j]>=1)
               p1[1][j]=0.99999999;   
         }   
      }
   }
}

//***********************************************************************************************************

// Avance de fase para mapa S
double mapa_s(pos1,pos2,pos3,state_in)
double state_in;
int pos1,pos2,pos3;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double bdir[NREAC],binv[NREAC],p1[2][NREAC],km[2*NREAC][NSUB];  
   // extern double a,p2;
   // extern int eps[NREAC][NSUB],n[NSUB];
   // extern int ***estado;
   // extern int **act;
   double state_out=0.0,dado,umbral;
   int flag1,flag2,m,j,gof=0,gob=0,state;

   // pos1 = reaccion a ser actualizada
   // pos2 = enzima dentro del grupo que lleva a cabo la reaccion
   // pos3 = subunidad dentro del oligomero

   umbral=p1[0][pos1]/(p1[0][pos1]+p1[1][pos1]);
   dado=drand48();

   // La enzima esta en reposo
   if(estado[pos1][pos2][pos3]==0)
   {
      gof=verifdir(pos1);			// sustratos reaccion directa disponibles?
      gob=verifinv(pos1);			// sustratos reaccion inversa disponibles?

      // Dos posibilidades:
      // a) La enzima es un oligomero y ninguna unidad de la enzima ha sido ocupada 
      // b) La enzima es un monomero 
      if(act[pos1][pos2]==0)
      {
         // Hay sustrato para la reaccion directa y para la reaccion inversa
         // La direccion de la reaccion se determina de forma probabilistica
         if((gof==1)&(gob==1))
         {
            if(dado<umbral)		
            {
               state_out=avandir(state_in,pos1,pos2,pos3,0);
               estado[pos1][pos2][pos3]=1;
            } 

            if(dado>umbral)
            {
               state_out=avaninv(state_in,pos1,pos2,pos3,0);
               estado[pos1][pos2][pos3]=-1;
            }
         }

         // Solo hay sustrato para reaccion directa
         if((gof==1)&(gob==0))
         {
            state_out=avandir(state_in,pos1,pos2,pos3,0);
            estado[pos1][pos2][pos3]=1;
         } 

         // Solo hay sustrato para reaccion inversa
         if((gof==0)&(gob==1))
         {
            state_out=avaninv(state_in,pos1,pos2,pos3,0);
            estado[pos1][pos2][pos3]=-1;
         }

         // No hay sustratos disponibles
         if((gof==0)&(gob==0))
            state_out=state_in;    

         // Revisa estados para tomar sustrato o liberar producto 
         sustprod(pos1,pos2,pos3,state_in,state_out);				

         return state_out;
      }

      // La enzima es un oligomero y tiene unidades que ya han sido ocupadas
      // Es necesario forzar el resto de las unidades en una misma direccion! 
      if(act[pos1][pos2]!=0)
      {
         // Hay sustrato para la reaccion directa y para la reaccion inversa
         if((gof==1)&(gob==1))
         {
            // Itera con mapa reaccion directa si las unidades ocupadas van en esa direccion	
            if(act[pos1][pos2]>0)		// La enzima va hacia adelante		
            {
               state_out=avandir(state_in,pos1,pos2,pos3,0);
               estado[pos1][pos2][pos3]=1;
            } 

            // Itera con mapa reaccion inversa si las unidades ocupadas van en esa direccion
            if(act[pos1][pos2]<0)		// La enzima va hacia atras
            {
               state_out=avaninv(state_in,pos1,pos2,pos3,0);
               estado[pos1][pos2][pos3]=-1;
            }
         }

         // Solo hay sustrato para reaccion directa: la reaccion procede si las unidades
         // ocupadas van en esa misma direccion
         if((gof==1)&(gob==0)&(act[pos1][pos2]>0))
         {
            state_out=avandir(state_in,pos1,pos2,pos3,0);
            estado[pos1][pos2][pos3]=1;
         } 

         // Solo hay sustrato para reaccion inversa: la reaccion procede si las unidades ocupadas
         // van en esa misma direccion
         if((gof==0)&(gob==1)&(act[pos1][pos2]<0))
         {
            state_out=avaninv(state_in,pos1,pos2,pos3,0);
            estado[pos1][pos2][pos3]=-1;
         }

         // No hay sustratos disponibles
         if((gof==0)&(gob==0))
            state_out=state_in;    

         // Revisa estados para tomar sustrato o liberar producto 
         sustprod(pos1,pos2,pos3,state_in,state_out);				

         return state_out;
      } 
   }

   // La enzima ya arranco
   // Avance de fase de la enzima, reaccion directa
   state=estado[pos1][pos2][pos3];
   if((state<10)&(abs(state)>1))
   {
      state=estado[pos1][pos2][pos3];
      if(state>0)		
         state_out=avandir(state_in,pos1,pos2,pos3,1);

      // Avance de fase de la enzima, reaccion inversa
      if(state<0)		
         state_out=avaninv(state_in,pos1,pos2,pos3,1);
   }

   // Revisa estados para tomar sustrato o liberar producto 
   sustprod(pos1,pos2,pos3,state_in,state_out);				

   // La enzima completa su ciclo, queda disponible
   if(estado[pos1][pos2][pos3]==3)	
   {
      if(fabs(state_in-state_out)>=1.0+bdir[pos1])			// Reaccion directa
         estado[pos1][pos2][pos3]=0;					// Enzima queda disponible
   }
   if(estado[pos1][pos2][pos3]==-3)	
   {
      if(fabs(state_in-state_out)>=1.0+binv[pos1])			// Reaccion inversa
         estado[pos1][pos2][pos3]=0;					// Enzima queda disponible
   }

   // La enzima esta inhibida
   if((estado[pos1][pos2][pos3]==10)||(estado[pos1][pos2][pos3]==11)||(estado[pos1][pos2][pos3]==12))
      state_out=state_in;

   return state_out;
}

//***********************************************************************************************************

// Avance de fase para mapa I
double mapa_i(pos1,pos2,pos3,state_in)
double state_in;
int pos1,pos2,pos3;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double bdir[NREAC],binv[NREAC],p1[2][NREAC],km[2*NREAC][NSUB];  
   // extern double a,p2;
   // extern int eps[NREAC][NSUB],n[NSUB];
   // extern int ***estado;
   // extern int **act;
   double state_out=0.0,pin,pout;
   int flag1,flag2,m,j,gof=0;

   // pos1 = reaccion a ser actualizada
   // pos2 = enzima dentro del grupo que lleva a cabo la reaccion
   // pos3 = subunidad dentro del oligomero

   // Determina si hay inhibidor, calcula probabilidades 
   gof=0;
   if(inh[pos1][0]!=-1) 			// R. directa, I. competitiva
   {
      if(n[inh[pos1][0]]>0)
      {
         gof=1;
         pin=7.0*n[inh[pos1][0]]*vmax[pos1][0]/ki[inh[pos1][0]];	
         pout=.02; 
      } 
   }
   else if(inh[pos1][1]!=-1) 			// R. inversa, I. competitiva
   { 
      if(n[inh[pos1][1]]>0)
      {
         gof=1; 
         pin=7.0*n[inh[pos1][1]]*vmax[pos1][1]/ki[inh[pos1][1]];	
         pout=.02; 
      }
   }
   else if(inh[pos1][2]!=-1) 			// R. directa, I. no competitiva
   {
      if(n[inh[pos1][2]]>0)
      { 
         gof=1; 
         pin=0.75*n[inh[pos1][2]]*vmax[pos1][0]/ki[inh[pos1][2]];	
         pout=.02; 
      }
   }
   else if(inh[pos1][3]!=-1) 			// R. inversa, I. no competitiva
   {
      if(n[inh[pos1][3]]>0)
      {  
         gof=1; 
         pin=0.75*n[inh[pos1][3]]*vmax[pos1][1]/ki[inh[pos1][3]];	
         pout=.02; 
      }
   }
       
   // INHIBICION COMPETITIVA
   if((inh[pos1][0]!=-1)||(inh[pos1][1]!=-1))
   { 
      // Hay inhibidor, el sustrato no ha entrado: el inhibidor comienza entrada
      if((gof==1)&(abs(estado[pos1][pos2][pos3-1])<2)&(estado[pos1][pos2][pos3]==0))
      {
         state_out=avaninh(state_in,pin);
         estado[pos1][pos2][pos3]=5;   
      }
      // Hay inhibidor, el sustrato no ha entrado: el inhibidor prosigue entrada 
      else if((gof==1)&(abs(estado[pos1][pos2][pos3-1])<2)&(estado[pos1][pos2][pos3]==5))  
      {
         state_out=avaninh(state_in,pin);

         // La enzima toma inhibidor
         if(state_out>1.0)
         {
            // Se forma complejo EI, el mapa S queda inhabilitado
            estado[pos1][pos2][pos3-1]=10;
                 
            // Actualiza concentracion de inhibidor
            if(inh[pos1][0]!=-1) 			// R. directa, I. competitiva
               n[inh[pos1][0]]--;
            else if(inh[pos1][1]!=-1) 		// R. inversa, I. competitiva
               n[inh[pos1][1]]--;
            else if(inh[pos1][2]!=-1) 		// R. directa, I. no competitiva
               n[inh[pos1][2]]--;
            else if(inh[pos1][3]!=-1) 		// R. inversa, I. no competitiva
               n[inh[pos1][3]]--;
                
            // Inhibidor comienza salida
            estado[pos1][pos2][pos3]=-5;
            state_out=0.005;   
         }  
      } 
      // SALIDA DEL INHIBIDOR
      else if(estado[pos1][pos2][pos3]==-5)	// El inhibidor prosigue salida
      {
         state_out=avaninh(state_in,pout);

         // La enzima libera inhibidor
         if(state_out>1.0)
         {
            // Mapa S queda habilitado
            estado[pos1][pos2][pos3-1]=0;
                 
            // Actualiza concentracion de inhibidor
            if(inh[pos1][0]!=-1) 			// R. directa, I. competitiva
               n[inh[pos1][0]]++;
            else if(inh[pos1][1]!=-1) 		// R. inversa, I. competitiva
               n[inh[pos1][1]]++;
            else if(inh[pos1][2]!=-1) 		// R. directa, I. no competitiva
               n[inh[pos1][2]]++;
            else if(inh[pos1][3]!=-1) 		// R. inversa, I. no competitiva
               n[inh[pos1][3]]++;
                
            // Mapa I queda en reposo
            estado[pos1][pos2][pos3]=0;
            state_out=drand48();   
         }  
      }   
      else 
         gof=0;					// La fase no avanza  
   }

   // INHIBICION NO COMPETITIVA
   if((inh[pos1][2]!=-1)||(inh[pos1][3]!=-1))
   {
      // Hay inhibidor: el inhibidor comienza entrada
      if((gof==1)&(estado[pos1][pos2][pos3]==0))
      {
         state_out=avaninh(state_in,pin);
         estado[pos1][pos2][pos3]=5;   
      }
      // Hay inhibidor: el inhibidor prosigue entrada 
      else if((gof==1)&(estado[pos1][pos2][pos3]==5))  
      {
         state_out=avaninh(state_in,pin);

         // La enzima toma inhibidor
         if(state_out>1.0)
         {
            // Se forma complejo EI, avance de la enzima queda detenido 
            if(abs(estado[pos1][pos2][pos3-1])<2)
               estado[pos1][pos2][pos3-1]=10;
            // Se forma complejo ESI, avance de la enzima queda detenido
            else if(estado[pos1][pos2][pos3-1]==2)		// Enzima va hacia adelante
               estado[pos1][pos2][pos3-1]=11;
            else if(estado[pos1][pos2][pos3-1]==-2)		// Enzima va hacia atras
               estado[pos1][pos2][pos3-1]=12;

            // Actualiza concentracion de inhibidor
            if(inh[pos1][0]!=-1) 		// R. directa, I. competitiva
               n[inh[pos1][0]]--;
            else if(inh[pos1][1]!=-1) 		// R. inversa, I. competitiva
               n[inh[pos1][1]]--;
            else if(inh[pos1][2]!=-1) 		// R. directa, I. no competitiva
               n[inh[pos1][2]]--;
            else if(inh[pos1][3]!=-1) 		// R. inversa, I. no competitiva
               n[inh[pos1][3]]--;
                
            // Inhibidor comienza salida
            estado[pos1][pos2][pos3]=-5;
            state_out=0.005;   
         }  
      } 
      // SALIDA DEL INHIBIDOR
      else if(estado[pos1][pos2][pos3]==-5)	// El inhibidor prosigue salida
      {
         state_out=avaninh(state_in,pout);

         // La enzima libera inhibidor
         if(state_out>1.0)
         {
            // Se libera inhibidor del complejo EI 
            if(estado[pos1][pos2][pos3-1]==10)
               estado[pos1][pos2][pos3-1]=0;
            // Se forma complejo ESI, avance de la enzima queda detenido
            else if(abs(estado[pos1][pos2][pos3-1])==11)	// Hacia adelante
               estado[pos1][pos2][pos3-1]=2;
            else if(estado[pos1][pos2][pos3-1]==12)		// Hacia atras
               estado[pos1][pos2][pos3-1]=-2;

                 
            // Actualiza concentracion de inhibidor
            if(inh[pos1][0]!=-1) 			// R. directa, I. competitiva
               n[inh[pos1][0]]++;
            else if(inh[pos1][1]!=-1) 		// R. inversa, I. competitiva
               n[inh[pos1][1]]++;
            else if(inh[pos1][2]!=-1) 		// R. directa, I. no competitiva
               n[inh[pos1][2]]++;
            else if(inh[pos1][3]!=-1) 		// R. inversa, I. no competitiva
               n[inh[pos1][3]]++;
                
            // Mapa I queda en reposo
            estado[pos1][pos2][pos3]=0;
            state_out=drand48();   
            gof=1; 
         }  
      }
      else 
         gof=0;					// La fase no avanza  
   }

   // No hubo avance de fase para el inhibidor
   if(gof==0)
      state_out=state_in;
  
   return state_out;
}

//***********************************************************************************************************

// Actualizacion de concentraciones 
void sustprod(pos1,pos2,pos3,state_in,state_out)
double state_in,state_out;
int pos1,pos2,pos3;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double phicdir[NREAC],phicinv[NREAC],km[2*NREAC][NSUB];
   // extern int eps[NREAC][NSUB],n[NSUB];
   // extern int ***estado; 
   int gof,gob,m,j;

   gof=verifdir(pos1);			// sustratos reaccion directa disponibles?
   gob=verifinv(pos1);			// sustratos reaccion inversa disponibles?

   // La enzima toma sustratos para la reaccion directa
   if(estado[pos1][pos2][pos3]==1)			
   {
      if((state_in<=1.0)&(state_out>1.0)&(gof==1))	
      {
         for(m=0;m<NSUB;m++)
         { 
            if(eps[pos1][m]<0)			// La sustancia es sustrato de la reaccion directa
               n[m]+=eps[pos1][m];		// Actualiza concentraciones
         }
         estado[pos1][pos2][pos3]=2;		// La enzima queda ocupada, region laminar directa = 2
         probab();				// Actualiza probabilidades
 
         if((inh[pos1][0]!=-1)&(inh[pos1][1]!=-1))	// Habia un inhibidor por entrar, I. competitiva
            estado[pos1][pos2][pos3+1]=0;
      }

      // La enzima sigue en la region caotica
      if((state_in<=1.0)&(state_out<=1.0))	
      {
         estado[pos1][pos2][pos3]=0;  
         return;
      }
   }        

   // La enzima toma sustratos para la reaccion inversa
   if(estado[pos1][pos2][pos3]==-1)			
   {
      if((state_in<=1.0)&(state_out>1.0)&(gob==1))	
      {  
         for(m=0;m<NSUB;m++)
         {  
            if(eps[pos1][m]>0) 			// La sustancia es sustrato de la reaccion inversa
            { 
               if(km[2*pos1+1][m]!=0.0)		// La sustancia participa en la reaccion inversa
                  n[m]-=eps[pos1][m];		// Actualiza concentraciones
            }
         } 
         estado[pos1][pos2][pos3]=-2;		// La enzima queda ocupada, region laminar inversa = -2
         probab();				// Actualiza probabilidades

         if((inh[pos1][0]!=-1)&(inh[pos1][1]!=-1))	// Habia un inhibidor por entrar, I. competitiva
            estado[pos1][pos2][pos3+1]=0;
      }

      // La enzima sigue en la region caotica
      if((state_in<=1.0)&(state_out<=1.0))	
      {
         estado[pos1][pos2][pos3]=0;  
         return;
      } 
   }
 
   // Liberacion de producto en la reaccion directa
   if(estado[pos1][pos2][pos3]==2)
   {
      // Reaccion directa libera producto, la enzima queda en recuperacion 
      if(state_out>=phicdir[pos1])	
      {
         for(j=0;j<NSUB;j++)
         {
            if(eps[pos1][j]>0)				// Libera producto
               n[j]+=eps[pos1][j];			// Actualiza concentracion
         }
      	 estado[pos1][pos2][pos3]=3;			// Enzima en recuperacion, reaccion directa
         probab();					// Actualiza probabilidades
      }
   }

   // Liberacion de producto en la reaccion inversa
   if(estado[pos1][pos2][pos3]==-2)
   {
      // Reaccion inversa libera producto, la enzima queda en recuperacion 
      if(state_out>=phicinv[pos1])	
      {
         for(j=0;j<NSUB;j++)
         {
            if(eps[pos1][j]<0)		 		// Libera producto
               n[j]-=eps[pos1][j];			// Actualiza concentracion
         }
         estado[pos1][pos2][pos3]=-3;			// Enzima en recuperacion, reaccion inversa
         probab();					// Actualiza probabilidades
      }
   }  
   return;
}



//***********************************************************************************************************

// Avance de fase reaccion directa
double avandir(state_in,pos1,pos2,pos3,arr)
double state_in;
int pos1,pos2,pos3,arr;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double bdir[NREAC],binv[NREAC],p1[2][NREAC];  
   // extern double a,p2;
   // extern int ***estado;
   double state_out,p;
   int flag1,flag2;

   p=p1[0][pos1];

   // La enzima es un dimero
   if(tipenz[pos1]==2)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
   }

   // La enzima es un trimero	
   if(tipenz[pos1]==3)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
      if(abs(act[pos1][pos2])==2)
         p*=pow(coop[pos1],2.00);  
   }

   // La enzima es un tetramero
   if(tipenz[pos1]==4)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
      if(abs(act[pos1][pos2])==2)
         p*=pow(coop[pos1],2.00);  
      if(abs(act[pos1][pos2])==3)
         p*=pow(coop[pos1],3.00);  
   }

   flag1=0;
   flag2=0;
   while(flag1==0)			// Toma en cuenta puntos fijos
   {
      while(flag2==0)
      {
         if(state_in==0.5)		// Fase = 0.5
            state_in=drand48();  
         else if(state_in!=0.5)
            flag2=1;
      } 
      state_out=enzima47(state_in,a,bdir[pos1],p,p2);               
      if(state_out!=state_in)
         flag1=1;
      else				// Punto fijo
      { 
         state_in=drand48();
         flag2=0;  
      }
   }
   if(arr==0)
      estado[pos1][pos2][pos3]=1;		// Iteracion mapa reaccion directa 
    
   return state_out;
}

//***********************************************************************************************************

// Avance de fase reaccion inversa
double avaninv(state_in,pos1,pos2,pos3,arr)
double state_in;
int pos1,pos2,pos3,arr;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double bdir[NREAC],binv[NREAC],p1[2][NREAC];  
   // extern double a,p2;
   // extern int ***estado;
   double state_out,p;
   int flag1,flag2;

   p=p1[1][pos1];

   // La enzima es un dimero
   if(tipenz[pos1]==2)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
   }

   // La enzima es un trimero	
   if(tipenz[pos1]==3)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
      if(abs(act[pos1][pos2])==2)
         p*=pow(coop[pos1],2.00);  
   }

   // La enzima es un tetramero
   if(tipenz[pos1]==4)
   {
      if(abs(act[pos1][pos2])==1)
         p*=coop[pos1];  
      if(abs(act[pos1][pos2])==2)
         p*=pow(coop[pos1],2.00);  
      if(abs(act[pos1][pos2])==3)
         p*=pow(coop[pos1],3.00);  
   }

   flag1=0;
   flag2=0;
   while(flag1==0)			// Toma en cuenta puntos fijos
   {
      while(flag2==0)
      {
         if(state_in==0.5)		// Fase = 0.5
            state_in=drand48();  
         else if(state_in!=0.5)
            flag2=1;
      } 
      state_out=enzima47(state_in,a,binv[pos1],p1[1][pos1],p2);               
      if(state_out!=state_in)
         flag1=1;
      else				// Punto fijo
      { 
         state_in=drand48();
         flag2=0;  
      }
   }
   
   if(arr==0)
      estado[pos1][pos2][pos3]=-1;		// Iteracion mapa reaccion inversa
      
   return state_out;
}

//***********************************************************************************************************

// Avance de fase para inhibidores
double avaninh(state_in,pr)
double state_in,pr;
{
   // Esta rutina emplea las siguientes variables externas:
   // extern double a,p2;
   double state_out;
   int flag1,flag2;

   flag1=0;
   flag2=0;
   while(flag1==0)			// Toma en cuenta puntos fijos
   {
      while(flag2==0)
      {
         if(state_in==0.5)		// Fase = 0.5
            state_in=drand48();  
         else if(state_in!=0.5)
            flag2=1;
      } 
      state_out=enzima47(state_in,a,5.00,pr,p2);               
      if(state_out!=state_in)
         flag1=1;
      else				// Punto fijo
      { 
         state_in=drand48();
         flag2=0;  
      }
   }
   return state_out;
}

//***********************************************************************************************************

// Enzima47
double enzima47(x,a,b,p1,p2)
double x,a,b,p1,p2;
{
   double out,m1,m2,m3,m4; 

   if(p1>=1)
      p1=0.9999;

   m1 = 2/(1-p1);
   m2 = 2*fabs(a)/p1;
   m3 = 2*(1-fabs(a))/(1-p2);
   m4 = 2*fabs(a)/p2;

   if(a>=0.0)
   { 
      if(x<0.0) 
         out = a+b+2+2*x;
      else if((x>=0)&(x<0.5-0.5*p1)) 
         out = m1*x;
      else if((x>=0.5-0.5*p1)&(x<0.5)) 
         out = 1+a-0.5*m2+m2*x;
      else if((x>=0.5)&(x<0.5+0.5*p1)) 
         out = 0.5*m2*(1+p1)+1-m2*x;
      else if((x>=0.5+0.5*p1)&(x<=1))  
         out = m1*(1-x);   
      else if((x>1)&(x<(1+b))) 
         out = a+x;
      else if((x>=(b+1))&(x<(b+1.5-p2/2))) 
         out = m3*(x-1.5-b+p2/2+(b+2)/m3);
      else if((x>=(b+1.5-p2/2))&(x<(b+1.5))) 
         out = m4*(x-1.5-b+(b+2+a)/m4);
      else if((x>=(b+1.5))&(x<(1.5+b+p2/2))) 
         out = b+2+a+m4*(1.5+b)-m4*x;
      else if((x>=(b+1.5+p2/2))&(x<(2+b))) 
         out = b+2+m3*(1.5+b+p2/2)-m3*x;
      else if(x>=(2+b)) 
         out = -1/a*(b+2)+1/a*x;
   }
   return out;
}

//***********************************************************************************************************

// Genera un numero aleatorio entero entre 0 y PHIMAX 
int drnd(phimax)
int phimax;
{
	int j;
        float dado;
        dado=drand48()*phimax;  
        for(j=0;j<phimax;j++)
        {
           if((j<=dado)&(dado<j+1))
              return j;
        } 
}

//***********************************************************************************************************

// Inicializa generador de numeros aleatorios 
void rand_init()
{
	int seed;
	time_t timenow;
	time(&timenow);
	seed=timenow % 100000;
        //printf("\nLa semilla es %d",seed);
	srand48(seed);
}

//***********************************************************************************************************

// Verifica presencia de todos los sustratos para reaccion directa
int verifdir(pos1)
int pos1;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern int eps[NREAC][NSUB],n[NSUB];
   int gof,m;

   gof=1;			       		// Verifica los sustratos esten disponibles (gof=1)
   for(m=0;m<NSUB;m++)
   {
      if(eps[pos1][m]<0)			// Si substancia es sustrato de la reaccion directa...
      {
         if(n[m]>=abs(eps[pos1][m]))  		// Si TODOS los sustratos estan disponibles...
            gof*=1;				// ...la reaccion procede
         else
            gof*=0;                        
      }  
   }
   return gof;
}

//***********************************************************************************************************

// Verifica presencia de todos los sustratos para reaccion inversa
int verifinv(pos1)
int pos1;
{
   // Esta rutina emplea las siguientes varaibles externas:
   // extern double km[2*NREAC][NSUB];
   // extern int eps[NREAC][NSUB],n[NSUB];
   int gob,m;

   gob=1;					// Verifica que los sustratos esten disponibles (gob=1)
   for(m=0;m<NSUB;m++)
   {
      if(eps[pos1][m]>0)			// Si sustancia es sustrato de la reaccion inversa...
      {
         if(km[2*pos1+1][m]!=0.0)		// La sustancia participa en la reaccion inversa
         {  
            if(n[m]>=eps[pos1][m]) 		// Si TODOS los sustratos estan disponibles...
               gob*=1;				// ...la reaccion procede
            else
               gob*=0;                        
         }
	 else
            gob*=0;		 
      }  
   }
   return gob;
}
//***********************************************************************************************************

/* Histograma */
void hist(tope)
double tope;
{
   int j,k;
   double delta,inf,sup,pimp,cuenta;

   delta=tope/(double)RES;
   cuenta=0.0;

   for(j=0;j<RES;j++)
   {
      inf=(double)(j*delta);
      sup=(double)((j+1)*delta);
      for(k=0;k<NENZ;k++)
      {
         pimp=fase[0][k][0];
         //if((inf<pimp)&(pimp<sup)&(estado[0][k][0]>=0))
         //if((inf<pimp)&(pimp<sup))
         if((inf<fase[0][k][0])&(fase[0][k][0]<sup))
            hh[j]+=1.0;
      }
      fas[j]=cuenta;
      cuenta+=delta;
   }

   // Normalizacion de las fases
   for(j=0;j<RES;j++)
      fas[j]/=cuenta;
}

//**********************************************************************************************************************

