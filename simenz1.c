/* 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	

   simenz2.c
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
   *
   * 02/06/2008: Problema detectado: la reaccion 5 no camina cuando el volumen es relativamente grande o las concentraciones son pequeñas!
   * 
   * Verificacion de secuencia de reacciones en glicolisis: OK!
   *
   * Ajusta el numero de puntos a escribir para que sea siempre 1e6
     
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<time.h>
//#include "cpgplot.h"

#define NSUB 17		// Numero de substancias
#define NREAC 10	// Numero de reacciones
#define NENZ 50		// Numero de enzimas
//#define NENZ 50		// Numero de enzimas
//#define DT 100e-6	// Incremento de tiempo correspondiente a una iteracion (segundos)
#define RES 10		// Resolucion del histograma
#define NA 6.022e23	// Numero de Avogadro
//#define VSIM 10.3e-18 	// Volumen de simulacion (litros)
#define VSIM 1.03e-19 	// Volumen de simulacion (litros)
#define NINJ 200.0 	// Concentracion final de sustrato (microMolar)
#define TTOP 100.0		// tiempo de simulacion, en segundos
//#define INTERV 10	// Intervalo para escritura de datos

// Variables globales
FILE *dptr;
double p1[2][NREAC]={0.0},vmax[NREAC][2]={0.0},km[2*NREAC][NSUB]={0.0},bdir[NREAC]={0.0},binv[NREAC]={0.0};
double phicdir[NREAC]={0.0},phicinv[NREAC]={0.0},a,b,p2,ki[NSUB]={0.0},coop[NENZ]={0.0},ccc[NREAC][2]={0.0};
float fas[RES]={0.0},hh[RES]={0.0};
int eps[NREAC][NSUB]={0},n[NSUB]={0},tenz,zz;  
int nenz[NREAC]={0},tipenz[NREAC]={0},inh[NREAC][5]={-1};
int *cont;
int ***estado;
int **act;
double ***fase;
double DT;
int INTERV;

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
   double fa1,fa2,fa3,fb1,fb2,fb3;
   double prods,B,prodp,tiempo=0.0,deltat,dado,tau,sigma,Aext_prev,Aext,tau_min,pimp;
   double mm[NREAC] = {0.0};
   int n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,TOP,lineas=0;
   
   rand_init();

   // Abre archivos de salida
   //if((dptr=fopen("discrete.dat","w"))==NULL)
      //printf("\nEl archivo de salida no puede ser abierto\n");

   //if((fhist=fopen("histo.dat","w"))==NULL)
   //   printf("\nEl archivo de salida no puede ser abierto\n");

   // Inicia pgplot
   //cpgopen("/XSERVE");	
   //cpgask(0);

   // ****************************** DEFINICION DE PARAMETROS *******************************************
 
   a = 0.1;
   p2 = 0.999999;

   // Velocidades y afinidades
   // IMPORTANTE: Las afinidades de los sustratos deben ser ingresadas en el mismo orden en el que los
   // sustratos aparecen en la lista de las sustancias!
   
   // PESOS MOLARES: en Daltons
   mm[0] = 53000.;
   mm[1] = 63000.;
   mm[2] = 64000.;
   mm[3] = 39000.;
   mm[4] = 31000.;
   mm[5] = 39000.;
   mm[6] = 47000.;
   mm[7] = 64000.;
   mm[8] = 50000.;
   mm[9] = 100000.;
   
   /*mm[0] = 53000.;
   mm[1] = 63000.;
   mm[2] = 64000.;
   mm[3] = 39000.;
   mm[4] = 31000.;
   mm[5] = 39000.;
   mm[6] = 47000.;
   mm[7] = 64000.;
   mm[8] = 50000.;
   mm[9] = 100000.;*/
   
   // VELOCIDADES: micromol/(mg-min)
   vmax[0][0]=105.0;		// Vmax reaccion directa
   vmax[0][1]=-1.0;		// Vmax reaccion inversa

   vmax[1][0]=541.0;
   vmax[1][1]=284.0;
   //vmax[1][1]=-1.0;

   vmax[2][0]=298.0;
   vmax[2][1]=392.0;
   //vmax[2][1]=-1.0;

   vmax[3][0]=24.0;
   vmax[3][1]=29.0;
   //vmax[3][1]=-1.0;

   vmax[4][0]=284.0;
   vmax[4][1]=1697.0;
   //vmax[4][1]=-1.0;

   vmax[5][0]=27.0;
   vmax[5][1]=36.0;
   //vmax[5][1]=-1.0;

   vmax[6][0]=628.0;
   vmax[6][1]=87.0;
   //vmax[6][1]=-1.0;

   vmax[7][0]=42.0;
   vmax[7][1]=13.0;
   //vmax[7][1]=-1.0;

   vmax[8][0]=103.0;
   vmax[8][1]=26.0;
   //vmax[8][1]=-1.0;

   vmax[9][0]=8.0;
   vmax[9][1]=2.3;
   //vmax[9][1]=-1.0;
   
   // Las velocidades se expresan en sustrato/(segundo-enzima)
   for(j=0;j<NREAC;j++)
   {
      vmax[j][0] *= 1.e-3*mm[j]/60.0;		
      vmax[j][1] *= 1.e-3*mm[j]/60.0;		
      printf("\nLa velocidad de la reaccion %d directa es %lf [sustratos/segundo-enzima]", j, vmax[j][0]);
      printf("\nLa velocidad de la reaccion %d inversa es %lf [sustratos/segundo-enzima]", j, vmax[j][1]);
   }

   tau_min = 1.e6;
   for(j=0;j<NREAC;j++)
   {
      if((1./vmax[j][0]<tau_min) && (vmax[j][0]>0.0))
         tau_min = 1./vmax[j][0];	      
      if((1./vmax[j][1]<tau_min) && (vmax[j][1]>0.0))
         tau_min = 1./vmax[j][1];	      
   }	   
   printf("\n");
   printf("\nEl valor de tau minimo es %lf microsegundos",tau_min/1.e-6);

   //DT = tau_min/10.0;
   DT = 1.e-6;
   printf("\nEl paso de tiempo es %lf microsegundos",DT/1.e-6);

   // Tiempo de simulacion
   TOP = (int)(TTOP/DT);

   // Ajusta el numero de puntos a escribir
   if(TOP>1e6)
      INTERV=TOP/1e6;
   else 
      INTERV=1;	   

   // Verifica que no se escriban mas de un millon de puntos
   if(TOP/INTERV>1000000)
   {	   
      printf("\nCUIDADO! SE ESCRIBEN %d PUNTOS\n",TOP); 	   
      sleep(3);
   }
   else
      printf("\nEl numero de puntos a calcular es %d\n",TOP);	   
      printf("\nEl numero de puntos a escribir es %d\n",TOP/INTERV);	   

   printf("\n");

   // AFINIDADES: microMolar	
   km[0][0]=33.0;		
   km[0][4]=84.0;

   km[2][7]=750;
   km[3][8]=130.0;

   km[4][4]=50.0;
   km[4][8]=455;
   km[5][2]=1440.0;
   km[5][9]=124.0;

   km[6][9]=4.0;
   km[7][10]=105.0;
   km[7][11]=108.0;

   km[8][10]=1400.0;
   km[9][11]=830.0;

   km[10][1]=33.0;
   km[10][3]=59.0;
   km[10][11]=33.0;
   km[11][5]=10.0;
   km[11][12]=246.0;

   km[12][2]=3400;
   km[12][12]=127.0;
   km[13][4]=3300.0;
   km[13][13]=570.0;

   km[14][13]=500.0;
   km[15][14]=66.0;

   km[16][14]=55.0;
   //km[17][6]=0.0;
   km[17][6]=20.0;
   km[17][15]=63.0;

   km[18][2]=20.0;
   km[18][15]=24.0;
   km[19][4]=68.0;
   km[19][16]=68.0;

   // Las afinidades se expresan en numero de sustratos
   for(j=0;j<2*NREAC;j++)
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
   ccc[0][0] = 0.01;		// Reaccion directa
   ccc[0][1] = 0.01;		// Reaccion inversa
   ccc[1][0] = 0.01;		// Reaccion directa
   ccc[1][1] = 0.01;		// Reaccion inversa
   ccc[2][0] = 0.01;		// Reaccion directa
   ccc[2][1] = 0.01;		// Reaccion directa
   ccc[3][0] = 0.01;		// Reaccion inversa
   ccc[3][1] = 0.01;		// Reaccion inversa
   ccc[4][0] = 0.01;		// Reaccion directa
   ccc[4][1] = 0.01;		// Reaccion inversa
   ccc[5][0] = 0.01;		// Reaccion directa
   ccc[5][1] = 0.01;		// Reaccion inversa
   ccc[6][0] = 0.01;		// Reaccion directa
   ccc[6][1] = 0.01;		// Reaccion inversa
   ccc[7][0] = 0.01;		// Reaccion directa
   ccc[7][1] = 0.01;		// Reaccion inversa
   ccc[8][0] = 0.01;		// Reaccion directa
   ccc[8][1] = 0.01;		// Reaccion inversa
   ccc[9][0] = 0.01;		// Reaccion directa
   ccc[9][1] = 0.01;		// Reaccion inversa

   // Numero de enzimas que catalizan cada reaccion
   nenz[0] = NENZ;
   nenz[1] = NENZ;
   nenz[2] = NENZ;
   nenz[3] = NENZ;
   nenz[4] = NENZ;
   nenz[5] = NENZ;
   nenz[6] = NENZ;
   nenz[7] = NENZ;
   nenz[8] = NENZ;
   nenz[9] = NENZ;

   // Tipo de enzima que cataliza cada reaccion
   tipenz[0] = 2;
   tipenz[1] = 2;
   tipenz[2] = 2;
   tipenz[3] = 4;
   tipenz[4] = 2;
   tipenz[5] = 4;
   tipenz[6] = 2;
   tipenz[7] = 1;
   tipenz[8] = 4;
   tipenz[9] = 4;
   
   /*tipenz[0] = 1;
   tipenz[1] = 1;
   tipenz[2] = 1;
   tipenz[3] = 1;
   tipenz[4] = 1;
   tipenz[5] = 1;
   tipenz[6] = 1;
   tipenz[7] = 1;
   tipenz[8] = 1;
   tipenz[9] = 1;*/

   // Reparte el peso molecular de la enzima entre las subunidades
   for(j=0;j<NREAC;j++)
   {
      if(tipenz[j]==2)
         mm[j]/=2.0;	      
      if(tipenz[j]==3)
         mm[j]/=4.0;	      
      if(tipenz[j]==4)
         mm[j]/=4.0;	      
   }	   

   // Cooperatividad de los oligomeros
   coop[0] = 10.0;
   coop[1] = 10.0;
   coop[2] = 10.0;
   coop[3] = 10.0;
   coop[4] = 10.0;
   coop[5] = 10.0;
   coop[6] = 10.0;
   coop[7] = 10.0;
   coop[8] = 10.0;
   coop[9] = 10.0;

   // Definicion de inhibiciones
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
   
   inh[2][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[2][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[2][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[2][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[2][4] = -1;		// Total

   inh[3][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[3][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[3][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[3][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[3][4] = -1;		// Total

   inh[4][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[4][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[4][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[4][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[4][4] = -1;		// Total

   inh[5][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[5][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[5][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[5][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[5][4] = -1;		// Total

   inh[6][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[6][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[6][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[6][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[6][4] = -1;		// Total

   inh[7][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[7][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[7][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[7][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[7][4] = -1;		// Total

   inh[8][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[8][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[8][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[8][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[8][4] = -1;		// Total

   inh[9][0] = -1;		// Substancia inhibidora, reaccion directa, ic
   inh[9][1] = -1;		// Substancia inhibida, reaccion inversa, ic	
   inh[9][2] = -1;		// Substancia inhibidora, reaccion directa, inc
   inh[9][3] = -1;		// Substancia inhibida, reaccion inversa,inc
   inh[9][4] = -1;		// Total
   
   // La reaccion es inhibida?
   for(j=0;j<NREAC;j++)
   {
      inh[j][4] = 0;   
      for(k=0; k<4; k++)
         inh[j][4] += inh[j][k];
   }

   // Afinidades para inhibidores
   ki[0] = 100.0;
   ki[1] = 100.0;
   ki[2] = 100.0;
   ki[3] = 100.0;
   ki[4] = 100.0;
   ki[5] = 100.0;
   ki[6] = 100.0;
   ki[7] = 100.0;
   ki[8] = 100.0;
   ki[9] = 100.0;
   ki[10] = 100.0;
   ki[11] = 100.0;
   ki[12] = 100.0;
   ki[13] = 100.0;
   ki[14] = 100.0;
   ki[15] = 100.0;
   ki[16] = 100.0;
   
   // Las afinidades para los inhibidores se expresan en numero de sustratos
   for(j=0;j<NSUB;j++)
      ki[j]*=NA*1.e-6*VSIM;	      

   // Halla numero total de enzimas
   tenz = 0; 
   for(j=0;j<NREAC;j++)
      tenz+=nenz[j];

   // Cantidades iniciales de sustrato 
   
   //n[0] = (int)NINJ*1.e-6*NA*VSIM;		// S0	// Cantidades minimas para completar la reaccion
   n[0]=1;		// Glucosa	// Cantidades minimas para completar la reaccion
   n[1]=2;		// Pi
   n[2]=2;		// ADP
   n[3]=2;		// NADox
   n[4]=2;		// ATP
   n[5]=0;		// NADred	
   n[6]=0;		// H2O
   n[7]=0;   		// Glucosa-6-fosfato
   n[8]=0;		// Fructosa 6-fosfato
   n[9]=0;		// Fructosa 1,6 bifosfato
   n[10]=0;		// Glicerona fosfato
   n[11]=0;    		// Gliceraldehido 3-fosfato
   n[12]=0;		// 3-fosfo-D-glicerol fosfato
   n[13]=0;		// 3-fosfo-D-glicerato
   n[14]=0;		// 2-fosfo-D-glicerato
   n[15]=0;   		// Fosfenolpiruvato
   n[16]=0;   		// Piruvato

   // Escala cantidades iniciales de sustratos
   for(j=0;j<NSUB;j++)
      n[j]*=50000;
      //n[j]*=10000;

   printf("\nLas concentraciones iniciales son:");
   printf("\nGLUCOSA: %lf microMolar",(double)n[0]*1.e6/(NA*VSIM));
   printf("\nPi: %lf microMolar",(double)n[1]*1.e6/(NA*VSIM));
   printf("\nADP: %lf microMolar",(double)n[2]*1.e6/(NA*VSIM));
   printf("\nATP: %lf microMolar",(double)n[4]*1.e6/(NA*VSIM));
   printf("\nNADox: %lf microMolar",(double)n[3]*1.e6/(NA*VSIM));
   printf("\n");

   // Se definen reacciones
   eps[0][0]=-1;  
   eps[0][1]=0;  
   eps[0][2]=1;  
   eps[0][3]=0;  
   eps[0][4]=-1;  
   eps[0][5]=0;  
   eps[0][6]=0;  
   eps[0][7]=1;  
   eps[0][8]=0;  
   eps[0][9]=0;  
   eps[0][10]=0;  
   eps[0][11]=0;  
   eps[0][12]=0;  
   eps[0][13]=0;  
   eps[0][14]=0;  
   eps[0][15]=0;  
   eps[0][16]=0;  
 
   eps[1][0]=0;  
   eps[1][1]=0;  
   eps[1][2]=0;  
   eps[1][3]=0;  
   eps[1][4]=0;  
   eps[1][5]=0;  
   eps[1][6]=0;  
   eps[1][7]=-1;  
   eps[1][8]=1;  
   eps[1][9]=0;  
   eps[1][10]=0;  
   eps[1][11]=0;  
   eps[1][12]=0;  
   eps[1][13]=0;  
   eps[1][14]=0;  
   eps[1][15]=0;  
   eps[1][16]=0;  

   eps[2][0]=0;  
   eps[2][1]=0;  
   eps[2][2]=1;  
   eps[2][3]=0;  
   eps[2][4]=-1;  
   eps[2][5]=0;  
   eps[2][6]=0;  
   eps[2][7]=0;  
   eps[2][8]=-1;  
   eps[2][9]=1;  
   eps[2][10]=0;  
   eps[2][11]=0;  
   eps[2][12]=0;  
   eps[2][13]=0;  
   eps[2][14]=0;  
   eps[2][15]=0;  
   eps[2][16]=0;  

   eps[3][0]=0;  
   eps[3][1]=0;  
   eps[3][2]=0;  
   eps[3][3]=0;  
   eps[3][4]=0;  
   eps[3][5]=0;  
   eps[3][6]=0;  
   eps[3][7]=0;  
   eps[3][8]=0;  
   eps[3][9]=-1;  
   eps[3][10]=1;  
   eps[3][11]=1;  
   eps[3][12]=0;  
   eps[3][13]=0;  
   eps[3][14]=0;  
   eps[3][15]=0;  
   eps[3][16]=0;  

   eps[4][0]=0;  
   eps[4][1]=0;  
   eps[4][2]=0;  
   eps[4][3]=0;  
   eps[4][4]=0;  
   eps[4][5]=0;  
   eps[4][6]=0;  
   eps[4][7]=0;  
   eps[4][8]=0;  
   eps[4][9]=0;  
   eps[4][10]=-1;  
   eps[4][11]=1;  
   eps[4][12]=0;  
   eps[4][13]=0;  
   eps[4][14]=0;  
   eps[4][15]=0;  
   eps[4][16]=0;  

   eps[5][0]=0;  
   eps[5][1]=-1;  
   eps[5][2]=0;  
   eps[5][3]=-1;  
   eps[5][4]=0;  
   eps[5][5]=1;  
   eps[5][6]=0;  
   eps[5][7]=0;  
   eps[5][8]=0;  
   eps[5][9]=0;  
   eps[5][10]=0;  
   eps[5][11]=-1;  
   eps[5][12]=1;  
   eps[5][13]=0;  
   eps[5][14]=0;  
   eps[5][15]=0;  
   eps[5][16]=0;  

   eps[6][0]=0;  
   eps[6][1]=0;  
   eps[6][2]=-1;  
   eps[6][3]=0;  
   eps[6][4]=1;  
   eps[6][5]=0;  
   eps[6][6]=0;  
   eps[6][7]=0;  
   eps[6][8]=0;  
   eps[6][9]=0;  
   eps[6][10]=0;  
   eps[6][11]=0;  
   eps[6][12]=-1;  
   eps[6][13]=1;  
   eps[6][14]=0;  
   eps[6][15]=0;  
   eps[6][16]=0;  

   eps[7][0]=0;  
   eps[7][1]=0;  
   eps[7][2]=0;  
   eps[7][3]=0;  
   eps[7][4]=0;  
   eps[7][5]=0;  
   eps[7][6]=0;  
   eps[7][7]=0;  
   eps[7][8]=0;  
   eps[7][9]=0;  
   eps[7][10]=0;  
   eps[7][11]=0;  
   eps[7][12]=0;  
   eps[7][13]=-1;  
   eps[7][14]=1;  
   eps[7][15]=0;  
   eps[7][16]=0;  

   eps[8][0]=0;  
   eps[8][1]=0;  
   eps[8][2]=0;  
   eps[8][3]=0;  
   eps[8][4]=0;  
   eps[8][5]=0;  
   eps[8][6]=1;  
   eps[8][7]=0;  
   eps[8][8]=0;  
   eps[8][9]=0;  
   eps[8][10]=0;  
   eps[8][11]=0;  
   eps[8][12]=0;  
   eps[8][13]=0;  
   eps[8][14]=-1;  
   eps[8][15]=1;  
   eps[8][16]=0;  

   eps[9][0]=0;  
   eps[9][1]=0;  
   eps[9][2]=-1;  
   eps[9][3]=0;  
   eps[9][4]=1;  
   eps[9][5]=0;  
   eps[9][6]=0;  
   eps[9][7]=0;  
   eps[9][8]=0;  
   eps[9][9]=0;  
   eps[9][10]=0;  
   eps[9][11]=0;  
   eps[9][12]=0;  
   eps[9][13]=0;  
   eps[9][14]=0;  
   eps[9][15]=-1;  
   eps[9][16]=1;  

   // ************************* RESERVA DE MEMORIA *********************************************

   // Reserva memoria para vector "cont"
   if((cont=(int*)calloc(tenz,sizeof(int)))==NULL)
   {
      printf("\nERROR: No hay memoria para el vector 'cont'!");
      exit(1);
   }

   // Reserva memoria para la matriz "estado"
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
   
   // Modificacion de las afinidades de los oligomeros
   for(j=0;j<2*NREAC;j++)
   {
      if(j%2==0) 		// Numero de la reaccion
         cdir=j/2; 

      // Modificacion de km 
      for(k=0;k<NSUB;k++)
      {
         if(tipenz[cdir] == 2)			// Dimero
            km[j][k]*=sqrt(coop[cdir]); 

         if(tipenz[cdir]==3)			// Trimero 
            km[j][k]*=1.05*coop[cdir];

         if(tipenz[cdir]==4) 
            km[j][k]*=pow(coop[cdir],1.53);	// Tetramero
      }
   }

   // Modificacion de las afinidades de los inhibidores
   for(j=0;j<NREAC;j++)
   {
      // Modificacion de ki, inhibicion competitiva
      if(inh[j][0]!=-1)				// La reaccion directa es inhibida
      {
         if(tipenz[cdir]==2)			// Dimero
            ki[inh[j][0]]*=sqrt(coop[cdir]); 

         if(tipenz[cdir]==3)			// Trimero 
            ki[inh[j][0]]*=1.05 * coop[cdir];

         if(tipenz[cdir]==4) 			// Tetramero
            ki[inh[j][0]]*=pow(coop[cdir],1.53);
      }
      if(inh[j][1]!=-1)				// La reaccion inversa es inhibida
      {
         if(tipenz[cdir]==2)			// Dimero
            ki[inh[j][1]]*=sqrt(coop[cdir]); 

         if(tipenz[cdir]==3)			// Trimero 
            ki[inh[j][1]]*=1.05 * coop[cdir];

         if(tipenz[cdir]==4) 			// Tetramero
            ki[inh[j][1]]*=pow(coop[cdir],1.53);
      }

      // Modificacion de ki, inhibicion no competitiva
      if(inh[j][2]!=-1)				// Reaccion directa es inhibida
      {
         if(tipenz[cdir]==2)			// Dimero
            ki[inh[j][2]]*=sqrt(coop[cdir]); 

         if(tipenz[cdir]==3)			// Trimero 
            ki[inh[j][2]]*=1.05 * coop[cdir];

         if(tipenz[cdir]==4) 			// Tetramero
            ki[inh[j][2]]*=pow(coop[cdir],1.53);
      }
      if(inh[j][3]!=-1)				// Reaccion inversa es inhibida
      {
         if(tipenz[cdir]==2)			// Dimero
            ki[inh[j][3]]*=sqrt(coop[cdir]); 

         if(tipenz[cdir]==3)			// Trimero 
            ki[inh[j][3]]*=1.05*coop[cdir];

         if(tipenz[cdir] == 4) 			// Tetramero
            ki[inh[j][3]]*=pow(coop[cdir],1.53);
       } 
   }  

   // Calcula los parametros "b" y "phic" para cada reaccion, ambas direcciones
   for(j=0;j<NREAC;j++)
   {
      bdir[j]=a/(vmax[j][0]*DT)-3.0*a;
      printf("\nEl numero de iteraciones para la r. %d dir. es %f, el No. de ciclos/iteracion es %lf",j,(bdir[j]/a)+3,(double)TOP/((bdir[j]/a)+3)); 
      //phicdir[j]=1.0+bdir[j]*ccc[j][0];
      phicdir[j]=1.0+bdir[j]*0.01;
      printf("\nphicdir[%d] es %lf",j,phicdir[j]); 
       
      if(vmax[j][1]<0.0)			// Reaccion inversa que no se da!
         binv[j]=0.0;
      else if(vmax[j][1]>0.0)   		// Si la reaccion inversa se da...	
         binv[j]=a/(vmax[j][1]*DT)-3.0*a;
      //phicinv[j]=1.0+binv[j]*ccc[j][1];
      phicinv[j]=1.0+binv[j]*0.01;
      if(vmax[j][1]>0.0)   		// Si la reaccion inversa se da...	
         printf("\nEl numero de iteraciones para la r. %d inv. es %f, el No. de ciclos/iteracion es %lf",j,(binv[j]/a)+3,(double)TOP/((binv[j]/a)+3)); 
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
   //   hh[j]=0.0;
   //hist(bdir[0]+1);
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
   //fprintf(dptr,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo,(double)n[0]*1.e6/(NA*VSIM),(double)n[1]*1.e6/(NA*VSIM),(double)n[2]*1.e6/(NA*VSIM),(double)n[3]*1.e6/(NA*VSIM),(double)n[4]*1.e6/(NA*VSIM),(double)n[5]*1.e6/(NA*VSIM),(double)n[6]*1.e6/(NA*VSIM),(double)n[7]*1.e6/(NA*VSIM),(double)n[8]*1.e6/(NA*VSIM),(double)n[9]*1.e6/(NA*VSIM),(double)n[10]*1.e6/(NA*VSIM),(double)n[11]*1.e6/(NA*VSIM),(double)n[12]*1.e6/(NA*VSIM),(double)n[13]*1.e6/(NA*VSIM),(double)n[14]*1.e6/(NA*VSIM),(double)n[15]*1.e6/(NA*VSIM),(double)n[16]*1.e6/(NA*VSIM)); 
   //
   if((dptr=fopen("discrete.dat","w"))==NULL)
      printf("\nEl archivo de salida no puede ser abierto\n");
   
   fprintf(dptr,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo,(double)n[0]*1.e6/(NA*VSIM*1.e3),(double)n[1]*1.e6/(NA*VSIM*1.e3),(double)n[2]*1.e6/(NA*VSIM*1.e3),(double)n[3]*1.e6/(NA*VSIM*1.e3),(double)n[4]*1.e6/(NA*VSIM*1.e3),(double)n[5]*1.e6/(NA*VSIM*1.e3),(double)n[6]*1.e6/(NA*VSIM*1.e3),(double)n[7]*1.e6/(NA*VSIM*1.e3),(double)n[8]*1.e6/(NA*VSIM*1.e3),(double)n[9]*1.e6/(NA*VSIM*1.e3),(double)n[10]*1.e6/(NA*VSIM*1.e3),(double)n[11]*1.e6/(NA*VSIM*1.e3),(double)n[12]*1.e6/(NA*VSIM*1.e3),(double)n[13]*1.e6/(NA*VSIM*1.e3),(double)n[14]*1.e6/(NA*VSIM*1.e3),(double)n[15]*1.e6/(NA*VSIM*1.e3),(double)n[16]*1.e6/(NA*VSIM*1.e3)); 
   
   fclose(dptr);

   lineas++;
   
//   fprintf(dptr,"%lf \t %d \t %d \t %lf \t %lf \n",tiempo/tau, n[0], n[1], (double)n[0]*1.e6/(NA*VSIM), (double)n[1]*1.e6/(NA*VSIM)); 
   tiempo+=DT;
  
   //***************************************** SIMULACION ************************************
   
   // Se agrega sustrato gradualmente
   sigma=13.0*pow(tau,2.0);
   Aext_prev = 0.0;
   A_acc = 0;
   flag = 0;
   flag2=0;
   //n[0] = 0;
   ntop = TOP*tau/DT;

   //printf("\nEL VALOR DE NTOP ES %d\n",ntop);

   // Ciclo de iteracion de las enzimas
   for(zz=0;zz<TOP;zz++)
   {
      /*if(flag==0)
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
      }*/

      /*if((tiempo/tau > 3.319) && (tiempo/tau < 3.321))
      {
         printf("\nCalculando histograma :)");     
         // Calcula histograma 
         for(j=0;j<RES;j++)
            hh[j]=0.0;
         hist(bdir[1]+1);
         flag2 = 1;
      } */	 

      // Baraja vector de recorrido de las enzimas
      //baraja(zz); 

      pimp = n[0];
      pimp = n[1];
      pimp = n[2];
      pimp = n[3];
      pimp = n[4];
      pimp = n[5];
      pimp = n[6];
      pimp = n[7];
      pimp = n[8];
      pimp = n[9];
      pimp = n[10];
      pimp = n[10];
      pimp = n[12];
      pimp = n[12];
      pimp = n[14];
      pimp = n[15];
      pimp = n[16];
  
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

	 if(reac==0)
	    pimp = 0.0;	 
	 if(reac==1)
	    pimp = 0.0;	 
	 if(reac==2)
	    pimp = 0.0;	 
	 if((reac==3)&&(n[9]!=0))
	    pimp = 0.0;	 
	 if(reac==4)
	    pimp = 0.0;	 
	 if(reac==5)
	    pimp = 0.0;	 
	 if(reac==6)
	    pimp = 0.0;	 
	 if(reac==7)
	    pimp = 0.0;	 
	 if(reac==8)
	    pimp = 0.0;	 
	 if(reac==9)
	    pimp = 0.0;	 

         // Calcula probabilidad de acuerdo al sustrato disponible
         probab(); 
	 
         n0 = n[0];
         n1 = n[1];
         n2 = n[2];
         n3 = n[3];
         n4 = n[4];
         n5 = n[5];
         n6 = n[6];
         n7 = n[7];
         n8 = n[8];
         n9 = n[9];
         n10 = n[10];
         n11 = n[11];
         n12 = n[12];
         n13 = n[13];
         n14 = n[14];
         n15 = n[15];
         n16 = n[16];

         // Monomero
         if(tipenz[reac]==1)
         {
            if(inh[reac][4]==-4)					// Sin inhibicion
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)||estado[reac][enz][0]>=2) 		// La reaccion se da
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
                  //if(abs(estado[reac][enz][k])>0)  
                  if(abs(estado[reac][enz][k])>=2)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)||(estado[reac][enz][j]>=2)) 		// La reaccion se da
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
                  //if(abs(estado[reac][enz][k])>0)  
                  if(abs(estado[reac][enz][k])>=2)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)||(estado[reac][enz][j]>=2)) 		// La reaccion se da
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
                  //if(abs(estado[reac][enz][k])>0)  
                  if(abs(estado[reac][enz][k])>=2)  
                     act[reac][enz]++;
               }
               if((p1[0][reac]!=0.0)||(p1[1][reac]!=0.0)||(estado[reac][enz][j]>=2)) 		// La reaccion se da
                  fase[reac][enz][j]=mapa_s(reac,enz,j,fase[reac][enz][j]);    	
            }
         }  
      }

      // Registra proporciones de enzimas en reposo, procesando y en recuperacion
      /*fa1 = 0.;
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
      }
      
      fa1 /= (double)NENZ;
      fa2 /= (double)NENZ;
      fa3 /= (double)NENZ;
      fb1 /= (double)NENZ;
      fb2 /= (double)NENZ;
      fb3 /= (double)NENZ;*/

      // Registra concentraciones de sustancias
      if(zz%INTERV==0)
      {	      
         if((dptr=fopen("discrete.dat","a"))==NULL)
            printf("\nEl archivo de salida no puede ser abierto\n");

         fprintf(dptr,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo,(double)n[0]*1.e6/(NA*VSIM*1.e3),(double)n[1]*1.e6/(NA*VSIM*1.e3),(double)n[2]*1.e6/(NA*VSIM*1.e3),(double)n[3]*1.e6/(NA*VSIM*1.e3),(double)n[4]*1.e6/(NA*VSIM*1.e3),(double)n[5]*1.e6/(NA*VSIM*1.e3),(double)n[6]*1.e6/(NA*VSIM*1.e3),(double)n[7]*1.e6/(NA*VSIM*1.e3),(double)n[8]*1.e6/(NA*VSIM*1.e3),(double)n[9]*1.e6/(NA*VSIM*1.e3),(double)n[10]*1.e6/(NA*VSIM*1.e3),(double)n[11]*1.e6/(NA*VSIM*1.e3),(double)n[12]*1.e6/(NA*VSIM*1.e3),(double)n[13]*1.e6/(NA*VSIM*1.e3),(double)n[14]*1.e6/(NA*VSIM*1.e3),(double)n[15]*1.e6/(NA*VSIM*1.e3),(double)n[16]*1.e6/(NA*VSIM*1.e3)); 

         fclose(dptr);
	 lineas++;
      }	 
         //fprintf(dptr,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tiempo , (double)n[0]*1.e6/(NA*VSIM), (double)n[1]*1.e6/(NA*VSIM), (double)n[2]*1.e6/(NA*VSIM), (double)n[3]*1.e6/(NA*VSIM), (double)n[4]*1.e6/(NA*VSIM), (double)n[5]*1.e6/(NA*VSIM), (double)n[6]*1.e6/(NA*VSIM), (double)n[7]*1.e6/(NA*VSIM), (double)n[8]*1.e6/(NA*VSIM), (double)n[9]*1.e6/(NA*VSIM), (double)n[10]*1.e6/(NA*VSIM), (double)n[11]*1.e6/(NA*VSIM), (double)n[12]*1.e6/(NA*VSIM), (double)n[13]*1.e6/(NA*VSIM), (double)n[14]*1.e6/(NA*VSIM), (double)n[15]*1.e6/(NA*VSIM), (double)n[16]*1.e6/(NA*VSIM)); 

      // Retira piruvato
      if(n[16]>0)
         n[16]--;

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
   //fclose(dptr);

   // Conserva histograma final
   //for(j=0;j<RES;j++)
   //   fprintf(fhist,"%f \t %f \n",fas[j],hh[j]);
   //fclose(fhist);

   // Termina pgplot
   //cpgend();
   //return(0);
   printf("\nSe escribieron %d lineas",lineas);
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
   double pimp;

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
      pimp = phicdir[pos1];
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

   /*if(x<0.0) 
      out = a+b+2+2*x;
   else if((x>=0)&(x<0.5-0.5*p1)) 
      out = m1*x;
   else if((x>=0.5-0.5*p1)&(x<0.5)) 
      out = 1+a-0.5*m2+m2*x;
   else if((x>=0.5)&(x<0.5+0.5*p1)) 
      out = 0.5*m2*(1+p1)+1-m2*x;
   else if((x>=0.5+0.5*p1)&(x<=1))  
      out = m1*(1-x);   
   else if((x>=1)&(x<(1+b))) 
      out = a+x;
   else if(x>=(1+b)) 
      out = (x-2.)/a;*/
   
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
   int gob,m,flag=1;

   gob=1;					// Verifica que los sustratos esten disponibles (gob=1)
   for(m=0;m<NSUB;m++)
   {
      if(eps[pos1][m]>0)			// Si sustancia es sustrato de la reaccion inversa...
      {
         if(km[2*pos1+1][m]!=0.0)		// La sustancia participa en la reaccion inversa
         {  
            if(n[m]>=eps[pos1][m]) 		// Si TODOS los sustratos estan disponibles...
            {		    
               gob*=1;				// ...la reaccion procede
	       flag=0;
	    }   
            else
               gob*=0;                        
         }
      }  
   }
   if((gob==1)&&(flag==1))			// La reaccion no se da
      gob=0;	   
   
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

