/* 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	 

   glicobis2.c
   Una red en las que hay varios grupos de enzimas trabajando en paralelo.
   * Las enzimas son recorridas aleatoriamente para su actualizacion; de esta 
   manera no se favorece a ninguna reaccion. 
   * El recorrido de las enzimas es 
   realizado leyendo secuencialmente las posiciones (reaccion-enzima) contenidas 
   en un vector previamente generado. 
   * Cada substancia tiene un TMAX distinto (afinidad) dependiendo de la direccion 
   de la reaccion en la que esta involucrada. 
   * Cada reaccion se define una vez, el codigo automaticamente toma en cuenta la reaccion inversa.
   * Se representa la glicolis: la matriz de reacciones es leida desde un archivo
   * La secuencia de las reacciones fue comprobada
   * Se limita a la glicolisis, el ciclo de Krebs se excluye
   * El valor x = 0.5 fue excluido de los mapas debido a problemas a bajas probabilidades
   (ver histograma del mapa tienda)
   * Los parametros (Vmax, K) de las enzimas se leen desde un archivo
   * La estructura del programa es modular: La actualizacion de las enzimas y el proceso de toma de sustrato
     y liberacion de producto se realizan desde una funcion (8/12/2005)
   * La cantidad de enzimas varia de una enzima a otra (9/12/2005)
   * Se sigue idea de A Parravano para el arranque de una enzima en direccion directa o reversa
   * Se modifico rutina "monomer"
   * Se incluyen oligomeros
   * Se incluye re-definicion de las afinidades para oligomeros
   * Una unidad de substancia = 1 microMolar
   * Una iteracion = 0.6 milisegundos (scal = 1e-6) (8/6/2006)
   * Se incluye re-definicion de las velocidades para oligomeros (13/02/06)
   * Se incluye inhibicion competitiva (14/02/06)
   * Se cambio la estructura de datos de las afinidades (22/02/06)

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<time.h>

#define NSUB 17
#define NREAC 10
#define TOP 300000
#define NENZ 10

// Variables globales
FILE *dptr;
double p1[2][NREAC]={0.0},vmax[NREAC][2]={0.0},km[2*NREAC][NSUB]={0.0},bdir[NREAC]={0.0},binv[NREAC]={0.0};
double phicdir[NREAC]={0.0},phicinv[NREAC]={0.0},a,b,p2,ki[NSUB]={0.0};
int eps[NREAC][NSUB]={0},n[NSUB]={0},tenz;  
int nenz[NREAC]={0},tipenz[NREAC]={0},inhib[2*NREAC][2]={-1};
int *cont;
int ***ocup;
int **act;

void rand_init();
void baraja();
void probab();
void sustprod();
double monomer();
double enzima47();
double avandir();
double avaninv();
int drnd();

main()
{
   // Declaraciones 
   FILE *dptr;
   //FILE *dptt;
   //FILE *dptq;
   int i,j=0,k=0,z,reac,enz,prod=0,m,cdir,cinv,index1,index2;
   double scal,ak,bk,ck,time=0.0;
   double ***xaux;
   
   //rand_init();
   //srand48(43878);
   srand48(34200);

   // Abre archivos de salida
   if((dptr=fopen("glico.dat","w"))==NULL)
      printf("\nEl archivo de salida no puede ser abierto\n");

   // Abre para lectura archivo de reacciones
   //if((dptt=fopen("glico2.rec","r"))==NULL)
   //   printf("\nEl archivo de reacciones no puede ser abierto\n");

   // Abre para lectura archivo con parametros de las enzimas
   //if((dptq=fopen("glico2.par","r"))==NULL)
   //   printf("\nEl archivo de salida no puede ser abierto\n");

   // ****************************** DEFINICION DE PARAMETROS *******************************************
 
   a=0.1;
   p2=0.99;

   // Cooperatividad de los multimeros
   ak=10.0;
   bk=ak;
   ck=ak;

   // Lee archivo con afinidades de las enzimas
   //for(j=0;j<2*NREAC;j++)
   //{  
   //   for(k=0;k<NSUB;k++)
   //      fscanf(dptq,"%lf",&km[j][k]);    
   //}  
   //fclose(dptq);

   //for(j=0;j<NREAC;j++)
   //   printf("\nPara la reaccion %d, phicdir es %lf, phicinv es %lf",j,phicdir[j],phicinv[j]);

   // Lee archivo de reacciones
   //for(j=0;j<NREAC;j++)
   //{
   //   for(k=0;k<NSUB;k++)
   //      fscanf(dptt,"%d",&eps[j][k]);    
   //}
   //fclose(dptt);

   // Velocidades y afinidades
   // IMPORTANTE: Las afinidades de los sustratos deben ser ingresadas en el mismo orden en el que los
   // sustratos aparecen en la lista de las sustancias!

   // VELOCIDADES
   vmax[0][0]=105.0;		// Vmax reaccion directa
   vmax[0][1]=-1.0;		// Vmax reaccion inversa

   vmax[1][0]=541.1;
   vmax[1][1]=284.0;

   vmax[2][0]=298.0;
   vmax[2][1]=392.0;

   vmax[3][0]=28.0;
   vmax[3][1]=29.0;

   vmax[4][0]=284.0;
   vmax[4][1]=1697.0;

   vmax[5][0]=27.0;
   vmax[5][1]=36.0;

   vmax[6][0]=628.0;
   vmax[6][1]=87.0;

   vmax[7][0]=42.0;
   vmax[7][1]=13.0;

   vmax[8][0]=103.0;
   vmax[8][1]=26.0;

   vmax[9][0]=8.0;
   vmax[9][1]=2.3;

   // AFINIDADES	
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
   km[17][6]=0.0;
   km[17][15]=63.0;

   km[18][2]=20.0;
   km[18][15]=24.0;
   km[19][4]=68.0;
   km[19][16]=68.0;

   // Numero de enzimas que catalizan cada reaccion
   nenz[0]=NENZ;
   nenz[1]=NENZ;
   nenz[2]=NENZ;
   nenz[3]=NENZ;
   nenz[4]=NENZ;
   nenz[5]=NENZ;
   nenz[6]=NENZ;
   nenz[7]=NENZ;
   nenz[8]=NENZ;
   nenz[9]=NENZ;

   // Tipo de enzima que cataliza cada reaccion
   tipenz[0]=2;
   tipenz[1]=2;
   tipenz[2]=2;
   tipenz[3]=4;
   tipenz[4]=2;
   tipenz[5]=4;
   tipenz[6]=2;			
   tipenz[7]=1;
   tipenz[8]=4;
   tipenz[9]=4;

   // Inhibidores
   //inhib[0][0]=-1;	// Substancia inhibidora
   //inhib[0][1]=-1;	// Substancia inhibida	
   inhib[0][0]=2;	// Substancia inhibidora
   inhib[0][1]=4;	// Substancia inhibida	

   inhib[1][0]=-1;
   inhib[1][1]=-1;

   inhib[2][0]=-1;
   inhib[2][1]=-1;

   inhib[3][0]=-1;
   inhib[3][1]=-1;

   inhib[4][0]=-1;
   inhib[4][1]=-1;

   inhib[5][0]=-1;
   inhib[5][1]=-1;

   inhib[6][0]=-1;
   inhib[6][1]=-1;

   inhib[7][0]=-1;
   inhib[7][1]=-1;

   inhib[8][0]=-1;
   inhib[8][1]=-1;

   inhib[9][0]=-1;
   inhib[9][1]=-1;

   inhib[10][0]=-1;
   inhib[10][1]=-1;

   inhib[11][0]=-1;
   inhib[11][1]=-1;

   inhib[12][0]=-1;
   inhib[12][1]=-1;

   inhib[13][0]=-1;
   inhib[13][1]=-1;

   inhib[14][0]=-1;
   inhib[14][1]=-1;

   inhib[15][0]=-1;
   inhib[15][1]=-1;

   inhib[16][0]=-1;
   inhib[16][1]=-1;

   inhib[17][0]=-1;
   inhib[17][1]=-1;

   inhib[18][0]=-1;
   inhib[18][1]=-1;

   inhib[19][0]=-1;
   inhib[19][1]=-1;

   // Afinidades para inhibidores
   ki[0]=0.0;
   ki[1]=0.0;
   ki[2]=60000.0;
   ki[3]=0.0;
   ki[4]=0.0;
   ki[5]=0.0;
   ki[6]=0.0;
   ki[7]=0.0;
   ki[8]=0.0;
   ki[9]=0.0;
   ki[10]=0.0;
   ki[11]=0.0;
   ki[12]=0.0;
   ki[13]=0.0;
   ki[14]=0.0;
   ki[15]=0.0;
   ki[16]=0.0;

   // Halla numero total de enzimas
   tenz=0; 
   for(j=0;j<NREAC;j++)
      tenz+=nenz[j];

   // Cantidades iniciales de sustrato 
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
      n[j]*=200;

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

   // Reserva memoria para la matriz "ocup"
   if((ocup=(int***)calloc(NREAC,sizeof(int**)))==NULL)
   {
     printf("\nNo hay memoria para la matriz 'ocup'");
     exit(1);
   }
   for(j=0;j<NREAC;j++)
   {
     if((ocup[j]=(int**)calloc(nenz[j],sizeof(int*))) == NULL)
     {
       printf("\nNo hay memoria para la matriz 'ocup'");
       exit(1);
     }
   }
   for(j=0;j<NREAC;j++)
   {
      for(k=0;k<nenz[j];k++)
      {
         if((ocup[j][k]=(int*)calloc(tipenz[j],sizeof(int))) == NULL)
         {
            printf("\nNo hay memoria para la matriz 'ocup'");
            exit(1);
         }
      } 
   }

   // Reserva memoria para la matriz "xaux"
   if((xaux=(double***)calloc(NREAC,sizeof(double**)))==NULL)
   {
     printf("\nNo hay memoria para la matriz 'xaux'");
     exit(1);
   }
   for(j=0;j<NREAC;j++)
   {
     if((xaux[j]=(double**)calloc(nenz[j],sizeof(double*))) == NULL)
     {
       printf("\nNo hay memoria para la matriz 'xaux'");
       exit(1);
     }
   }
   for(j=0;j<NREAC;j++)
   {
      for(k=0;k<nenz[j];k++)
      {
         if((xaux[j][k]=(double*)calloc(tipenz[j],sizeof(double))) == NULL)
         {
            printf("\nNo hay memoria para la matriz 'ocup'");
            exit(1);
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

   // Modificacion de las velocidades de los oligomeros
   for(j=0;j<NREAC;j++)
   {
      if(tipenz[j]==2)			// Dimero
      { 
         vmax[j][0]/=2.0; 
         vmax[j][1]/=2.0; 
      }  

      if(tipenz[j]==3)			// Trimero 
      { 
         vmax[j][0]/=3.0;
         vmax[j][1]/=3.0;
      }

      if(tipenz[j]==4) 			// Tetramero
      {
         vmax[j][0]/=4.0;		
         vmax[j][1]/=4.0;		
      }
   }

   // Modificacion de las afinidades de los oligomeros
   for(j=0;j<2*NREAC;j++)
   {
      if(fmod((double)j,2.0)==0.0) 		// Numero de la reaccion
         cdir=j/2; 

      // Modificacion de km 
      for(k=0;k<NSUB;k++)
      {
         if(tipenz[cdir]==2)			// Dimero
            km[j][k]*=sqrt(ak); 

         if(tipenz[cdir]==3)			// Trimero 
            km[j][k]*=1.05*ak;

         if(tipenz[cdir]==4) 
            km[j][k]*=pow(ak,1.53);		// Tetramero
      }

      // Modificacion de ki
      if(inhib[j][0]>=0)			// La reaccion es inhibida
      {
         if(tipenz[cdir]==2)			// Dimero
            ki[inhib[j][0]]*=sqrt(ak); 

         if(tipenz[cdir]==3)			// Trimero 
            ki[inhib[j][0]]*=1.05*ak;

         if(tipenz[cdir]==4) 			// Tetramero
            ki[inhib[j][0]]*=pow(ak,1.53);
      }
   }  

   // Calcula los parametros "b" y "phic" para cada reaccion, ambas direcciones
   //scal=2.0e-6;
   scal=25.0e-6;
   for(j=0;j<NREAC;j++)
   {
      bdir[j]=a/(vmax[j][0]*scal)-3.0*a;
      phicdir[j]=bdir[j]-0.001;
       
      if(vmax[j][1]<0.0)			// Reaccion inversa que no se da!
         binv[j]=0.0;
      else if(vmax[j][1]>0.0)   		// Si la reaccion inversa se da...	
         binv[j]=a/(vmax[j][1]*scal)-3.0*a;
      phicinv[j]=binv[j]-0.001;
   }

   //******************************************* INICIALIZACIONES *************************************	

   // Fija condiciones iniciales de las enzimas, region caotica
   for(i=0;i<NREAC;i++)
   {
      for(j=0;j<nenz[i];j++)
      {
         for(k=0;k<tipenz[i];k++)
         {
            xaux[i][j][k]=drand48();
            ocup[i][j][k]=0;
         }
      }
   }

   // Inicializa registro de control de actualizacion de enzimas
   for(j=0;j<tenz;j++)
      cont[j]=j;

   // Inicializa matriz 'act'. Si se trata de un monomero 'act' siempre vale cero!
   for(j=0;j<NREAC;j++)
   {
      for(k=0;k<nenz[j];k++)
         act[j][k]=0;
   }

   // Registra concentraciones de sustancias
   fprintf(dptr,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",time,n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11],n[12],n[13],n[14],n[15],n[16]); 
   time+=0.6e-3;
  
   //**************************** COMIENZO DE SIMULACION ************************************

   // Ciclo de las enzimas
   for(z=0;z<TOP;z++)
   {
      // Baraja vector de recorrido de las enzimas
      baraja(z); 
     
      // Modificacion de las afinidades por inhibicion
      for(j=0;j<2*NREAC;j++)
      {
         if(inhib[j][0]>=0)			// Si hay inhibicion...
            km[j][inhib[j][1]]*=1.0+n[inhib[j][0]]/ki[inhib[j][0]];
      }

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
            xaux[reac][enz][0]=monomer(reac,enz,0,xaux[reac][enz][0]);    	
          
         // Dimero 
         if(tipenz[reac]==2)
         {
            // Encuentra numero de unidades ocupadas en la enzima
            act[reac][enz]=0;
            for(j=0;j<2;j++)
            {
               if(ocup[reac][enz][j]>0)  
                  act[reac][enz]++;
               if(ocup[reac][enz][j]<0)
                  act[reac][enz]--;  
            }

            // Modifica probabilidades de acuerdo a numero de unidades ocupadas
            if(abs(act[reac][enz])==1)			// Una unidad ocupada
            {
               p1[0][reac]*=ak;  
               p1[1][reac]*=ak;  
            }

            // Actualiza fases para unidades en la enzima
            for(j=0;j<2;j++)
               xaux[reac][enz][j]=monomer(reac,enz,j,xaux[reac][enz][j]);    	
         }  

         // Trimero
         if(tipenz[reac]==3)
         {
            // Encuentra numero de unidades ocupadas en la enzima
            act[reac][enz]=0;
            for(j=0;j<3;j++)
            {
               if(ocup[reac][enz][j]>0)  
                  act[reac][enz]++;
               if(ocup[reac][enz][j]<0)
                  act[reac][enz]--;  
            }

            // Modifica probabilidades de acuerdo a numero de unidades ocupadas
            if(abs(act[reac][enz])==1)			// Una unidad ocupada
            {
               p1[0][reac]*=ak;  
               p1[1][reac]*=ak;  
            }
            else if(abs(act[reac][enz])==2)		// Dos unidades ocupadas
            {
               p1[0][reac]*=ak*bk;  
               p1[1][reac]*=ak*bk;  
            }

            // Actualiza fases para unidades en la enzima
            for(j=0;j<3;j++)
               xaux[reac][enz][j]=monomer(reac,enz,j,xaux[reac][enz][j]);    	
         }  

         // Tetramero 
         if(tipenz[reac]==4)
         {
            // Encuentra numero de unidades ocupadas en la enzima
            act[reac][enz]=0;
            for(j=0;j<4;j++)
            {
               if(ocup[reac][enz][j]>0)  
                  act[reac][enz]++;
               if(ocup[reac][enz][j]<0)
                  act[reac][enz]--;  
            }

            // Modifica probabilidades de acuerdo a numero de unidades ocupadas
            if(abs(act[reac][enz])==1)			// Una unidad ocupada
            {
               p1[0][reac]*=ak;  
               p1[1][reac]*=ak;  
            }
            else if(abs(act[reac][enz])==2)		// Dos unidades ocupadas
            {
               p1[0][reac]*=ak*bk;  
               p1[1][reac]*=ak*bk;  
            }
            else if(abs(act[reac][enz])==3)		// Tres unidades ocupadas
            {
               p1[0][reac]*=ak*bk*ck;  
               p1[1][reac]*=ak*bk*ck;  
            }
          
            // Actualiza fases para unidades en la enzima
            for(j=0;j<4;j++)
               xaux[reac][enz][j]=monomer(reac,enz,j,xaux[reac][enz][j]);    	
         }  
      }

      // Se retira piruvato
      /*if(fmod((double)z,200)==0) 
      {
         if(n[16]>0)
            n[16]--;
      }*/

      // Registra concentraciones de sustancias
      fprintf(dptr,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",time,n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11],n[12],n[13],n[14],n[15],n[16]); 
      time+=0.6e-3;
   }

   fclose(dptr);
}

//******************************************************************************************************
//***************************************** SUBRUTINAS *************************************************
//******************************************************************************************************

// Barajado del orden en que las enzimas son recorridas
void baraja(z)
int z;
{
   extern int *cont;
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
   else if(fmod((double)z,20.0)==0)		// Baraja cada 20 iteraciones
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
   extern double p1[2][NREAC],km[2*NREAC][NSUB],vmax[NREAC][2];
   extern int eps[NREAC][NSUB],n[NSUB];
   int j,k,cdir,cinv;
 
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
            p1[0][j]*=pow((double)n[k]*(vmax[j][0]/km[2*j][k]),(double)abs(eps[j][k]));   
         }

         if(eps[j][k]>0) 			// Reaccion inversa
         {
            if(vmax[j][1]<0.0)			// La reaccion no se da
               p1[1][j]*=0.0;   
            //if((km[2*j+1][k]==0.0)&(vmax[j][1]>0.0))	// Substancia no participa en la reaccion
            if(km[2*j+1][k]==0.0)		// Substancia no participa en la reaccion
               p1[1][j]*=1.0;   
            if(km[2*j+1][k]!=0.0)
               p1[1][j]*=pow((double)n[k]*(vmax[j][1]/km[2*j+1][k]),(double)abs(eps[j][k]));   
         }   
      }
   }
}

//***********************************************************************************************************

// Avance de fase para monomeros (o subunidades de oligomeros)
double monomer(pos1,pos2,pos3,state_in)
double state_in;
int pos1,pos2,pos3;
{
   extern double bdir[NREAC],binv[NREAC],p1[2][NREAC];  
   extern double a,p2;
   extern int eps[NREAC][NSUB],n[NSUB];
   extern int ***ocup;
   extern int **act;
   double state_out=0.0,umbral=0,dado;
   int flag1,flag2,m,j,gof=0,gob=0;

   // La enzima esta en reposo
   if(ocup[pos1][pos2][pos3]==0)
   {
      gof=1;			       		// Verifica los sustratos esten disponibles (gof=1)
      for(m=0;m<NSUB;m++)
      {
         if(eps[pos1][m]<0)			// Si substancia es sustrato de la reaccion directa...
         {
            if(n[m]>=abs(eps[pos1][m]))  	// Si TODOS los sustratos estan disponibles...
               gof*=1;				// ...la reaccion procede
            else
               gof*=0;                        
         }  
      }

      gob=1;					// Verifica que los sustratos esten disponibles (gob=1)
      for(m=0;m<NSUB;m++)
      {
         if(eps[pos1][m]>0)			// Si sustancia es sustrato de la reaccion inversa...
         {
            if(n[m]>eps[pos1][m]) 		// Si TODOS los sustratos estan disponibles...
               gob*=1;				// ...la reaccion procede
            else
               gob*=0;                        
         }  
      }

      // Determina si la enzima cataliza la reaccion directa o la inversa
      dado=drand48();
      umbral=p1[0][pos1]/(p1[0][pos1]+p1[1][pos1]);		// Umbral reaccion directa

      // Dos posibilidades:
      // a) La enzima es un oligomero y ninguna unidad de la enzima ha sido ocupada 
      // b) La enzima es un monomero 
      if(act[pos1][pos2]==0)
      {
         // Hay sustrato para la reaccion directa y para la reaccion inversa
         if((gof==1)&(gob==1))
         {
            // Itera con mapa reaccion directa	
            if(dado<=umbral)		
               state_out=avandir(state_in,pos1,pos2,pos3,0);

            // Itera con mapa reaccion inversa
            if(dado>umbral)
               state_out=avaninv(state_in,pos1,pos2,pos3,0);
         }

         // Solo hay sustrato para reaccion directa
         if((gof==1)&(gob==0))
            state_out=avandir(state_in,pos1,pos2,pos3,0);

         // Solo hay sustrato para reaccion inversa
         if((gof==0)&(gob==1))
            state_out=avaninv(state_in,pos1,pos2,pos3,0);

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
            if((dado<=umbral)&(act[pos1][pos2]>0))		
               state_out=avandir(state_in,pos1,pos2,pos3,0);

            // Itera con mapa reaccion inversa si las unidades ocupadas van en esa direccion
            if((dado>umbral)&(act[pos1][pos2]<0))
               state_out=avaninv(state_in,pos1,pos2,pos3,0);
         }

         // Solo hay sustrato para reaccion directa: la reaccion procede si las unidades
         // ocupadas van en esa misma direccion
         if((gof==1)&(gob==0)&(act[pos1][pos2]>0))
            state_out=avandir(state_in,pos1,pos2,pos3,0);

         // Solo hay sustrato para reaccion inversa: la reaccion procede si las unidades ocupadas
         // van en esa misma direccion
         if((gof==0)&(gob==1)&(act[pos1][pos2]<0))
            state_out=avaninv(state_in,pos1,pos2,pos3,0);

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
   if((ocup[pos1][pos2][pos3]>0))			
      state_out=avandir(state_in,pos1,pos2,pos3,1);

   // Avance de fase de la enzima, reaccion inversa
   if((ocup[pos1][pos2][pos3]<0))
      state_out=avaninv(state_in,pos1,pos2,pos3,1);

   // Revisa estados para tomar sustrato o liberar producto 
   sustprod(pos1,pos2,pos3,state_in,state_out);				

   // La enzima completa su ciclo, queda disponible
   if(ocup[pos1][pos2][pos3]==3)	
   {
      if(fabs(state_in-state_out)>=1.0+bdir[pos1])	
         ocup[pos1][pos2][pos3]=0;					// Enzima queda disponible
   }
   if(ocup[pos1][pos2][pos3]==-3)	
   {
      if(fabs(state_in-state_out)>=1.0+binv[pos1])	
         ocup[pos1][pos2][pos3]=0;					// Enzima queda disponible
   }
   return state_out;
}

//***********************************************************************************************************

// Actualizacion de concentraciones 
void sustprod(pos1,pos2,pos3,state_in,state_out)
double state_in,state_out;
int pos1,pos2,pos3;
{
   extern double phicdir[NREAC],phicinv[NREAC],km[2*NREAC][NSUB];
   extern int eps[NREAC][NSUB],n[NSUB];
   extern int ***ocup; 
   int gof,gob,m,j,count;

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

   gob=1;					// Verifica que los sustratos esten disponibles (gob=1)
   for(m=0;m<NSUB;m++)
   {
      if(eps[pos1][m]>0)			// Si sustancia es sustrato de la reaccion inversa...
      {
         if(km[2*pos1+1][m]!=0.0)		// La sustancia participa en la reaccion inversa
         {  
            if(n[m]>eps[pos1][m]) 		// Si TODOS los sustratos estan disponibles...
               gob*=1;				// ...la reaccion procede
            else
               gob*=0;                        
         }
      }  
   }

   // La enzima toma sustratos para la reaccion directa
   if(ocup[pos1][pos2][pos3]==1)			
   {
      if((state_in<=1.0)&(state_out>1.0)&(gof==1))	
      {
         for(m=0;m<NSUB;m++)
         { 
            if(eps[pos1][m]<0)			// La sustancia es sustrato de la reaccion directa
               n[m]+=eps[pos1][m];		// Actualiza concentraciones
         }
         ocup[pos1][pos2][pos3]=2;		// La enzima queda ocupada, region laminar directa = 2
         probab();				// Actualiza probabilidades
      }

      // La enzima sigue en la region caotica
      if((state_in<=1.0)&(state_out<=1.0))	
      {
         ocup[pos1][pos2][pos3]=0;  
         return;
      }
   }        

   // La enzima toma sustratos para la reaccion inversa
   if(ocup[pos1][pos2][pos3]==-1)			
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
         ocup[pos1][pos2][pos3]=-2;		// La enzima queda ocupada, region laminar inversa = -2
         probab();				// Actualiza probabilidades
      }

      // La enzima sigue en la region caotica
      if((state_in<=1.0)&(state_out<=1.0))	
      {
         ocup[pos1][pos2][pos3]=0;  
         return;
      } 
   }
 
   // Liberacion de producto en la reaccion directa
   if(ocup[pos1][pos2][pos3]==2)
   {
      // Reaccion directa libera producto, la enzima queda en recuperacion 
      if(state_out>=phicdir[pos1])	
      {
         for(j=0;j<NSUB;j++)
         {
            if(eps[pos1][j]>0)				// Libera producto
               n[j]+=eps[pos1][j];			// Actualiza concentracion
         }
      	 ocup[pos1][pos2][pos3]=3;			// Enzima en recuperacion, reaccion directa
         probab();					// Actualiza probabilidades
      }
   }

   // Liberacion de producto en la reaccion inversa
   if(ocup[pos1][pos2][pos3]==-2)
   {
      // Reaccion inversa libera producto, la enzima queda en recuperacion 
      if(state_out>=phicinv[pos1])	
      {
         for(j=0;j<NSUB;j++)
         {
            if(eps[pos1][j]<0)		 		// Libera producto
               n[j]-=eps[pos1][j];			// Actualiza concentracion
         }
         ocup[pos1][pos2][pos3]=-3;			// Enzima en recuperacion, reaccion inversa
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
   extern double bdir[NREAC],binv[NREAC],p1[2][NREAC];  
   extern double a,p2;
   extern int ***ocup;
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
      state_out=enzima47(state_in,a,bdir[pos1],p1[0][pos1],p2);               
      if(state_out!=state_in)
         flag1=1;
      else				// Punto fijo
      { 
         state_in=drand48();
         flag2=0;  
      }
   }
   if(arr==0)
      ocup[pos1][pos2][pos3]=1;		// Iteracion mapa reaccion directa 
    
   return state_out;
}

//***********************************************************************************************************

// Avance de fase reaccion inversa
double avaninv(state_in,pos1,pos2,pos3,arr)
double state_in;
int pos1,pos2,pos3,arr;
{
   extern double bdir[NREAC],binv[NREAC],p1[2][NREAC];  
   extern double a,p2;
   extern int ***ocup;
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
      ocup[pos1][pos2][pos3]=-1;		// Iteracion mapa reaccion inversa
      
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
        printf("\nLa semilla es %d",seed);
	srand48(seed);
}

        
