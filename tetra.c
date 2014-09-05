/* tetra.c 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	

   Modelo simple de una enzima alosterica (tetramero). Se utiliza un mapa caotico 
   modificado para representar cada enzima, y se emplea el modelo de cooperatividad
   simple de Adair-Pauling.

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


#define DT 0.1e-6	// Incremento de tiempo correspondiente a una iteracion (segundos)
#define VMAX 8000.	// Velocidad maxima, micromol/(min-mg)
#define S_05 25		// [S]0.5 en micromolar
#define VSIM 10.3e-18	// Volumen de simulacion
#define NA 6.022e23	// Numero de Avogadro


#define N 10
#define TOP 5000000
#define VECES 1
#define SMAX 2000

// variables globales
int sort[N]={0};
double phia[2][N]={0.0},phib[2][N]={0.0},phic[2][N]={0.0},phid[2][N]={0.0};
int libera[N]={0},liberb[N]={0},liberc[N]={0},liberd[N]={0},act[N]={0};
int acta[N]={0},actb[N]={0},actc[N]={0},actd[N]={0};
int total1;
double phil;

double enzima47();
void avanfase();
void rand_init();
void baraja();
int drnd();

main()
{
	// Declaraciones 
	FILE *dptr;    
	int i,j,k,m=0,n,z,flag,pos,conc;
        int sg[SMAX]={0};
        double v[VECES]={0.0},vprom[SMAX-1]={0.0};
	double n0[SMAX]={0},n1[SMAX]={0},n2[SMAX]={0},n3[SMAX]={0},n4[SMAX]={0},p[N]={0.0},hill[SMAX]={0.0},tetram[SMAX]={0.0}; 
	double a,b,p2,sum,count,pa,pb,pc,pd,ak,bk,ck,n0c,n1c,n2c,n3c,n4c,sum0,sum1,sum2,sum3,sum4,stot,smax,km0;
        double vmax,km_tetra,km,pimp,s,pmolar,trec;  

        //srand48(time(NULL));
        rand_init();

	if((dptr=fopen("tetra.dat","w"))==NULL)
	   printf("\nEl archivo de salida no puede ser abierto\n");

        // ********************** Parametros de las enzimas *******************************
        // Cooperatividad de oligomeros 
        ak=10.0;
        bk=ak;
        ck=ak; 

        a=0.1;
        p2=0.999999;

        // Parametros de los monomeros
	
        // PESOS MOLARES: en Daltons
        pmolar = 40000.;

	// Fraccion del tiempo de procesamiento
        trec = 0.01;		

	// Calcula velocidad para una sola enzima
        vmax=1.02*VMAX*1.e-3*pmolar/60.0;		
	vmax /= 4.0;
   
        // Calcula los parametros "b" y "phic"
        b=a/(vmax*DT)-3.0*a;
        printf("\nEl numero de iteraciones para la r. dir. es %f, el No. de ciclos/iteracion es %lf",(b/a)+3,(double)TOP/((b/a)+3)); 
        phil=1.0+b;

        // Parametros del tetramero
        //km=S_05*pow(ak,1.53);			// ak = 10
        km=S_05*pow(ak,1.50);			// ak = 10
	km0=km;
        km=km*NA*VSIM/1.e6;	
        //km=km_tetra*pow(ak,1.63);		// ak = 15
	printf("\nEl valor de Km es %lf",km);
	
        smax=km/vmax;

        // Inicializa vector de recorrido de las enzimas
        for(j=0;j<N;j++)
           sort[j]=j;

	// Ciclo de las enzimas 
	conc=0;    			// Cantidad inicial de producto 
	for(conc=1;conc<SMAX;conc++)
	{
           sg[conc-1]=conc;  			// Guarda valores de cantidad de sustrato
           s=(double)conc*NA*VSIM/1.e6;	

           // Fija afinidades
           pa=DT*s/smax;
           pb=ak*pa;
           pc=ak*bk*pa;
           pd=ak*bk*ck*pa;


           // Inicializa vector de velocidades
           for(k=0;k<VECES;k++)
              v[k]=0.0; 

           for(n=0;n<VECES;n++)
           {

              // Inicializa vectores de activacion 
              for(k=0;k<N;k++)
              {
                 acta[k]=0;
                 actb[k]=0;
                 actc[k]=0;
                 actd[k]=0;
                 libera[k]=0; 
                 liberb[k]=0; 
                 liberc[k]=0; 
                 liberd[k]=0; 
              }

	      // La enzimas arrancan en la region caotica 	
	      for(i=0;i<N;i++)
              {
	         phia[0][i]=drand48();
	         phia[1][i]=0;
	         phib[0][i]=drand48();
	         phib[1][i]=0;
	         phic[0][i]=drand48();
	         phic[1][i]=0;
	         phid[0][i]=drand48();
	         phid[1][i]=0;
              }

              // Inicializa control de activacion
              for(i=0;i<N;i++)
                 act[i]=0;

	      n0c=0;
              n1c=0;
              n2c=0;
              n3c=0;
              n4c=0; 
              total1=0; 
              for(z=0;z<TOP;z++)
              {
                 // Aleatoriza recorrido de las enzimas  
                 baraja(z);

	         // Actualizacion de las enzimas
	         for(k=0;k<N;k++)  
	         {	
                    pos=sort[k];
                    //pos=k;

                    // Control de activacion
                    act[pos]=acta[pos]+actb[pos]+actc[pos]+actd[pos];

                    if(act[pos]==0)		// Ninguna unidad ocupada
		    {
		       n0c++;	    
                       p[pos]=pa;  
		    }   
                    else if(act[pos]==1)	// Una unidad ocupada
                    {
                       n1c++;
                       p[pos]=pb;  
                    }
                    else if(act[pos]==2)	// Dos unidades ocupadas
                    {
                       n2c++;
                       p[pos]=pc;
                    } 
                    else if(act[pos]==3)	// Tres unidades ocupadas
                    {
                       n3c++; 
                       p[pos]=pd;   
                    } 
                    else if(act[pos]==4)	// Cuatro unidades ocupadas
                       n4c++; 

                    ///////////////////// ENZIMA A /////////////////////		
                    avanfase(phia[0][pos],a,b,p[pos],p2,pos,0);               
                   
                    ///////////////////// ENZIMA B /////////////////////		
                    avanfase(phib[0][pos],a,b,p[pos],p2,pos,1);               

                    ///////////////////// ENZIMA C /////////////////////		
                    avanfase(phic[0][pos],a,b,p[pos],p2,pos,2);               
                 
                    ///////////////////// ENZIMA D /////////////////////		
                    avanfase(phid[0][pos],a,b,p[pos],p2,pos,3);               
                 }
              }  
              v[n]=(double)total1*1.e6*60./(NA*VSIM*DT*TOP);

	      sum0+=n0c;
              sum1+=n1c; 
              sum2+=n2c; 
              sum3+=n3c; 
              sum4+=n4c; 
           }  
	   printf("\n[S] = %d, pa = %lf, pb = %lf, pc = %lf, pd = %lf",conc,pa,pb,pc,pd);
           sum=0;
           for(k=0;k<VECES;k++)
              sum+=v[k];

           vprom[conc-1]=sum/(double)VECES;
	   stot=(sum0+sum1+sum2+sum3+sum4);
	   
           n0[conc-1]=sum0/stot;
           n1[conc-1]=sum1/stot;
           n2[conc-1]=sum2/stot;
           n3[conc-1]=sum3/stot;
           n4[conc-1]=sum4/stot;

           // Ecuacion de Hill
           //pimp=km_tetra*1.25e5;
           //pimp=4*vmax*1e-6*pow(sg[conc-1],4.0)/(km+pow(sg[conc-1],4.0));
           //hill[conc-1]=vmax*1e-6*pow(sg[conc-1],4.0)/(km_tetra*1.25e5+pow(sg[conc-1],4.0));
           //hill[conc-1]=VMAX*(N*pmolar*1.e3/(NA*VSIM))*pow(sg[conc-1],4.0)/(S_05*1.25e5+pow(sg[conc-1],4.0));
           hill[conc-1]=VMAX*(N*pmolar*1.e3/(NA*VSIM))*pow(sg[conc-1],4.0)/(pow(S_05,4.0)+pow(sg[conc-1],4.0));
           tetram[conc-1]=VMAX*(N*pmolar*1.e3/(NA*VSIM))*(sg[conc-1]/km0 + 3.0*ak*pow(sg[conc-1]/km0,2.0) + 3.0*pow(ak,2.0)*bk*pow(sg[conc-1]/km0,3.0) 
			   + pow(ak,3.0)*pow(bk,2.0)*ck*pow(sg[conc-1]/km0,4.0))/(1.0 + 4.0*sg[conc-1]/km0 + 6.0*ak*pow(sg[conc-1]/km0,2.0)
			   + 4.0*pow(ak,2.0)*bk*pow(sg[conc-1]/km0,3.0) + pow(ak,3.0)*pow(bk,2.0)*ck*pow(sg[conc-1]/km0,4.0));

	}
           
        for(k=0;k<SMAX-1;k++)
	   fprintf(dptr,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",(double)sg[k]/km0,vprom[k],tetram[k],hill[k],n0[k],n1[k],n2[k],n3[k],n4[k]);

	fclose(dptr);
	return(0);
}

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
      else if((x>1)&(x<=(1+b))) 
         out = a+x;
      else if((x>(b+1))&(x<=(b+1.5-p2/2))) 
         out = m3*(x-1.5-b+p2/2+(b+2)/m3);
      else if((x>(b+1.5-p2/2))&(x<=(b+1.5))) 
         out = m4*(x-1.5-b+(b+2+a)/m4);
      else if((x>(b+1.5))&(x<=(1.5+b+p2/2))) 
         out = b+2+a+m4*(1.5+b)-m4*x;
      else if((x>(b+1.5+p2/2))&(x<=(2+b))) 
         out = b+2+m3*(1.5+b+p2/2)-m3*x;
      else if(x>=(2+b)) 
         out = -1/a*(b+2)+1/a*x;
   }
   return out;
}

void baraja(z)
int z;
{
   int i,j,k,posic;
   extern int sort[N];

   // Baraja vector de recorrido de las enzimas
   if(z==0)
   {
      for(k=0;k<10;k++)
      {   
         for(j=0;j<N;j++)
         {
            posic=drnd(N); 
            i=sort[posic];
            sort[posic]=sort[j];
            sort[j]=i;
         }
      }        
   }
   else if(z!=0)
   {
      for(j=0;j<N;j++)
      {
         posic=drnd(N); 
         i=sort[posic];
         sort[posic]=sort[j];
         sort[j]=i;
      }
   }
}

// Avance de fase reaccion directa
void avanfase(state_in,a,b,p1,p2,pos,num)
double state_in,a,b,p1,p2;
int num,pos;
{
   extern double phia[2][N],phib[2][N],phic[2][N],phid[2][N];
   extern int acta[N],actb[N],actc[N],actd[N];
   extern int libera[N],liberb[N],liberc[N],liberd[N],act[N];
   extern int total1;
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
      state_out=enzima47(state_in,a,b,p1,p2);               
      if(state_out!=state_in)
         flag1=1;
      else				// Punto fijo
      { 
         state_in=drand48();
         flag2=0;  
      }
   }

   // Primera subunidad 
   if(num==0)
   {
      if((state_out>=1.0)&(state_in<1.0))	// Enzima "a" toma sustrato
      {
         libera[pos]=1;
         acta[pos]=1;				// Enzima "a" esta ocupada
      }

      if((state_out>=phil)&(libera[pos]==1))	// Enzima "a" libera producto
      {
         libera[pos]=-1; 
         total1++;  
      }

      if(fabs(state_out-state_in)>b+1)		// Enzima "a" completa ciclo
      { 
         libera[pos]=0;				// Enzima "a" queda disponible
         acta[pos]=0;				// Cesa activacion
         state_out=drand48();			// Re-inicializa en region laminar
      } 
      phia[0][pos]=state_out;
   }

   // Segunda subunidad 
   if(num==1)
   {
      if((state_out>=1.0)&(state_in<1.0))	// Enzima "b" toma sustrato
      {
         liberb[pos]=1;
         actb[pos]=1;				// Enzima "b" esta ocupada
      }

      if((state_out>=phil)&(liberb[pos]==1))	// Enzima "b" libera producto
      {
         liberb[pos]=-1; 
         total1++;  
      }

      if(fabs(state_out-state_in)>b+1)		// Enzima "b" completa ciclo
      { 
         liberb[pos]=0;				// Enzima "b" queda disponible
         actb[pos]=0;				// Cesa activacion
         state_out=drand48();			// Re-inicializa en region laminar
      } 
      phib[0][pos]=state_out;
   }

   // Tercera subunidad 
   if(num==2)
   {
      if((state_out>=1.0)&(state_in<1.0))	// Enzima "c" toma sustrato
      {
         liberc[pos]=1;
         actc[pos]=1;				// Enzima "c" esta ocupada
      }

      if((state_out>=phil)&(liberc[pos]==1))	// Enzima "c" libera producto
      {
         liberc[pos]=-1; 
         total1++;  
      }

      if(fabs(state_out-state_in)>b+1)		// Enzima "c" completa ciclo
      { 
         liberc[pos]=0;				// Enzima "c" queda disponible
         actc[pos]=0;				// Cesa activacion
         state_out=drand48();			// Re-inicializa en region laminar
      } 
      phic[0][pos]=state_out;
   }

   // Cuarta subunidad 
   if(num==3)
   {
      if((state_out>=1.0)&(state_in)<1.0)	// Enzima "d" toma sustrato
      {
         liberd[pos]=1;
         actd[pos]=1;				// Enzima "d" esta ocupada
      }

      if((state_out>=phil)&(liberd[pos]==1))	// Enzima "d" libera producto
      {
         liberd[pos]=-1; 
         total1++;  
      }

      if(fabs(state_out-state_in)>b+1)		// Enzima "d" completa ciclo
      { 
         liberd[pos]=0;				// Enzima "d" queda disponible
         actd[pos]=0;				// Cesa activacion
         state_out=drand48();			// Re-inicializa en region laminar
      } 
      phid[0][pos]=state_out;
   }
}


// Inicializa generador de numeros aleatorios 
void rand_init()
{
	int seed;
	time_t timenow;
	time(&timenow);
	seed=timenow % 100000;
	srand48(seed);
}

// Genera un numero aleatorio entero entre 0 y PHIMAX 
int drnd(phimax)
int phimax;
{
	int j;
        float dado;
        dado=drand48()*phimax;  
        for(j=0;j<phimax;j++)
        {
           if((j<=dado)&(dado<=j+1))
              return j;
        } 
}


