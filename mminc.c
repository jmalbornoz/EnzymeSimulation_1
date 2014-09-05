/* mminc_a.c 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	 

   Modelo simple de una enzima Michaelis-Menten. Se utiliza un mapa caotico 
   modificado. La enzima toma un solo sustrato
   La velocidad maxima y la afinidad se ajustaron. Con un factor scal = 10e-6 cada "sustrato" representa
   un miliMol; de este modo la velocidad maxima queda multiplicada por 10. En consecuencia, las velocidades
   registradas en la simulacion son divididas por 10. Una "unidad" de sustrato representa un micromolar; de 
   esta manera una iteracion corresponde a 6e-4 segundos
   * Se incluye inhibicion no competitiva a traves de un segundo mapa

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define DT 10.0e-6	// Incremento de tiempo correspondiente a una iteracion (segundos)
#define VMAX 100.	// Velocidad maxima, micromol/(min-mg)
#define Km 50.		// Afinidad por el sustrato, micromolar
#define Ki 25.		// Afinidad por el inhibidor, micromolar
#define I 25.		// Concentracion de inhibidor
#define N 10		// Numero de enzimas
#define VECES 1		// Numero de veces para calcular promedio
#define SS 500		// Maxima cantidad de sustrato
#define TOPE 25000
#define VSIM 10.3e-18
#define NA 6.022e23

double enzima47();
void rand_init();
int drnd();

int main()
{
	// Declaraciones 
	FILE *dptr;    
	int i,j,k,m,n,z,P,flag1,flag2,pos,KM,KI,I0,Ic,conc;
	int estado[N][2]={0},sg[SS]={0},sort[N];
	double fase[N][2]={0.0},v[VECES]={0.0},vprom[SS-1]={0.0},mm[SS-1]={0.0};
	double a,b,p1d,p2,pinh=0.0,pinh1=0.0,pinh2=0.0,sum,phic,smax,vtop=0;
        double fase_ins,fase_ini,fase_outs,fase_outi,kmapp,pmolar,trec,vmax,s;

	if((dptr=fopen("mminc.dat","w"))==NULL)
	   printf("\nEl archivo de salida no puede ser abierto\n");

        rand_init();

        a=0.1;
        p2=0.999999;

        // PESOS MOLARES: en Daltons
        pmolar = 60000.;

	// Fraccion del tiempo de procesamiento
        trec = 0.01;		

	// Calcula velocidad para una sola enzima
        vmax=VMAX*1.e-3*pmolar/60.0;		
   
        // Calcula los parametros "b" y "phic"
        b=a/(vmax*DT)-3.0*a;
        printf("\nEl numero de iteraciones para la r. dir. es %f, el No. de ciclos/iteracion es %lf",(b/a)+3,(double)TOPE/((b/a)+3)); 
        phic=1.0+b;

        // Afinidad por el sustrato
        KM=(int)Km*NA*VSIM/1.e6;	

	// Afinidad por el inhibidor
        KI=(int)Ki*NA*VSIM/1.e6;

	// Concentracion de inhibidor
        I0=I*NA*VSIM/1.e6;	

        // Parametros de la enzima
        smax=(double)KM/vmax;

        printf("\nsmax = %f, b = %f",smax,b);

        // Inicializa vector de recorrido de las enzimas
        for(j=0;j<N;j++)
           sort[j]=j;

	// Ciclo de las enzimas 
	for(conc=1;conc<SS;conc++)
	{
           sg[conc-1]=conc;  			// Guarda valores de concentracion de sustrato
           s=(double)conc*NA*VSIM/1.e6;	

           p1d=DT*s/smax;		// Define probabilidad reaccion directa

           for(k=0;k<VECES;k++)
              v[k]=0.0; 

           for(n=0;n<VECES;n++)
           {
	      // Inicializacion: 
              // Las enzimas arrancan en la region caotica 	
              rand_init(); 
	      for(i=0;i<N;i++)
              {
	         fase[i][0]=drand48();
	         fase[i][1]=drand48();
              }

              // Las enzimas estan en reposo
	      for(i=0;i<N;i++)
              {
	         estado[i][0]=0;
	         estado[i][1]=1;
              } 

              P=0;  
              Ic=I0;
              for(z=0;z<TOPE;z++)
              {
                 // Baraja vector de recorrido de las enzimas
                 if(z==0)
                 {
                    for(k=0;k<10;k++)
                    {   
                       for(j=0;j<N;j++)
                       {
                          pos=drnd(N); 
                          i=sort[pos];
                          sort[pos]=sort[j];
                          sort[j]=i;
                       }
                    }        
                 }
                 else if(z%20==0)
                 {
                    for(j=0;j<N;j++)
                    {
                       pos=drnd(N); 
                       i=sort[pos];
                       sort[pos]=sort[j];
                       sort[j]=i;
                   }
                 }

	         // Actualizacion de las fases de las enzimas
	         for(k=0;k<N;k++)  
	         {	
                    pos=sort[k];
 
                    // Actualiza fase para sustrato
                    if((estado[pos][0]==0)||(estado[pos][0]==1)||(estado[pos][0]==-1))	// No hay inhibidor adherido
                    { 
                       flag1=0;
                       flag2=0;
                       fase_ins=fase[pos][0];
                       while(flag1==0)			// Toma en cuenta puntos fijos
                       {  
                          while(flag2==0)
                          { 
                             if(fase_ins==0.5)
                                fase_ins=drand48();
                             else if(fase_ins!=0.5)
                                flag2=1;
                          }  
                          fase_outs=enzima47(fase_ins,a,b,p1d,p2);
                          if(fase_outs!=fase_ins)
                             flag1=1;
                          else
                             fase_ins=drand48();
                       }    
                       fase[pos][0]=fase_outs; 
                    } 

                    // Actualiza fase para inhibidor
                    if((I0!=0)&(Ic>=0))			// Hay inhibidor...
                    {  
                       pinh1=(double)Ic*vmax/(10*KI);	// Probabilidad de entrada del inhibidor
                       pinh2=(double)Ic*vmax/230;  	// Probabilidad de salida del inhibidor

                       if(estado[pos][1]==1)		// El inhibidor entra
                          pinh=pinh1;
                       if(estado[pos][1]==-1)		// El inhibidor sale
                          pinh=pinh2;

                       flag1=0;
                       flag2=0;
                       fase_ini=fase[pos][1];
                       while(flag1==0)			// Toma en cuenta puntos fijos
                       {  
                          while(flag2==0)
                          { 
                             if(fase_ini==0.5)
                                fase_ini=drand48();
                             else if(fase_ini!=0.5)
                                flag2=1;
                          }  
  	                  fase_outi=enzima47(fase_ini,a,b,pinh,p2);
                          if(fase_outi!=fase_ini)
                             flag1=1;
                          else
                             fase_ini=drand48();
                       }
                       fase[pos][1]=fase_outi; 
                    }

                    // Actualiza concentraciones

                    // Formacion de complejo ES: el sustrato entra primero y pasa a region laminar
                    if((estado[pos][0]==0)&(estado[pos][1]==1)&(fase_outs>1)&(fase_ins<=1))	
                       estado[pos][0]=1;

                    // Formacion de complejo ESI inactivo a partir de complejo ES
                    if((estado[pos][0]==1)&(estado[pos][1]==1)&(fase_outi>1)&(fase_ini<=1)&(Ic>0)) 
                    {
                       estado[pos][0]=11;		// Detiene avance de la enzima
                       estado[pos][1]=-1;		// Inhibidor entra
                       fase[pos][1]=drand48();
                       fase_outi=0.0;
                       fase_ini=0.0;
                       Ic--;  
                    }

                    // Formacion de complejo ES a partir de complejo ESI
                    if((estado[pos][0]==11)&(estado[pos][1]==-1)&(fase_outi>1)&(fase_ini<=1))
                    {
                       estado[pos][0]=1;		// Reanuda avance de la enzima
                       estado[pos][1]=1;		// Inhibidor salio
                       fase[pos][1]=drand48();
                       fase_outi=0.0;
                       fase_ini=0.0;
                       Ic++; 
                    }

                    // Complejo ES libera producto
                    if((estado[pos][0]==1)&(estado[pos][1]==1)&(fase_outs>phic))  
                    {
                       estado[pos][0]=-1;		// Enzima queda en recuperacion
                       P++;
                    }

                    // Enzima completa ciclo y queda en reposo
                    if((estado[pos][0]==-1)&(fabs(fase_outs-fase_ins)>=phic)) 
                       estado[pos][0]=0;

                    // Formacion de complejo EI: el inhibidor entro primero
                    if((estado[pos][0]==0)&(estado[pos][1]==1)&(fase_outi>1)&(fase_ini<=1)&(Ic>0))	
                    {
                       estado[pos][0]=10; 		// Detiene avance de la enzima
                       estado[pos][1]=-1;		// Inhibidor sale
                       fase[pos][1]=drand48();
                       fase_outi=0.0;
                       fase_ini=0.0;
                       Ic--;
                    }

                    // Complejo EI libera inhibidor
                    if((estado[pos][0]==10)&(estado[pos][1]==-1)&(fase_outi>1)&(fase_ini<=1))
                    {
                       estado[pos][0]=0;
                       estado[pos][1]=1;
                       fase[pos][1]=drand48();  
                       fase_outi=0.0;
                       fase_ini=0.0;
                       Ic++;
                    }   

                    // Formacion de complejo EI a partir de complejo EIS
                    //if((estado[pos][0]==11)&(estado[pos][1]==-1)&(fase_outs>1)&(fase_ins<=1))	
                    //   estado[pos][0]=11; 

                 }
              }  
              v[n]=(double)P*1.e6*60./(NA*VSIM*DT*TOPE);
           }  
           sum=0.0;
           for(k=0;k<VECES;k++)
              sum+=v[k];

           vprom[conc-1]=sum/(double)VECES; 

           printf("\nEl valor de la concentracion es %d",conc);
           
	}

        for(k=0;k<SS-1;k++)
           mm[k]=VMAX*(N*pmolar*1.e3/(NA*VSIM))*sg[k]/(Km+sg[k]); 

        for(k=0;k<SS-1;k++)
	   fprintf(dptr,"%d \t %f \t %f \n",sg[k],vprom[k],mm[k]);

        for(k=0;k<SS-1;k++)
        {
           if(vprom[k]>vtop)
              vtop=vprom[k];  
        }

        printf("\nVmax es %lf, Vmax/2 es %lf",vtop,vtop/2);

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


