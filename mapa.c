/* mapa.c 

   Jose M Albornoz
   Grupo de Caos y Sistemas Complejos
   Universidad de Los Andes	 

   Muestra el mapa caotico modificado

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TOP 5000

double enzima46();
double enzima47();
double enzima47bis();
double enzima48();
double enzima49();
void rand_init();

main()
{
   // Declaraciones 
   FILE *dptr;    
   int i,j,k,m,n;
   double x[TOP]={0.0},y[TOP]={0.0},xlin[TOP]={0.0},mapa[TOP]={0.0};
   double a,b,p1,p2,paso,pinh;

   rand_init();  
 
   if((dptr=fopen("mapa.dat","w"))==NULL)
      printf("\nEl archivo de salida no puede ser abierto\n");

   a=0.1;
   b=1.0;
   p1 = 0.25;
   p2 = 0.999999;
   pinh=0.00625;

   if(a>=0.0) 
      //x[0]=drand48();
      x[0]=0.265;
   else if(a<0.0) 
      x[0]=b+1.0+drand48();
  
   y[0]=0;

   // Iteraciones
   for(n=1;n<34;n++)
   {
      if(fmod((double)n,2.0)!=0.0) 
      {
         x[n]=x[n-1];
         y[n]=enzima49(x[n-1],a,b,p1,p2);
         //y[n]=enzima47(x[n-1],a,b,p1,p2);
         //y[n]=enzima48(x[n-1],a,b,p1);
         //y[n]=enzima47bis(x[n-1],a,b,p1,p2,pinh);
      }
      else
      {
         x[n]=y[n-1];
         y[n]=y[n-1];
      }
   }

   // linea x=y
   paso=-1.0;
   for(i=0;i<TOP;i++)
   {
      paso+=(b+4.0)/TOP;
      xlin[i]=paso;
   } 
  
   // mapa 
   for(n=0;n<TOP;n++)
      mapa[n]=enzima49(xlin[n],a,b,p1,p2);
      //mapa[n]=enzima47(xlin[n],a,b,p1,p2);
      //mapa[n]=enzima47bis(xlin[n],a,b,p1,p2,pinh);
      //mapa[]=enzima48(xlin[n],a,b,p1);
   
   for(i=0;i<TOP;i++)
      fprintf(dptr,"\n%f\t%f\t%f\t%f",xlin[i],x[i],y[i],mapa[i]);

   fclose(dptr);

}

// Enzima48
double enzima48(x,a,b,p)
double x,a,b,p;
{
   double out,m1,m2;

   m1 = 2.0/(1-p);
   m2 = 2.0*fabs(a)/p;

   if(a>=0.0)
   {
      if(x<0.0) 
         out = a+b+2.0+2.0*x;
      else if((x>=0.0)&(x<0.5-0.5*p)) 
         out = m1*x;
      else if((x>=0.5-0.5*p)&(x<0.5)) 
         out = 1.0+a-0.5*m2+m2*x;
      else if((x>=0.5)&(x<0.5+0.5*p)) 
         out = 0.5*m2*(1.0+p)+1.0-m2*x;
      else if((x>=0.5+0.5*p)&(x<1.0))  
         out = m1*(1.0-x);   
      else if((x>=1.0)&(x<(1.0+b)))
        out = a+x;
      else if((x>=(b+1))&(x<(b+1+a)))
        out = 2+2*b+1+a-x;
      else if((x>=(b+1+a))&(x<(2+b-a)))
        out = 2+b;
      else if((x>=(2+b-a))&(x<(2+b)))
        out = x+a;
      else if(x>=(2+b))
        out = -1/a*(b+2)+1/a*x;
  }
  else if(a<0.0)
  {
     if(x<0.0) 
        out = a+b+2+2*x;
     else if((x>=0)&(x<0.5-0.5*p)) 
        out = m1*x;
     else if((x>=0.5-0.5*p)&(x<0.5)) 
        out = 1+a-0.5*m2+m2*x;
     else if((x>=0.5)&(x<0.5+0.5*p)) 
        out = 0.5*m2*(1+p)+1-m2*x;
     else if((x>=0.5+0.5*p)&(x<1))  
        out = m1*(1-x);   
     else if((x>=1)&(x<(1+b))) 
        out = a+x;
     else if((x>=(1+b))&(x<(b+1.5-p/2)))
       out = 2+b+m1*(1+b)-m1*x;
     else if((x>=(b+1.5-p/2))&(x<(b+1.5)))
       out = 1+b+a+m2*(1+b+0.5)-m2*x;
     else if((x>=(b+1.5))&(x<(b+1.5+p/2)))
       out = m2*(x-b-1.5+(b+1-abs(a))/m2);
     else if((x>=(b+1.5+p/2))&(x<(2+b)))
       out = m1*(x-b-1.5-p/2+(b+1)/m1);
     else if(x>=(2+b))
       out = -1/abs(a)*(b+2)+1/abs(a)*x;
  }
  return out;
}   

// Enzima47
double enzima47(x,a,b,p1,p2)
double x,a,b,p1,p2;
{
   double out,m1,m2,m3,m4; 

   if(p1>=1)
      p1=0.99999;
   
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
      else if((x>=1)&(x<(1+b))) 
         out = a+x;
      else if((x>=(b+1))&(x<(1+b+a))) 
         out = b+x;
      else if((x>=(1+b+a))&(x<(b+2))) 
         out = m4*((b+2+a)/m4);
      else if(x>=(2+b)) 
         out = -(b+2)/a+x/a;
   }
/*   else if(a<0.0) 
   {
   if(x<0.0) 
      out = b+2+1/abs(a)*x;
   else if((x>=0.0) & (x<0.5-p2/2)) 
      out = 1-abs(a)-m3*x;
   else if((x>=0.5-p2/2) & (x<0.5)) 
      out = m4*(0.5-p2/2)-m4*x;
   else if((x>=0.5) & (x<0.5+p2/2)) 
      out = m4*(x-0.5-p2/2);
   else if((x>=0.5+p2/2) & (x<1.0)) 
      out = m3*(x-0.5-p2/2);
   else if((x>=1.0) & (x<(1+b))) 
      out = x+a;
   else if((x>=(1+b)) & (x<(b+1.5-p1/2))) 
      out = b+2+m1*(b+1)-m1*x;
   else if((x>=(b+1.5-p1/2)) & (x<(b+1.5))) 
      out = b+1+m2*(b+1.5-p1/2)-m2*x;
   else if((x>=(b+1.5)) & (x<(b+1.5+p1/2))) 
      out = m2*(x-b-1.5+(b+1-abs(a))/m2);
   else if((x>=(b+1.5+p1/2)) & (x<(2+b)))
      out = m1*(x-b-1.5-p1/2+(b+1)/m1);
   else if(x>=(2+b)) 
      out = -1/abs(a)*(b+2)+1/abs(a)*x;
   }*/
   return out;
}

// Enzima47bis
// pinh = probabilidad de inhibicion
double enzima47bis(x,a,b,p1,p2,pinh)
double x,a,b,p1,p2,pinh;
{
   double out,m1,m2,m3,m4,mp,bp,p,ap,ki,I;
   ki=10;
   I=10;

   // Inhibicion competitiva
   if(pinh!=0.0)
   {
      p=p1+pinh;
      if(p<=1.0)
         p=0.99999;  
      m2 = 2*fabs(a)/p;
      bp=0.5*m2*(1+p)+1-m2*0.5*(1+p1);
      bp-=1.0;
      ap=bp*ki*p1/(2*I);
      p1+=pinh; 
   }
   else if(pinh==0.0)
      ap=a;

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
   else if((x>1)&(x<=1+bp)) 
      out = ap+x;
   else if((x>1+bp)&(x<=1+b+bp)) 
      out = a+x;
   else if((x>1+b+bp)&(x<=b+bp+1.5-p2/2)) 
      out = m3*(x-1.5-b+p2/2+(b+2)/m3);
   else if((x>b+bp+1.5-p2/2)&(x<=b+bp+1.5)) 
      out = m4*(x-1.5-b+(b+2+a)/m4);
   else if((x>b+bp+1.5)&(x<=1.5+b+bp+p2/2)) 
      out = b+2+a+m4*(1.5+b)-m4*x;
   else if((x>b+bp+1.5+p2/2)&(x<=2+b+bp)) 
      out = b+2+m3*(1.5+b+p2/2)-m3*x;
   else if(x>=2+b+bp) 
      out = -1/a*(b+2)+1/a*x;

   return out;
}

// Enzima46
double enzima46(x,a,b,w0,w1,s,p2)
double x,a,b,w0,w1,p2;
int s;
{
   double out,m1,m2,m3,m4; 
   
   m1 = 2/(1-w0-w1*s);
   m2 = 2*fabs(a)/(w1*s);
   m3 = 2*(1-fabs(a))/(1-p2);
   m4 = 2*fabs(a)/p2;

   if(w0+w1*s>=1)
      w1=(0.99999-w0)/(float)s; 
   
   if(a>=0.0)
   { 
      if(x<0.0) 
         out = a+b+2+2*x;
      else if((x>=0)&(x<0.5-0.5*(w0+w1*s))) 
         out = m1*x;
      else if((x>=0.5-0.5*(w0+w1*s))&(x<0.5-0.5*w0)) 
         out = 1+(w0+w1*s-1)*a/(w1*s)+m2*x;
      else if((x>=0.5-0.5*w0)&(x<0.5+0.5*w0)) 
         out = 1+a;
      else if((x>=0.5+0.5*w0)&(x<0.5+0.5*(w0+w1*s))) 
         out = 1+(w0+w1*s+1)*a/(w1*s)-m2*x;
      else if((x>=0.5+0.5*(w0+w1*s))&(x<1))  
         out = m1*(1-x);   
      else if((x>=1)&(x<(1+b))) 
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
/*   else if(a<0.0) 
   {
   if(x<0.0) 
      out = b+2+1/abs(a)*x;
   else if((x>=0.0) & (x<0.5-p2/2)) 
      out = 1-abs(a)-m3*x;
   else if((x>=0.5-p2/2) & (x<0.5)) 
      out = m4*(0.5-p2/2)-m4*x;
   else if((x>=0.5) & (x<0.5+p2/2)) 
      out = m4*(x-0.5-p2/2);
   else if((x>=0.5+p2/2) & (x<1.0)) 
      out = m3*(x-0.5-p2/2);
   else if((x>=1.0) & (x<(1+b))) 
      out = x+a;
   else if((x>=(1+b)) & (x<(b+1.5-p1/2))) 
      out = b+2+m1*(b+1)-m1*x;
   else if((x>=(b+1.5-p1/2)) & (x<(b+1.5))) 
      out = b+1+m2*(b+1.5-p1/2)-m2*x;
   else if((x>=(b+1.5)) & (x<(b+1.5+p1/2))) 
      out = m2*(x-b-1.5+(b+1-abs(a))/m2);
   else if((x>=(b+1.5+p1/2)) & (x<(2+b)))
      out = m1*(x-b-1.5-p1/2+(b+1)/m1);
   else if(x>=(2+b)) 
      out = -1/abs(a)*(b+2)+1/abs(a)*x;
   }*/
   return out;
}

// Enzima49
double enzima49(x,a,b,p1,p2)
double x,a,b,p1,p2;
{
   double out,m1,m2,m3,m4; 

   if(p1>=1)
      p1=0.99999;
   
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
      else if((x>=1)&(x<(1+b))) 
         out = a+x;
      else if(x>=(1+b)) 
         out = (x-2.)/a;
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

