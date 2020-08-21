//
//  main.c
//  Shredinger
//
//  Created by Никита Пухляков on 18/04/2020.
//  Copyright © 2020 Никита Пухляков. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI = 3.1415926535897932384626433832795;
const double alpha = 10;
const double beta = 4;

void scheme(int M, int N, int L, double h_x, double h_y, double tau, double *v, double *w,
            double *G, double *R, double *c_v, double *c_w, double *c_G, double *c_R,
            double *f, double *d, double *bas);
double func(double x, double y);
void print_file(FILE *f1, int M, int N, int L, double h_x, double h_y, double tau,
                double *v, double *w);
void print_screen(int M, int N, int L, double h_x, double h_y, double tau,
                  double *v, double *w);
double scalar_x(double h_x, int M, double *f, double *g);
double scalar_y(double h_y, int N, double *f, double *g);
double scalar_classic(int N, double *f, double *g);
void basis_x(int M, int m, double *f);
void basis_y(int N, int n, double *f);
void basis_2_x(int M, double x, double *f);
void basis_2_y(int N, double y, double *f);
void FtoC( int M, int N, double *F, double*c, double *f, double *d, double *bas);
void CtoF(int l, int M, int N, double *pr, double*c, double *f, double *d, double *bas);
void zero(int n, double *x);
double lambda(int m, int n, double h_x, double h_y);
void print_vector(int n, double *x);
double f_test(double x, double y);
void Print_matrix(int M, int N, double *x);
double eps(int K, int M, int N, int L, double h_x, double h_y, double tau, double *v, double *w,
           double *G, double *R, double *c_v, double *c_w, double *c_G, double *c_R,
           double *f, double *d, double *bas);
void animation(FILE* f0,FILE* f1,FILE* f2,FILE* f3,FILE* f4,FILE* f5,FILE* f6,FILE* f7,
                 int M, int N,int L, double h_x, double h_y, 
                double tau, double *v, double *w);

int main(void)
{
    FILE* f0; FILE* f1; FILE* f2; FILE* f3; FILE* f4; FILE* f5; FILE* f6; FILE* f7;
    FILE* f8;FILE* f9;FILE* f10;
    int N,M,L,max,K;
    double T,h_x,h_y,tau,*v,*w, *d, *bas, *f, *c_v,*c_w,*c_G,*c_R,*G,*R;
    f0=fopen("out0.txt","w");
    f1=fopen("out1.txt","w");
    f2=fopen("out2.txt","w");
    f3=fopen("out3.txt","w");
    f4=fopen("out4.txt","w");
    f5=fopen("out5.txt","w");
    f6=fopen("out6.txt","w");
    f7=fopen("out7.txt","w");
    printf("Введите T: ");
    scanf("%lf",&T);
    printf("Введите M: ");
    scanf("%d",&M);
    printf("Введите N: ");
    scanf("%d",&N);
    printf("Введите L: ");
    scanf("%d",&L);
    h_x = 1/(double)M;
    h_y = 1/(double)N;
    tau = T/(double)L;
    max = M+1;
    if (max<N+1) {max = N+1;}
    
    v = (double*)malloc(sizeof(double)*(M+1)*(N+1)*(L+1));
    w = (double*)malloc(sizeof(double)*(M+1)*(N+1)*(L+1));
    d = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    bas = (double*)malloc(sizeof(double)*max);
    f = (double*)malloc(sizeof(double)*max);
    c_G = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    c_R = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    c_v = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    c_w = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    G = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    R = (double*)malloc(sizeof(double)*(M+1)*(N+1));
    
    zero((M+1)*(N+1)*(L+1),v);
    zero((M+1)*(N+1)*(L+1),w);
    zero((M+1)*(N+1),d);
    zero(max,bas);
    zero(max,f);
    zero((M+1)*(N+1),c_G);
    zero((M+1)*(N+1),c_R);
    zero((M+1)*(N+1),c_v);
    zero((M+1)*(N+1),c_w);
    zero((M+1)*(N+1),G);
    zero((M+1)*(N+1),R);
    
    scheme(M,N,L,h_x,h_y,tau,v,w,G,R,c_v,c_w,c_G,c_R,f,d,bas);
    

    for(K=2;K<=16;K*=2)
    { printf(" Сходимость по сетке при K = %d ",K);
      printf("  %e \n", eps(K,M,N,L,h_x,h_y,tau,v,w,G,R,c_v,c_w,c_G,c_R,f,d,bas)
               /(tau+pow(h_x,2)+pow(h_y,2)));
    }


    //print_file(f1,M,N,L,h_x,h_y,tau,v,w);
    //print_screen(M,N,L,h_x,h_y,tau,v,w);
    
    //animation(f0,f1,f2,f3,f4,f5,f6,f7,M,N,L,h_x,h_y,tau,v,w);

    free(v);
    free(w);
    free(d);
    free(bas);
    free(f);
    free(c_G);
    free(c_R);
    free(c_v);
    free(c_w);
    free(G);
    free(R);
    fclose(f0);
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);

    return 0;
}

void animation(FILE* f0,FILE* f1,FILE* f2,FILE* f3,FILE* f4,FILE* f5,FILE* f6,FILE* f7,
                 int M, int N,int L, double h_x, double h_y, 
                double tau, double *v, double *w)
{
        int i,j,l; 

    //Печатаем в файл f_i поверхность в момент t=i*const

      for(l=0;l<=7*50;l+=50)      
       {

        switch ( l ) {
        case 0:   
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f0,"%e   ", i*h_x);
            //fprintf(f0,"%e   ", j*h_y);
            //fprintf(f0,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);  
            //fprintf(f0,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);  
            fprintf(f0,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f0," \n");
           }
            break;    
        case 50*1:   
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f1,"%e   ", i*h_x);
            //fprintf(f1,"%e   ", j*h_y);
            //fprintf(f1,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]); 
            //fprintf(f1,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);    
            fprintf(f1,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f1," \n");
           }
            break;
        case 50*2:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f2,"%e   ", i*h_x);
            //fprintf(f2,"%e   ", j*h_y);
            //fprintf(f2,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);    
            //fprintf(f2,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]); 
            fprintf(f2,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f2," \n");
           }
            break;
        case 50*3:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f3,"%e   ", i*h_x);
            //fprintf(f3,"%e   ", j*h_y);
            //fprintf(f3,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]); 
            //fprintf(f3,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);    
            fprintf(f3,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f3," \n");
           }
            break;
        case 50*4:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f4,"%e   ", i*h_x);
            //fprintf(f4,"%e   ", j*h_y);
            //fprintf(f4,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);   
            //fprintf(f4,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);  
            fprintf(f4,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            } 
            fprintf(f4," \n"); 
           }
            break;
        case 50*5:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f5,"%e   ", i*h_x);
            //fprintf(f5,"%e   ", j*h_y);
           // fprintf(f5,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);   
            //fprintf(f5,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);  
            fprintf(f5,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            } 
            fprintf(f5," \n"); 
           }
            break; 
        case 50*6:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f6,"%e   ", i*h_x);
            //fprintf(f6,"%e   ", j*h_y);
            //fprintf(f6,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]); 
            //fprintf(f6,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);    
            fprintf(f6,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f6," \n");
           }
            break;
        case 50*7:
            for(i=0;i<=M;i++)
           {
             for(j=0;j<=N;j++)
            {       
            //fprintf(f7,"%e   ", i*h_x);
            //fprintf(f7,"%e   ", j*h_y);
            //fprintf(f7,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);    
            //fprintf(f7,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]); 
            fprintf(f7,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            }  
            fprintf(f7," \n");
           }
           
            break;           
        default:
            printf( "Неправильный ввод.\n" );

        }
    }
}

void print_file(FILE *f1, int M, int N, int L, double h_x, double h_y, double tau,
                double *v, double *w)
{
    int i,j,l; 

    //Фиксируем t
    l = L;
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            
            fprintf(f1,"%e   ", i*h_x);
            fprintf(f1,"%e   ", j*h_y);
            //Re(u)
            //fprintf(f1,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);
            //Im(u)
            //fprintf(f1,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);
            //|u|
            fprintf(f1,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            fprintf(f1," \n");
        }
    }
    

    //Фиксируем x и печатаем в файл
  /*  i = M/2;
    for(l=0;l<=L;l++)
    {
        for(j=0;j<=N;j++)
        {
            fprintf(f1,"%e   ", l*tau);
            fprintf(f1,"%e   ", j*h_y);
            //Re(u)
            //fprintf(f1,"%e   ", v[l*(M+1)*(N+1)+i*(N+1)+j]);
            //Im(u)
            fprintf(f1,"%e   ", w[l*(M+1)*(N+1)+i*(N+1)+j]);
            //|u|
            //fprintf(f1,"%e   ", sqrt(pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2)));
            
            fprintf(f1," \n");
            
        }
    }
   */ 
    
}

void print_screen(int M, int N, int L, double h_x, double h_y, double tau,
                  double *v, double *w)
{
    
    int i,j,l;
    //Фиксируем t
    l = L/2;
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            
            printf("%e   ", i*h_x);
            printf("%e   ", j*h_y);
            printf("%e   ", pow(v[l*(M+1)*(N+1)+i*(N+1)+j],2)+
                   pow(w[l*(M+1)*(N+1)+i*(N+1)+j],2));
            printf(" \n");
        }
    }
    
    printf(" \n");
    
}
void print_vector(int n, double *x)
{
    int i;
    for(i=0;i<n;i++)
    {
        printf("  %e",x[i]);
    }
    printf("\n");
}

void Print_matrix(int M, int N, double *x)
{
    int i,j;
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            printf("%10.3e ",x[i*(N+1)+j]);
        }
        printf("\n");
    }
    
    printf("\n");
}
void scheme(int M, int N, int L, double h_x, double h_y, double tau, double *v, double *w,
            double *G, double *R, double *c_v, double *c_w, double *c_G, double *c_R,
            double *f, double *d, double *bas)
{
    int i,j,m,n,l,new,old;
    //заполняем краевые и начальные условия
    for(l=0;l<=L;l++)
    {
        new = l*(M+1)*(N+1);
        
        for(i=0;i<=M;i++)
        {
            for(j=0;j<=N;j++)
            {
                if (l==0)
                {
                    v[i*(N+1)+j] = func(i*h_x,j*h_y);
                    w[i*(N+1)+j] = 0;
                }
                
                if ((i==0)||(i==M)||(j==0)||(j==N))
                {
                    v[new+i*(N+1)+j] =  0;
                    w[new+i*(N+1)+j] =  0;
                }
            }
        }
    }
    
    //находим u на остальных слоях методом Фурье
    for(l=1;l<=L;l++)
    {
        new = l*(M+1)*(N+1);
        old = (l-1)*(M+1)*(N+1);
        
        //заполним массивы правых частей уравнений
        for(i=0;i<=M;i++)
        {
            for(j=0;j<=N;j++)
            {
                G[i*(N+1)+j] = v[old+i*(N+1)+j]/tau+alpha*w[old+i*(N+1)+j]*
                (pow(v[old+i*(N+1)+j],2)+pow(w[old+i*(N+1)+j],2));
                R[i*(N+1)+j] = w[old+i*(N+1)+j]/tau-alpha*v[old+i*(N+1)+j]*
                (pow(v[old+i*(N+1)+j],2)+pow(w[old+i*(N+1)+j],2));
                
            }
        }
        
        //находим коэф-ты разложения по базису правой части первого уравнения G
        FtoC(M,N,G,c_G,f,d,bas);
        //находим коэф-ты разложения по базису правой части второго уравнения R
        FtoC(M,N,R,c_R,f,d,bas);
        

        //находим коэф-ты разложения по базису cv и cw соответсвенно действ и мнимой части u
        for(m=0;m<=M;m++)
        {
            for(n=0;n<=N;n++)
            {
                if ((m==0)||(m==M)||(n==0)||(n==N))
                {
                    c_v[m*(N+1)+n] = 0;
                    c_w[m*(N+1)+n] = 0;
                }
                else
                {
                    
                    
                    c_v[m*(N+1)+n] = -(c_R[m*(N+1)+n]-c_G[m*(N+1)+n]/(tau*lambda(m,n,h_x,h_y)))
                     /(lambda(m,n,h_x,h_y)+1/(lambda(m,n,h_x,h_y)*pow(tau,2)));
                    c_w[m*(N+1)+n] = (c_G[m*(N+1)+n]-c_v[m*(N+1)+n]/tau)/lambda(m,n,h_x,h_y);
                    
                }

            }
        }
        


        //вычисляем действительную часть u, т.е. v
        CtoF(l,M,N,v,c_v,f,d,bas);
        //вычисляем мнимую часть u, т.е. w
        CtoF(l,M,N,w,c_w,f,d,bas);
        
    }
    

    zero((M+1)*(N+1),d);
    zero((N+1),bas);
    zero((N+1),f);
    zero((M+1)*(N+1),c_G);
    zero((M+1)*(N+1),c_R);
    zero((M+1)*(N+1),c_v);
    zero((M+1)*(N+1),c_w);
    zero((M+1)*(N+1),G);
    zero((M+1)*(N+1),R);
}

double eps(int K, int M, int N, int L, double h_x, double h_y, double tau, double *v, double *w,
           double *G, double *R, double *c_v, double *c_w, double *c_G, double *c_R,
           double *f, double *d, double *bas)
{
    int max,i,j,new;
    double norm,tau_2,h_x_2,h_y_2,*v_2,*w_2, *d_2, *bas_2, *f_2,
    *c_v_2,*c_w_2,*c_G_2,*c_R_2,*G_2,*R_2;
    
    max = M+1;
    if (max<N+1) {max = N+1;}
    v_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*(L+1)*K*K*K);
    w_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*(L+1)*K*K*K);
    d_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    bas_2 = (double*)malloc(sizeof(double)*max*K);
    f_2 = (double*)malloc(sizeof(double)*max*K);
    c_G_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    c_R_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    c_v_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    c_w_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    G_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    R_2 = (double*)malloc(sizeof(double)*(M+1)*(N+1)*K*K);
    
    tau_2=tau/(double)K;
    h_x_2=h_x/(double)K;
    h_y_2=h_y/(double)K;
    
    scheme(M,N,L,h_x_2,h_y_2,tau,v,w,G,R,c_v,c_w,c_G,c_R,f,d,bas);
    //L*K
    scheme(M*K,N*K,L*K,h_x_2,h_y_2,tau_2,v_2,w_2,G_2,R_2,c_v_2,c_w_2,c_G_2,c_R_2,f_2,d_2,bas_2);
    
    new = L*(M+1)*(N+1);
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            v[new+i*(N+1)+j] = fabs(v[new+i*(N+1)+j]-v_2[K*K*K*new+i*(N+1)*K+j*K]);
            w[new+i*(N+1)+j] = fabs(w[new+i*(N+1)+j]-w_2[K*K*K*new+i*(N+1)*K+j*K]);
        }
    }
    
    norm = 0;
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
        {
            norm+= pow(v[new+i*(N+1)+j],2) + pow(w[new+i*(N+1)+j],2);
        }
    }
  
    norm *= h_x*h_y;
    norm = sqrt(norm);

    free(v_2);
    free(w_2);
    free(d_2);
    free(bas_2);
    free(f_2);
    free(c_G_2);
    free(c_R_2);
    free(c_v_2);
    free(c_w_2);
    free(G_2);
    free(R_2);
    return norm;
}



double func(double x, double y)
{
    return exp(-beta*(pow((x-0.5),2)+pow((y-0.5),2)));
}

double f_test(double x, double y)
{
    return 2*sin(3*PI*x)*sin(4*PI*y)+2*sin(1*PI*x)*sin(1*PI*y);
}

double lambda(int m, int n, double h_x, double h_y)
{
    return 4*pow(sin(PI*m*h_x/2.0)/h_x,2)+4*pow(sin(PI*n*h_y/2.0)/h_y,2);
}

double scalar_x(double h_x, int M, double *f, double *g)
{
    int i;
    double res;
    res=0;
    for (i=1;i<=M-1;i++)
    {
        res += g[i]*f[i];
    }
    res *= h_x;
    return res;
}

double scalar_y(double h_y, int N, double *f, double *g)
{
    int i;
    double res;
    res=0;
    //print_vector(N,f);
    //print_vector(N,g);
    for (i=1;i<=N-1;i++)
    {
        res += g[i]*f[i];
        //printf("f[i] %e\n", f[i]);
        //printf("g[i] %e\n", g[i]);
    }
    res *= h_y;
    //printf("h_y %e\n", h_y);
    return res;
}

double scalar_classic(int N, double *f, double *g)
{
    int i;
    double res,temp;
    res=0;
    temp = 0;
    for (i=0;i<=N;i++)
    { temp = g[i]*f[i];
        res += temp;
    }
    return res;
}

void basis_x(int M, int m, double *f)
{   // f  это вектор значений в узлах m-й базисной функции
    int j;
    double h_x = 1/(double)M;
    
    if (m==0)
    {
        for (j=0;j<=M;j++)
            f[j] = 0;
    }
    if (m==M)
    {
        for (j=0;j<=M;j++)
            f[j] = 0;
    }
    if ((0<m)&(m<M))
    {   for (j=0;j<=M;j++)
        f[j] = sqrt(2)*sin(PI*m*j*h_x);
    }
}

void basis_y(int N, int n, double *f)
{   // f  это вектор значений в узлах n-й базисной функции
    int j;
    double h_y = 1/(double)N;
    if (n==0)
    {
        for (j=0;j<=N;j++)
            f[j] = 0;
    }
    if (n==N)
    {
        for (j=0;j<=N;j++)
            f[j] = 0;
    }
    if ((0<n)&(n<N))
    {   for (j=0;j<=N;j++)
        f[j] = sqrt(2)*sin(PI*n*j*h_y);
    }
    
}

void basis_2_x(int M, double x, double *f)
{   //f  это вектор значений базисных функций в точке x
    int m;
    f[0] = 0;
    f[M] = 0;
    for (m=1;m<=M-1;m++)
        f[m] = sqrt(2)*sin(PI*m*x);
    
}

void basis_2_y(int N, double y, double *f)
{   //f  это вектор значений базисных функций в y
    int n;
    f[0] = 0;
    f[N] = 0;
    for (n=1;n<=N-1;n++)
        f[n] = sqrt(2)*sin(PI*n*y);
    
}

void FtoC( int M, int N, double *F, double*c, double *f, double *d, double *bas)
{
    int n,i,j;
    double h;
    
    //находим матрицу d промежуточных коэффициентов
    for(i=0;i<=M;i++)
    {
        for(j=0;j<=N;j++)
            f[j] = F[i*(N+1)+j];
        
        
        
        for(n=1;n<=N-1;n++)
        {
            basis_y(N,n,bas);
            h = 1/(double)N;
            d[i*(N+1)+n] = scalar_y(h,N,bas,f);
        }
    }
    
    
    
    //находим матрицу с коэффициентов Фурье разложения нашей функции F по базису
    for (n=1;n<=N-1;n++)
    {
        for(i=0;i<=M;i++)
        {
            basis_x(M,i,bas);
            h = 1/(double)M;
            for(j=0;j<=M;j++)
                f[j] = d[j*(M+1)+n];
            
            c[i*(N+1)+n] = scalar_x(h,M,bas,f);
            
        }
    }
    
    zero((M+1)*(N+1),d);
    
}

void CtoF(int l, int M, int N, double *pr, double*c, double *f, double *d, double *bas)
{
    int n,i,j,k,new;
    double h;
    new = l*(M+1)*(N+1);
    
    //ищем приближенные значения
    //сначала опять матрицу d
    h = 1/(double)M;
    for(n=1;n<=N-1;n++)
        for(i=0;i<=M;i++)
        {
            for(j=0;j<=M;j++)
                f[j] = c[j*(N+1)+n];
            basis_2_x(M,i*h,bas);
            d[i*(N+1)+n] = scalar_classic(N,bas,f);
        }
    // окончательно находим приближенные значения
    for(i=0;i<=M;i++)
        for(j=0;j<=N;j++)
        {
            f[0] = 0;
            f[N] = 0;
            for(k=1;k<=N-1;k++)
                f[k] = d[i*(N+1)+k];
            
            basis_2_y(N,j*h,bas);
            pr[new+i*(N+1)+j] = scalar_classic(N,f,bas);
            
        }
    
    zero((M+1)*(N+1),d);
}


void zero(int n, double *x)
{
    int i;
    for(i=0;i<n;i++)
    {
        x[i] = 0;
    }
    
}


