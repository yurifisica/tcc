#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define length(a) sizeof(a)/sizeof(a[0])
const float pi = 3.1415926;
const float mu=4*3.1415926E-7;
const float eps= 8.85E-12;

// double _Complex ctanh(double _Complex x)
// {
//     if (creal(x)<=86)
//     {
//         return (cexp(x)-cexp(-x))/(cexp(x)+cexp(-x));
//     } else
//     {
//         return 1;
//     }
// }

double _Complex mult_complex(double _Complex Z1, double _Complex Z2)
{
  double _Complex Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2))+(creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I;
  return Z;
}

double _Complex div_complex(double _Complex Z1, double _Complex Z2)
{
  double _Complex Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2)+(-creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I)/(pow(creal(Z2),2)-pow(cimag(Z2),2));
  return Z;
}

//cálculo do campo primário
void campo_prim(int n,double sigma[],int n_h, double h[], double frq, double zp[], int np,double _Complex Hp[])
{
  FILE *file;
  file=fopen("a.dat","w");
  int cont,ic;
  if (file==NULL)
  {
    printf("caminho inválido");
  }
  double w=2*pi*frq;
  double _Complex Z=I*w*mu;
  double _Complex H[n];
  H[0]=1;
  double z[n_h+1];
  z[0]=0;
  for (cont=0;cont<n_h;cont++)
  {
    z[cont+1]=z[cont]+h[cont];
  }
  double _Complex Y[n+1];
  double _Complex u[n+1];
  double _Complex Zint[n+1];
  Y[0]=I*w*eps;
  u[0]=I*csqrt(-Z*Y[0]);
  Zint[0]=u[0]/Y[0];
  for (cont=1;cont<=n;cont++)
  {
    Y[cont]=sigma[cont-1]+I*w*eps;
    u[cont]=I*csqrt(-Z*Y[cont]);
    Zint[cont]=u[cont]/Y[cont];
  }
  double _Complex Zap[n];
  Zap[n-1]=Zint[n];
  for (cont=n-2;cont>=0;cont--)
  {
    double _Complex aux=u[cont+1]*h[cont];
    Zap[cont]=Zint[cont+1]*(Zap[cont+1]+Zint[cont+1]*ctanh(aux))/(Zint[cont+1]+Zap[cont+1]*ctanh(aux));
  }
  double _Complex rtm[n-1];
  rtm[0]=(Zint[0]-Zap[0])/(Zint[0]+Zap[0]);
  for (cont=1;cont<n;cont++)
  {
    double _Complex aux=u[cont]*h[cont-1];
    rtm[cont]=(Zint[cont]-Zap[cont])/(Zint[cont]+Zap[cont]);
    H[cont]=H[cont-1]*(1+rtm[cont-1])*cexp(-aux)/(1+rtm[cont]*cexp(-2*aux));
  }
  H[n]=H[n-1]*(1+rtm[n-1]);
  int cmd=0;
  for (cont=0;cont<np;cont++)
  {
    int flag=0;
    for (ic=0;ic<n;ic++)
    {
      if (zp[cont]<z[ic] && flag==0)
      {
        cmd=ic;
        flag=1;
      }
    }
    if (zp[cont]>z[0])
    {
      if (zp[cont]<z[n_h])
      {
        double _Complex aux=u[cmd]*(zp[cont]-z[cmd]);
        Hp[cont]=H[cmd]*(cexp(-aux)+rtm[cmd]*cexp(aux));
      }
      else
      {
        double _Complex aux=u[n]*(zp[cont]-z[cmd]);
        Hp[cont]=H[n]*cexp(-aux);
      }
    }
    else
    {
      double _Complex aux=u[0]*zp[cont];
      Hp[cont]=H[0]*(cexp(-aux)+rtm[0]*cexp(aux));
    }
    fprintf(file,"%i %f %f\n",cont,creal(Hp[cont]),cimag(Hp[cont]));
  }
  fclose(file);
  return;
}
double _Complex a(double x, double z){return 1;}
double _Complex b(double x, double z){return 0;}
double _Complex c(double x, double z){return 1;}
double _Complex d(double x, double z){return 0;}
double _Complex e(double x, double z){return 0;}
double _Complex f(double _Complex Z, double _Complex Y){return -Z*Y;}
double _Complex g(double _Complex Z, double _Complex Yj, double _Complex Yp, double _Complex Hp){
  // printf("%lf %lf\n", creal(Z),cimag(Z));
  // printf("%lf %lf\n", creal(Yj),cimag(Yj));
  // printf("%lf %lf\n", creal(Yp),cimag(Yp));
  // printf("%lf %lf\n", creal(Hp),cimag(Hp));
  // printf("dentro %lf %lf\n",creal(-Z*(Yj-Yp)*Hp),cimag(-Z*(Yj-Yp)*Hp));
  return -Z*(Yj-Yp)*Hp;}
  //programa principal
  int main()
  {
    int n_cam=2,cont,j,i;//numero de camadas
    double frq=5*pow(10,3);//frequência
    double espessura[1]={200};//espessura das camadas intermediárias
    size_t n_h=length(espessura);
    double sigma[2];
    sigma[0]=pow(10,-3);sigma[1]=0.5;
    size_t n_s=length(sigma);
    double z0=0;
    double zn=espessura[0]*2;//profundidade da malha
    int nz=10;//numero de pontos em z
    double dz1=abs(zn-z0)/(nz-1);
    double z[nz];
    for (cont=0;cont<nz;cont++)
    {
      z[cont]=cont*dz1;
      // printf("%lf\n",z[cont] );
    }
    nz=length(z);
    double x0=0;
    double xn=10000;//largura da malha
    int nx=15;//numero de pontos em x, na malha reduzida
    double x[nx+4],xt[nx];
    x[0]=-xn+xn/4;x[1]=-xn+xn/2;
    double dx1=abs(xn-x0)/(nx-1);
    for (cont=0;cont<nx;cont++)
    {
      xt[cont]=cont*dx1;
      x[cont+2]=cont*dx1;
      // printf("%lf\n", x[cont+2]);
    }
    x[nx+2]=xn+xn/4;x[nx+3]=xn+xn/2;
    nx=length(x);//numero de pontos em x, na malha completa
    int nx2=length(xt);
    int np=nx*nz;
    int ne=(nx-1)*(nz-1)*2;
    double xv[np],zv[np];
    double _Complex Y[np],Yj[np],sig[np];
    double w =2*pi*frq;
    Y[0]=I*w*eps;
    double _Complex Z=I*w*mu;
    // printf("%lf %lf\n",creal(Z),cimag(Z));
    double _Complex K[np][np];
    double _Complex k[3][3],m[3];
    int np2=(nx-2)*(nz-2);
    double dt;
    double dx[3],dz[3];
    int i1,j1,l;
    int jc=0;
    double sigma1=0.5;//solo embaixo dos prédios
    double sigma2=4.8;//água salgada
    double sigma3=2.8;//proporção de 1/2 de areia e cimento (prédios), com fibras de aço carbono
    double _Complex Hp[np];
    double _Complex al[3];
    double _Complex bl[3];
    double _Complex cl[3];
    double _Complex dl[3];
    double _Complex El[3];
    double _Complex fl[3];
    double _Complex gl[3];
    double xv2[nz*nx2];
    double _Complex Hs[np];
    int vc[np],qt_front;
    int cc;
    int el[(nz-1)*(nx-1)*2][4];
    double _Complex M[np];
    int flag;
    int ic=0;
    for (j=0;j<nz;j++)
    {
      for(i=0;i<nx;i++)
      {
        if (z[j]>=espessura[0])
        {
          sig[ic]=sigma[1];
          Yj[ic]=sigma[1]+I*w*eps;
          Y[ic]=sigma[1]+I*w*eps;
        }
        else
        {
          sig[ic]=sigma[0];
          Yj[ic]=sigma[0]+I*w*eps;
          Y[ic]=sigma[0]+I*w*eps;
        }
        xv[ic]=x[i];zv[ic]=z[j];
        ic++;
      }
    }


    ic=0;
    campo_prim(n_s,sigma,n_h,espessura,frq,zv,np,Hp);
    for (j=0;j<nz;j++)
    {
      for (i=0;i<nx2;i++)
      {
        xv2[ic]=x[i];
        ic++;
      }
    }
    //definição de coeficientes da edo
    cc=-1;
    ic=0;
    FILE *fileaux;
    fileaux=fopen("el.dat","w");

    for (j=0;j<(nz-1);j++)
    {
      for (i=0;i<nx;i++)
      {
        if(i<nx-1)
        {
          cc++;
          el[cc][0]=ic;
          el[cc][1]=ic+1;
          el[cc][2]=ic+nx;
          el[cc][3]=cc;
          fprintf(fileaux, "%i\t",el[cc][0]);
          fprintf(fileaux, "%i\t",el[cc][1]);
          fprintf(fileaux, "%i\t",el[cc][2]);
          fprintf(fileaux, "%i\n",el[cc][3]);
          cc++;
          el[cc][0]=ic+1;
          el[cc][1]=ic+nx+1;
          el[cc][2]=ic+nx;
          el[cc][3]=cc;
          fprintf(fileaux, "%i\t",el[cc][0]);
          fprintf(fileaux, "%i\t",el[cc][1]);
          fprintf(fileaux, "%i\t",el[cc][2]);
          fprintf(fileaux, "%i\n",el[cc][3]);
        } else
        {
          cont++;
        }
        ic++;
      }
    }
    // printf("%i\n",cc);
    // printf("%i\n",(nz-1)*(nx-1)*2);
    // printf("%i %i\n",nz,nx);
    //condições de fronteira
    qt_front=0;
    for (i=0;i<np;i++)
    {
      vc[i]=0;
      if (xv[i]==xv[0] || xv[i]==xv[np-1] || zv[i]==zv[0]|| zv[i]==zv[np-1])
      {
        Hs[i]=0; vc[i]=1;
        qt_front++;
      }
      // printf("%d\n",vc[i] );
    }
    double _Complex Kr[np-qt_front][np-qt_front];
    double _Complex Mr[np-qt_front];
    //definição dos meios secundários
    //meio secundário 1 (solo)
    // printf("%d\n",qt_front );

    for (i=0;i<np;i++)
    {
      // printf("%lf %lf\n",xv[i],zv[i]);
      if (xv[i]>=3000 && xv[i]<=10000 && zv[i]>espessura[0]/2 && zv[i]<=(3*espessura[0]/2))
      {
        Yj[i]=sigma1+I*w*eps; sig[i]=sigma1;
        // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      }
      M[i]=0;
    }

    // printf("asd\n");
    //meio secundário 2 (mar)
    for (i=0;i<np;i++)
    {
      if (xv[i]>=0 && xv[i]<=3000)
      {
        if ((zv[i]>=espessura[0]) && (zv[i]<(5*espessura[0])/4))
        {
          Yj[i]=sigma2+I*w*eps; sig[i]=sigma2;
          // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
        }
      }
    }
    ic=0;
    flag=0;
    for (j=1;j<(nz-1);j++)
    {
      for (i=1;i<(nx-1);i++)
      {
        ic++;
        if (xv[i]>=2000 && x[i]<=5000)
        {
          if (z[j]>=(((z[j]-z[j-1])/(x[i]-x[i-1]))*x[i]) && z[j]<=(3*espessura[0])/2 && z[j]>(13*espessura[0])/12)
          {
            Yj[i]=sigma2+I*w*eps; sig[i]=sigma2;
            // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
          }
        }
      }
    }
    //meio secundário 3 (prédios)
    for (i=0;i<np;i++)
    {
      // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      if (xv[i]>=6500 && xv[i]<=7000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
        Yj[i]=sigma3+I*w*eps;
        sig[i]=sigma3;
        // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      }
      if (xv[i]>=7500 && xv[i]<=8000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
        Yj[i]=sigma3+I*w*eps;
        sig[i]=sigma3;
        // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      }
      if (xv[i]>=8500 && xv[i]<=9500 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
        Yj[i]=sigma3+I*w*eps;
        sig[i]=sigma3;
        // printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      }
    }
    for (i=0;i<np;i++)
    {
      for (j=0;j<np;j++)
      {
          K[i][j]=0;
      }
      M[i]=0;
    }
    // construção da matriz e vetor da direita
    for (i=0;i<ne;i++){
      dx[0]=zv[el[i][1]]-zv[el[i][2]];
      dx[1]=zv[el[i][2]]-zv[el[i][0]];
      dx[2]=zv[el[i][0]]-zv[el[i][1]];
      dz[0]=xv[el[i][2]]-xv[el[i][1]];
      dz[1]=xv[el[i][0]]-xv[el[i][2]];
      dz[2]=xv[el[i][1]]-xv[el[i][0]];
      dt=abs(xv[el[i][1]]*zv[el[i][2]]+xv[el[i][0]]*zv[el[i][1]]+xv[el[i][2]]*zv[el[i][0]]-xv[el[i][1]]*zv[el[i][0]]-xv[el[i][2]]*zv[el[i][1]]-xv[el[i][0]]*zv[el[i][2]])/2;//delta
      //printf("%f\n",dt);
      for (cont=0;cont<3;cont++){
        al[cont]=a(xv[el[i][cont]],zv[el[i][cont]]);
        bl[cont]=b(xv[el[i][cont]],zv[el[i][cont]]);
        cl[cont]=c(xv[el[i][cont]],zv[el[i][cont]]);
        dl[cont]=d(xv[el[i][cont]],zv[el[i][cont]]);
        El[cont]=e(xv[el[i][cont]],zv[el[i][cont]]);
        fl[cont]=f(Z,Yj[el[i][cont]]);
        gl[cont]=g(Z,Yj[el[i][cont]],Y[el[i][cont]],Hp[el[i][cont]]);
        printf("fora %lf %lf\n", creal(gl[cont]),cimag(gl[cont]));
        //fprintf(fileaux,"%f\n",creal(fl[cont]));
        // if (cabs(gl[cont])!= 0)
        // {
        // printf("\n%i %f %f",el[i][cont],creal(gl[cont]),cimag(gl[cont]));
        // }
        // printf("\n%i %f %f",el[i][cont],creal(gl[cont]),cimag(gl[cont]));
      }
      //matriz de cada elemento
      for(cc=0;cc<3;cc++)
      {
        for(j=0;j<3;j++)
        {
          k[cc][j]=0;
        }
      }
      for (cc=0;cc<3;cc++)
      {
        for(j=0;j<3;j++)
        {
          k[cc][j]=(-1/(4*dt))*(dx[cc]*dx[j]*al[j]+dx[cc]*dz[j]*bl[j]+dz[cc]*dz[j]*cl[j])+(1/6)*(dx[j]*dl[j]+dz[j]*El[j]);
          if (j==cc)
          {
            k[cc][j]+=2*fl[j];
          } else
          {
            k[cc][j]+=fl[j];
          }
          // printf("%d %d >> %lf %lf\n",cc,j, creal(k[cc][j]),cimag(k[cc][j]));
        }
      }

      //vetor da direita em cada elemento
      m[0]=(dt/12)*(2*gl[0]+gl[1]+gl[2]);
      m[1]=(dt/12)*(gl[0]+2*gl[1]+gl[2]);
      m[2]=(dt/12)*(gl[0]+gl[1]+2*gl[2]);
      // printf("%lf\n",creal(m[0]));
      // printf("%lf\n",creal(m[1]));
      // printf("%lf\n",creal(m[2]));
      //Matriz global vetor da direita
      for (j1=0;j1<3;j1++){
        for (i1=0;i1<3;i1++){
          K[el[i][j1]][el[i][i1]]=K[el[i][j1]][el[i][i1]]+k[j1][i1];
          // printf("\n%d %d %5.15f %5.15f",el[i][j1],el[i][i1],creal(K[el[i][j1]][el[i][i1]]),cimag(K[el[i][j1]][el[i][i1]]));
          // printf("\n%5.15f %5.15f",creal(M[el[i][j1]]),cimag(M[el[i][j1]]));
        }
        M[el[i][j1]]=M[el[i][j1]]+m[j1];
      }
    }
    for (i=0;i<np;i++)
    {
      for (j=0;j<np;j++)
      {
        // printf("%2.2lf ",creal(K[i][j]));
      }
      // printf("%d >> %lf %lf\n",i,creal(M[i]),cimag(M[i]));
      // printf("\n");
    }
    //redução do sistema - Inclusão das condições de fronteira

    for (i=0;i<np-qt_front;i++)
    {
      for (j=0;j<np-qt_front;j++)
      {
        Kr[i][j]=0;
      }
      Mr[i]=0;
    }
    ic=0;
    jc=0;
    for (i=0;i<np;i++)
    {
      if (vc[i]!=1)
      {
        Mr[ic]=M[i];
        jc=0;
        for (j=0;j<np;j++)
        {
          if (vc[j]!=1)
          {
            Kr[ic][jc]=K[i][j];
            jc++;
          } else
          {
            Mr[ic]-=K[i][j]*Hs[j];
          }
        }
        ic++;
      }
    }
    // printf("%d %d\n",ic,jc);
    // printf("%lf %lf\n",creal(Kr[0][0]),cimag(Kr[0][0]));
    for (i=0;i<np-qt_front;i++)
    {
      for (j=0;j<np-qt_front;j++)
      {
        // printf("%2.2lf ",creal(Kr[i][j]));
      }
      // printf("%2.2lf\n",creal(Mr[i]));
      // printf("\n");
    }
    // printf("%2.2lf\n",creal(Kr[1][1]));
    //solução por eliminação gaussiana
    //1-triangularização
    double _Complex Ur[np-qt_front];
    double _Complex temp;
    for (i=0;i<(np-qt_front-1);i++)
    {
      for (j=i+1;j<np-qt_front;j++)
      {
        temp=Kr[j][i]/Kr[i][i];
        //printf(">> %lf <> %lf\n",creal(Kr[i][i]),cimag(Kr[i][i]));
        for (l=i;l<np-qt_front;l++)
        {
          Kr[j][l]-=temp*Kr[i][l];
        }
        Mr[j]-=temp*Mr[i];
      }
    }
    //2-retrossubstituição
    Ur[np-qt_front]=Mr[np-qt_front]/Kr[np-qt_front][np-qt_front];
    for (i=(np-qt_front-1);i>=0;i--)
    {
      Ur[i]=Mr[i];
      for (j=i+1;j<np-qt_front;j++)
      {
        Ur[i]-=Kr[i][j]*Ur[j];
      }
      Ur[i]=Ur[i]/Kr[i][i];
    }
    //Incorporando as condições de fronteira
    double _Complex U[np];
    ic=0;
    for (i=0;i<np;i++)
    {
      if (vc[i]!=1)
      {
        U[i]=Ur[ic];
        ic++;
      }
      else
      {
        U[i]=Hs[i];
      }
    }
    double _Complex Htotal[np];
    FILE *file;
    file=fopen("H_sec.dat","w");
    for (i=0;i<np;i++)
    {
      Htotal[i]=U[i]+Hp[i];
      //printf("%i %f %f\n",i,creal(Htotal[i]),cimag(Htotal[i]));
      fprintf(file,"%i %f %f\n",i,creal(Htotal[i]),cimag(Htotal[i]));
    }

    return 0;
  }
