#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <cmath>
#include <math.h>
#include <unistd.h>
#define length(a) sizeof(a)/sizeof(a[0])

//variaveis
const float pi = acos(-1.);
const float mu=4*3.1415926E-7;
const float eps= 8.85E-12;
std::complex<double> *Hp,*M, **Kr, *Mr,*Y,*Yj,**K,k[3][3],m[3],*Hs,*Htotalc,*Htotal,*Ur,temp,*U;
double *xv,*zv,*xv2,*zv2,*sigv2;
double x0=0;
double xn=10000;//largura da malha
int nx=30;//numero de pontos em x, na malha reduzida
int n_cam=2,cont,j,i,aux,cc,ic;//numero de camadas
double frq=5E3;//frequência
double z0=0;
int nz=30;//numero de pontos em z
int **elc,**el,*vc,qt_front;

void campo_prim(int n,double sigma[],int n_h, double h[], double zp[], int np)
{
  int cont,ic;
  FILE *file;
  file=fopen("a.dat","w");

  if (file==NULL)
  {
    printf("caminho inválido");
  }
  double w=2*pi*frq;
  std::complex<double> Z (0,w*mu);
  std::complex<double> H[n];
  std::complex<double> tgh;
  std::complex<double> aux;
  H[0]=1;
  double z[n_h+1];
  z[0]=0;
  for (cont=0;cont<n_h;cont++)
  {
    z[cont+1]=z[cont]+h[cont];
  }
  std::complex<double> Yp[n+1];
  std::complex<double> u[n+1];
  std::complex<double> Zint[n+1];
  Yp[0]={0,w*eps};
  // u[0]={0,sqrt(pow(w,2)*mu*eps)};
  u[0]=sqrt(Z*Yp[0]);
  Zint[0]=u[0]/Yp[0];
  for (cont=1;cont<n+1;cont++)
  {
    Yp[cont]={sigma[cont-1],w*eps};
    u[cont]=sqrt(Z*Yp[cont]);
    Zint[cont]=u[cont]/Yp[cont];
  }
  std::complex<double> Zap[n];
  Zap[n-1]=Zint[n];
  for (cont=n-2;cont>=0;cont--)
  {
    aux=u[cont+1]*h[cont];
    if (real(aux)<=86)
    {
      tgh=(exp(aux)-exp(-aux))/(exp(aux)+exp(-aux));
    } else {
      tgh=1;
    }
    Zap[cont]=Zint[cont+1]*(Zap[cont+1]+Zint[cont+1]*tgh)/(Zint[cont+1]+Zap[cont+1]*tgh);
  }
  std::complex<double> rtm[n-1];
  rtm[0]=(Zint[0]-Zap[0])/(Zint[0]+Zap[0]);
  for (cont=1;cont<n;cont++)
  {
  aux=u[cont]*h[cont-1];
    rtm[cont]=(Zint[cont]-Zap[cont])/(Zint[cont]+Zap[cont]);
    H[cont]=H[cont-1]*(1.+rtm[cont-1])*exp(-aux)/(1.+rtm[cont]*exp(-2.*aux));
  }
  H[n]=H[n-1]+H[n-1]*rtm[n-1];
  int cmd;
  int flag,i;
  for (ic=0;ic<np;ic++)
  {
    flag=0;
    for (i=0;i<n;i++)
    {
      if (zp[ic]<z[i] && flag==0)
      {
        cmd=i;
        flag=1;
      }
    }
    if (zp[ic]>z[0])
    {
      if (zp[ic]<z[n-1])
      {
      aux=u[cmd]*(zp[ic]-z[cmd]);
        Hp[ic]=H[cmd]*(exp(-aux)+rtm[cmd]*exp(aux));
      }
      else
      {
      aux=u[n]*(zp[ic]-z[n-1]);
        Hp[ic]=H[n]*exp(-aux);
      }
    }
    else
    {
    aux=u[0]*zp[ic];
      Hp[ic]=H[0]*(exp(-aux)+rtm[0]*exp(aux));
    }
    fprintf(file,"%f %f %f\n",zp[ic],real(Hp[ic]),imag(Hp[ic]));
  }
  fclose(file);
  return;
}
std::complex<double> a(double x, double z){return 1;}
std::complex<double> b(double x, double z){return 0;}
std::complex<double> c(double x, double z){return 1;}
std::complex<double> d(double x, double z){return 0;}
std::complex<double> e(double x, double z){return 0;}
std::complex<double> f(std::complex<double> Z, std::complex<double> Y){return -Z*Y;}
std::complex<double> g(std::complex<double> Z, std::complex<double> Yj, std::complex<double> Yp, std::complex<double> Hp){
  return -Z*(Yj-Yp)*Hp;}
  //programa principal
  int main()
  {
    double w=2*pi*frq;
    double espessura[n_cam-1]={400};//espessura das camadas intermediárias
    int n_h=length(espessura);
    double sigma[n_cam];
    sigma[0]=pow(10,-3);sigma[1]=0.5;
    int n_s=length(sigma);
    double zn=espessura[0]*2;//profundidade da malha
    double dz1=abs(zn-z0)/(nz-1);
    double z[nz];
    for (cont=0;cont<nz;cont++)
    {
      z[cont]=cont*dz1;
    }
    nz=length(z);
    double x[nx+4],xt[nx];
    x[0]=-xn+xn/4;x[1]=-xn+xn/2;
    double dx1=abs(xn-x0)/(nx-1);
    for (cont=0;cont<nx;cont++)
    {
      xt[cont]=cont*dx1;
      x[cont+2]=cont*dx1;
    }
    x[nx+2]=xn+xn/4;x[nx+3]=xn+xn/2;
    nx=length(x);//numero de pontos em x, na malha completa
    int nx2=length(xt);
    int np=nx*nz;
    int ne=(nx-1)*(nz-1)*2;
    double sig[np];
    std::complex<double> Z={0,w*mu};
    int np2=(nx-4)*(nz);
    double dt;
    double dx[3],dz[3];
    int i1,j1,l;
    int jc=0;
    double sigma2=4.8;//água salgada
    double sigma3=2.8;//proporção de 1/2 de areia e cimento (prédios), com fibras de aço carbono
    std::complex<double> a1;
    std::complex<double> a2;
    std::complex<double> a3;
    std::complex<double> b1;
    std::complex<double> b2;
    std::complex<double> b3;
    std::complex<double> c1;
    std::complex<double> c2;
    std::complex<double> c3;
    std::complex<double> d1;
    std::complex<double> d2;
    std::complex<double> d3;
    std::complex<double> e1;
    std::complex<double> e2;
    std::complex<double> e3;
    std::complex<double> f1;
    std::complex<double> f2;
    std::complex<double> f3;
    std::complex<double> g1;
    std::complex<double> g2;
    std::complex<double> g3;
    // int nec=(nx2-1)*(nz-1)*2;//numeros de elementos da malha reduzida
    Hp=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    Y=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    Yj=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    xv=(double*)malloc(np*sizeof(double));
    zv=(double*)malloc(np*sizeof(double));
    for (j=0;j<nz;j++)
    {
      for(i=0;i<nx;i++)
      {
        if (z[j]>=espessura[0])
        {
          sig[ic]=sigma[1];
          Yj[ic]={sigma[1],w*eps};
          Y[ic]={sigma[1],w*eps};
        }
        else
        {
          sig[ic]=sigma[0];
          Yj[ic]={sigma[0],w*eps};
          Y[ic]={sigma[0],w*eps};
        }
        xv[ic]=x[i];zv[ic]=z[j];
        ic++;
      }
    }

    campo_prim(n_s,sigma,n_h,espessura,zv,np);
    xv2=(double*)malloc(np2*sizeof(double));
    zv2=(double*)malloc(np2*sizeof(double));
    ic=0;
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
    el=(int**)malloc(ne*sizeof(int*));
    for (i=0;i<ne;i++)
    {
      el[i]=(int*)malloc(4*sizeof(int));
    }
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
    fclose(fileaux);
    //condições de fronteira
    qt_front=0;
    Hs=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    vc=(int*)malloc(np*sizeof(int));
    for (i=0;i<np;i++)
    {
      vc[i]=0;
      if (xv[i]==xv[0] || xv[i]==xv[np-1] || zv[i]==zv[0]|| zv[i]==zv[np-1])
      {
        Hs[i]=0; vc[i]=1;
        qt_front++;
      }
    }


    //definição dos meios secundários
    //meio secundário 2 (mar)
    for (i=0;i<np;i++)
    {
      if (xv[i]>=0 && xv[i]<=3000)
      {
        if ((zv[i]>=espessura[0]) && (zv[i]<(5*espessura[0])/4))
        {
          Yj[i]={sigma2,w*eps}; sig[i]=sigma2;
          // printf("\n%i %f",i,real(Yj[i]-Y[i]));
        }
      }
    }
    //meio secundário 3 (prédios)
    for (i=0;i<np;i++)
    {
      // printf("\n%i %f",i,real(Yj[i]-Y[i]));
      if (xv[i]>=6500 && xv[i]<=7000 && zv[i]>=150 && zv[i]<espessura[0]) {
        Yj[i]={sigma3,w*eps};
        sig[i]=sigma3;
      }
      if (xv[i]>=7500 && xv[i]<=8000 && zv[i]>=150 && zv[i]<espessura[0]) {
        Yj[i]={sigma3,w*eps};
        sig[i]=sigma3;
      }
      if (xv[i]>=8500 && xv[i]<=9500 && zv[i]>=150   && zv[i]<espessura[0]) {
        Yj[i]={sigma3,w*eps};
        sig[i]=sigma3;
      }
    }
    K=(std::complex<double>**)malloc(np*sizeof(std::complex<double>*));
    M=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    for (i=0;i<np;i++)
    {
      K[i]=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
      for (j=0;j<np;j++)
      {
        K[i][j]=0;
      }
      M[i]=0;
    }
    // construção da matriz e vetor da direita
    for (i=0;i<ne;i++)
    {
      dx[0]=zv[el[i][1]]-zv[el[i][2]];
      dx[1]=zv[el[i][2]]-zv[el[i][0]];
      dx[2]=zv[el[i][0]]-zv[el[i][1]];
      dz[0]=xv[el[i][2]]-xv[el[i][1]];
      dz[1]=xv[el[i][0]]-xv[el[i][2]];
      dz[2]=xv[el[i][1]]-xv[el[i][0]];
      dt=abs(xv[el[i][1]]*zv[el[i][2]]+xv[el[i][0]]*zv[el[i][1]]+xv[el[i][2]]*zv[el[i][0]]-xv[el[i][1]]*zv[el[i][0]]-xv[el[i][2]]*zv[el[i][1]]-xv[el[i][0]]*zv[el[i][2]])/2;//delta
      //printf("%f\n",dt);
      a1=a(xv[el[i][0]],zv[el[i][0]]);
      a2=a(xv[el[i][1]],zv[el[i][1]]);
      a3=a(xv[el[i][2]],zv[el[i][2]]);

      b1=b(xv[el[i][0]],zv[el[i][0]]);
      b2=b(xv[el[i][1]],zv[el[i][1]]);
      b3=b(xv[el[i][2]],zv[el[i][2]]);

      c1=c(xv[el[i][0]],zv[el[i][0]]);
      c2=c(xv[el[i][1]],zv[el[i][1]]);
      c3=c(xv[el[i][2]],zv[el[i][2]]);

      d1=d(xv[el[i][0]],zv[el[i][0]]);
      d2=d(xv[el[i][1]],zv[el[i][1]]);
      d3=d(xv[el[i][2]],zv[el[i][2]]);

      e1=e(xv[el[i][0]],zv[el[i][0]]);
      e2=e(xv[el[i][1]],zv[el[i][1]]);
      e3=e(xv[el[i][2]],zv[el[i][2]]);

      f1=f(Z,Yj[el[i][0]]);
      f2=f(Z,Yj[el[i][1]]);
      f3=f(Z,Yj[el[i][2]]);

      g1=g(Z,Yj[el[i][0]],Y[el[i][0]],Hp[el[i][0]]);
      g2=g(Z,Yj[el[i][1]],Y[el[i][1]],Hp[el[i][1]]);
      g3=g(Z,Yj[el[i][2]],Y[el[i][2]],Hp[el[i][2]]);

      // g1=(-Z*(Yj[el[i][0]]-Y[el[i][0]])*Hp[el[i][0]]);
      // g2=(-Z*(Yj[el[i][1]]-Y[el[i][1]])*Hp[el[i][1]]);
      // g3=(-Z*(Yj[el[i][2]]-Y[el[i][2]])*Hp[el[i][2]]);


      //matriz de cada elemento
      k[0][0]=(-1/(4*dt))*dx[0]*dx[0]*a1
      + (-1/(4*dt))*dx[0]*dz[0]*b1
      + (-1/(4*dt))*dz[0]*dz[0]*c1
      + (1/6)*dx[0]*d1
      + (1/6)*dz[0]*e1
      + (dt/12)*2*f1;
      k[0][1]=(-1/(4*dt))*dx[0]*dx[1]*a2
      + (-1/(4*dt))*dx[0]*dz[1]*b2
      + (-1/(4*dt))*dz[0]*dz[1]*c2
      + (1/6)*dx[1]*d2
      + (1/6)*dz[1]*e2
      + (dt/12)*f2;
      k[0][2]=(-1/(4*dt))*dx[0]*dx[2]*a3
      + (-1/(4*dt))*dx[0]*dz[2]*b3
      + (-1/(4*dt))*dz[0]*dz[2]*c3
      + (1/6)*dx[2]*d3
      + (1/6)*dz[2]*e3
      + (dt/12)*f3;

      k[1][0]=(-1/(4*dt))*dx[1]*dx[0]*a1
      + (-1/(4*dt))*dx[1]*dz[0]*b1
      + (-1/(4*dt))*dz[1]*dz[0]*c1
      + (1/6)*dx[0]*d1
      + (1/6)*dz[0]*e1
      + (dt/12)*f1;
      k[1][1]=(-1/(4*dt))*dx[1]*dx[1]*a2
      + (-1/(4*dt))*dx[1]*dz[1]*b2
      + (-1/(4*dt))*dz[1]*dz[1]*c2
      + (1/6)*dx[1]*d2
      + (1/6)*dz[1]*e2
      + (dt/12)*2*f2;
      k[1][2]=(-1/(4*dt))*dx[1]*dx[2]*a3
      + (-1/(4*dt))*dx[1]*dz[2]*b3
      + (-1/(4*dt))*dz[1]*dz[2]*c3
      + (1/6)*dx[2]*d3
      + (1/6)*dz[2]*e3
      + (dt/12)*f3;

      k[2][0]=(-1/(4*dt))*dx[2]*dx[0]*a1
      + (-1/(4*dt))*dx[2]*dz[0]*b1
      + (-1/(4*dt))*dz[2]*dz[0]*c1
      + (1/6)*dx[0]*d1
      + (1/6)*dz[0]*e1
      + (dt/12)*f1;
      k[2][1]=(-1/(4*dt))*dx[2]*dx[1]*a2
      + (-1/(4*dt))*dx[2]*dz[1]*b2
      + (-1/(4*dt))*dz[2]*dz[1]*c2
      + (1/6)*dx[1]*d2
      + (1/6)*dz[1]*e2
      + (dt/12)*f2;
      k[2][2]=(-1/(4*dt))*dx[2]*dx[2]*a3
      + (-1/(4*dt))*dx[2]*dz[2]*b3
      + (-1/(4*dt))*dz[2]*dz[2]*c3
      + (1/6)*dx[2]*d3
      + (1/6)*dz[2]*e3
      + (dt/12)*2*f3;

      //vetor da direita em cada elemento
      m[0]=(dt/12)*(2.*g1+g2+g3);
      m[1]=(dt/12)*(g1+2.*g2+g3);
      m[2]=(dt/12)*(g1+g2+2.*g3);

      //Matriz global vetor da direita
      for (j1=0;j1<3;j1++){
        for (i1=0;i1<3;i1++){
          K[el[i][j1]][el[i][i1]]=K[el[i][j1]][el[i][i1]]
          +k[j1][i1];
        }
        M[el[i][j1]]=M[el[i][j1]]+m[j1];
      }
      printf("elemento %d\n",i);
    }
    free(Y);
    free(Yj);
    Kr=(std::complex<double>**)malloc((np-qt_front)*sizeof(std::complex<double>*));
    Mr=(std::complex<double>*)malloc((np-qt_front)*sizeof(std::complex<double>));
    //redução do sistema - Inclusão das condições de fronteira

    for (i=0;i<np-qt_front;i++)
    {
      Kr[i]=(std::complex<double>*)malloc((np-qt_front)*sizeof(std::complex<double>));
      for (j=0;j<np-qt_front;j++)
      {
        Kr[i][j]=0;
      }
      Mr[i]=0;
    }
    ic=-1;
    for (i=0;i<np;i++)
    {
      if (vc[i]!=1)
      {
        jc=-1;
        ic++;
        Mr[ic]=M[i];
        for (j=0;j<np;j++)
        {
          if (vc[j]!=1)
          {
            jc++;
            Kr[ic][jc]=K[i][j];
          } else
          {
            Mr[ic]=Mr[ic]-K[i][j]*Hs[j];
          }
        }
      }
    }
free(K);
free(M);
fileaux=fopen("lixo.dat","w");
for (ic=0;ic<np-qt_front;ic++)
{
  fprintf(fileaux, "%f %f\n",real(Mr[ic]),imag(Mr[ic]));
}
fclose(fileaux);
    //solução por eliminação gaussiana
    //1-triangularização

    for (i=0;i<(np-qt_front-1);i++)
    {
      for (j=i+1;j<np-qt_front;j++)
      {
        temp=Kr[j][i]/Kr[i][i];
        for (l=i;l<np-qt_front;l++)
        {
          Kr[j][l]-=temp*Kr[i][l];
        }
        Mr[j]-=temp*Mr[i];
      }
    }
    Ur=(std::complex<double>*)malloc((np-qt_front)*sizeof(std::complex<double>));
    //2-retrossubstituição
    Ur[np-qt_front-1]=Mr[np-qt_front-1]/Kr[np-qt_front-1][np-qt_front-1];
    for (i=(np-qt_front-1);i>=0;i--)
    {
      Ur[i]=Mr[i];
      for (j=i+1;j<np-qt_front;j++)
      {
        Ur[i]-=Kr[i][j]*Ur[j];
      }
      Ur[i]=Ur[i]/Kr[i][i];
    }
free(Kr);
free(Mr);
U=(std::complex<double>*)malloc((np)*sizeof(std::complex<double>));

    //Incorporando as condições de fronteira
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
    free(Ur);
    Htotal=(std::complex<double>*)malloc(np*sizeof(std::complex<double>));
    FILE *file;
    // file=fopen("H_sec.dat","w");
    file=fopen("b.dat","w");
    for (i=0;i<np;i++)
    {
       Htotal[i]=U[i]+Hp[i];
    }
    ic=0;
    for (i=0;i<nz;i++)
    {
      for (j=0;j<nx;j++)
      {
            fprintf(file,"%f %f %f %f\n",xv2[ic],zv2[ic],real(U[ic]),imag(U[ic]));
          ic++;
      }
      fprintf(file,"\n");
    }
    fclose(file);
    free(U);
    free(Hp);

    fclose(file);
    //volta para a malha reduzida, para retirar os erros das bordas
    Htotalc=(std::complex<double>*)malloc(np2*sizeof(std::complex<double>));
    sigv2=(double*)malloc(np2*sizeof(double));
    ic=0;
    aux=7000;
    for (i=0;i<np;i++)
    {
      if((zv[i]<=aux))
      {
        if((xv[i]>=x0) & (xv[i]<=xn))
        {
          xv2[ic]=xv[i];
          zv2[ic]=zv[i];
          sigv2[ic]=sig[i];
          Htotalc[ic]=Htotal[i];
          ic=ic+1;
        }
      }
    }
    free(xv);
    free(zv);
    free(Htotal);
    int nz2=ic/nx2;
    file=fopen("c.dat","w");
    ic=0;
    for (i=0;i<nz2;i++)
    {
      for (j=0;j<nx2;j++)
      {
        fprintf(file,"%f %f %f %f\n",xv2[ic],zv2[ic],(sigv2[ic]),(sigv2[ic]));
        // fprintf(file,"%f %f %f %f\n",xv[ic],zv[ic],real(U[ic]),imag(U[ic]));
        ic++;
      }
      fprintf(file,"\n");
    }
    free(sigv2);
    //@@@
    // elc=(int**)malloc(nec*sizeof(int*));
    // for (i=0;i<nec;i++)
    // {
    //   elc[i]=(int*)malloc(4*sizeof(int));
    // }
    // cont=1;cc=0;ic=0;
    // for (i=0;i<nz2-1;i++)
    // {
    //   for (j=0;j<nx2;j++)
    //   {
    //     if (j<nx2)
    //     {
    //       cc=cc+1;
    //       elc[cc][0]=ic+1;
    //       elc[cc][1]=ic+2;
    //       elc[cc][2]=ic+nx2+1;
    //       elc[cc][3]=cc;
    //       cc=cc+1;
    //       elc[cc][0]=ic+2;
    //       elc[cc][1]=ic+nx2+2;
    //       elc[cc][2]=ic+nx2+1;
    //       elc[cc][3]=cc;
    //     }else
    //     {
    //       cont=cont+1;
    //     }
    //     ic=ic+1;
    //   }
    // }
    // free(elc);
    file=fopen("H_total_.dat","w");
    // for (i=0;i<np2;i++)
    // {
    //   fprintf(file,"%f %f %f %f\n",xv2[i],zv2[i],real(Htotalc[i]),imag(Htotalc[i]));
    // }
    ic=0;
    for (i=0;i<nz2;i++)
    {
      for (j=0;j<nx2;j++)
      {
        fprintf(file, "%f %f %f %f\n",xv2[ic],zv2[ic],real(Htotalc[ic]),imag(Htotalc[ic]));
        ic++;
      }
      fprintf(file, "\n");
    }
    fclose(file);
    file=fopen("H_plot.gnu","w");
    fprintf(file, "reset; \n");
    fprintf(file, "set title 'Parte real do Campo total'; \n");
    fprintf(file, "set xlabel 'X(m)'; \n");
    fprintf(file, "set ylabel 'Z(m)'; \n");
    fprintf(file, "set zlabel 'H'; \n");
    fprintf(file, "set contour; \n");
    fprintf(file, "splot 'H_total_.dat' using 1:2:3 with pm3d;\n");
    fprintf(file, "pause -1");
    fclose(file);
    system("gnuplot H_plot.gnu");
    // file=fopen("H_total.dat","w");
    // std::complex<double>Mtotal[nx2][nz2];
    // ic=0;
    // for (i=0;i<nx2;i++)
    // {
    //   for (j=0;j<nz;j++)
    //   {
    //     Mtotal[i][j]=Htotalc[ic];
    //     fprintf(file, "%f %f \t",real(Mtotal[i][j]),imag(Mtotal[i][j]));
    //     ic++;
    //   }
    //   fprintf(file, "\n");
    // }
    // fclose(file);
    return 0;
  }
