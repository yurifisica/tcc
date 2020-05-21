#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <cmath>
#include <math.h>
#include <unistd.h>
#include <math.h>
#include <locale>
#ifdef __APPLE__
#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif
#endif

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

#include <string.h>
#include <GL/glui.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define ENABLED_PRIM 101
#define ENABLED_SEC 201
#define DISABLED 100
#define CONFIRMED 111
#define CALCULATE 222
#define CLEAR 0
#define length(a) sizeof(a)/sizeof(a[0])

//variaveis interface
GLUI *glui, *glui_prim,*glui_sec,*glui_alert;
GLUI_Spinner    *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel,*medium_panel_prim,*medium_panel_sec;
GLUI_EditText *condutivity_prim,*condutivity_sec;
int main_window;

char texto[30];
GLfloat win, xf,zf;
GLfloat coord[10][2];
float tw1,th1;
GLint view_w, view_h, primitiva;
int nx2=20,np2;
float sigma[1];
float sigma_aux;
float sigma_sec[30];
double dx1,dz1;
double *xt,*zt,*xp,*zp;
double xvp,zvp;//coordenada do vertice mais proximo ao clique do usuario
double dist,aux_dist;
int *click_xv, *click_zv;
int qt_prim,qt_sec;
// float medium_prim[10][2][2];
float medium_sec[30][2][2];
int sec,prim,click=0;
int menu;
int tx, ty, tw, th;

//variaveis
const float pi = acos(-1.);
const float mu=4*3.1415926E-7;
const float eps= 8.85E-12;
std::complex<double> *Hp,*M, **Kr, *Mr,*Y,*Yj,**K,k[3][3],m[3],*Hs,*Htotalc,*Htotal,*Ur,temp,*U;
double *xv,*zv,*xv2,*zv2,*sigv2;
float x0=0;
float xn=10000;//largura da malha
int nx=30;//numero de pontos em x, na malha reduzida
int n_cam=1,cont,j,i,aux,cc,ic;//numero de camadas
float frq=5E3;//frequência
float z0=0;
float zn=800;
int nz=30;//numero de pontos em z
int **elc,**el,*vc,qt_front;


/*---------------rotina de calculo de campo primario------------------*/
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


/*-------------rotina de elementos finitos-----------------*/
void elementos_finitos()
{
  double w=2*pi*frq;
  double espessura[n_cam-1]={zn/2};//espessura das camadas intermediárias
  int n_h=length(espessura);
  double sigma[n_cam];
  sigma[0]=pow(10,-3);//sigma[1]=0.5;
  int n_s=length(sigma);
  // double zn=espessura[0]*2;//profundidade da malha
  double dz1=abs(zn-z0)/(nz-1);
  double z[nz];
  // FILE *fp = fopen("asd.dat","w");
  for (cont=0;cont<nz;cont++)
  {
    z[cont]=cont*dz1;
  }
  // nz=length(z);
  double x[nx2+4];
  x[0]=-xn+xn/4;x[1]=-xn+xn/2;
  double dx1=abs(xn-x0)/(nx2-1);
  for (cont=0;cont<nx2;cont++)
  {
    x[cont+2]=cont*dx1;
  }
  x[nx2+2]=xn+xn/4;x[nx2+3]=xn+xn/2;
  nx=length(x);//numero de pontos em x, na malha completa
  // int nx2=length(xt);
  int np=nx*nz;
  int ne=(nx-1)*(nz-1)*2;
  double sig[np];
  std::complex<double> Z={0,w*mu};
  double dt;
  double dx[3],dz[3];
  int i1,j1,l;
  int jc=0;
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
  ic=0;
  // FILE *fp=fopen("asd.dat","w");
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
  // fclose(fp);

  campo_prim(n_s,sigma,n_h,espessura,zv,np);
  // xv2=(double*)malloc(np2*sizeof(double));
  // zv2=(double*)malloc(np2*sizeof(double));
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
  float aux_menor[2],aux_maior[2];
for (i=0;i<qt_sec;i++)
{
  aux_menor[0]=medium_sec[i][1][0];
  aux_maior[0]=medium_sec[i][0][0];
  aux_menor[1]=medium_sec[i][1][1];
  aux_maior[1]=medium_sec[i][0][1];
  if (medium_sec[i][0][0]<medium_sec[i][1][0])
  {
    aux_menor[0]=medium_sec[i][0][0];
    aux_maior[0]=medium_sec[i][1][0];
  }
  if (medium_sec[i][0][1]<medium_sec[i][1][1])
  {
    aux_menor[1]=medium_sec[i][0][1];
    aux_maior[1]=medium_sec[i][1][1];
  }
  // printf("maior %f %f\n", aux_maior[0]*xn/tw1,aux_maior[1]*zn/th1);
  // printf("menor %f %f\n", aux_menor[0]*xn/tw1,aux_menor[1]*zn/th1);

  for (j=0;j<np;j++)
  {
    if (xv[j]>=aux_menor[0]*xn/tw1 && xv[j]<=aux_maior[0]*xn/tw1)
    {
      if (zv[j]>=aux_menor[1]*zn/th1 && zv[j]<=aux_maior[1]*zn/th1)
      {
        sig[j]=sigma_sec[i];
        Yj[j]={sigma_sec[i],w*eps};
      }
    }
  }
}
FILE *fp =fopen("asd.dat","w");
for (i=0;i<np;i++)
{
  fprintf(fp,"%f %f %f\n",xv[i],zv[i],sig[i]);
}
fclose(fp);


  //meio secundário 2 (mar)
  // for (i=0;i<np;i++)
  // {
  //   if (xv[i]>=0 && xv[i]<=3000)
  //   {
  //     if ((zv[i]>=zn/2) && (zv[i]<(5*zn/2)/4))
  //     {
  //       Yj[i]={sigma2,w*eps}; sig[i]=sigma2;
  //       // printf("\n%i %f",i,real(Yj[i]-Y[i]));
  //     }
  //   }
  // }
  // //meio secundário 3 (prédios)
  // for (i=0;i<np;i++)
  // {
  //   // printf("\n%i %f",i,real(Yj[i]-Y[i]));
  //   if (xv[i]>=6500 && xv[i]<=7000 && zv[i]>=150 && zv[i]<zn/2) {
  //     Yj[i]={sigma3,w*eps};
  //     sig[i]=sigma3;
  //   }
  //   if (xv[i]>=7500 && xv[i]<=8000 && zv[i]>=150 && zv[i]<zn/2) {
  //     Yj[i]={sigma3,w*eps};
  //     sig[i]=sigma3;
  //   }
  //   if (xv[i]>=8500 && xv[i]<=9500 && zv[i]>=150   && zv[i]<espessura[0]) {
  //     Yj[i]={sigma3,w*eps};
  //     sig[i]=sigma3;
  //   }
  // }
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
     Htotal[i]=U[i]+Hp[i];  }
  ic=0;
  for (i=0;i<nz;i++)
  {
    for (j=0;j<nx;j++)
    {
          fprintf(file,"%f %f %f %f\n",xv[ic],zv[ic],real(Htotal[ic]),imag(Htotal[ic]));
        ic++;
    }
    fprintf(file,"\n");
  }
  fclose(file);
  free(U);
  free(Hp);
  //volta para a malha reduzida, para retirar os erros das bordas
  Htotalc=(std::complex<double>*)malloc(np2*sizeof(std::complex<double>));
  sigv2=(double*)malloc(np2*sizeof(double));
  ic=0;
  aux=7000;
  for (i=0;i<np;i++)
  {
    if((zv[i]<=aux))
    {
      if((xv[i]>=x0) & (xv[i]<=xn+1))
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
  fclose(file);
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
  fprintf(file, "set view 40,130,1,1; \n");
  fprintf(file, "set multiplot; \n");
  fprintf(file, "set size 1,0.5; \n");
  fprintf(file, "set origin 0,0.5; \n");
  fprintf(file, "set title 'Parte real do Campo total'; \n");
  fprintf(file, "set xlabel 'X'; \n");
  fprintf(file, "set ylabel 'Z'; \n");
  fprintf(file, "set zlabel 'H'; \n");
  fprintf(file, "set xrange [%f:%f] reverse; \n",xn,x0);
  fprintf(file, "set contour; \n");
  fprintf(file, "splot 'H_total_.dat' using 1:2:3 with pm3d;\n");
  fprintf(file, "set size 1,0.6; \n");
  fprintf(file, "set origin 0,0; \n");
  fprintf(file, "set title 'Parte imaginaria do Campo total'; \n");
  fprintf(file, "set xlabel 'X'; \n");
  fprintf(file, "set ylabel 'Z'; \n");
  fprintf(file, "set zlabel 'H'; \n");
  fprintf(file, "set xrange [%f:%f] reverse; \n",xn,x0);
  fprintf(file, "set contour; \n");
  fprintf(file, "splot 'H_total_.dat' using 1:2:4 with pm3d;\n");
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
  return;
}


/*------------rotinas graficas-------------------*/
void myGlutIdle()
{
  /* According to the GLUT specification, the current window is
  undefined during an idle callback.  So we need to explicitly change
  it if necessary */
  if ( glutGetWindow() != main_window )
  glutSetWindow(main_window);

  /*  GLUI_Master.sync_live_all();  -- not needed - nothing to sync in this
  application  */

  glutPostRedisplay();
}

void control_cb(int control)
{
  // if (control==ENABLED_PRIM)
  // {
  //   prim=1;
  //   sec=0;
  //   click=0;
  //   glui_prim->show();
  // }
  if (control==ENABLED_SEC)
  {
    //prim=0;
    sec=1;
    click=0;
    glui_sec->show();
    // glui_sec->enable();
  }
  if (control==DISABLED)
  {
    prim=0;
    sec=0;
    click=0;
    glui_prim->hide();
    glui_sec->hide();
  }
  if (control==CONFIRMED)
  {
    if (prim==1)
    {
      qt_prim++;
      glui_prim->hide();
    }else
    if (sec==1) {
      sigma_sec[qt_sec]=sigma_aux;
      qt_sec++;
      glui_sec->hide();
    }
    obj_panel->disable();


  }
  if (control==CALCULATE)
  {
    elementos_finitos();
  }
  if (control==CLEAR)
  {
    qt_prim=0;
    qt_sec=0;
    prim=0;
    sec=0;
    click=0;
    obj_panel->enable();
    for (i=0;i<np2;i++)
    {
      click_xv[i]=0;
      click_zv[i]=0;
    }
    glutPostRedisplay();
  }
}

void DesenhaPontos(void)
{
  np2=nx2*nz;
  dx1=abs(xn-x0)/(nx2-1);
  dz1=abs(zn-z0)/(nz-1);
  xt=(double*)malloc(nx2*sizeof(double));
  zt=(double*)malloc(nz*sizeof(double));
  xv2=(double*)malloc(np2*sizeof(double));
  zv2=(double*)malloc(np2*sizeof(double));
  xp=(double*)malloc(np2*sizeof(double));
  zp=(double*)malloc(np2*sizeof(double));
  click_xv=(int*)malloc(np2*sizeof(int));
  click_zv=(int*)malloc(np2*sizeof(int));
  i=0;
  for (cont=0;cont<nx2;cont++)
  {
    xt[cont]=cont*dx1;
    for (ic=0;ic<nz;ic++)
    {
      zt[ic]=ic*dz1;
      xv2[i]=xt[cont];
      zv2[i]=zt[ic];
      i++;
    }
  }
  tw1=tw-20;
  th1=th-20;
  for (i=0;i<np2;i++)
  {
    xp[i]=xv2[i]*tw1/xn;
    zp[i]=zv2[i]*th1/zn;
  }
  glPointSize(3);
  glBegin(GL_POINTS);
  for (cont=0;cont<np2;cont++)
  {
    for (ic=0;ic<np2;ic++)
    {
      glColor3f(0,0,0);
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np2;cont++)
  {
    for (ic=0;ic<np2;ic++)
    {
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
      glVertex2f((xp[cont+1])-tw1/2,zp[ic]-th1/2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np2;cont++)
  {
    for (ic=0;ic<np2;ic++)
    {
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
      glVertex2f((xp[cont])-tw1/2,zp[ic+1]-th1/2);
    }
  }
  glEnd();
}
// Desenha um texto na janela GLUT
void DesenhaTexto(char *string, int x_m, int z_m)
{
  glPushMatrix();
  // Posição no universo onde o texto será colocado
  glRasterPos2f(x_m,z_m);
  glColor3f(1,0,0);
  // Exibe caracter a caracter
  while(*string)
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*string++);
  glPopMatrix();
}
//desenha Meios
void DesenhaMeio(void)
{
    if (qt_sec!=0)
    {
      // glColor3f(1,0,0);
      // glBegin(GL_QUADS);
      // for (i=0;i<qt_prim;i++)
      // {
      //   if ((i<qt_prim-1) | (click==2) | (sec==1))
      //   {
      //     glVertex2f(medium_prim[i][0][0]-tw1/2,-medium_prim[i][0][1]+th1/2);
      //     glVertex2f(medium_prim[i][1][0]-tw1/2,-medium_prim[i][0][1]+th1/2);
      //     glVertex2f(medium_prim[i][1][0]-tw1/2,-medium_prim[i][1][1]+th1/2);
      //     glVertex2f(medium_prim[i][0][0]-tw1/2,-medium_prim[i][1][1]+th1/2);
      //   }
      // }
      // glEnd();
      // printf("%d\n",qt_sec );
      for (i=0;i<qt_sec;i++)
      {
        if ((i<qt_sec-1) | (click==2))
        {
        glColor3f(0,1,0);
        glBegin(GL_QUADS);
          glVertex2f(medium_sec[i][0][0]-tw1/2,-medium_sec[i][0][1]+th1/2);
          glVertex2f(medium_sec[i][1][0]-tw1/2,-medium_sec[i][0][1]+th1/2);
          glVertex2f(medium_sec[i][1][0]-tw1/2,-medium_sec[i][1][1]+th1/2);
          glVertex2f(medium_sec[i][0][0]-tw1/2,-medium_sec[i][1][1]+th1/2);
        glEnd();
          sprintf(texto,"%.4f",sigma_sec[i]);
          // printf("%d -> %f\n",i,sigma_sec[i]);
          glColor3f(1,0,0);
          DesenhaTexto(texto,(medium_sec[i][1][0]+medium_sec[i][0][0])/2-tw1/2,(-medium_sec[i][1][1]-medium_sec[i][0][1])/2+th1/2);
        }
      }
      }
    // glutPostRedisplay();
  }
// Função callback chamada para fazer o desenho
void Desenha(void)
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT);
  DesenhaMeio();
  DesenhaPontos();
  // DesenhaTexto(texto,0,0);
  // Exibe a posição do mouse na janela
  glColor3f(1.0f,1.0f,1.0f);
  glutSwapBuffers();
  // glutPostRedisplay();
}
void Inicializa (void)
{
  // Define a cor de fundo da janela de visualização como branca
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  for (i=0;i<np2;i++)
  {
    click_xv[i]=0;
    click_zv[i]=0;
  }
  win=50;
}

// Função callback chamada quando o tamanho da janela é alterado
void AlteraTamanhoJanela(GLsizei w, GLsizei h)
{
  // Especifica as dimensões da Viewport
  GLUI_Master.get_viewport_area ( &tx, &ty, &tw, &th );
  glViewport( tx, ty, tw, th );
  // printf("%i %i\n",tw, th );
  view_w = w;
  view_h = h;
  // Inicializa o sistema de coordenadas
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D (-tw/2, tw/2, -th/2, th/2);
}

void GerenciaMouse(int button, int state, int x_m, int z_m)
{
  if (button==GLUT_LEFT_BUTTON)
  {
    if (state==GLUT_DOWN)
    {
      // if (prim==1)
      // {
      //   if (click<2)
      //   {
      //     aux_dist=10000;
      //     for (i=0;i<np2;i++)
      //     {
      //       dist=sqrt(pow(x_m-xp[i]-xn/2,2)+pow(z_m-zp[i]-zn/2,2));
      //       if (dist<aux_dist)
      //       {
      //         xvp=xp[i];
      //         zvp=zp[i];
      //         click_xv[i]=1;
      //         click_zv[i]=1;
      //         aux_dist=dist;
      //       }
      //     }
      //     medium_prim[qt_prim-1][click][0]=xvp;
      //     medium_prim[qt_prim-1][click][1]=zvp;
      //     click++;
      //   }
      // }
      if (sec==1)
      {
        if (click<2)
        {
          aux_dist=10000;
          for (i=0;i<np2;i++)
          {
            dist=sqrt(pow(x_m-xp[i]-10,2)+pow(z_m-zp[i]-10,2));
            if (dist<aux_dist)
            {
              xvp=xp[i];
              zvp=zp[i];
              click_xv[i]=1;
              click_zv[i]=1;
              aux_dist=dist;
            }
          }
          medium_sec[qt_sec-1][click][0]=xvp;
          medium_sec[qt_sec-1][click][1]=zvp;
          click++;
        }
      }
    }
  }
  glutPostRedisplay();
}
void MoveMouse(int x_m, int z_m)
{
  sprintf(texto, "(%d , %d)", x_m, z_m);
  glutPostRedisplay();
}



  /*------------programa principal---------------------*/
  int main(int argc, char *argv[])
  {
    // setlocale(LC_ALL, "Portuguese");
    qt_prim=0;
    qt_sec=0;
    /****************************************/
    /*   Initialize GLUT and create window  */
    /****************************************/

    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
    glutInitWindowPosition( 50, 50 );
    glutInitWindowSize( 1100, 700 );
    main_window = glutCreateWindow( "Elementos Finitos - Ondas planas eletromagneticas" );
    GLUI_Master.set_glutReshapeFunc( AlteraTamanhoJanela);
    GLUI_Master.set_glutMouseFunc( GerenciaMouse );
    glutDisplayFunc( Desenha );
    //GLUI_Master.set_glutSpecialFunc( TeclasEspeciais);
    // glutPassiveMotionFunc(MoveMouse);

    printf("Simulação da propagação de ondas eletromagnéticas planas");
    /*** Create the side subwindow ***/
    glui = GLUI_Master.create_glui_subwindow( main_window,GLUI_SUBWINDOW_RIGHT );
    obj_panel = new GLUI_Rollout(glui, "Propriedades da malha", true );
    GLUI_Spinner *x0_spinner = new GLUI_Spinner(obj_panel,"xi",&x0);
    x0_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    x0_spinner->set_float_limits(0,100);
    GLUI_Spinner *xn_spinner = new GLUI_Spinner(obj_panel,"xn",&xn);
    xn_spinner->set_float_limits(0,10000);
    GLUI_Spinner *nx_spinner = new GLUI_Spinner(obj_panel,"nx",&nx2);
    nx_spinner->set_int_limits(10,50);
    nx_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    xn_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    GLUI_Spinner *z0_spinner = new GLUI_Spinner(obj_panel,"zi",&z0);
    z0_spinner->set_float_limits(0,100);
    z0_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    GLUI_Spinner *zn_spinner = new GLUI_Spinner(obj_panel,"zn",&zn);
    zn_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    zn_spinner->set_float_limits(0,10000);
    GLUI_Spinner *nz_spinner = new GLUI_Spinner(obj_panel,"nz",&nz);
    nz_spinner->set_alignment( GLUI_ALIGN_RIGHT );
    nz_spinner->set_int_limits(10,50);
    new GLUI_StaticText( glui, "" );
    GLUI_EditText *freq = new GLUI_EditText(glui,"Frequencia da onda(Hz):",&frq);
    medium_panel_prim = new GLUI_Panel(glui,"Meio primario");
    condutivity_prim = new GLUI_EditText(medium_panel_prim,"Condutividade (S/m)",&sigma[0]);
    condutivity_prim->set_float_limits(0.000000001,100000000);
    new GLUI_StaticText( glui, "" );
    medium_panel_sec =  new GLUI_Panel(glui,"Meio Secundario");
    new GLUI_Button(medium_panel_sec,"Adicionar",ENABLED_SEC,control_cb);
    new GLUI_Button(medium_panel_sec,"Limpar Malha",CLEAR,control_cb);
    new GLUI_StaticText(glui,"");
    GLUI_Button *calcular = new GLUI_Button( glui, "Simular",CALCULATE,control_cb);

    new GLUI_StaticText(glui,"");
    /****** A 'quit' button *****/
    new GLUI_Button( glui, "Sair", 0,(GLUI_Update_CB)exit );
    new GLUI_StaticText(glui,"");
    new GLUI_StaticText(glui,"Instruçao:");
    new GLUI_StaticText(glui,"");
    new GLUI_StaticText(glui,"Para desenhar o meio secundario, clique");
    new GLUI_StaticText(glui,"nos nos correspondentes a as duas");
    new GLUI_StaticText(glui,"extremidades do retangulo desejado");

    glui_prim = GLUI_Master.create_glui("Meio primario");
    condutivity_prim = new GLUI_EditText(glui_prim,"Condutividade (S/m): ",&sigma[0]);
    GLUI_Panel *panel1 = new GLUI_Panel(glui_prim,"",false);
    new GLUI_Button(panel1,"Confirmar",CONFIRMED,control_cb);
    new GLUI_Column(panel1,false);
    new GLUI_Button(panel1,"Cancelar",DISABLED,control_cb);
    glui_prim->hide();
    glui_sec = GLUI_Master.create_glui("Meio secundario");
    condutivity_sec = new GLUI_EditText(glui_sec,"Condutividade (S/m)",&sigma_aux);
    GLUI_Panel *panel2 = new GLUI_Panel(glui_sec,"",false);
    new GLUI_Button(panel2,"Confirmar",CONFIRMED,control_cb);
    new GLUI_Column(panel2,false);
    new GLUI_Button(panel2,"Cancelar",DISABLED,control_cb);
    glui_sec->hide();

    /**** Link windows to GLUI, and register idle callback ******/

    glui->set_main_gfx_window( main_window );
    #if 0
    /**** We register the idle callback with GLUI, *not* with GLUT ****/
    GLUI_Master.set_glutIdleFunc( myGlutIdle );
    #endif
    /**** Regular GLUT main loop ****/
    Inicializa();
    glutMainLoop();

    return EXIT_SUCCESS;
  }
