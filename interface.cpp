/****************************************************************************
Programa de simulação computacional utilizando a tecnica de elementos finitos
Autor: Yuri Ferreira dos Santos
****************************************************************************/
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
#define CLEAR 0
/** Pointers to the windows and some of the controls we'll create **/
GLUI *glui, *glui_prim,*glui_sec,*glui_alert;
GLUI_Spinner    *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel,*medium_panel_prim,*medium_panel_sec;
GLUI_EditText *condutivity;
int main_window;

char texto[30];
GLfloat win, r, g, b, xf,zf;
GLfloat coord[10][2];
float tw1,th1;
GLint view_w, view_h, primitiva;
int nx=50,nz=40,np;
float x0=0;
float xn=100;
float z0=0;
float zn=80;
float sigma,sigma_aux,sigma_sec[30];
double dx,dz;
double *x,*z,*xv,*zv,*xp,*zp;
double xvp,zvp;//coordenada do vertice mais proximo ao clique do usuario
double dist,aux_dist;
int *click_xv, *click_zv;
int cont,ic,i,j,qt_prim,qt_sec;
// float medium_prim[10][2][2];
float medium_sec[30][2][2];
int sec,prim,click=0;
int menu;
int tx, ty, tw, th;
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
    // if (prim==1)
    // {
    //   qt_prim++;
    //   glui_prim->hide();
    // }else
    if (sec==1) {
      sigma_sec[qt_sec]=sigma_aux;
      qt_sec++;
      glui_sec->hide();
    }
  }
  if (control==CLEAR)
  {
    qt_prim=0;
    qt_sec=0;
    prim=0;
    sec=0;
    click=0;
    glutPostRedisplay();
    for (i=0;i<np;i++)
    {
      click_xv[i]=0;
      click_zv[i]=0;
    }
  }
}

void DesenhaPontos(void)
{
  np=nx*nz;
  dx=abs(xn-x0)/(nx-1);
  dz=abs(zn-z0)/(nz-1);
  x=(double*)malloc(nx*sizeof(double));
  z=(double*)malloc(nz*sizeof(double));
  xv=(double*)malloc(np*sizeof(double));
  zv=(double*)malloc(np*sizeof(double));
  xp=(double*)malloc(np*sizeof(double));
  zp=(double*)malloc(np*sizeof(double));
  click_xv=(int*)malloc(np*sizeof(int));
  click_zv=(int*)malloc(np*sizeof(int));
  i=0;
  for (cont=0;cont<nx;cont++)
  {
    x[cont]=cont*dx;
    for (ic=0;ic<nz;ic++)
    {
      z[ic]=ic*dz;
      xv[i]=x[cont];
      zv[i]=z[ic];
      i++;
    }
  }
  tw1=tw-xn;
  th1=th-zn;
  for (i=0;i<np;i++)
  {
    xp[i]=xv[i]*tw1/xn;
    zp[i]=zv[i]*th1/zn;
  }
  glPointSize(3);
  glBegin(GL_POINTS);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
      glColor3f(0,0,0);
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
      glVertex2f((xp[cont+1])-tw1/2,zp[ic]-th1/2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
      glVertex2f((xp[cont])-tw1/2,zp[ic]-th1/2);
      glVertex2f((xp[cont])-tw1/2,zp[ic+1]-th1/2);
    }
  }
  glEnd();
}
// Desenha um texto na janela GLUT
void DesenhaTexto(char *string, int x, int z)
{
  glPushMatrix();
  // Posição no universo onde o texto será colocado
  glRasterPos2f(x,z);
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
  for (i=0;i<np;i++)
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
      //     for (i=0;i<np;i++)
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
          for (i=0;i<np;i++)
          {
            dist=sqrt(pow(x_m-xp[i]-xn/2,2)+pow(z_m-zp[i]-zn/2,2));
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
/**************************************** main() ********************/

int main(int argc, char* argv[])
{
  setlocale(LC_ALL, "");
  qt_prim=0;
  qt_sec=0;
  /****************************************/
  /*   Initialize GLUT and create window  */
  /****************************************/

  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 1100, 700 );
  main_window = glutCreateWindow( "Elementos Finitos - ondas planas" );
  GLUI_Master.set_glutReshapeFunc( AlteraTamanhoJanela);
  GLUI_Master.set_glutMouseFunc( GerenciaMouse );
  glutDisplayFunc( Desenha );
  //GLUI_Master.set_glutSpecialFunc( TeclasEspeciais);
  // glutPassiveMotionFunc(MoveMouse);
  printf("GLUI version: %3.2f\n", GLUI_Master.get_version());

  /*** Create the side subwindow ***/
  glui = GLUI_Master.create_glui_subwindow( main_window,GLUI_SUBWINDOW_RIGHT );
  obj_panel = new GLUI_Rollout(glui, "Propriedades da malha", true );
  GLUI_Spinner *x0_spinner = new GLUI_Spinner(obj_panel,"xi",&x0);
  x0_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  x0_spinner->set_float_limits(0,100);
  GLUI_Spinner *xn_spinner = new GLUI_Spinner(obj_panel,"xn",&xn);
  xn_spinner->set_float_limits(0,100);
  xn_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  GLUI_Spinner *z0_spinner = new GLUI_Spinner(obj_panel,"zi",&z0);
  z0_spinner->set_float_limits(0,100);
  z0_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  GLUI_Spinner *zn_spinner = new GLUI_Spinner(obj_panel,"zn",&zn);
  zn_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  zn_spinner->set_float_limits(0,100);
  GLUI_Spinner *nx_spinner = new GLUI_Spinner(obj_panel,"nx",&nx);
  nx_spinner->set_int_limits(10,50);
  nx_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  GLUI_Spinner *nz_spinner = new GLUI_Spinner(obj_panel,"nz",&nz);
  nz_spinner->set_alignment( GLUI_ALIGN_RIGHT );
  nz_spinner->set_float_limits(10,50);
  new GLUI_StaticText( glui, "" );
  medium_panel_prim = new GLUI_Panel(glui,"Meio primario");
  condutivity = new GLUI_EditText(medium_panel_prim,"Condutividade (siemens/m): ",&sigma);
  new GLUI_StaticText( glui, "" );
  medium_panel_sec =  new GLUI_Panel(glui,"Meio Secundario");
  new GLUI_Button(medium_panel_sec,"Adicionar",ENABLED_SEC,control_cb);
  new GLUI_Button(medium_panel_sec,"Limpar Malha",CLEAR,control_cb);
  new GLUI_StaticText(glui,"");
  /****** A 'quit' button *****/
  new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );
  new GLUI_StaticText(glui,"");
  new GLUI_StaticText(glui,"Instrução:");
  new GLUI_StaticText(glui,"");
  new GLUI_StaticText(glui,"Para desenhar o meio secundario, clique");
  new GLUI_StaticText(glui,"no no correspondente ao canto inferior");
  new GLUI_StaticText(glui,"esquerdo e apos no no´ correspondente");
  new GLUI_StaticText(glui,"ao canto superior  direito do retangulo");
  new GLUI_StaticText(glui,"desejado ");

  // glui_prim = GLUI_Master.create_glui("Meio primario");
  // condutivity = new GLUI_EditText(glui_prim,"Condutividade: ",&sigma);
  // GLUI_Panel *panel1 = new GLUI_Panel(glui_prim,"",false);
  // new GLUI_Button(panel1,"Confirmar",CONFIRMED,control_cb);
  // new GLUI_Column(panel1,false);
  // new GLUI_Button(panel1,"Cancelar",DISABLED,control_cb);
  // glui_prim->hide();
  glui_sec = GLUI_Master.create_glui("Meio secundario");
  condutivity = new GLUI_EditText(glui_sec,"Condutividade(siemens/m): ",&sigma_aux);
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
