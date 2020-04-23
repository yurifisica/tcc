/****************************************************************************

  example5.cpp

  A GLUI program demonstrating subwindows, rotation controls,
	translation controls, and listboxes

  -----------------------------------------------------------------------

  7/10/98 Paul Rademacher (rademach@cs.unc.edu)

****************************************************************************/
#include <math.h>
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

/** Pointers to the windows and some of the controls we'll create **/
GLUI *glui, *glui2;
GLUI_Spinner    *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel;
int main_window;

char texto[30];
GLfloat win, r, g, b, xf,zf;
GLfloat coord[10][2];
GLint view_w, view_h, primitiva;
int nx=50,nz=40,np;
float x0=0;
float xn=100;
float z0=0;
float zn=80;
double dx,dz;
double *x,*z,*xv,*zv;
double xvp,zvp;//coordenada do vertice mais proximo ao clique do usuario
double dist,aux_dist;
int *click_xv, *click_zv;
int cont,ic,i,j;
double *polygon[2];
int prim,click=0;
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
void DesenhaPontos(void)
{
  np=nx*nz;
  dx=(xn-x0)/(nx+1);
  dz=(zn-z0)/(nz+1);
  x=(double*)malloc(nx*sizeof(double));
  z=(double*)malloc(nz*sizeof(double));
  xv=(double*)malloc(np*sizeof(double));
  zv=(double*)malloc(np*sizeof(double));
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
  i=0;
  float tw1=tw-xn,th1=th-zn;
  float lixo=tw1/2;//(x[nx-1]-x[0])*tw1/200;
  float lixo2=th1/2;//(z[nz-1]-z[0])*tw1/160;
  glPointSize(3);
  glBegin(GL_POINTS);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
      if (click_xv[cont]==1 && click_zv[ic]==1)
      {
          glColor3f(1,0,0);
      }
      else
      {
          glColor3f(0,0,0);
      }
        glVertex2f((xv[cont])*tw1/x[nx-1]-lixo,zv[ic]*th1/z[nz-1]-lixo2);
        // printf("%f %f\n",xv[cont]-tw1/2,zv[ic]+th1/2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
        glVertex2f((xv[cont])*tw1/x[nx-1]-lixo,zv[ic]*th1/z[nz-1]-lixo2);
        glVertex2f((xv[cont+1])*tw1/x[nx-1]-lixo,zv[ic]*th1/z[nz-1]-lixo2);
    }
  }
  glEnd();
  glBegin(GL_LINES);
  for (cont=0;cont<np;cont++)
  {
    for (ic=0;ic<np;ic++)
    {
      glVertex2f((xv[cont])*tw1/x[nx-1]-lixo,zv[ic]*th1/z[nz-1]-lixo2);
      glVertex2f((xv[cont])*tw1/x[nx-1]-lixo,zv[ic+1]*th1/z[nz-1]-lixo2);
        // printf("%f %f\n",xv[cont]-tw/2,zv[ic]+th/2);
    }
  }
  glEnd();
}
// Desenha um texto na janela GLUT
void DesenhaTexto(char *string)
{
  glPushMatrix();
  // Posição no universo onde o texto será colocado
  glRasterPos2f(-win,win-2);
  glColor3f(1,0,0);
  // Exibe caracter a caracter
  while(*string)
  glutBitmapCharacter(GLUT_BITMAP_9_BY_15,*string++);
  glPopMatrix();
}
//desenha Meios
void DesenhaMeio(void)
{
  glColor3f(1,0,0);
  if (prim==1)
  {
    printf("%d\n",click );
      if (click>1)
      {
        for (ic=2;ic<=click;ic++)
        {
          printf("%f\n",polygon[0][ic-2]);
          printf("%f\n",polygon[0][ic-1]);
        }
        glBegin(GL_LINES);
        glVertex2f(polygon[0][click-2],polygon[1][click-2]);
        glVertex2f(polygon[0][click-1],polygon[1][click-1]);
        glEnd();
      }
      glBegin(GL_LINE_STRIP);
        for (cont=0;cont<click;cont++)
        {
          glVertex2f(polygon[0][cont],polygon[1][cont]);
        }
      glEnd();
      // if (polygon[0][click-1]==polygon[0][0] && polygon[1][click-1]==polygon[1][0])
      //   {
      //     glBegin(GL_POLYGON);
      //     for (cont=0;cont<click;cont++)
      //     {
      //       glVertex2f(polygon[0][cont],polygon[1][cont]);
      //     }
      //     glEnd();
      //     prim=0;
      //   }
  }
}

// Função callback chamada para fazer o desenho
void Desenha(void)
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT);
  DesenhaMeio();
  DesenhaPontos();
  DesenhaTexto(texto);
  // Exibe a posição do mouse na janela
  glColor3f(1.0f,1.0f,1.0f);
  glutSwapBuffers();
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
  printf("%i %i\n",tw, th );
  view_w = w;
  view_h = h;
  // Inicializa o sistema de coordenadas
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D (-tw/2, tw/2, -th/2, th/2);
}

// Função callback chamada para gerenciar eventos do teclado
// para teclas especiais, tais como F1, PgDn e Home
void TeclasEspeciais(int key, int x, int y)
{
  if(key == GLUT_KEY_UP) {
    win -= 10;
    if (win < 10) win = 10;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D (-win, win, -win, win);
  }
  if(key == GLUT_KEY_DOWN) {
    win += 10;
    if (win > 500) win = 500;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D (-win, win, -win, win);
  }
  glutPostRedisplay();
}

void GerenciaMouse(int button, int state, int x, int z)
{
if (button==GLUT_LEFT_BUTTON)
  {
    if (state==GLUT_DOWN)
    {
        if (click<2)
        {
          sprintf(texto,"Digite coordenada %i",click+1);
          coord[click][0]=x;
          coord[click][1]=z;
          printf("%f , %f\n",coord[click][0],coord[click][1]);
          aux_dist=10000;
          for (i=0;i<np;i++)
          {
            printf("%i-> %f , %f\n",i, xv[i], zv[i]);
            dist=sqrt(pow(coord[click][0]-xv[i],2)+pow(coord[click][1]-zv[i],2));
            if (dist<aux_dist)
              {
                xf=xv[i];
                zf=zv[i];
                aux_dist=dist;
              }
            }
          printf("%f , %f\n", xf, zf );
          printf("%f , %f\n",coord[click][0],coord[click][1]);
           // printf("%f %f\n",xn/2,zn/2 );
          click++;
        }
      // if (prim==1)
      // {
      //   xf = ( (2 * win * x) / view_w) - win;
      //   zf = ( ( (2 * win) * (z-view_h) ) / -view_h) - win;
      //       polygon[0][click]=xf;
      //       polygon[1][click]=zf;
      //       click++;
      // }
    }
  }

  glutPostRedisplay();
}
void MoveMouse(int x, int z)
{
        sprintf(texto, "(%d , %d)", x, z);
        glutPostRedisplay();
}
/**************************************** main() ********************/

int main(int argc, char* argv[])
{
  /****************************************/
  /*   Initialize GLUT and create window  */
  /****************************************/

  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 1000, 700 );
  main_window = glutCreateWindow( "GLUI Example 5" );
  GLUI_Master.set_glutReshapeFunc( AlteraTamanhoJanela);
  glutDisplayFunc( Desenha );
  //GLUI_Master.set_glutSpecialFunc( TeclasEspeciais);
  GLUI_Master.set_glutMouseFunc( GerenciaMouse );
  glutPassiveMotionFunc(MoveMouse);
  printf("GLUI version: %3.2f\n", GLUI_Master.get_version());

  /*** Create the side subwindow ***/
  glui = GLUI_Master.create_glui_subwindow( main_window,GLUI_SUBWINDOW_RIGHT );

  obj_panel = new GLUI_Rollout(glui, "Properties", true );
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

  /****** A 'quit' button *****/
  new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );


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
