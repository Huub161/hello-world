/*
 * Silo.c
 *
 * Multiple particle dynamics simulation
 *
 * (c) 2016 Joris Remmers TU/e
 * 
 * Edited by S.A. de Milliano
 * Edited by H. Wezenbeek
 */


#include <stdlib.h>
#include "mylib.h"

void plot

  ( char*           name  ,
    Plist*          pl     )

{
  FILE *of;
  int iPar;
  
  char colors[5][10] = {"red" , "green" , "blue" , "grey" , "black"};
  
  of=fopen(name,"w");

  fprintf(of,"<?xml version='1.0' standalone='no'?>\n");
  fprintf(of,"<!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.1//EN'\n");
  fprintf(of,"'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>\n");
  
  fprintf(of,"<svg width='600px' height='800px' version='1.1'\n");
  fprintf(of,"xmlns='http://www.w3.org/2000/svg'>\n");
  fprintf(of,"<g transform='translate(300,0)'>\n");
  fprintf(of,"<g transform='scale(1.5)'>\n");
  
  for ( iPar = 0 ; iPar < pl->ntot ; iPar++ )
  {
    fprintf(of,"<circle cx='%f' cy='%f' r='%f' fill='%s'/>\n",
            100*pl->p[iPar].r.x,
            500-100*pl->p[iPar].r.y,
            100*pl->p[iPar].radius,
            colors[pl->p[iPar].type]);
  }
  
  fprintf(of,"</g>\n</g>\n</svg>\n");
  fclose(of);
}

//----------------------------------------------------------
//
//----------------------------------------------------------

void reset_list

  ( LINKED_LIST *linked_list )
{
  for ( int i = 0; i < MAX_CELLS; i++ )
  {
	linked_list -> head[i] = -1;
  }
  for ( int i = 0; i < MAX_PARTICLES; i++ )
  {
	linked_list -> next[i] = -1;
  }
}

//----------------------------------------------------------
//
//----------------------------------------------------------

void fill_linkedlist

  ( LINKED_LIST     *linked_list ,
    Plist           *pl )

{
  int cell;
  int row;
  int colomn;

  Vec2 position;

  reset_list( linked_list );

  for( int iPar = 0; iPar < pl->ntot; iPar++)
  {
    position = pl -> p[iPar].r;		

    colomn = (position.x + 1.3) * NCOLUMNS/2.6;

    row = position.y * NROWS/5;
    cell = colomn + NCOLUMNS * row;

    if ( linked_list -> head[cell] == -1 )
    {
	  linked_list -> head[cell] = iPar;
    }
    else
    {
      int a = linked_list -> head[cell];
 
      while ( linked_list -> next[a] != -1 )
      {
        a = linked_list -> next[a];
      }

      linked_list -> next[a] = iPar;
    }

    pl -> p[iPar].cellnumber = cell;
  }
}



//----------------------------------------------------------
//
//----------------------------------------------------------


void read_input
  
  ( char*           name  ,
    Plist*          pl    )
  
{
  FILE         *fp;
  float        x,y,r;
  int          iPar;
  int          nFix,nGap,i;
  
  if ( ( fp=fopen(name,"r") ) == NULL) 
  {
    printf("Cannot open file.\n");
  }
  
  i = fscanf(fp, "%d %d", &nFix , &nGap );
  
  pl->nwall = nFix;
  pl->ndoor = nGap;

  pl->ntot  = nFix + nGap;
  
  for ( iPar = 0 ; iPar < pl->ntot ; iPar++ )
  {    
    i = fscanf(fp,"%e %e %e",&x,&y,&r);
    
    pl->p[iPar].r.x = x;
    pl->p[iPar].r.y = y;
    
    pl->p[iPar].radius = r;
    
    pl->p[iPar].mass = pl->p[iPar].radius*pl->p[iPar].radius;
   
    if ( iPar < nFix )
    {  
      pl->p[iPar].type = 0;
    }
    else
    {
      pl->p[iPar].type = 1;
    }
  }

  fclose( fp );
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void calc_interaction

  ( Plist         *pl          ,
    LINKED_LIST   *linked_list )
 
{
  int      iPar,jPar,i;
  int      pcp[MAX_PARTICLES_9CEL]; // ADDED possible contact particles
  int      amount_pcp=0;
 
  const int nPar = pl->ntot;
  
  for( iPar = 0 ; iPar < nPar ; iPar++ )
   {

amount_pcp = possible_contact   ( pl , linked_list , iPar , pcp);

    for ( i=0 ; i< amount_pcp; i++)
    {
// Als regels 195, 196 en 199 aan staan, lijkt de berekening meer op de originele, maar dit geeft een gekke melding op het eind. 
//wellicht zit hier ook het verschil in van de kinetic energy. (hij zou er ook sneller van moeten worden)
//			if (pcp[i]>iPar) 
//			{
               jPar = pcp[i];
               int_force( &pl->p[iPar] , &pl->p[jPar] );
//			}    
    }
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void int_force

  ( Particle *pi ,
    Particle *pj )

{
  Vec2   dr,f;
  double dist,eps,force;
  
  dr.x = pj->r.x - pi->r.x;
  dr.y = pj->r.y - pi->r.y;
      
  dist = sqrt( dr.x*dr.x + dr.y*dr.y );

  eps = pi->radius + pj->radius - dist;
   
  if ( eps > 0. )
  {
    force = eps * K_CONST / dist;
    
    f.x = -force*dr.x + B_CONST * ( pj->v.x - pi->v.x );
    f.y = -force*dr.y + B_CONST * ( pj->v.y - pi->v.y);
    
    pi->f.x += f.x;
    pi->f.y += f.y;
    
    pj->f.x += -f.x;
    pj->f.y += -f.y;
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void add_gravity

  ( Plist   * pl )

{
  int iPar;
  
  for ( iPar = 0 ; iPar < pl->  ntot ; iPar++ )
  {        
    pl->p[iPar].f.y += -GRAVITY*pl->p[iPar].mass;
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void add_particle

  ( Plist*         pl )

{
  double xpos;
  double rad;
  
  int iPar = pl->ntot;	

  init_particle( &pl->p[iPar] );

  xpos = -0.25+rand()%11*0.05;
      
  pl->p[iPar].r.x = xpos;
  pl->p[iPar].r.y = 4.0;
  
  pl->p[iPar].v.y = -2.0;
    
  if( iPar%2 == 0 )
  {
    rad = 0.04;
    pl->p[iPar].type = 2;
  }
  else
  {
    rad = 0.03;
    pl->p[iPar].type = 3;
  }
  
  pl->p[iPar].radius  = rad;
  
  pl->p[iPar].mass    = 3.1415*rad*rad;
       
  pl->ntot++; 

  if ( pl->ntot == MAX_PARTICLES )
  {
    printf("Too many particles\n");
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void init_particle

  ( Particle*      p )
  
{
  p->r.x = 0.;
  p->r.y = 0.;
    
  p->radius = 0.;

  p->v.x = 0.;
  p->v.y = 0.;
    
  p->a.x = 0.;
  p->a.y = 0.;

  p->f.x = 0.;
  p->f.y = 0.;
  
  p->mass = 0.;
   
  p->type = 0;
}


//----------------------------------------------------------
//
//----------------------------------------------------------


double solve

  ( LINKED_LIST *linked_list ,
	Plist       *pl          )

{
  const double dt2 = DT * DT;

  double ekin = 0.;
  
  int iPar;

  const int nPar = pl->ntot;
  
  for ( iPar = pl->nwall+pl->ndoor ; iPar < nPar ; iPar++ )
  {
    pl->p[iPar].r.x += DT * pl->p[iPar].v.x + 0.5*dt2*pl->p[iPar].a.x;
    pl->p[iPar].r.y += DT * pl->p[iPar].v.y + 0.5*dt2*pl->p[iPar].a.y;
    
    pl->p[iPar].v.x += 0.5*DT*pl->p[iPar].a.x;
    pl->p[iPar].v.y += 0.5*DT*pl->p[iPar].a.y;
    
    pl->p[iPar].f.x = 0.;
    pl->p[iPar].f.y = 0.;
  }


  fill_linkedlist( linked_list , pl );  

  calc_interaction( pl, linked_list );
  
  add_gravity( pl );
     
  for ( iPar= pl->nwall + pl->ndoor ; iPar < nPar ; iPar++ )
  { 
    pl->p[iPar].a.x = pl->p[iPar].f.x/pl->p[iPar].mass;
    pl->p[iPar].a.y = pl->p[iPar].f.y/pl->p[iPar].mass;
    
    pl->p[iPar].v.x += 0.5*DT*pl->p[iPar].a.x;
    pl->p[iPar].v.y += 0.5*DT*pl->p[iPar].a.y;
    
    ekin += pl->p[iPar].mass * ( pl->p[iPar].v.x*pl->p[iPar].v.x +
                                 pl->p[iPar].v.y*pl->p[iPar].v.y );
  }
  
  return 0.5*ekin;
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void check_particles
  
  ( Plist       *pl )

{
  int i;
  
  for ( i = pl->ndoor+pl->nwall ; i < pl->ntot ; i++ )
  {
    if ( pl->p[i].r.y < -1.0 )
    {
      remove_particle( pl , i );
    }
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void open_door
  
  ( Plist       *pl )

{
  int i,j;
  
  for ( i = 0 ; i < pl->ndoor ; i++ )
  {
    pl->p[pl->nwall+i]=pl->p[pl->ntot - pl->ndoor+i];
  }

  pl->ntot  = pl->ntot - pl->ndoor;
  pl->ndoor = 0;

  for ( i = pl->nwall ; i < pl->ntot ; i++ )
  {
    j = (int)(pl->p[i].r.y / 0.3 );

    if ( j % 2 == 0 )
    {
      pl->p[i].type = 2;
    }
    else
    {
      pl->p[i].type = 3;
    }
  }
}

//----------------------------------------------------------
//
//----------------------------------------------------------


void remove_particle

  ( Plist        *pl    ,
    int          iPar   )

{
  if ( iPar >= pl-> ntot  || iPar < pl->ndoor+pl->nwall )
  {
    printf("Error\n");
  }
 
  pl->ntot = pl->ntot-1;

  pl->p[iPar] = pl->p[pl->ntot];
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void get_filename

  ( char   *names ,
    int    k      )

{
  int m;
  char cc[4];
  
  strcpy(names,"output000.svg");
  
  m = sprintf(cc,"%d",k);
    
  if ( k < 10 )
  {
    names[8]=cc[0];
  }
  else  if ( k <100 )
  { 
    names[7]=cc[0];
    names[8]=cc[1];
  }
  else
  {
    names[6]=cc[0];
    names[7]=cc[1];
    names[8]=cc[2];
  }
}


//----------------------------------------------------------
//
//----------------------------------------------------------


void show_info

  ( char*     svgfile ,
    double    ekin    ,
    int       ntot    )

{
  printf("Writing file        : %s\n",svgfile);
  printf("Kinetic energy      : %f\n",ekin);
  printf("Number of particles : %d\n\n",ntot);
}


//----------------------------------------------------------
// find out which particles may possibly be in contact 
// return those in a array possible_contact_particles
//----------------------------------------------------------

int possible_contact 

  ( Plist       *pl          , 
    LINKED_LIST *linked_list ,
    int         partnr      , 
    int*        possible_contact_particles ) //Als je van deze waarde een pointer kan maken, kan dat de snelheid aanzienlijk tengoede komen.

{
  int i,j,k,n;
//  int inside_cell = pl -> p[partnr].cellnumber;
  int cell;
  
  j=0;
  n=0;

  for( i=-1 ; i<2 ; i++)
	{

	for( k=-NCOLUMNS ; k < NCOLUMNS*2 ; k=k+NCOLUMNS)
		{
		cell = pl -> p[partnr].cellnumber + i + k;
			if (cell >= 0 && cell <= MAX_CELLS)
			  {
			    if (linked_list->head[cell] != -1) // 
			    {
				    if (linked_list->head[cell] == partnr)
				    {
					    n=linked_list->next[linked_list->head[cell]];
				    }
				    else
				    {
				        possible_contact_particles[j] = linked_list->head[cell];
				        j++;
				        n=linked_list->next[linked_list->head[cell]];
				    }
					
				    while (n != -1)
					    {
						    if ( n==partnr)
						    {
							    n=linked_list->next[n];
						    }
						    else
						    {
							    possible_contact_particles[j] = n;
							    j++;
							    n=linked_list->next[n];
						    }

					    }

			    }
			  }
		}
	}
  return j;

}
//----------------------------------------------------------
