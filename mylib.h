#ifndef MYLIB_H
#define MYLIB_H

/*
 * Header file for Library mylib.h for silo.c
 *
 * Multiple particle dynamics simulation
 *
 * (c) 2016 Joris Remmers TU/e
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GRAVITY 9.81
#define MAX_PARTICLES 2200

#define B_CONST  0.05
#define K_CONST  100.0
#define DT       0.001
#define PI_CONST 3.1415

#define MAX_CELLS 1233
#define NROWS 44
#define NCOLUMNS 28
#define MAX_PARTICLES_9CEL 50

//
// struct for the linked list to store particles in for the calculations
//

typedef struct
{
  int head[MAX_CELLS], next[MAX_PARTICLES];
} LINKED_LIST;


//
// Definition of new structure Vec2 (vector with two components
//

typedef struct
{
  double   x,y;
} Vec2;

//
// Definition of new type Particle
//

typedef struct
{
  Vec2     r,v,a,f;  // position, velocity, accelleration, force
  double   radius,mass;
  int      type;
  int	   cellnumber;
} Particle;



//
// Definition of new type Plist (list of particles, unordered)
//

typedef struct
{
  Particle p[MAX_PARTICLES];
  int      nwall;
  int      ndoor;
  int      ntot;
} Plist;


//-----------------------------------------------------------------------------
//  Plot the silo and particle in SVG format
//  pre:    correct filename (ends with .svg) 
//  post:   writes a file to the current directory (.svg)
//  return: -

void plot

  ( char            *name  ,
    Plist           *plist );


//-----------------------------------------------------------------------------
//  Reads the input file
//  pre:    valid input txt file and empty particle list
//  post:   Plist is filled
//  return: -

void read_input

  ( char            *name  ,
    Plist           *p     );

//-----------------------------------------------------------------------------
//  Adds a particle to the list from a point on top of the silo
//  pre:    -
//  post:   a particle is added with a velocity v.y=-2.0
//  return: -

void add_particle

  ( Plist           *pl    );


//-----------------------------------------------------------------------------
//  Initialises particle
//  pre:    -
//  post:   Initialises particle
//  return: -

void init_particle

  ( Particle        *p     );


//-----------------------------------------------------------------------------
//  Opens the silo door
//  pre:    -
//  post:   removes the particles that represent the door and change
//          the colors of the dynamic particles to make a banded 
//          structures (for visualisation)
//  return: -

void open_door

  ( Plist           *p     );


//-----------------------------------------------------------------------------
//  Solve a single step in the verlet algorithm
//  pre:    -
//  post:   Updated particle positions and velocities
//  return: kinetic energy

double solve

  ( LINKED_LIST		*linked_list ,
	Plist           *p           );


//-----------------------------------------------------------------------------
//  Calculate all particle interaction forces
//  pre:    -
//  post:   update force vectors in particle list
//  return: -

void calc_interaction
  
  ( Plist           *plist        , 
    LINKED_LIST     *linked_list  );
 
//-----------------------------------------------------------------------------
//  Calculate the interaction force between two partcles
//  pre:    two valid particles
//  post:   if in contact, updated forces, if not inte forces are zero.
//  return: -

void int_force
  
  ( Particle        *pi ,
    Particle        *pj );
 
  
//-----------------------------------------------------------------------------
//  Add gravity force to all particles
//  pre:    -
//  post:   updated force
//  return: -

void add_gravity

  ( Plist           *plist );


//-----------------------------------------------------------------------------
//  Removes particle from the list
//  pre:    non-empty particle list
//  post:   removes the particle iPar, updates position
//  return: -

void remove_particle

  ( Plist           *plist ,
    int             iPar  );

//-----------------------------------------------------------------------------
//  Get filename for svg output file
//  pre:    -
//  post:   returns filename
//  return: -

void get_filename

  ( char           *names ,
    int            k      );

//-----------------------------------------------------------------------------
//  pre:    -
//  post:   writes number of particles, current ouptut file and kinetic energy to the screen
//  return: -

void show_info

  ( char*    svgfile ,
    double   ekin    ,
    int      ntot    );


// ALSO ADDED THESE, was missing for some reason
void check_particles

  ( Plist 		*p );


void reset_list

  ( LINKED_LIST *linked_list );


void fill_linkedlist

	( LINKED_LIST	*linked_list ,
	  Plist         *plist       );
// STOPPED ADDING THESE

int possible_contact 

  ( Plist       *plist                     , 
    LINKED_LIST *linked_list               ,
    int         ipar                       , 
    int*        pcp                        );

#endif
