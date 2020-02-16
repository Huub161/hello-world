/*
 * Silo.c
 *
 * Multiple particle dynamics simulation
 *
 * (c) 2016 Joris Remmers TU/e
 */


#include "mylib.h"
#include <time.h>

int main( void )

{  
  int          iCyc  = 0;            // Cycle counter
  int          iPlot = 0;            // Plot counter
  char         svgfile[20];          // File name for output
  double       ekin  = 0.0;          // Kinetic Energy
 
 clock_t tic = clock();
 
  Plist        plist;
     
  read_input( "silo.dat" , &plist );
    
  while( iCyc < 100 || ekin > 1.0e-8 )
  {
    iCyc++;
        
    if ( iCyc%50 == 0 && plist.ntot < plist.nwall + plist.ndoor + 2000 && plist.ndoor > 0 )
    {	
      add_particle( &plist );
    }
    
    ekin = solve( &plist );
     
    check_particles( &plist );
 
    if ( iCyc%1000 == 0 )
    {
      iPlot++;

      get_filename( svgfile , iPlot );
      
      plot( svgfile , &plist );

      show_info( svgfile , ekin , plist.ntot );
    }

    if ( iCyc > 100 && ekin < 1.0e-4 && plist.ndoor > 0 )
    {
      open_door( &plist );
    }   
  }

  clock_t toc = clock();

  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  return 0;
}

