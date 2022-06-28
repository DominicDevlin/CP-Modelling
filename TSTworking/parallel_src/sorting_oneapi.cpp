/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "oneapi/tbb.h"

// #ifdef QTGRAPHICS
// #include "qtgraph.h"
// #else
#include "x11graph.h"
// #endif


using namespace std;
using namespace oneapi::tbb;

INIT {

  try {

    // Define initial distribution of cells
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->ConstructInitCells(*this);
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file using 
    // parameters n_init_cells, size_init_cells and divisions
    for (int i=0;i<par.divisions;i++) {
      CPM->DivideCells();
    }
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();
    
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }

}

TIMESTEP { 
 
  try {

    static int i=0;
    // cout << " PRINT SOMETHING " << endl;
  
    static Dish *dish=new Dish();
    static Info *info=new Info(*dish, *this);
    
    dish->CPM->AmoebaeMove(dish->PDEfield);
    
    //cerr << "Done\n";
    if (par.graphics && !(i%par.storage_stride)) {
      
      
      BeginScene();
      ClearImage();
      dish->Plot(this);

      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
     
    }
  
    if (par.store && !(i%par.storage_stride)) {
      char fname[200];
      sprintf(fname,"%s/extend%07d.png",par.datadir,i);
    
      BeginScene();
      ClearImage();    
      dish->Plot(this);
      
      EndScene();
    
      Write(fname);
        
    }

    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

void timestep()
{
  try {

    static int i=0;
    // cout << " PRINT SOMETHING " << endl;
  
    static Dish *dish=new Dish();
    // static Info *info=new Info(*dish, *this);
    
    dish->CPM->AmoebaeMove(dish->PDEfield);
    
    //cerr << "Done\n";

    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
  
}




int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}



int main(int argc, char *argv[]) {
  
	
  try {

    std::vector<int> sims;
    for (int i=0;i<8;++i)
    {
      sims.push_back(i);
    }
    // Read parameters (single sim)
    // par.Read(argv[1]);
    
    // Seed(par.rseed);
    
    // X11Graphics g(par.sizex*2,par.sizey*2);
    // int t;

    // for (t=0;t<par.mcs;t++) {

    //   g.TimeStep();
    
    // }
    // XInitThreads();


    tbb::parallel_for( tbb::blocked_range<int>(0,sims.size()),
                      [&](tbb::blocked_range<int> r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            
            par.inject();
            
            /// seed is based on time so need to change this so they have different seeds. 
            
            if (i == 0)
            {
              Seed(i);
              X11Graphics g1(par.sizex*2, par.sizey*2);
              int t;
              for (t=0;t<par.mcs;t++) 
              {
                std::cout << "got here with sim n " << i <<  std::endl;
                g1.TimeStep();                           
              }
              std::cout << "simulation number " << i << " is done!" << std::endl;            
              
            }
            if (i == 2)
            {
              Seed(i);
              X11Graphics g2(par.sizex*2, par.sizey*2);
              int t;
              for (t=0;t<par.mcs;t++) 
              {
                std::cout << "got here with sim n " << i <<  std::endl;
                g2.TimeStep();                           
              }
              std::cout << "simulation number " << i << " is done!" << std::endl;            
              
            }
            if (i == 3)
            {
              Seed(i);
              X11Graphics g3(par.sizex*2, par.sizey*2);
              int t;
              for (t=0;t<par.mcs;t++) 
              {
                std::cout << "got here with sim n " << i <<  std::endl;
                g3.TimeStep();                           
              }
              std::cout << "simulation number " << i << " is done!" << std::endl;            
              
            }
        }
    });







    
  } catch(const char* error) {
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}