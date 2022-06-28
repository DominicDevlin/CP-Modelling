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
#include "omp.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

using namespace std;

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->ConstructInitCells(*this);
    CPM->SetRandomTypes();
    CPM->start_network(par.start_matrix);
    
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }
}

TIMESTEP 
{
  cerr << "Error" << endl;
}




int main(int argc, char *argv[]) {

  
  par.graphics=false;
  Parameter();
  const int n_orgs = 4;
  
  Dish* dishes = new Dish[n_orgs];

  // omp_set_num_threads(4);
  // #pragma omp parallel for 
  for (int i=0; i < n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i+1);
    // Define initial distribution of cells
    dishes[i].Init();

    int t=0;

    for (t=0;t<par.mcs;t++) 
    {
      if (t < 1100)
      {
        if (t%200==0 && t > 0)
          dishes[i].CPM->Programmed_Division();
        
        if (t > 999 && t % 50 == 0)
          dishes[i].CPM->update_network();
      }
      else
      {
        if (t % 50 == 0)
          dishes[i].CPM->update_network();

        dishes[i].CPM->CellGrowthAndDivision();
      }
      dishes[i].CPM->AmoebaeMove(t);

      if (t % 1000 == 0)
      {
        dishes[i].CPM->print_random_cell();
        dishes[i].CPM->get_ntypes();
        cout << t << " TIME STEPS HAVE PASSED." << endl;
      }
    }
  }
  delete[] dishes;

  return 0;
}
