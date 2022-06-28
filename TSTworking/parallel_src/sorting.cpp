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
#include "parameter.h"
#include "sqr.h"
#include "omp.h"
#include <ctime>


using namespace std;


int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}



int main(int argc, char *argv[]) 
{
  
  // const bool parallel_on{0};
	
  par.inject();
  // Seed(seed_val);

  const int n_orgs = 1;
  par.n_orgs = n_orgs;
  std::cerr << "here" << '\n';


  Dish* dishes = new Dish[n_orgs];

  // for (int i = 0; i < n_orgs; ++i)
  // {
  //   dishes[i].CPM->set_seed(i*3);
  // }


  // omp_set_num_threads(20);
  // #pragma omp parallel for  
  for (int i=0; i < n_orgs; ++i)
  {
    int t;
    dishes[i].CPM->start_network();

    for (t=0;t<par.mcs;t++) 
    {
      dishes[i].CPM->AmoebaeMove();
      dishes[i].CPM->CellGrowthAndDivision();
      // std::cout << t << " TIME STEPS HAVE PASSED" << std::endl;
      // std::cout << "The number of cells is: " << dishes[i].CountCells() << std::endl;
      if (t % 20 == 0 && t > 0)
      	dishes[i].CPM->update_network();
      if (t % 1000 == 0)
        dishes[i].CPM->print_random_cell();
    }
    std::cout << "The number of cells is: " << dishes[i].CountCells() << std::endl;
    std::cout << "simulation number " << i+1 << " is done!" << std::endl;
  }
  return 0;
}
