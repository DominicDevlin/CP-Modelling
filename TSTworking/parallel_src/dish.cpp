/* 

Copyright 1996-2006 Roeland Merks, Paulien Hogeweg

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
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "dish.h"
#include "sticky.h"
#include "parameter.h"
#include "info.h"
#include "crash.h"
#include "pde.h"

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

// void Dish::Start_Dish(void) {

//   ConstructorBody();
//   CPM = new CellularPotts(&cell, par.sizex, par.sizey);
//   if (par.n_chem)
//     PDEfield=new PDE(par.n_chem,par.sizex, par.sizey);
  
   
Dish::Dish(void) {
  
  
  ConstructorBody();
  
  CPM=new CellularPotts(&cell, par.sizex, par.sizey);
  if (par.n_chem)
    PDEfield=new PDE(par.n_chem,par.sizex, par.sizey);

  CPM->set_seed();
  par.seeder = par.seeder + 1;
  
  // Initial cell distribution is defined by user in INIT {} block
    
  // if (par.target_area>0)
  //   for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
  //     c->SetTargetArea(par.target_area);
  //   } 

  CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
  
  CPM->ConstructInitCells(*this);

  for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) 
  {
    std::cout << c->Area() << "-" << c->TargetArea() << "   ";
  }
  std::cout << "<- before || after -> " << std::endl;

  for (int j=0;j<par.divisions;j++) {
    CPM->DivideCells();
  }

  for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) 
  {
    std::cout << c->Area() << "-" << c->TargetArea() << "   ";
  }
  std::cout << "<- before || after -> " << std::endl;  



  // Assign a random type to each of the cells
  CPM->SetRandomTypes();
  
  
}

Dish::~Dish() 
{
    cell.clear();
    delete CPM;
}

void Dish::set_amount(void)
{
  std::cout << amount << std::endl;
  std::cout << maxsigma << std::endl;
  amount = 0;
  maxsigma = 0;
}

int Dish::get_amount()
{
  return amount;
}


void Dish::add_amount()
{
  amount += 1;
}


void Dish::sub_amount(void)
{
  --amount;
}

void Dish::add_maxsigma(void)
{
  ++maxsigma;
}

int Dish::get_maxsigma(void)
{
  return maxsigma;
}

void Dish::set_maxsigma(int max)
{
  maxsigma = max;
}
 


// void Dish::kill_dish(void)
// {
//   delete CPM;
//   cell.clear();
// }

void Dish::start_CPM(void)
{
  CPM = new CellularPotts(&cell, par.sizex, par.sizey);
  if (par.n_chem)
    PDEfield=new PDE(par.n_chem,par.sizex, par.sizey);
  else 
    PDEfield=0;
}

void Dish::kill_CPM(void)
{
  delete CPM;
}

void Dish::start_cells(void)
{
  if (par.target_area>0)
    for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
      c->SetTargetArea(par.target_area);
    }
}

void Dish::kill_cells(void)
{
  cell.clear();
}

// void Dish::Plot(Graphics *g) {
//     if (CPM)
//       CPM.Plot(g);
//  }


void Dish::ConstructorBody() {
  
  maxsigma = 0;
  
  // Allocate the first "cell": this is the medium (tau=0)
  cell.push_back(*(new Cell(*this,0)));
  
  // indicate that the first cell is the medium
  cell.front().sigma=0; 
  cell.front().tau=0;
  
  CPM=0;
  PDEfield=0;

}


bool Dish::CellLonelyP(const Cell &c, int **neighbours) const {

  int i;

  for (i=0;i<(int)cell.size();i++) {
    if (neighbours[c.sigma][i]==EMPTY) 
      break;
    else
      if (neighbours[c.sigma][i]>0)
	return false;
  }
  
  return true;
  
}



// Based on code by Paulien Hogeweg.
void Dish::CellGrowthAndDivision(void) {
  // std::cout << "G1" << std::endl;
  vector<bool> which_cells(cell.size());

 
  // if called for the first time: calculate mem_area
  // if (!mem_area) {
  //   mem_area=TargetArea()/CountCells();
  // }
  
  int cell_division=0;
  
  vector<Cell>::iterator c;
  for ( (c=cell.begin(), c++);
	c!=cell.end();
	c++) {
    // std::cout << c->target_area << "  ";
    if ( (c->Area()-c->TargetArea())>c->GrowthThreshold() ) {
      c->IncrementTargetArea();
      
    }
    if ( (c->TargetArea() > par.div_threshold ) ) {
      which_cells[c->Sigma()]=true;
      cell_division++;
    }

    
  }
  
  for (int i = 0; i < which_cells.size(); ++i)
  {
    std::cout << which_cells.at(i) << " ";
  }
  std::cout << std::endl;
   
  // Divide scheduled cells
  if (cell_division) {
    CPM->DivideCells(which_cells);
  }
  // std::cout << "G2" << std::endl;  

}


int Dish::CountCells(void) const {
  
  int amount=0;
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++); i!=cell.end(); i++) {
    if (i->AliveP()) {
      amount++;
    } else {
      // cerr << "Dead cell\n";
    }
  }
  return amount;
}

 

int Dish::Area(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    total_area+=i->Area();
    
  }
  return total_area;
}

int Dish::TargetArea(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    if (i->AliveP()) 
      total_area+=i->TargetArea();
    
  }
  return total_area;
}



void Dish::SetCellOwner(Cell &which_cell) {
  which_cell.owner=this;
}



void Dish::ClearGrads(void) {

  vector<Cell>::iterator i;
  for ( (i=cell.begin(), i++); i!=cell.end(); i++) {
    i->ClearGrad();
  }
}


int Dish::ZygoteArea(void) const {
    return CPM->ZygoteArea();
}

int Dish::Time(void) const {
    return CPM->Time();
}


void Dish::MeasureChemConcentrations(void) {
 
  // clear chemical concentrations
  for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
    for (int ch=0;ch<par.n_chem;ch++) 
      c->chem[ch]=0.;
  }

  // calculate current ones
  for (int ch=0;ch<par.n_chem;ch++)
    for (int i=0;i<SizeX()*SizeY();i++) {
      
      int cn=CPM->Sigma(0,i);
      if (cn>=0) 
	cell[cn].chem[ch]+=PDEfield->Sigma(ch,0,i);
	
    }

    for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
      for (int ch=0;ch<par.n_chem;ch++) 
	c->chem[ch]/=(double)c->Area();
    }

}

int Dish::SizeX(void) { return CPM->SizeX(); }
int Dish::SizeY(void) { return CPM->SizeY(); }	
