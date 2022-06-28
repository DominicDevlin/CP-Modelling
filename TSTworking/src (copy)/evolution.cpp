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
#include <chrono>
#include <random>

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

using namespace std;

auto start = chrono::steady_clock::now();
//rng for making random networks
auto mseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(mseed) );
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> genes_dist(0, par.n_genes-1);
std::uniform_int_distribution<> activ_dist(0, par.n_activators-1);


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



void swap(int *xp, int *yp)
{
  int temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapv(vector<vector<int>> *xp, vector<vector<int>> *yp)
{
  vector<vector<int>> temp = *xp;
  *xp = *yp;
  *yp = temp;

}

void sorter(vector<vector<vector<int>>> &networks, vector<int> &fitlist)
{
  int i, j, max_idx;
  int n = par.n_orgs;

  for (i = 0; i < n-1; i++)
  {
    max_idx = i;
    for (j=i+1;j<n;j++)
    {
      if (fitlist.at(j) > fitlist.at(max_idx))
        max_idx = j;
    }
    // swap largest element with first element
    swap(&fitlist.at(max_idx), &fitlist.at(i));
    swapv(&networks.at(max_idx), &networks.at(i));
  }
}


// randomise a new network.  
vector<vector<int>> get_random_network()
{
  vector<vector<int>> matrix;
  matrix.resize(par.n_genes);
  for (int i=0; i < par.n_genes;++i)
  {
    matrix.at(i).resize(par.n_activators);
  }

  for (int i = 0; i < par.n_genes; ++i)
  {
    for (int j = 0; j < par.n_activators; ++j)
    {
      double val = double_num(mersenne);
      // slight ON bias for random networks. This is due to theta = -0.3. 
      if (val < 0.01)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.18)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.8)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.98)
      {
        matrix[i][j] = 1;
      }
      else
      {
        matrix[i][j] = 2;
      }
    }
  }
  return matrix;
}


// mutate a network. 
void mutate(vector<vector<int>> &network)
{
  for (int k=0; k < par.n_mutations; ++k)
  {
    double val = double_num(mersenne);
    int i = genes_dist(mersenne);
    int j = activ_dist(mersenne);
    if (val < 0.05)
    {
      network[i][j] = -2;
    }
    else if (val < 0.2)
    {
      network[i][j] = -1;
    }
    else if (val < 0.8)
    {
      network[i][j] = 0;
    }
    else if (val < 0.95)
    {
      network[i][j] = 1;
    }
    else
    {
      network[i][j] = 2;
    }
  }
}



// main function that runs the population for a single evolutionary step. 
vector<int> process_population(vector<vector<vector<int>>>& network_list)
{
  vector<int> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  // run organisms in parallel. 
  omp_set_num_threads(32);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();

    vector<int> intra_org_fitness{};
    int t;

    dishes[i].CPM->start_network(network_list.at(i));

    // run simulation for single org
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

      if (t > par.mcs * 0.8 && t % 100)
      {
        intra_org_fitness.push_back(dishes[i].CPM->get_ntypes());
      }

    }
    // calculate org fitness as minimum over last 20% of time steps. 
    int min_fitness{1000};
    for (int j=0; j < static_cast<int>(intra_org_fitness.size()); ++j)
    {
      if (intra_org_fitness.at(j) < min_fitness)
        min_fitness = intra_org_fitness.at(j);
    }
    inter_org_fitness.at(i) = min_fitness;
    
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }
  delete[] dishes;
  
  // do sorting algorithm and return fitness
  sorter(network_list, inter_org_fitness);

  vector<vector<vector<int>>> nextgen{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    // Currently no random networks are added if largest fitness > 4
    if (inter_org_fitness.front() > 4)
    {
      nextgen.push_back(network_list.at(j));

      //mutate network with probability = 0.5
      double mu = double_num(mersenne);
      if (mu > 0.5)
      {
        mutate(nextgen.back());;
      }
      if (j > par.n_orgs / 4)
        j=0;
      else
        ++j;
    }
    else
    {
      // the last 1/4 are random networks
      if (i > (par.n_orgs * 3)/4)
        nextgen.push_back(get_random_network());
      else 
        nextgen.push_back(network_list.at(j));

      //mutate network with probability = 0.5
      double mu = double_num(mersenne);
      if (mu > 0.5)
      {
        mutate(nextgen.back());
      }
      if (j > par.n_orgs / 4)
        j=0;
      else
        ++j;
    }
  }
  // set nextgen = current pop. 
  for (int i=0;i<par.n_orgs;++i)
    network_list.at(i) = nextgen.at(i);

  return inter_org_fitness;
}


void printn(vector<vector<int>> netw, vector<int> fitn)
{
  // create and open file
  std::string var_name = "gene_networks.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

  for (int i=0;i<par.n_genes;++i)
  {
    if (i == 0)
      outfile << "{ ";
    for (int j=0;j<par.n_activators;++j)
    {
      if (j==0)
        outfile << "{ " << netw[i][j] << ", ";
      else if (j==par.n_activators-1)
        outfile << netw[i][j] << " }, ";
      else 
        outfile << netw[i][j] << ", ";
    }
    if (i == par.n_genes -1)
      outfile << "}" << endl;
  }
  outfile.close();

  // max fitness 
  int max_fit = fitn.front();

  //average fitness
  double avgfit = 0;
  for (int i : fitn)
  {
    avgfit += i;
  }
  avgfit = avgfit / par.n_orgs;

  //calculate time since begin
  auto end = chrono::steady_clock::now();
  auto diff = end - start;

  //output fitness and time since beginning simulation. 
  var_name = "fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  Parameter();

  // make initial random networks. 
  vector<vector<vector<int>>> networks{};
  for (int i=0;i<par.n_orgs;++i)
  {
    networks.push_back(get_random_network());
  }


  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<int> fit = process_population(networks);

    // output every x evolution steps. 
    if (t%1==0)
    {
      printn(networks.front(), fit);
    }
  }
  // finished

  return 0;
}
