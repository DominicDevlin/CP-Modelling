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
#include <stdlib.h>
#include <sys/timeb.h>
#include <iostream>
#include "random.h"
#include <cmath>
#include "parameter.h"

extern Parameter par;


// thread safe xooshiro256**
uint64_t rotl(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

// static uint64_t s[4]; // need to be set to non zero, different one for each simulation in private CellularPotts

uint64_t next(uint64_t s[4]) 
{
	const uint64_t result = rotl(s[1] * 5, 7) * 9;

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}


void jump(uint64_t s[4]) 
{
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 1;
	uint64_t s1 = 1;
	uint64_t s2 = 1;
	uint64_t s3 = 1;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			next(s);	
		}
		
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

double to_double(uint64_t x) 
{
   const union { uint64_t i; double d; } u = {.i = UINT64_C(0x3FF) << 52 | x >> 12 };
   return u.d - 1.0;
}

double RANDOM(uint64_t s[4])
{
  double val = to_double(next(s));
  // std::cout << val << std::endl;
  return val;
}

/*! Returns a random integer value between 1 and 'max'
  \param The maximum value (long)
  \return A random integer (long)
**/
long RandomNumber(long max, uint64_t s[4])
{
  double val = ((long)(RANDOM(s)*max+1));
  // std::cout << val << std::endl;
  return val;
}


void Seed(uint64_t s[4])
{
  int x = par.seeder;
  s = new uint64_t[4];
  s[0] = x;
  s[1] = x + 1;
  s[2] = x + 2;
  s[3] = x + 3;
  std::cout << s[0] << " " << s[1] << std::endl;
}


