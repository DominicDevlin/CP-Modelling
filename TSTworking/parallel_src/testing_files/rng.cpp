#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <iostream>
#include <random>
#include <stdint.h>

uint64_t x=123456789, y=362436069, z=521288629;

/// this one was from stackexchange
uint64_t xorshf96(void) 
{          //period 2^96-1
    uint64_t t = 1;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;
    return z;
}

///non thread-safe one from wikipedia
uint64_t xorshift64star(void) {
    static uint64_t x = 1; /* initial seed must be nonzero, don't use a static variable for the state if multithreaded */
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
}



// thread safe xooshiro256**
uint64_t rotl(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}


static uint64_t s[4]; // need to be set to non zero

uint64_t next(void) {
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


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void jump(void) 
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
			next();	
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



int main()
{
    s[0] = 4;
    s[1] = 5;
    s[2] = 6;
    s[3] = 7;
    jump();

    int it = 1000;
    for (int i=0; i < it; ++i)
    {
        std::cout << next() << std::endl;
        std::cout << to_double(next()) << std::endl;

    }
    return 0;
}
