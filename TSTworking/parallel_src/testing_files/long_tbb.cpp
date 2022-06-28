#include "oneapi/tbb.h"
#include <iostream>
#include <vector>
#include <random>
#include <stdint.h>




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
// void jump(void) 
// {
// 	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

// 	uint64_t s0 = 1;
// 	uint64_t s1 = 1;
// 	uint64_t s2 = 1;
// 	uint64_t s3 = 1;
// 	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
// 		for(int b = 0; b < 64; b++) {
// 			if (JUMP[i] & UINT64_C(1) << b) {
// 				s0 ^= s[0];
// 				s1 ^= s[1];
// 				s2 ^= s[2];
// 				s3 ^= s[3];
// 			}
// 			next();	
// 		}
		
// 	s[0] = s0;
// 	s[1] = s1;
// 	s[2] = s2;
// 	s[3] = s3;
// }

double to_double(uint64_t x) 
{
   const union { uint64_t i; double d; } u = {.i = UINT64_C(0x3FF) << 52 | x >> 12 };
   return u.d - 1.0;
}


class mytask {
private:
  size_t m_n{};
public: 
  mytask(size_t n)
    : m_n{n}
  {}

  void operator()(int val) 
  {
    for (int i=0;i<1000000;++i) 
    {}  // Deliberately run slow
    std::cerr << "[" << m_n << "]";
    std::cout << "val is " << val << std::endl;
  }
};






class executor
{
private:
  std::vector<mytask>& m_task;
public:
  executor(std::vector<mytask>& t): m_task{t}
  {}

  void operator()(const tbb::blocked_range<size_t>& r) const {
    for (size_t i=r.begin();i!=r.end();++i)
    {
      int val = 2;
      m_task[i](val);
    }
  }

  //std::vector<mytask>& _tasks;
};

int main(int,char**) {

  // tbb::task_scheduler_init init;  // Automatic number of threads
  // tbb::task_scheduler_init init(2);  // Explicit number of threads

  std::vector<mytask> tasks;
  for (int i=0;i<1000;++i)
    tasks.push_back(mytask(i));


  executor exec(tasks);
  tbb::parallel_for(tbb::blocked_range<size_t>(0,tasks.size()),exec);
  std::cerr << std::endl;

  return 0;
}