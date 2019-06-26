/// Xoroshift reimplemented in C
///
///	Algoritm Written in 2016-2018 by David Blackman and 
/// Sebastiano Vigna (vigna@acm.org)
///
/// http://vigna.di.unimi.it/xorshift/xoroshiro128plus.c///
/// This implementation written in 2019 by Jef Wagner
/// (wagnerj@union.edu)
///
/// The uniform random floats is implemented as described
/// by Andy Gainey in 2016:
///
/// https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers#cast-method

#include <stdint.h>

/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128**. It passes all tests we are aware of except for the four
   lower bits, which might fail linearity tests (and just those), so if
   low linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 64-bit outputs, too; moreover, this
   generator has a very mild Hamming-weight dependency making our test
   (http://prng.di.unimi.it/hwd.php) fail after 5 TB of output; we believe
   this slight bias cannot affect any application. If you are concerned,
   use xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. 

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static inline uint64_t splitmix64( uint64_t x) {
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

/// Sets the interal state of the random number generator.
static inline void rng_seed( uint64_t *s, uint64_t seed){
	s[0] = splitmix64(seed);
	s[1] = splitmix64(s[0]);
}

/// Gives a random long integer between 0 and 2^(64)-1
static inline uint64_t rng_next( uint64_t *s) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c

	return result;
}

/// Gives a random double uniformly distributed in the
/// range (0,1]:  exclusive zero, inclusive one
static inline double rng_uniform( uint64_t *s){
	uint64_t rand_int = rng_next(s);
	uint64_t bits = 0x3ff0000000000000 | (rand_int >> 12);
	return 2.0 - *((double *) &bits);
}
