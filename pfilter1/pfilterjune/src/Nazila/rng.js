var MersenneTwister = {}
   
  /* initializes mt[N] with a seed */
  MersenneTwister.init_genrand = function(s) {
    MersenneTwister.mt[0] = s >>> 0;
    for (MersenneTwister.mti=1; MersenneTwister.mti<MersenneTwister.N; MersenneTwister.mti++) {
        var s = MersenneTwister.mt[MersenneTwister.mti-1] ^ (MersenneTwister.mt[MersenneTwister.mti-1] >>> 30);
     MersenneTwister.mt[MersenneTwister.mti] = (((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253)
    + MersenneTwister.mti;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        MersenneTwister.mt[MersenneTwister.mti] >>>= 0;
        /* for >32 bit machines */
    }
  }
   
  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  MersenneTwister.init_by_array = function(init_key, key_length) {
    var i, j, k;
    MersenneTwister.init_genrand(19650218);
    i=1; j=0;
    k = (MersenneTwister.N>key_length ? MersenneTwister.N : key_length);
    for (; k; k--) {
      var s = MersenneTwister.mt[i-1] ^ (MersenneTwister.mt[i-1] >>> 30)
      MersenneTwister.mt[i] = (MersenneTwister.mt[i] ^ (((((s & 0xffff0000) >>> 16) * 1664525) << 16) + ((s & 0x0000ffff) * 1664525)))
        + init_key[j] + j; /* non linear */
      MersenneTwister.mt[i] >>>= 0; /* for WORDSIZE > 32 machines */
      i++; j++;
      if (i>=MersenneTwister.N) { MersenneTwister.mt[0] = MersenneTwister.mt[MersenneTwister.N-1]; i=1; }
      if (j>=key_length) j=0;
    }
    for (k=MersenneTwister.N-1; k; k--) {
      var s = MersenneTwister.mt[i-1] ^ (MersenneTwister.mt[i-1] >>> 30);
      MersenneTwister.mt[i] = (MersenneTwister.mt[i] ^ (((((s & 0xffff0000) >>> 16) * 1566083941) << 16) + (s & 0x0000ffff) * 1566083941))
        - i; /* non linear */
      MersenneTwister.mt[i] >>>= 0; /* for WORDSIZE > 32 machines */
      i++;
      if (i>=MersenneTwister.N) { MersenneTwister.mt[0] = MersenneTwister.mt[MersenneTwister.N-1]; i=1; }
    }
  
    MersenneTwister.mt[0] = 0x80000000; /* MSB is 1; assuring non-zero initial array */ 
  }
   
  /* generates a random number on [0,0xffffffff]-interval */
  MersenneTwister.genrand_int32 = function() {
    var y;
    var mag01 = new Array(0x0, MersenneTwister.MATRIX_A);
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
    if (MersenneTwister.mti >= MersenneTwister.N) { /* generate N words at one time */
      var kk;
  
      if (MersenneTwister.mti == MersenneTwister.N+1)   /* if init_genrand() has not been called, */
        MersenneTwister.init_genrand(5489); /* a default initial seed is used */
  
      for (kk=0;kk<MersenneTwister.N-MersenneTwister.M;kk++) {
        y = (MersenneTwister.mt[kk]&MersenneTwister.UPPER_MASK)|(MersenneTwister.mt[kk+1]&MersenneTwister.LOWER_MASK);
        MersenneTwister.mt[kk] = MersenneTwister.mt[kk+MersenneTwister.M] ^ (y >>> 1) ^ mag01[y & 0x1];
      }
      for (;kk<MersenneTwister.N-1;kk++) {
        y = (MersenneTwister.mt[kk]&MersenneTwister.UPPER_MASK)|(MersenneTwister.mt[kk+1]&MersenneTwister.LOWER_MASK);
        MersenneTwister.mt[kk] = MersenneTwister.mt[kk+(MersenneTwister.M-MersenneTwister.N)] ^ (y >>> 1) ^ mag01[y & 0x1];
      }
      y = (MersenneTwister.mt[MersenneTwister.N-1]&MersenneTwister.UPPER_MASK)|(MersenneTwister.mt[0]&MersenneTwister.LOWER_MASK);
      MersenneTwister.mt[MersenneTwister.N-1] = MersenneTwister.mt[MersenneTwister.M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
  
      MersenneTwister.mti = 0;
    }
  
    y = MersenneTwister.mt[MersenneTwister.mti++];
  
    /* Tempering */
    y ^= (y >>> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >>> 18);
  
    return y >>> 0;
  }
   
  /* generates a random number on [0,0x7fffffff]-interval */
  MersenneTwister.genrand_int31 = function() {
    return (MersenneTwister.genrand_int32()>>>1);
  }
   
  /* generates a random number on [0,1]-real-interval */
  MersenneTwister.genrand_real1 = function() {
    return MersenneTwister.genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
  }
  
  /* generates a random number on [0,1)-real-interval */
  MersenneTwister.random = function() {
    return MersenneTwister.genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
  }
   
  /* generates a random number on (0,1)-real-interval */
  MersenneTwister.genrand_real3 = function() {
    return (MersenneTwister.genrand_int32() + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
  }
   
  /* generates a random number on [0,1) with 53-bit resolution*/
  MersenneTwister.genrand_res53 = function() { 
    var a=MersenneTwister.genrand_int32()>>>5, b=MersenneTwister.genrand_int32()>>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
  } 

  MersenneTwister.N = 624;
    MersenneTwister.M = 397;
    MersenneTwister.MATRIX_A = 0x9908b0df;   /* constant vector a */
    MersenneTwister.UPPER_MASK = 0x80000000; /* most significant w-r bits */
    MersenneTwister.LOWER_MASK = 0x7fffffff; /* least significant r bits */
   
    MersenneTwister.mt = new Array(MersenneTwister.N); /* the array for the state vector */

  MersenneTwister.init = function(seed) {
    if (seed == undefined) {
      seed = new Date().getTime();
    } 
    /* Period parameters */  
    MersenneTwister.N = 624;
    MersenneTwister.M = 397;
    MersenneTwister.MATRIX_A = 0x9908b0df;   /* constant vector a */
    MersenneTwister.UPPER_MASK = 0x80000000; /* most significant w-r bits */
    MersenneTwister.LOWER_MASK = 0x7fffffff; /* least significant r bits */
   
    MersenneTwister.mt = new Array(MersenneTwister.N); /* the array for the state vector */
    MersenneTwister.mti=MersenneTwister.N+1; /* mti==N+1 means mt[N] is not initialized */
  
    MersenneTwister.init_genrand(seed);
  }  
  module.exports = MersenneTwister