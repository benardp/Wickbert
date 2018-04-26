// Include the STL hash_set implementation
#ifndef MDSIMHASH_H
#define MDSIMHASH_H


#ifdef __GNUC__

#ifdef __ICC
#include <ext/hash_set>
#include <ext/hash_map>
using namespace stlport;
#define HASH_VERSION stlport
#else

#include <ext/hash_set>
#include <ext/hash_map>
using namespace __gnu_cxx;
#define HASH_VERSION __gnu_cxx

#endif

#elif WIN32
#include <hash_set>
#include <hash_map>
using namespace stdext;
#define HASH_VERSION stdext
#endif


#endif
