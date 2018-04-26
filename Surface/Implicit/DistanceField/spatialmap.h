#ifndef SPATIALMAP_H
#define SPATIALMAP_H

#include "mdsimhash.h"
#include "vector3d.h"
#include "scalar.h"
#include "cell.h"

#ifdef WIN32
typedef hash_set<unsigned int> VoxelHashSet;
#else
typedef HASH_VERSION::hash_set<unsigned int> VoxelHashSet;
#endif


template<class T> 
class VoxelCache{
 public:
  VoxelHashSet hset;
	std::vector<T *> cache;
  bool dirty;

  VoxelCache(){
    dirty = false;
  }

  ~VoxelCache(){
    cache.clear();
    hset.clear();
  }

	std::vector< T * > * getContents(){
    if(dirty){
      dirty = false;
      
      unsigned int size = hset.size();
      
      cache.resize(size);
      
      VoxelHashSet::const_iterator i   = hset.begin();
      for(unsigned int n = 0; n < size; n++){
	cache[n] = (T *) (*i);
	i++;
      }
    }
    //    cout << cache.size() << endl;
    return &cache;    
  }
  

  bool empty() const{
    return hset.empty();
  }  

  void add(const T * e){
    dirty = true;
    hset.insert((unsigned int)e);      
  }

  void remove(const T * e){
    dirty = true;
    hset.erase((unsigned int)e);
  }

    
  unsigned int size() const{
    return hset.size();
  }

  bool contains(const T * e){
    if(hset.count((unsigned int)e) == 0){
      return false;
    } else {
      return true;
    }
  } 
};




template <class T > 
class VoxelNeighborhoodMap{
 private:
  scalar cellsize;
  scalar cellsize_inverse;

#ifdef WIN32
  typedef  HASH_VERSION::hash_map< Cell, VoxelCache<T>, cell_hash_compare> voxelcachemap;
#else
  typedef  HASH_VERSION::hash_map< Cell, VoxelCache<T> , hash_cell, equal_cell >  voxelcachemap;
#endif
  voxelcachemap map;

	std::vector< T *> emptyvector;

 public:
  
  VoxelNeighborhoodMap(scalar cell_size){ assert(cell_size > 0.0); cellsize = cell_size, cellsize_inverse = (scalar) 1.0 / cellsize; }
  
  void insert(const Cell& c, const T * tp){
    for(int i = -1; i <= 1; i++){
      for(int j = -1; j <= 1; j++){
	for(int k = -1; k <= 1; k++){
	  insertSingle(c + Cell(i,j,k),tp);
	}
      }
    }
  }

  void erase(const Cell& c, const T * tp){
    for(int i = -1; i <= 1; i++){
      for(int j = -1; j <= 1; j++){
	for(int k = -1; k <= 1; k++){
	  eraseSingle(c + Cell(i,j,k),tp);
	}
      }
    }
  }

	std::vector< T * > * neighbors(const Cell& c){
    typename voxelcachemap::iterator i = map.find(c); 
    
    if(i != map.end()){ 
      return i->second.getContents(); 
    } else { 
      return &emptyvector; 
    } 
  } 

  Cell toCell(const Vector3d &v) const { return toCellInverse(v,cellsize_inverse); }

 protected:
  void insertSingle(const Cell& c, const T * tp){
    map[c].add(tp);
  }

  void eraseSingle(const Cell& c, const T * tp){
    map[c].remove(tp);
    if(map[c].size() == 0){
      map.erase(c);
    }
  }
     
};



template <class T> 
class VoxelContents : public std::vector< std::pair< T, T* > > {
  public:
    void remove(T * tp){
      for(typename VoxelContents::iterator i = VoxelContents::begin(); i != VoxelContents::end(); i++){
	if (i->second == tp){
	  (*i) = VoxelContents::back();
	  VoxelContents::pop_back();
	  return;
	}
      }
    }

    void insert(T * tp){
      push_back(std::pair<T,T*>(*tp, tp));
    }

    void contains(T * tp) const{
      for(typename VoxelContents::iterator i = VoxelContents::begin(); i != VoxelContents::end(); i++){
	if (i->second == tp){
	  return true;
	}
      }
      return false;
    }
};


#ifdef WIN32
template <class T > 
class VoxelMap : public HASH_VERSION::hash_map< Cell, VoxelContents<T>, cell_hash_compare> {
#else
template <class T > 
class VoxelMap : public HASH_VERSION::hash_map< Cell, VoxelContents<T> , hash_cell, equal_cell > {
#endif
 private:
  scalar cellsize;
  scalar cellsize_inverse;

  VoxelContents<T> empty_voxel;
 public:
  
  VoxelMap(scalar cell_size){ assert(cell_size > 0.0); cellsize = cell_size, cellsize_inverse = (scalar) 1.0 / cellsize; }
  
  void insert(const Cell& c, const T * tp){
  }

  void erase(const Cell& c, const T * tp){
  }
  

  Cell toCell(const Vector3d &v) const { return toCellInverse(v,cellsize_inverse); }

};

#ifdef WIN32
template <class T > 
class CellHashMap : public HASH_VERSION::hash_map< Cell, T, cell_hash_compare> {
  
};
#else
template <class T > 
class CellHashMap : public HASH_VERSION::hash_map< Cell, T, hash_cell, equal_cell > {
  // public:
  //CellHashMap<T> (size_t  n) : HASH_VERSION::hash_map< Cell, T, hash_cell, equal_cell >(n){}
};
#endif

#endif
