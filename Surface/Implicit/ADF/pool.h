#ifndef SPACEPOOL_INCLUDED // -*- C++ -*-
#define SPACEPOOL_INCLUDED

#ifdef _MSC_VER
#pragma once
#pragma warning(disable:4786)
#endif

/************************************************************************
  Built on top of the array class, to provide a memory pool interface.
  Maintains a list of free nodes, providing constant-time allocation and
  de-allocation.  Allocated nodes must be explicitly returned to the
  memory pool, of course.
  Added features to allow the pool to be resized.  In this case, the pool
  doubles in size when it runs out of items.

  $Id: pool.h,v 1.1.1.1 2004/02/04 03:19:51 wensu Exp $
 ************************************************************************/

// #include <gfx/gfx.h>
// #include <gfx/array.h>

#include <vector>

/**
 * Blocks of pool items are allocated in chunks, to decrease allocation
 * overhead.  Chunks are simply arrays with a pointer to the next chunk
 * (used for freeing a pool).
 */
template<class T>
class chunk 
{
public:
    chunk(int n) { items = new T[n]; }
    ~chunk() { delete[] items; }

public:
    chunk<T> *  nextp;
    T *  items;
};


template<class T>
class pool
{
public: pool(int cs, int nc, bool grow = true) : free_list(cs*nc)
    { 
        chunk_size = cs;
        is_growable = grow;
        pool_size = 0;
        allocated_count = 0;

        chunkp = NULL;
        for( int i=0; i<nc; i++ )
            alloc_chunk();
    }

    ~pool()
    {
        while( chunkp != NULL )
        {
            chunk<T> *  oldp = chunkp;
            chunkp = chunkp->nextp;
            delete oldp;
        }
    }
    
    T *  alloc()
    {
        if ( free_list.size() == 0 )
        {
            if ( !is_growable )
                return NULL;

            alloc_chunk();
        }

        allocated_count++;
		T *element=free_list.back();
        free_list.pop_back();
		return element;
    }

    void dealloc(T *  t)
    {
        allocated_count--;
        free_list.push_back(t);
    }

    void reset()
    {
        free_list.reset();
        chunk<T> *  cp = chunkp;
        while( cp != NULL )
        {
            for( int i=0; i<chunk_size; i++ )
                free_list.add(&cp->items[i]);
        }
        allocated_count = 0;
    }

    int get_allocated_count()
    {
        return allocated_count;
    }


protected:
    pool() { }     // Prevent no-arg constructor.

    void alloc_chunk()
    {
        chunk<T> *  new_chunkp = new chunk<T>(chunk_size);
        pool_size += chunk_size;
        for( int i=0; i<chunk_size; i++ )
            free_list.push_back( &(new_chunkp->items[i]) );
        new_chunkp->nextp = chunkp;
        chunkp = new_chunkp;
    }


private:
    chunk<T> *  chunkp;
    int chunk_size;
    std::vector<T *>  free_list;
    bool is_growable;
    int pool_size;
    int allocated_count;
};

// SPACEPOOL_INCLUDED
#endif
