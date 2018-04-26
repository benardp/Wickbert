/************************************************************************
	written by Wen Su for ADF project 2001/11/26
	changed to a vector of maps
	hashmap <keyT, valueT, hasherT, mapT>
	stores a pair<keyT, valueT>
 ************************************************************************/

#ifndef hashmap_h
#define hashmap_h

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <vector>
#include <map>

template <class keyT, class dataT, class hasherT, class mapT>
class hashmap
{
	typedef std::map<keyT, dataT, mapT> stdmap;
	typedef std::vector< stdmap > hash_array;
public:

	// local iterator
	class iterator
	{
	public:
		iterator()
		{
			vbegin=(typename hash_array::iterator)NULL;
			vend=(typename hash_array::iterator)NULL;
			mbegin=(typename stdmap::iterator)NULL;
			mend=(typename stdmap::iterator)NULL;
		}

		iterator(typename hash_array::iterator vfirst, typename hash_array::iterator vlast, typename stdmap::iterator mfirst, typename stdmap::iterator mlast) {
			vbegin = vfirst;
			vend = vlast;
			mbegin = mfirst;
			mend = mlast;
		}

		bool operator==(const iterator& rhs)
		{
			return ((rhs.mbegin==mbegin) && (rhs.vbegin==vbegin));
		}

		bool operator!=(const iterator& rhs)
		{
			return !(*this == rhs);
		}

		/**
		 * increment
		 * if current list next element not null, use it
		 * else find next list and point to the first non-null element.
		 */
		iterator& operator++()
		{
			mbegin++;
			if (mbegin==mend)
			{
				vbegin++;
				// make sure map is not empty
				while ((vbegin!=vend)&&(vbegin->size()==0))
				{
					vbegin++;
				}
				if (vbegin==vend)
					mbegin=NULL;
				else
					mbegin=vbegin->begin();
				mend=vbegin->end();
			}
			// if to the end return the end;
			return *this;
		}

		std::pair<const keyT, dataT>& operator*() const
		{
			return (*mbegin);
		}
		std::pair<const keyT, dataT>* operator->() const
		{
			return (&(*mbegin));
		}

	private:
		/// current iterator
		typename hash_array::iterator vbegin, vend;
		/// current iterator
		typename stdmap::iterator mbegin, mend;
	};

	/* all public methods for hashmap */

	/// return the first element valid element
	iterator begin()
	{
		// find the first vector that is not size 0;
		for(typename hash_array::iterator hit=table.begin(); hit!=table.end(); hit++)
			if (hit->begin()!=hit->end())
				return iterator(hit, table.end(), hit->begin(), hit->end());
		// not found so return end
		return end();
	}

	/// return the end iterator
	iterator end()
	{
		return iterator(table.end(), table.end(), NULL, NULL);
	}

	/// default constructor initialize the table of size n=127
	hashmap(unsigned int n=127, unsigned int lf=10)
	{
		tableSize=0;
		vectorSize=n;
		table.resize(n);
		loadFactor=lf;
	}

	/// to the users, they do no know about the vector size they just care the number of elements
	unsigned int size() const
	{
		return tableSize;
	}

	/**
	 * insert an item, this may cause expand to be called
	 * @param key the key
	 * @param data the data
	 * @return iterator of the item just inserted
	 */
	iterator insert(keyT key, dataT data)
	{
		// test if need to rehash
		if (needToRehash())
			rehash();

		int i = hasherT()(key) % vectorSize;

		std::pair<stdmap::iterator, bool> pa=table[i].insert(std::pair<keyT, dataT>(key, data));

		if (pa.second)
		{
			tableSize++;
			return iterator(table.begin() + i, table.end(), table[i].begin(), table[i].end());
		}
		else
			return end();
	}

	/// insert an item, calls insert(key, data)
	iterator insert(std::pair<keyT, dataT> pa)
	{
		// call the other insert
		return insert(pa.first, pa.second);
	}

	/// operation []
	dataT & operator[](const keyT& key)
	{
		iterator i=find(key);
		// test if exist
		if (i!=end())
			return i->second;
		else
			i=insert(key, dataT());
		return i->second; 
	}

	/**
	 * return an iterator points to the key, if not found return end
	 * @param v key
	 * @return iterator points to the key
	 */
	iterator find(keyT v)
	{
		int i = hasherT()(v) % vectorSize;
		// find valid items
		typename stdmap::iterator it=table[i].find(v);
		if (it!=table[i].end())
			return iterator(table.begin()+i,table.end(),it,table[i].end());
		// not found
		return end();
	}

	/**
	 * remove a key
	 * @param x key
	 * @return if the key is found and deleted
	 */
	bool remove(keyT x)
	{
		int i = hasherT()(x) % vectorSize;
		// erase the item with key=x
		int j=table[i].erase(x);
		if (j==0)
			// nothing removed
			return false;
		else
		{
			// key found and remove it
			tableSize--;
			return true;
		}
	}
	

	/// erase a key, same as remove.
	bool erase(keyT x)
	{
		return remove(x);
	}

protected:
	/// array that hold the actual values
	hash_array table;
	/// size of the vector
	unsigned int vectorSize;
	/// number of elements
	unsigned int tableSize;
	/// load factor controls how the rehash condition will be set.
	unsigned int loadFactor;

	/// debug 
	unsigned int collision;

	/// check the condition to see if we need to rehash
	bool needToRehash()
	{
		return ((tableSize/vectorSize)>loadFactor);
		// we could also check if one of the map is too long.
	}

	/// expand the table and rehash all keys when the table is too small
	void rehash()
	{
		// create copy of table
		hash_array oldTable(table);
		// clear the current table
		tableSize=0;
		// double the size
		vectorSize=vectorSize+vectorSize;
		table.resize(vectorSize);
		// insert all items
		for (iterator i = begin(); i!= end(); ++i)
			insert(i->first,i->second);
		// debug
		// std::cout << "in hashmap expand : vector size: " << vectorSize << std::endl;
	}
};

#endif
