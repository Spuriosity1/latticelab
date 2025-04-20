#ifndef POINTMAP_CXX_H
#define POINTMAP_CXX_H

#include <stdexcept>
#include <vector>
#include "modulus.hpp"
/** A specialised "hashmap" implementation
 * exploiting the knowledge that operator[] takes in an integer vector of the 
 * form nx, ny, nz, entrywise 0 <= nx < Lx
 * 
 */
template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
class pointMap {
public:
	pointMap(const KeyType& L_, const size_t max_memory_=2000);
	/**
	 * Attempts to insert sublattice 'value' at position 'key'.
	 * Returns 'true' if the value was overwritten
	 */
	bool insert(const KeyType& key, ResType value);

	// Fast / dangerous access methods (bounds check only enabled buy DEBUG)
	ResType& operator[](const KeyType& key);
	ResType operator[](const KeyType& key) const;
	
	// always does bounds check
	ResType at(const KeyType& key) const;
	// method throws if key is not defined
	ResType get_sl(const KeyType& key);

	const KeyType L;
	const size_t max_memory;


	bool inbounds(const KeyType& x) const;

protected:
	std::vector<ResType> values;
/*	TODO reindex everything if there is a common factor of 2^n
 *	to save memory
 *	auto gcd(const KeyType& x){
		return std::gcd(std::gcd(x[0],x[1]),x[2]);
	}*/

};


////////////////////////////////////////
/////// IMPLEMENTATION ///////////////
/// has to be up here because it's templated
///
template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
pointMap<KeyType,ResType>::pointMap(const KeyType& L_, const size_t max_memory_):
	L(L_),
	max_memory(max_memory_)
{
	if (L[0] <= 0 || L[1] <= 0 || L[2] <= 0) {
		throw std::invalid_argument("Unit cell lengths given are negative");
	}
	this->values.resize(L[0]*L[1]*L[2]);
	for (auto& v : values){
		v = -1;
	}
	if ( static_cast<size_t>(L[0]*L[1]*L[2]) > max_memory ) {
		throw std::out_of_range("Unit cell is to large");
	}
}

/**
 * Attempts to insert sublattice 'value' at position 'key'.
 * Returns 'true' if the value was overwritten
 */
template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
bool pointMap<KeyType,ResType>::insert(const KeyType& key, ResType value){
	if (!inbounds(key))
		throw std::out_of_range("key out of bounds in insert");
	ResType& v = values[(key[2]*L[1]+key[1])*L[0]+key[0]];
	if (v == -1){
		v = value;
		return true;
	} else {
		// Something already there
		return false;
	}
}

template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
ResType& pointMap<KeyType,ResType>::operator[](const KeyType& key){
#ifdef DEBUG
	if (!inbounds(key))
		throw std::out_of_range("key out of bounds in pointMap::operator[]");
#endif
	return values[(key[2]*L[1]+key[1])*L[0]+key[0]];
}

template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
ResType pointMap<KeyType,ResType>::operator[](const KeyType& key) const {
#ifdef DEBUG
	if (!inbounds(key))
		throw std::out_of_range("key out of bounds in pointMap::operator[]");
#endif
	return values[(key[2]*L[1]+key[1])*L[0]+key[0]];
}


template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
ResType pointMap<KeyType,ResType>::at(const KeyType& key) const {
	if (!inbounds(key))
		throw std::out_of_range("key out of bounds in pointMap::at");
	return values[(key[2]*L[1]+key[1])*L[0]+key[0]];
}

template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
ResType pointMap<KeyType,ResType>::get_sl(const KeyType& key){
#ifdef DEBUG
	if (values[(key[2]*L[1]+key[1])*L[0]+key[0]] == -1)
		throw std::out_of_range("Sublattice is not initialised");
#endif
	return (*this)[key];
}

template<typename KeyType, typename ResType>
requires has_static_size<KeyType>
bool pointMap<KeyType,ResType>::inbounds(const KeyType& x) const {
	return (x[0] >= 0 && x[0] < L[0] && x[1] >= 0 
			&& x[1] < L[1] && x[2] >= 0 && x[2] < L[2]);
}



#endif
