//
//  simpleVector.h
//  Wrapper around a basic C-style vector (pointer to block of memory), but with
//  C++ allocation using new.
//
//  Created by Luc Janss on 23/10/2020.
//

#ifndef simpleVector_h
#define simpleVector_h

#include <Rcpp.h>
#include <algorithm>
#include "rbayzExceptions.h"

template <typename T>
class simpleVector {
public:

   simpleVector() = default;
   
   // explicit should mean ONLY accepting a size_t argument, not allowing
   // something else to be implicitly converted to size_t.
   // explicit is only possible for constructors with a single argument.
   explicit simpleVector(size_t n) {
      doalloc(n);
      std::fill_n(data, n, T{});
   }

   simpleVector(size_t n, T initvalue) {
      doalloc(n);
      std::fill_n(data, n, initvalue);

   }
   
   virtual ~simpleVector() {
      if (nelem > 0) {
         delete[] data;
      }
   }

   T& operator[](size_t i) { return data[i]; }   // retrieving and setting value using []
   const T& operator[](size_t i) const { return data[i]; }

   void initWith(size_t n, T initvalue) {
      doalloc(n);
      for (size_t i = 0; i < nelem; i++) {
         data[i] = initvalue;
      }
   }

   // These two following init were originally for use with Rcpp::NumericVector, will that
   // still work OK with these templates? Maybe it needs initWith<Rcpp::NumericVector> explicitly.
   // However, it does not mean that you can just throw any type at it; the type must have a size()
   // method and support indexing, and when compiling all actual uses will be checked if they comply.
   template <typename T2>
   void initWith(const T2& v, size_t useElem) {
      if (useElem > static_cast<size_t>(v.size())) {
         throw(generalRbayzError("useElem is larger than actual nelem in simpleVector"));
      }
      doalloc(useElem);
      for (size_t i = 0; i < useElem; i++) {
         data[i] = static_cast<T>(v[i]);
      }
   }

   template <typename T2>
   void initWith(const T2& v) {
      initWith(v, static_cast<size_t>(v.size()));
   }

   void initWith(const simpleVector<T>& X) {
      doalloc(X.nelem);
      for (size_t i = 0; i < nelem; i++) {
         data[i] = X.data[i];
      }
   }

   void swap(simpleVector<T>* other) {
      std::swap(data, other->data);
      std::swap(nelem, other->nelem);
   }

   void add(const simpleVector<T>& other) {
      if (nelem != other.nelem) {
         throw(generalRbayzError("Vector sizes do not match in simpleVector::add"));
      }
      for (size_t i = 0; i < nelem; i++) {
         data[i] += other.data[i];
      }
   }

   void sub(const simpleVector<T>& other) {
      if (nelem != other.nelem) {
         throw(generalRbayzError("Vector sizes do not match in simpleVector::sub"));
      }
      for (size_t i = 0; i < nelem; i++) {
         data[i] -= other.data[i];
      }
   }

   void fill(T value) {
      for (size_t i = 0; i < nelem; i++) {
         data[i] = value;
      }
   }

   T* data = 0;
   size_t nelem = 0;

private:
   void doalloc(size_t n) {
      if (n == 0) {
         throw(generalRbayzError("Zero size in initialisation in simpleVector"));
      }
      if (nelem > 0) {  // already allocated
         throw(generalRbayzError("Attempting realloc in simpleVector not supported"));
      }
      data = new T[n];
      nelem = n;
   }
};

using simpleIntVector = simpleVector<int>;
using simpleDblVector = simpleVector<double>;
using simpleSizeVector = simpleVector<size_t>;

#endif /* simpleVector_h */
