//GEM utilites header.
//Note: const_cast to allow situational ptr arithmetic.

#pragma once //Declares header variables just once.

template <typename T>
void Allocate2dPointerArray(const T* const startPtr, T**& ptrArray, const int& xdim, const int& ydim)
{
    ptrArray = new T*[ydim];
    for(auto j = 0; j < ydim; ++j)
    {
        ptrArray[j] = const_cast<T*>(startPtr) + j*xdim;
    }
};

template <typename T>
void Allocate3dPointerArray(const T* const startPtr, T***& ptrArray, const int& xdim, const int& ydim, const int& zdim)
{
   ptrArray = new T**[zdim];
   for(int k = 0; k < zdim; k++)
   {
      T** array3dptry_c = new T*[ydim];
      ptrArray[k] = array3dptry_c;
      for(int j = 0; j < ydim; j++)
      {
         ptrArray[k][j] = const_cast<T*>(startPtr)+ k*ydim*xdim + j*xdim;
      }
   }
};

template <typename T>
void Allocate4dPointerArray(const T* const startPtr, T****& ptrArray, const int& xdim, const int& ydim, const int& zdim, const int& tdim)
{
    ptrArray = new T***[tdim];
   for(int l = 0; l < tdim; l++)
   {
      T *** array4dptrz_c = new T**[zdim];
      ptrArray[l] = array4dptrz_c;
      for(int k = 0; k < zdim; k++)
      {
         int ** array4dptry_c = new T*[ydim];
         ptrArray[l][k] = array4dptry_c;
         for(int j = 0; j < ydim; j++)
         {
            ptrArray[l][k][j] = const_cast<T*>(startPtr) + l*zdim*ydim*xdim + k*ydim*xdim + j*xdim;
         }
      }
   }
};