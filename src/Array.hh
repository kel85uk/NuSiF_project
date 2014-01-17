#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <iomanip> //std::setw
#include <vector>
#include <algorithm>

#include "Types.hh"
#include "Debug.hh"

//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements are stored in a contiguous chunk of memory
*/
//*******************************************************************************************************************

template<typename T>
class Array
{

private:

   // For storing field data
   std::vector<T> Data;
   int Xsize,Ysize,Zsize;
   int noOfDim;

public:

   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );

   // Destructor
   ~Array();

   // Depending on your implementation you might need the following:
   // Array(const Array& s);
   // Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline T & operator () ( int i );
   inline T & operator () ( int i ,int j );
   inline T & operator () ( int i, int j, int k );

   // Access operators for const Arrays
   inline const T & operator () ( int i ) const;
   inline const T & operator () ( int i ,int j ) const;
   inline const T & operator () ( int i, int j, int k ) const;
   
   // Minimum and maximum operators for Arrays
   inline T minimum();
   inline T maximum();

   // Minimum and maximum operators for const Arrays
   inline const T minimum() const;
   inline const T maximum() const;

   // initialize the whole array with a constant value
   void fill( T value );


   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;


   // Print the whole array ( for debugging purposes )
   const void print() const;   

};


//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

// For 1D Array
template <typename T>
Array<T>::Array( int xSize ):Data(xSize), Xsize(xSize),  noOfDim(0)
{}

// For 2D Array
template <typename T>
Array<T>::Array( int xSize, int ySize ):Data(xSize*ySize), Xsize(xSize), Ysize(ySize),  noOfDim(1)
{}

// For 3D Array
template <typename T>
Array<T>::Array( int xSize, int ySize, int zSize ):Data(xSize*ySize*zSize), Xsize(xSize), Ysize(ySize),
                                                Zsize(zSize), noOfDim(2)
{}

// Destructor
template <typename T>
Array<T>::~Array()
{}

//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
template <typename T>
inline T& Array<T>::operator ()(int i)
{
   ASSERT_MSG(0 == noOfDim, "Array is not 1D");
   ASSERT_MSG(i<Xsize, "Array out of bounds");
   return Data[i];
}

// Operator() 2D
template <typename T>
inline T& Array<T>::operator ()(int i,int j)
{
   ASSERT_MSG(1 == noOfDim, "Array is not 2D");
   ASSERT_MSG(i<Xsize && j<Ysize, "Array out of bounds");
   return Data[i*(Ysize) + j];
}

// Operator() 3D
template <typename T>
inline T& Array<T>::operator ()(int i, int j, int k)
{
   ASSERT_MSG(2 == noOfDim, "Array is not 3D");
   ASSERT_MSG(i<Xsize && j<Ysize && k<Zsize, "Array out of bounds");
   return Data[i*Ysize*Zsize + j*Zsize + k];
}

// Operator() const 1D
template <typename T>
inline const T& Array<T>::operator ()(int i) const
{
   ASSERT_MSG(0 == noOfDim, "Array is not 1D");
   ASSERT_MSG(i<Xsize, "Array out of bounds");
   return Data[i];
}

// Operator() const 2D
template <typename T>
inline const T& Array<T>::operator ()(int i,int j) const
{
   ASSERT_MSG(1 == noOfDim, "Array is not 2D");
   ASSERT_MSG(i<Xsize && j<Ysize, "Array out of bounds");
   return Data[i*(Ysize) + j];
}

// Operator() const 3D
template <typename T>
inline const T& Array<T>::operator ()(int i, int j, int k) const
{
   ASSERT_MSG(2 == noOfDim, "Array is not 3D");
   ASSERT_MSG(i<Xsize && j<Ysize && k<Zsize, "Array out of bounds");
   return Data[i*Ysize*Zsize + j*Zsize + k];
}

// Minimum operators for Arrays
template <typename T>
inline T Array<T>::minimum()
{
   return (*min_element(Data.begin(), Data.end()));
}

// Maximum operators for Arrays
template <typename T>
inline T Array<T>::maximum()
{
   return (*max_element(Data.begin(), Data.end()));
}

// Minimum operators for const Arrays
template <typename T>
inline const T Array<T>::minimum() const
{
   return (*min_element(Data.begin(), Data.end()));
}

// Maximum operators for const Arrays
template <typename T>
inline const T Array<T>::maximum() const
{
   return (*max_element(Data.begin(), Data.end()));
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
template <typename T>
void Array<T>::fill( T value )
{
   Data.assign(Data.size(),value);
}

// Print the whole array (for debugging purposes)
template <typename T>
const void Array<T>::print() const
{
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value is printed first
   switch(noOfDim)
   {
      case 0: for (int i=0; i<Xsize; i++)
                 std::cout<< std::left<< std::setprecision(6)<< std::setw(12)<<Data[i];

              break;

      case 1: for (int j=Ysize-1; j>=0; j--) {
                 for (int i=0; i<Xsize; i++)
                     std::cout<< std::left<< std::setprecision(6)<< std::setw(12)<<Data[i*Ysize+j];
                 std::cout<<"\n";}
              break;
      case 2: for (int k=0; k<Zsize; k++){
                 std::cout<<"\nz = "<< k << "\n";
                 for (int j=Ysize-1; j>=0; j--) {
                    for (int i=0; i<Xsize; i++)
                        std::cout<< std::left<< std::setprecision(6)<<std::setw(12)<<Data[i*Ysize*Zsize + j*Zsize + k]<<" ";
                    std::cout<<"\n"; } }
              break;
      default:break;
   }
}

//returns required dimension of the array
template <typename T>
int Array<T>::getSize( int dimension ) const
{
   CHECK_MSG(dimension <= noOfDim, "Invalid request for dimension");
   switch(dimension)
   {
      case 0: return(Xsize);
      case 1: return(Ysize);
      case 2: return(Zsize);
      default:return 0;
   }
}

//returns total size of the array
template <typename T>
int Array<T>::getSize() const
{
   switch(noOfDim)
   {
      case 0: return(Xsize);
      case 1: return(Xsize*Ysize);
      case 2: return(Xsize*Ysize*Zsize);
      default:return 0;
   }
}

#endif //ARRAY_HH

