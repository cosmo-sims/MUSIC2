#ifndef __HDF_IO_HH
#define __HDF_IO_HH

/*
    HDF_IO.hh  --  templated C++ HDF5 front-end functions, v1.2b
      
    Copyright (C) 2006-7  Oliver Hahn  --  ojha@gmx.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <cstdlib>
#include "hdf5.h"

template<typename T>
hid_t GetDataType( void )
{
  if( typeid(T) == typeid(int) )
    return H5T_NATIVE_INT;

  if( typeid(T) == typeid(unsigned) )
    return H5T_NATIVE_UINT;

  if( typeid(T) == typeid(float) )
    return H5T_NATIVE_FLOAT;

  if( typeid(T) == typeid(double) )
    return H5T_NATIVE_DOUBLE;
  
	if( typeid(T) == typeid(long long) )
		return H5T_NATIVE_LLONG;
	
	if( typeid(T) == typeid(unsigned long long) )
		return H5T_NATIVE_ULLONG;
	
	if( typeid(T) == typeid(size_t) )
		return H5T_NATIVE_ULLONG;
	
  
  std::cerr << " - Error: [HDF_IO] trying to evaluate unsupported type in GetDataType\n\n";
  return -1;
}

#include <stdexcept>
 
class HDFException : public std::runtime_error {
	public:
	 HDFException( const std::string &errtxt ) : std::runtime_error(errtxt) { }
};


inline bool DoesFileExist( std::string Filename ){
        bool flag = false;
        std::fstream fin(Filename.c_str(),std::ios::in|std::ios::binary);
        if( fin.is_open() )
                flag=true;
        fin.close();
        return flag;
}

inline void AssertFileOpen( std::string Filename )
{
  if( !DoesFileExist( Filename ) ){
    std::fstream fout( Filename.c_str(), std::ios::out|std::ios::binary);
    fout.close();
  }
}

inline void HDFCreateFile( std::string Filename )
{
  hid_t HDF_FileID;
  HDF_FileID = H5Fcreate( Filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5Fclose( HDF_FileID );
}

template< typename T>
inline void HDFReadVector( const std::string Filename, const std::string ObjName, std::vector<T> &Data )
{
  HDFReadDataset( Filename, ObjName, Data );
}




inline void HDFGetDatasetExtent( const std::string Filename, const std::string ObjName, std::vector<int> &Extent )
{
  hid_t HDF_FileID, HDF_DatasetID, HDF_DataspaceID;

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  HDF_DatasetID = H5Dopen( HDF_FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( HDF_DatasetID < 0 ){
	  std::stringstream ss;
	  ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
	  H5Fclose( HDF_FileID );
	  throw HDFException(ss.str());
	  return;
  }

  //... get space associated with dataset and its extensions
  HDF_DataspaceID = H5Dget_space( HDF_DatasetID );

  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  
  hsize_t *dimsize = new hsize_t[ndims];

  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );

  Extent.clear();
  for(int i=0; i<ndims; ++i )
    Extent.push_back( dimsize[i] );

  delete[] dimsize;
	
  H5Sclose( HDF_DataspaceID );
  H5Dclose( HDF_DatasetID );
  H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFReadDataset( const std::string Filename, const std::string ObjName, std::vector<T> &Data )
{
 
  hid_t HDF_Type, HDF_FileID, HDF_DatasetID, HDF_DataspaceID;
  hsize_t HDF_StorageSize;

  HDF_Type = GetDataType<T>();

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );


  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  HDF_DatasetID = H5Dopen( HDF_FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( HDF_DatasetID < 0 ){
	  std::stringstream ss;
	  ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
	  Data.clear();
	  H5Fclose( HDF_FileID );
	  throw HDFException(ss.str());
	  return;
  }

  //... get space associated with dataset and its extensions
  HDF_DataspaceID = H5Dget_space( HDF_DatasetID );

  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  
  hsize_t dimsize[ndims];

  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );

  HDF_StorageSize = 1;
  for(int i=0; i<ndims; ++i )
    HDF_StorageSize *= dimsize[i];

  //... adjust the array size to hold the data
  Data.clear();
  Data.reserve( HDF_StorageSize );
  Data.assign( HDF_StorageSize, (T)1 );

  if( Data.capacity() < HDF_StorageSize ){
    std::cerr << " - Error: not enough memory to store all data in HDFReadDataset!\n\n";
    abort();
  }

  //... read the dataset
  H5Dread( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Data[0] );

  if( Data.size() != HDF_StorageSize ){
    std::cerr << " - Error: something went wrong while reading!\n\n";
    abort();
  }

  H5Sclose( HDF_DataspaceID );
  H5Dclose( HDF_DatasetID );
  H5Fclose( HDF_FileID );
}

template<typename T >
inline void HDFReadSelect( const std::string Filename, const std::string ObjName, const std::vector<unsigned>& ii, std::vector<T> &Data ){

  hid_t HDF_Type, HDF_FileID, HDF_DatasetID, HDF_DataspaceID, HDF_MemspaceID;
  hsize_t HDF_StorageSize;

  HDF_Type = GetDataType<T>();

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );


  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  HDF_DatasetID = H5Dopen( HDF_FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( HDF_DatasetID < 0 ){
    std::stringstream ss;
    ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
    Data.clear();
    H5Fclose( HDF_FileID );
    throw HDFException(ss.str());
  }

  //... get space associated with dataset and its extensions
  HDF_DataspaceID = H5Dget_space( HDF_DatasetID );
 
  hsize_t block[2];
  block[0]  = ii.size();
  block[1]  = 1;


  Data.clear();
  Data.reserve( block[0]*block[1] );
  Data.assign( block[0]*block[1], (T)1 );
  
  HDF_MemspaceID = H5Screate_simple( 2, block, NULL );
  //  H5Sselect_hyperslab( FilespaceID, H5S_SELECT_SET, offset, stride, count, block );
  H5Sselect_elements( HDF_DataspaceID, H5S_SELECT_SET, ii.size(), (const hsize_t *)&ii[0] );

  H5Dread( HDF_DatasetID, HDF_Type, HDF_MemspaceID, HDF_DataspaceID, H5P_DEFAULT, &Data[0] );
  
  H5Sclose( HDF_DataspaceID );
  H5Sclose( HDF_MemspaceID );
  H5Dclose( HDF_DatasetID );
  H5Fclose( HDF_FileID );

}

template<typename T >
inline void HDFReadVectorSelect( const std::string Filename, const std::string ObjName, const std::vector<unsigned>& ii, std::vector<T> &Data ){

  hid_t HDF_Type, HDF_FileID, HDF_DatasetID, HDF_DataspaceID, HDF_MemspaceID;
//  hsize_t HDF_StorageSize;

  HDF_Type = GetDataType<T>();

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );


  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  HDF_DatasetID = H5Dopen( HDF_FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( HDF_DatasetID < 0 ){
    std::stringstream ss;
    ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
    Data.clear();
    H5Fclose( HDF_FileID );
    throw HDFException(ss.str());
    return;
  }

  //... get space associated with dataset and its extensions
  HDF_DataspaceID = H5Dget_space( HDF_DatasetID );
  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  hsize_t dimsize[ndims];
  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );

  hsize_t block[2];
  block[0]  = ii.size();
  block[1]  = 3;
  
  std::vector<hsize_t> coord;
  for( unsigned i=0; i<ii.size(); ++i )
    for( unsigned k=0; k<3; ++k ){
      coord.push_back(ii[i]);
      coord.push_back(k);
    }
  //std::vector<unsigned>().swap(ii);
  

  

  if( ii.size() == 0 ){
    std::cerr << "attempted to read empty block. skipping....\n"; 
    return;
  }
  //std::cerr << "starting 2 read...\n";
  H5Sselect_none( HDF_DataspaceID );
  if( H5Sselect_elements( HDF_DataspaceID, H5S_SELECT_SET, coord.size()/2, (const hsize_t *)&coord[0] ) < 0 )//(const hsize_t**)&coord[0] ) < 0 )
    std::cerr << " - could not select elements properly\n";

  if(H5Sselect_valid( HDF_DataspaceID )<=0 ){
    std::cerr << "\n - sorry, invalid element selection in file \'"<< Filename.c_str() << "\'. \n - dumping 10 first indices...\n";

    /*for( unsigned i=0; i<10; ++i ){
      for( unsigned k=0; k<3; ++k ){
        std::cerr << coord[3*i+k] << " ";
      }
      std::cerr << "\n";
    }*/

    return;
  }

  std::vector<hsize_t>().swap(coord);
  Data.assign( block[0]*block[1], (T)0 );
  HDF_MemspaceID = H5Screate_simple( 2, &block[0], NULL );

  H5Dread( HDF_DatasetID, HDF_Type, HDF_MemspaceID, HDF_DataspaceID, H5P_DEFAULT, &Data[0] );
  

  H5Sclose( HDF_DataspaceID );
  H5Sclose( HDF_MemspaceID );
  H5Dclose( HDF_DatasetID );
  H5Fclose( HDF_FileID );

}

template< typename T >
inline void HDFReadVectorSlab( const std::string Filename, const std::string ObjName, unsigned nStart, unsigned nCount, std::vector<T> &Data )
{
  hsize_t 
      offset[2],
      stride[2],
      count[2],
      block[2];
    
  hid_t MemspaceID, FilespaceID, DatasetID, FileID;
  hid_t Type =  GetDataType<T>();
  
  FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );



 //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  DatasetID =  H5Dopen( FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( DatasetID < 0 ){
    std::stringstream ss;
    ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
    Data.clear();
    H5Fclose( FileID );
    throw HDFException(ss.str());
    return;
  }

  FilespaceID = H5Dget_space( DatasetID );
  
  offset[0] = nStart;
  offset[1] = 0;
  
  count[0]  = 1;
  count[1]  = 1;
  
  stride[0] = 1;
  stride[1] = 1;
  
  block[0]  = nCount;
  block[1]  = 3;


  Data.clear();
  Data.reserve( block[0]*block[1] );
  Data.assign( block[0]*block[1], (T)1 );
  
  MemspaceID = H5Screate_simple( 2, block, NULL );
  H5Sselect_hyperslab( FilespaceID, H5S_SELECT_SET, offset, stride, count, block );

  H5Dread( DatasetID, Type, MemspaceID, FilespaceID, H5P_DEFAULT, &Data[0] );

  H5Sclose( FilespaceID );
  H5Sclose( MemspaceID );
  H5Dclose( DatasetID );
  H5Fclose( FileID );
}

template< typename T >
inline void HDFReadDatasetSlab( const std::string Filename, const std::string ObjName, unsigned nStart, unsigned nCount, std::vector<T> &Data )
{
  hsize_t 
      offset[2],
      stride[2],
      count[2],
      block[2];
    
  hid_t MemspaceID, FilespaceID, DatasetID, FileID;
  hid_t Type =  GetDataType<T>();
  
  FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );



 //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... probe dataset opening
  DatasetID =  H5Dopen( FileID, ObjName.c_str() );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  //... dataset did not exist or was empty
  if( DatasetID < 0 ){
    std::stringstream ss;
    ss << " - Warning: dataset \'" << ObjName.c_str() << "\' does not exist or is empty.\n";
    Data.clear();
    H5Fclose( FileID );
    throw HDFException(ss.str());
    return;
  }

  FilespaceID = H5Dget_space( DatasetID );
  
  offset[0] = nStart;
  offset[1] = 0;
  
  count[0]  = 1;
  count[1]  = 1;
  
  stride[0] = 1;
  stride[1] = 1;
  
  block[0]  = nCount;
  block[1]  = 1;


  Data.clear();
  Data.reserve( block[0]*block[1] );
  Data.assign( block[0]*block[1], (T)1 );
  
  MemspaceID = H5Screate_simple( 2, block, NULL );
  H5Sselect_hyperslab( FilespaceID, H5S_SELECT_SET, offset, stride, count, block );

  H5Dread( DatasetID, Type, MemspaceID, FilespaceID, H5P_DEFAULT, &Data[0] );

  H5Sclose( FilespaceID );
  H5Sclose( MemspaceID );
  H5Dclose( DatasetID );
  H5Fclose( FileID );
}

template< typename T>
inline void HDFReadGroupAttribute( const std::string Filename, const std::string GroupName, const std::string ObjName, T &Data )
{

  hid_t HDF_Type, HDF_FileID, HDF_GroupID, HDF_AttributeID;
  //  hsize_t HDF_StorageSize;

  HDF_Type = GetDataType<T>();

  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... attempt to open attribute
  
  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, GroupName.c_str() );
  HDF_AttributeID = H5Aopen_name( HDF_GroupID, ObjName.c_str() );

  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_AttributeID < 0 ){
    std::stringstream ss;
    ss << " - Warning: attribute \'" << GroupName.c_str() << "/" << ObjName.c_str() << "\' does not exist or is empty.\n";
    H5Fclose( HDF_FileID );
    throw HDFException(ss.str());
    return;
  }


  H5Aread( HDF_AttributeID, HDF_Type, &Data );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);


  H5Aclose( HDF_AttributeID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );

}

template< typename T>
inline void HDFReadGroupAttribute( const std::string Filename, const std::string GroupName, const std::string ObjName, std::vector<T> &Data )
{

  hid_t HDF_Type, HDF_FileID, HDF_GroupID, HDF_AttributeID, HDF_DataspaceID;
  hsize_t HDF_StorageSize;

  HDF_Type = GetDataType<T>();

  //... save old error handler
  herr_t   (*old_func)(void*);
  void     *old_client_data;

  H5Eget_auto(&old_func, &old_client_data);

  //... turn off error handling by hdf5 library
  H5Eset_auto(NULL, NULL);

  //... attempt to open attribute

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, GroupName.c_str() );
  HDF_AttributeID = H5Aopen_name( HDF_GroupID, ObjName.c_str() );

  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_AttributeID < 0 ){
    std::stringstream ss;
    ss << " - Warning: attribute \'" << GroupName.c_str() << "/" << ObjName.c_str() << "\' does not exist or is empty.\n";
    H5Fclose( HDF_FileID );
    throw HDFException(ss.str());
    return;
  }

  //... get space associated with dataset and its extensions
  HDF_DataspaceID = H5Aget_space( HDF_AttributeID );

  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  
  hsize_t dimsize[ndims];

  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );

  HDF_StorageSize = 1;
  for(int i=0; i<ndims; ++i )
    HDF_StorageSize *= dimsize[i];

  //... adjust the array size to hold the data
  Data.clear();
  Data.reserve( HDF_StorageSize );
  Data.assign( HDF_StorageSize, (T)1 );

  H5Aread( HDF_AttributeID, HDF_Type, &Data[0] );

  //... restore previous error handler
  H5Eset_auto(old_func, old_client_data);

  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFWriteDataset( const std::string Filename, const std::string ObjName, const std::vector<T> &Data )
{

  hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;

  hsize_t HDF_Dims;

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

  HDF_Type                = GetDataType<T>();

  HDF_Dims                = Data.size();
  HDF_DataspaceID         = H5Screate_simple(1, &HDF_Dims, NULL);
  HDF_DatasetID           = H5Dcreate( HDF_FileID, ObjName.c_str(), HDF_Type,
                                       HDF_DataspaceID, H5P_DEFAULT );
  H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &Data[0] );
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );

  H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFWriteGroupDataset( const std::string Filename, const std::string GrpName, const std::string ObjName, const std::vector<T> &Data )
{

  hid_t
    HDF_FileID,
    HDF_GroupID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;

  hsize_t HDF_Dims;

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  
  HDF_GroupID = H5Gopen( HDF_FileID, GrpName.c_str() );

  HDF_Type                = GetDataType<T>();

  HDF_Dims                = Data.size();
  HDF_DataspaceID         = H5Screate_simple(1, &HDF_Dims, NULL);
  HDF_DatasetID           = H5Dcreate( HDF_GroupID, ObjName.c_str(), HDF_Type,
                                       HDF_DataspaceID, H5P_DEFAULT );
  H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &Data[0] );
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );
  
  H5Gclose( HDF_GroupID );

  H5Fclose( HDF_FileID );
}


template< typename T >
inline void HDFWriteDataset2D( const std::string Filename, const std::string ObjName, const std::vector<  std::vector<T> > &Data )
{

  hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;

  hsize_t HDF_Dims[2];

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

  HDF_Type                = GetDataType<T>();

  HDF_Dims[0]             = Data.size();
  HDF_Dims[1]             = Data[0].size();
  HDF_DataspaceID         = H5Screate_simple(2, HDF_Dims, NULL);
  HDF_DatasetID           = H5Dcreate( HDF_FileID, ObjName.c_str(), HDF_Type,
                                       HDF_DataspaceID, H5P_DEFAULT );

  T *tmp = new T[HDF_Dims[0]*HDF_Dims[1]];

  unsigned k=0;
  for(unsigned i=0; i<HDF_Dims[0]; ++i )
    for(unsigned j=0; j<HDF_Dims[1]; ++j ){
      tmp[k++] = (Data[i])[j];
    }

  H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, tmp );

  delete[] tmp;

  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );

  H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFWriteDataset3D( const std::string Filename, const std::string ObjName, unsigned nd[3], const std::vector< T > &Data )
{

  hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;

  hsize_t HDF_Dims[3];

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

  HDF_Type                = GetDataType<T>();

  HDF_Dims[0]             = nd[0];
  HDF_Dims[1]             = nd[1];
  HDF_Dims[2]             = nd[2];
  
  //std::cerr << nd[0]<<nd[1]<<nd[2]<<"\n";
  HDF_DataspaceID         = H5Screate_simple(3, HDF_Dims, NULL);
  HDF_DatasetID           = H5Dcreate( HDF_FileID, ObjName.c_str(), HDF_Type,
                                       HDF_DataspaceID, H5P_DEFAULT );

  H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &Data[0] );

  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );

  H5Fclose( HDF_FileID );
}


template< typename T >
struct HDFHyperslabWriter3Ds
{
	hid_t dset_id_, type_id_, file_id_;
	
	HDFHyperslabWriter3Ds( const std::string Filename, const std::string ObjName, size_t nd[3] )
	{
		hid_t filespace;
		
		hsize_t sizes[4] = { 1, nd[0], nd[1], nd[2] };
		
		type_id_	= GetDataType<T>();
		file_id_	= H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		
		//std::cerr << "creating filespace : 1 x " << nd[0] << " x " << nd[1] << " x " << nd[2] << std::endl;
		filespace	= H5Screate_simple( 4, sizes, NULL );
		dset_id_	= H5Dcreate( file_id_, ObjName.c_str(), type_id_, filespace, H5P_DEFAULT );
		
		H5Sclose(filespace);
	}
	
	~HDFHyperslabWriter3Ds()
	{
		H5Dclose( dset_id_ );
		H5Fclose( file_id_ );
	}
	
	void write_slab( T* data, size_t* count, size_t* offset )
	{
		
		hsize_t counts[4] = { 1, count[0], count[1], count[2] };
		hsize_t offsets[4] = { 0, offset[0], offset[1], offset[2] };
		
		hid_t filespace = H5Dget_space(dset_id_);
		
		//std::cerr << "creating memspace : 1 x " << count[0] << " x " << count[1] << " x " << count[2] << std::endl;
		hid_t memspace  = H5Screate_simple(4, counts, NULL);
		H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offsets, NULL, counts, NULL );
		
		herr_t status;
		status = H5Dwrite(dset_id_, type_id_, memspace, filespace, H5P_DEFAULT, reinterpret_cast<void*>(data));
		H5Sclose(filespace);
		H5Sclose(memspace);
	}
	
};


template< typename T >
inline void HDFWriteDataset3Ds( const std::string Filename, const std::string ObjName, unsigned nd[3], const std::vector< T > &Data )
{
	
	hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;
	
	hsize_t HDF_Dims[4];
	
	HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	
	HDF_Type                = GetDataType<T>();
	
	HDF_Dims[0]             = 1;
	HDF_Dims[1]             = nd[0];
	HDF_Dims[2]             = nd[1];
	HDF_Dims[3]             = nd[2];
	
	//std::cerr << nd[0]<<nd[1]<<nd[2]<<"\n";
	HDF_DataspaceID         = H5Screate_simple(4, HDF_Dims, NULL);
	HDF_DatasetID           = H5Dcreate( HDF_FileID, ObjName.c_str(), HDF_Type,
										HDF_DataspaceID, H5P_DEFAULT );
	
	H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, &Data[0] );
	
	H5Dclose( HDF_DatasetID );
	H5Sclose( HDF_DataspaceID );
	
	H5Fclose( HDF_FileID );
}


template< typename T >
inline void HDFWriteDatasetVector( const std::string Filename, const std::string ObjName, const std::vector<T> &Data )
{

  hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_DataspaceID,
    HDF_Type;

  hsize_t HDF_Dims[2];

  //  hsize_t HDF_Dims;

  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

  HDF_Type                = GetDataType<T>();

  HDF_Dims[0]             = (hsize_t)(Data.size()/3);
  HDF_Dims[1]             = 3;

  if( Data.size() % 3 != 0 ){
    std::cerr << " - Warning: Trying to write vector data in HDFWriteDatasetVector\n"
              << "            but array length not divisible by 3!\n\n";

  }

  HDF_DataspaceID         = H5Screate_simple(2, HDF_Dims, NULL);
  HDF_DatasetID           = H5Dcreate( HDF_FileID, ObjName.c_str(), H5T_NATIVE_FLOAT,
                                       HDF_DataspaceID, H5P_DEFAULT );
  H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &Data[0] );
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );

  H5Fclose( HDF_FileID );
}

inline void HDFCreateGroup( const std::string Filename, const std::string GroupName )
{
	hid_t HDF_FileID, HDF_GroupID;

	HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	HDF_GroupID = H5Gcreate( HDF_FileID, GroupName.c_str(), 0 );
	H5Gclose( HDF_GroupID );
	H5Fclose( HDF_FileID );

}

inline void HDFCreateSubGroup( const std::string Filename, const std::string SuperGroupName, const std::string GroupName )
{
	hid_t HDF_FileID, HDF_GroupID, HDF_SuperGroupID;

	HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  HDF_SuperGroupID = H5Gopen( HDF_FileID, SuperGroupName.c_str() );
	HDF_GroupID = H5Gcreate( HDF_SuperGroupID, GroupName.c_str(), 0 );
	H5Gclose( HDF_GroupID );
  H5Gclose( HDF_SuperGroupID );
	H5Fclose( HDF_FileID );

}

template< typename T >
inline void HDFWriteGroupAttribute( const std::string Filename, const std::string GroupName, const std::string ObjName, const std::vector< T > &Data )
{
 hid_t    HDF_FileID,    
	  HDF_GroupID,    
	  HDF_AttributeID,    
	  HDF_DataspaceID,    
	  HDF_DatatypeID;
 
 hsize_t  HDF_Dims;
 
 HDF_DatatypeID = GetDataType<T>();  
 
 HDF_Dims = (hsize_t)(Data.size());

 
 HDF_FileID      = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );  
 HDF_GroupID     = H5Gopen( HDF_FileID, GroupName.c_str() );  
 HDF_DataspaceID = H5Screate_simple(1, &HDF_Dims, NULL);

 HDF_AttributeID = H5Acreate(HDF_GroupID, ObjName.c_str(), HDF_DatatypeID, HDF_DataspaceID, H5P_DEFAULT);  
 H5Awrite( HDF_AttributeID, HDF_DatatypeID, &Data[0] );  
 H5Aclose( HDF_AttributeID );  
 H5Sclose( HDF_DataspaceID );  
 H5Gclose( HDF_GroupID );  
 H5Fclose( HDF_FileID ); 
}

template< typename T >
inline void HDFWriteDatasetAttribute( const std::string Filename, const std::string DatasetName, const std::string ObjName, const std::vector< T > &Data )
{
	hid_t    HDF_FileID,    
	HDF_DatasetID,    
	HDF_AttributeID,    
	HDF_DataspaceID,    
	HDF_DatatypeID;
	
	hsize_t  HDF_Dims;
	
	HDF_DatatypeID = GetDataType<T>();  
	
	HDF_Dims = (hsize_t)(Data.size());
	
	
	HDF_FileID      = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );  
	HDF_DatasetID   = H5Dopen( HDF_FileID, DatasetName.c_str() );  
	HDF_DataspaceID = H5Screate_simple(1, &HDF_Dims, NULL);
	
	HDF_AttributeID = H5Acreate(HDF_DatasetID, ObjName.c_str(), HDF_DatatypeID, HDF_DataspaceID, H5P_DEFAULT);  
	H5Awrite( HDF_AttributeID, HDF_DatatypeID, &Data[0] );  
	H5Aclose( HDF_AttributeID );  
	H5Sclose( HDF_DataspaceID );  
	H5Dclose( HDF_DatasetID );  
	H5Fclose( HDF_FileID ); 
}


template< typename T >
inline void HDFWriteGroupAttribute( const std::string Filename, const std::string GroupName, const std::string ObjName, T
Data )
{

  hid_t
    HDF_FileID,
    HDF_GroupID,
    HDF_AttributeID,
    HDF_DataspaceID,
    HDF_DatatypeID;

  HDF_DatatypeID = GetDataType<T>();

  
  
  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, GroupName.c_str() );
  HDF_DataspaceID         = H5Screate(H5S_SCALAR);
  HDF_AttributeID         = H5Acreate(HDF_GroupID, ObjName.c_str(), HDF_DatatypeID,
                                      HDF_DataspaceID, H5P_DEFAULT);
  H5Awrite( HDF_AttributeID, HDF_DatatypeID, &Data );
  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFWriteDatasetAttribute( const std::string Filename, const std::string DatasetName, const std::string ObjName, T Data )
{
	
	hid_t
    HDF_FileID,
    HDF_DatasetID,
    HDF_AttributeID,
    HDF_DataspaceID,
    HDF_DatatypeID;
	
	HDF_DatatypeID = GetDataType<T>();
	
	
	
	HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	HDF_DatasetID = H5Dopen( HDF_FileID, DatasetName.c_str() );
	HDF_DataspaceID         = H5Screate(H5S_SCALAR);
	HDF_AttributeID         = H5Acreate(HDF_DatasetID, ObjName.c_str(), HDF_DatatypeID,
										HDF_DataspaceID, H5P_DEFAULT);
	H5Awrite( HDF_AttributeID, HDF_DatatypeID, &Data );
	H5Aclose( HDF_AttributeID );
	H5Sclose( HDF_DataspaceID );
	H5Dclose( HDF_DatasetID );
	H5Fclose( HDF_FileID );
}

template< typename T >
inline void HDFWriteSubGroupAttribute( const std::string Filename, const std::string GroupName,  const std::string SubGroupName, const std::string ObjName, T
Data )
{

  hid_t
    HDF_FileID,
    HDF_GroupID,
    HDF_SubGroupID,
    HDF_AttributeID,
    HDF_DataspaceID,
    HDF_DatatypeID;

  HDF_DatatypeID = GetDataType<T>();

  
  
  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  std::cerr << "opening " << GroupName.c_str() << std::endl;
  HDF_GroupID = H5Gopen( HDF_FileID, GroupName.c_str() );
  std::cerr << "opening " << SubGroupName.c_str() << std::endl;
  HDF_SubGroupID = H5Gopen( HDF_GroupID, SubGroupName.c_str() );
  HDF_DataspaceID         = H5Screate(H5S_SCALAR);
  HDF_AttributeID         = H5Acreate(HDF_SubGroupID, ObjName.c_str(), HDF_DatatypeID,
                                      HDF_DataspaceID, H5P_DEFAULT);
  H5Awrite( HDF_AttributeID, HDF_DatatypeID, &Data );
  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_SubGroupID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

template<>
inline void HDFWriteGroupAttribute<std::string>( const std::string Filename, const std::string GroupName, const std::string ObjName, std::string Data )
{

  hid_t
    HDF_FileID,
    HDF_GroupID,
    HDF_AttributeID,
    HDF_DataspaceID,
    HDF_DatatypeID;

  HDF_DatatypeID = H5Tcopy( H5T_C_S1 );

  H5Tset_size( HDF_DatatypeID, Data.size() );
  H5Tset_strpad(HDF_DatatypeID, H5T_STR_NULLPAD);
  
  HDF_FileID = H5Fopen( Filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, GroupName.c_str() );
  HDF_DataspaceID         = H5Screate(H5S_SCALAR);
  HDF_AttributeID         = H5Acreate(HDF_GroupID, ObjName.c_str(), HDF_DatatypeID,
                                      HDF_DataspaceID, H5P_DEFAULT);
  H5Awrite( HDF_AttributeID, HDF_DatatypeID, Data.c_str() );
  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}
#endif
