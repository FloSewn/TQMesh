/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once


#include <vector>         
#include <fstream>         
#include <iostream>         
#include <iomanip>         
#include <functional>


namespace CppUtils {

/*********************************************************************
* Simple function to add whitespace to a file
*********************************************************************/
static inline void write_whitespaces(std::ofstream& of, size_t ns)
{
  std::string ws(ns, ' ');
  of << ws;
}


/*********************************************************************
* Implement a traits class to enable different output based on the
* data type
*
* Reference:
* ----------
* https://stackoverflow.com/questions/2033110/passing-a-string-\
* literal-as-a-type-argument-to-a-class-template
*
*********************************************************************/
template <class T>
struct VtkIOTypeTraits
{
  static const char* name;
};

template <>
struct  VtkIOTypeTraits<int32_t>
{ static const char* name; };
inline const char* VtkIOTypeTraits<int32_t>::name = "Int32";

template <>
struct  VtkIOTypeTraits<int64_t>
{ static const char* name; };
inline const char* VtkIOTypeTraits<int64_t>::name = "Int64";

template <>
struct  VtkIOTypeTraits<float>
{ static const char* name; };
inline const char* VtkIOTypeTraits<float>::name = "Float32";

template <>
struct  VtkIOTypeTraits<double>
{ static const char* name; };
inline const char* VtkIOTypeTraits<double>::name = "Float64";


/*********************************************************************
* Implement interface to store multiple data containers with 
* different types in a single vector.
*
* Reference:
* ----------
* https://dawnarc.com/2019/04/c-one-std-vector-containing-template-\
* class-of-multiple-types/
*********************************************************************/
class VtkIODataInterface
{
public:
  virtual const std::string& name() const = 0;
  virtual size_t dim() const = 0;
  virtual const char* type() const = 0;
  virtual void write_data(std::ofstream& of, size_t n_max) const = 0;
};

template<class T>
class VtkIOData : public VtkIODataInterface
{
public:
  VtkIOData(const std::vector<T>& data, 
            const std::string& name,
            size_t dim) 
  : data_ { data }
  , name_ { name }
  , dim_  { dim }
  {}

  const std::string& name() const { return name_; }
  size_t dim() const { return dim_; }
  const char* type() const { return VtkIOTypeTraits<T>::name; }

  void write_data(std::ofstream& outfile, size_t n_max_row) const
  { 

    write_whitespaces(outfile, 10);
    for (size_t j = 0; j < data_.size(); ++j)
    {
      // Add new line
      if ( j % n_max_row == 0 && j > 0)
      {
        outfile << std::endl;
        write_whitespaces(outfile, 10);
      }

      outfile << data_[j] << " ";
    }
    outfile << std::endl;

  } // VtkIOData::write_data()

private:
  std::vector<T> data_;
  std::string    name_;
  size_t         dim_;

};



/*********************************************************************
* This class handles the output of VTU ASCII files
*
* https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*
*********************************************************************/
class VtuWriter
{
public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  VtuWriter(const std::vector<double>& points,
            const std::vector<size_t>& connectivity,
            const std::vector<size_t>& offsets,
            const std::vector<size_t>& types) 
  : points_ { std::move( points ) }
  , connectivity_ { std::move( connectivity ) }
  , offsets_ { std::move( offsets ) }
  , types_  { std::move( types  ) }
  {}


  /*------------------------------------------------------------------
  | Add cell data
  ------------------------------------------------------------------*/
  template <class T>
  void add_cell_data(const std::vector<T>& data, 
                     const std::string& name,
                     size_t dim)
  {
    cell_data_.emplace_back( 
      new VtkIOData<T> { data, name, dim }
    );
  }

  /*------------------------------------------------------------------
  | Add point data
  ------------------------------------------------------------------*/
  template <class T>
  void add_point_data(const std::vector<T>& data,
                     const std::string& name,
                     size_t dim)
  {
    point_data_.emplace_back( 
      new VtkIOData<T> { data, name, dim }
    );
  }

  /*------------------------------------------------------------------
  | Write the VTU file
  ------------------------------------------------------------------*/
  void write(const std::string& file_name)
  {
    std::ofstream outfile;
    outfile.open(file_name);

    size_t n_points = points_.size() / 3;
    size_t n_cells  = offsets_.size();

    outfile << "<VTKFile type=\"UnstructuredGrid\" "
               "version=\"0.1\" "
               "byte_order=\"LittleEndian\">"
            << std::endl;

    write_whitespaces(outfile, 2);
    outfile << "<UnstructuredGrid>"
            << std::endl;

    write_whitespaces(outfile, 4);
    outfile << "<Piece NumberOfPoints=\""<< n_points << "\" "
            << "NumberOfCells=\"" << n_cells << "\">"
            << std::endl;

    write_point_data(outfile);
    write_cell_data(outfile);

    write_points(outfile);

    write_cells(outfile);

    write_whitespaces(outfile, 4);
    outfile << "</Piece>"
            << std::endl;

    write_whitespaces(outfile, 2);
    outfile << "</UnstructuredGrid>"
            << std::endl;

    outfile << "</VTKFile>" 
            << std::endl;

    outfile.close();

  } // VtuWriter::write()

private:

  /*------------------------------------------------------------------
  | Write point data to a vtu file
  ------------------------------------------------------------------*/
  void write_point_data(std::ofstream& outfile)
  {
    if (point_data_.size() < 1)
      return;

    write_whitespaces(outfile, 6);
    outfile << "<PointData Scalars=\"scalars\">" << std::endl;

    for ( size_t i = 0; i < point_data_.size(); ++i )
    {
      auto type = point_data_[i]->type();
      auto name = point_data_[i]->name();
      auto dim  = point_data_[i]->dim();

      write_whitespaces(outfile, 8);
      outfile << "<DataArray type=\"" << type << "\" "
                 "Name=\"" << name << "\" "
                 "NumberOfComponents=\"" << dim << "\" "
                 "Format=\"ascii\">"
              << std::endl;

      point_data_[i]->write_data( outfile, 10 );

      write_whitespaces(outfile, 8);
      outfile << "</DataArray>" << std::endl;
    }

    write_whitespaces(outfile, 6);
    outfile << "</PointData>" << std::endl;

  } // VtuWriter::write_point_data()

  /*------------------------------------------------------------------
  | Write cell data to a vtu file
  ------------------------------------------------------------------*/
  void write_cell_data(std::ofstream& outfile)
  {
    if (cell_data_.size() < 1)
      return;

    write_whitespaces(outfile, 6);
    outfile << "<CellData Scalars=\"scalars\">" << std::endl;

    for ( size_t i = 0; i < cell_data_.size(); ++i )
    {
      auto type = cell_data_[i]->type();
      auto name = cell_data_[i]->name();
      auto dim  = cell_data_[i]->dim();

      write_whitespaces(outfile, 8);
      outfile << "<DataArray type=\"" << type << "\" "
                 "Name=\"" << name << "\" "
                 "NumberOfComponents=\"" << dim << "\" "
                 "Format=\"ascii\">"
              << std::endl;

      cell_data_[i]->write_data( outfile, 10 );

      write_whitespaces(outfile, 8);
      outfile << "</DataArray>" << std::endl;
    }

    write_whitespaces(outfile, 6);
    outfile << "</CellData>" << std::endl;

  } // VtuWriter::write_cell_data()

  /*------------------------------------------------------------------
  | Write the points to file
  ------------------------------------------------------------------*/
  void write_points(std::ofstream& outfile)
  {
    write_whitespaces(outfile, 6);
    outfile << "<Points>"
            << std::endl;

    write_whitespaces(outfile, 8);
    outfile << "<DataArray type=\"Float32\" "
               "NumberOfComponents=\"3\" "
               "Format=\"ascii\">"
            << std::endl;

    write_whitespaces(outfile, 10);
    for (size_t i = 0; i < points_.size(); ++i)
    {
      if ( i % n_max_row_ == 0 && i > 0)
      {
        outfile << std::endl;
        write_whitespaces(outfile, 10);
      }

      outfile << std::setprecision(5) << std::fixed 
              << points_[i] << " ";
    }
    outfile << std::endl;

    write_whitespaces(outfile, 8);
    outfile << "</DataArray>" << std::endl;

    write_whitespaces(outfile, 6);
    outfile << "</Points>" << std::endl;

  } // VtuWriter::write_points()

  /*------------------------------------------------------------------
  | Write cells to file
  ------------------------------------------------------------------*/
  void write_cells(std::ofstream& outfile)
  {
    write_whitespaces(outfile, 6);
    outfile << "<Cells>" << std::endl;

    // Print connectivity
    write_connectivity( outfile );

    // Print offsets
    write_offsets( outfile );

    // Print types
    write_types( outfile );

    write_whitespaces(outfile, 6);
    outfile << "</Cells>" << std::endl;

  } // VtuWriter::write_cells()

  /*------------------------------------------------------------------
  | Write cell connectivities to file
  ------------------------------------------------------------------*/
  void write_connectivity(std::ofstream& outfile)
  {
    write_whitespaces(outfile, 8);
    outfile << "<DataArray type=\"Int32\" "
               "Name=\"connectivity\" "
               "Format=\"ascii\">"
            << std::endl;

    write_whitespaces(outfile, 10);
    for (size_t i = 0; i < connectivity_.size(); ++i)
    {
      if ( i % n_max_row_ == 0 && i > 0)
      {
        outfile << std::endl;
        write_whitespaces(outfile, 10);
      }

      outfile << std::setw(2) << connectivity_[i] << " ";
    }
    outfile << std::endl;

    write_whitespaces(outfile, 8);
    outfile << "</DataArray>" << std::endl;

  } // VtuWriter::write_connectivity()

  /*------------------------------------------------------------------
  | Write cell offsets to file
  ------------------------------------------------------------------*/
  void write_offsets(std::ofstream& outfile)
  {
    write_whitespaces(outfile, 8);
    outfile << "<DataArray type=\"Int32\" "
               "Name=\"offsets\" "
               "Format=\"ascii\">"
            << std::endl;

    write_whitespaces(outfile, 10);
    for (size_t i = 0; i < offsets_.size(); ++i)
    {
      if ( i % n_max_row_ == 0 && i > 0)
      {
        outfile << std::endl;
        write_whitespaces(outfile, 10);
      }

      outfile << std::setw(2) << offsets_[i] << " ";
    }

    outfile << std::endl;

    write_whitespaces(outfile, 8);
    outfile << "</DataArray>"
            << std::endl;

  } // VtuWriter::write_offsets()

  /*------------------------------------------------------------------
  | Write cell types to file
  ------------------------------------------------------------------*/
  void write_types(std::ofstream& outfile)
  {
    write_whitespaces(outfile, 8);
    outfile << "<DataArray type=\"Int32\" "
               "Name=\"types\" "
               "Format=\"ascii\">"
            << std::endl;

    write_whitespaces(outfile, 10);
    for (size_t i = 0; i < offsets_.size(); ++i)
    {
      if ( i % n_max_row_ == 0 && i > 0)
      {
        outfile << std::endl;
        write_whitespaces(outfile, 10);
      }

      outfile << std::setw(2) << types_[i] << " ";
    }

    outfile << std::endl;

    write_whitespaces(outfile, 8);
    outfile << "</DataArray>" << std::endl;

  } // VtuWriter::write_types()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  std::vector<double> points_;
  std::vector<size_t> connectivity_;
  std::vector<size_t> offsets_;
  std::vector<size_t> types_;

  size_t n_max_row_ { 10 };

  std::vector<std::unique_ptr<VtkIODataInterface>> cell_data_;
  std::vector<std::unique_ptr<VtkIODataInterface>> point_data_;

}; // VtuWriter 

} // namespace CppUtils
