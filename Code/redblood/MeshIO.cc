// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/MeshIO.h"

#include <filesystem>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "log/Logger.h"
#include "redblood/Cell.h"
#include "redblood/CellBase.h"
#include "redblood/CellEnergy.h"
#include "redblood/VTKError.h"
#include "util/Iterator.h"

namespace hemelb::redblood {

    auto MeshIO::readFile(std::string const &filename, bool fixFacetOrientation) const -> MeshPtr {
      return read(Storage::file, filename, fixFacetOrientation);
    }
    auto MeshIO::readString(std::string const &data, bool fixFacetOrientation) const -> MeshPtr {
      return read(Storage::string, data, fixFacetOrientation);
    }

    void MeshIO::writeFile(std::string const &filename,
			   MeshData const& mesh, util::UnitConverter const* c, PointScalarData const& data) const {
      writeFile(filename, mesh.vertices, mesh.facets, c, data);
    }

    void MeshIO::writeFile(std::string const &filename,
			   MeshData::Vertices const& vertices, MeshData::Facets const& facets,
			   util::UnitConverter const* c, PointScalarData const& data) const {
      write(Storage::file, filename, vertices, facets, c, data);
    }

    std::string MeshIO::writeString(MeshData const& mesh, util::UnitConverter const* c, PointScalarData const& data) const {
      return writeString(mesh.vertices, mesh.facets, c, data);
    }

    std::string MeshIO::writeString(MeshData::Vertices const& vertices, MeshData::Facets const& facets, util::UnitConverter const* c, PointScalarData const& data) const {
      return write(Storage::string, std::string{}, vertices, facets, c, data);
    }

    void MeshIO::writeFile(std::string const& filename, CellBase const& cell, util::UnitConverter const* c) const {
      writeFile(filename, cell.GetVertices(), cell.GetFacets(), c);
    }
    std::string MeshIO::writeString(CellBase const& cell, util::UnitConverter const* c) const {
      return writeString(cell.GetVertices(), cell.GetFacets(), c);
    }

    //! Write cell-mesh including individual forces
    static MeshIO::ScalarField make_force_field(Cell const& cell,
						const std::string& fieldName,
						const std::vector<LatticeForceVector>& forces)
    {
      std::vector<LatticeForce> force_magnitudes;
      force_magnitudes.reserve(cell.GetNumberOfNodes());
      for(auto force: forces) {
	force_magnitudes.push_back(force.GetMagnitude());
      }
	  
      return std::make_pair(fieldName, force_magnitudes);
    }

    static MeshIO::PointScalarData make_forces(Cell const& cell) {
      MeshIO::PointScalarData fields;

      auto&& moduli = cell.moduli;

      std::vector<LatticeForceVector> forces(cell.GetNumberOfNodes(), LatticeForceVector::Zero());
      cell.facetBending(forces);
      fields.push_back(make_force_field(cell, "bending", forces));

      volumeEnergy(cell.GetVertices(),
		   *cell.GetTemplateMesh().GetData(),
		   moduli.volume,
		   forces,
		   cell.GetScale());
      fields.push_back(make_force_field(cell, "volume", forces));

      surfaceEnergy(cell.GetVertices(),
		    *cell.GetTemplateMesh().GetData(),
		    moduli.surface,
		    forces,
		    cell.GetScale());
      fields.push_back(make_force_field(cell, "surface", forces));

      strainEnergy(cell.GetVertices(),
		   *cell.GetTemplateMesh().GetData(),
		   moduli.strain,
		   moduli.dilation,
		   forces,
		   cell.GetScale());
      fields.push_back(make_force_field(cell, "strain", forces));
      return fields;
    }

    void MeshIO::writeFile(std::string const& filename, Cell const& cell, util::UnitConverter const* c) const {
      writeFile(filename, cell.GetVertices(), cell.GetTemplateMesh().GetFacets(), c, make_forces(cell));
    }
    std::string MeshIO::writeString(Cell const& cell, util::UnitConverter const* c) const {
      return writeString(cell.GetVertices(), cell.GetTemplateMesh().GetFacets(), c, make_forces(cell));
    }

    // 
    // Timm Krueger's mesh format
    //
    
    // Rvalue ref OK cos ifstream has virtual d'tor
    static std::shared_ptr<MeshData> read_krueger_mesh(std::istream&& stream, bool fixFacetOrientation)
    {
        log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from stream");

        std::string line;

        // Drop header
        for (int i(0); i < 4; ++i)
            std::getline(stream, line);

        // Number of vertices
        unsigned int num_vertices;
        stream >> num_vertices;

        // Create Mesh data
        auto result = std::make_shared<MeshData>();
        result->vertices.resize(num_vertices);

        // Then read in first and subsequent lines
        MeshData::Facet::value_type offset;
        {
            auto &v = result->vertices[0];
            stream >> offset
                   >> v[0] >> v[1] >> v[2];
        }
        log::Logger::Log<log::Trace, log::Singleton>("Vertex 0 at %d, %d, %d",
                                                     result->vertices[0].x(),
                                                     result->vertices[0].y(),
                                                     result->vertices[0].z());

        for (unsigned i = 1; i < num_vertices; ++i)
        {
            unsigned index;
            auto& v = result->vertices[i];
            stream >> index >> v[0] >> v[1] >> v[2];
            // No gaps in vertex index list
            assert(index == (i + offset));
            log::Logger::Log<log::Trace, log::Singleton>("Vertex %i at %d, %d, %d",
                                                         i,
                                                         result->vertices[i].x(),
                                                         result->vertices[i].y(),
                                                         result->vertices[i].z());
        }

        // Drop mid-file headers
        for (int i = 0; i < 3; ++i)
            std::getline(stream, line);

        // Read facet indices
        unsigned num_facets;
        MeshData::Facet indices;
        stream >> num_facets;
        result->facets.resize(num_facets);

        for (unsigned i = 0; i < num_facets; ++i)
        {
            unsigned dummy;
            stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
                   >> indices[0] >> indices[1] >> indices[2];

            for (unsigned dim = 0; dim < 3; ++dim)
                result->facets[i][dim] = indices[dim] - offset;
            log::Logger::Log<log::Trace, log::Singleton>("Facet %i with %i, %i, %i",
                                                         i,
                                                         result->facets[i][0],
                                                         result->facets[i][1],
                                                         result->facets[i][2]);
        }

        log::Logger::Log<log::Debug, log::Singleton>("Read %i vertices and %i triangular facets",
                                                     num_vertices,
                                                     num_facets);

        if (fixFacetOrientation)
            orientFacets(*result);

        return result;
    }

    auto KruegerMeshIO::read(Storage mode, std::string const &filename_or_data, bool fixFacetOrientation) const -> MeshPtr {
        switch (mode) {
            case Storage::file:
                if (!std::filesystem::exists(filename_or_data.c_str()))
                    throw (Exception() << "Red-blood-cell mesh file '" << filename_or_data << "' does not exist");

                log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from %s",
                                                             filename_or_data.c_str());
                // Open file if it exists
                return read_krueger_mesh(std::ifstream{filename_or_data}, fixFacetOrientation);
            case Storage::string:
                return read_krueger_mesh(std::istringstream{filename_or_data}, fixFacetOrientation);
        }
    }

    static void write_krueger_mesh(std::ostream &stream, MeshData::Vertices const& vertices, MeshData::Facets const& facets, util::UnitConverter const* converter)
    {
      // Write Header
      stream << "$MeshFormat\n2 0 8\n$EndMeshFormat\n" << "$Nodes\n" << vertices.size()
	     << "\n";

      auto i_vertex = vertices.begin();
      auto const i_vertex_end = vertices.end();

      for (unsigned i = 1; i_vertex != i_vertex_end; ++i_vertex, ++i)
      {
          if (converter) {
              auto const vertex = converter->ConvertPositionToPhysicalUnits(*i_vertex);
              stream << i << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
          } else {
              auto&& vertex = *i_vertex;
              stream << i << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
          }
      }

      stream << "$EndNode\n" << "$Elements\n" << facets.size() << "\n";

      auto i_facet = facets.begin();
      auto const i_facet_end = facets.end();

      for (unsigned i(1); i_facet != i_facet_end; ++i_facet, ++i)
        stream << i << " 1 2 3 4 5 " << (*i_facet)[0] + 1 << " " << (*i_facet)[1] + 1 << " "
            << (*i_facet)[2] + 1 << "\n";

      stream << "$EndElement\n";
    }

    std::string KruegerMeshIO::write(Storage m, std::string const& filename,
                                     MeshData::Vertices const& vertices, MeshData::Facets const& facets,
                                     util::UnitConverter const* c, PointScalarData const& data) const {
        if (!data.empty()) {
            std::string msg{"Krueger mesh IO does not support point data. Omitting fields:"};
            for (auto&& pair: data) {
                msg += " " + pair.first;
            }
            log::Logger::Log<log::Warning, log::Singleton>(msg);
        }

        switch (m) {
            case Storage::file:
            {
                log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell to %s",
                                                             filename.c_str());
                std::ofstream file(filename.c_str());
                assert(file.is_open());
                write_krueger_mesh(file, vertices, facets, c);
                return {};
            }
            case Storage::string:
            {
                std::ostringstream ss;
                write_krueger_mesh(ss, vertices, facets, c);
                return ss.str();
            }

        }
    }

    //
    // VTK format
    // 

    auto VTKMeshIO::readUnoriented(Storage mode, std::string const &filename_or_data) const -> std::tuple<MeshPtr, PolyDataPtr>
    {
        // Read in VTK polydata object
        VtkErrorsThrow t;
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        switch (mode) {
            case Storage::file:
                log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from VTK polydata file");
                reader->ReadFromInputStringOff();
                reader->SetFileName(filename_or_data.c_str());
                break;
            case Storage::string:
                log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from VTK polydata string");
                reader->ReadFromInputStringOn();
                reader->SetInputString(filename_or_data);
                break;
        }

      reader->Update();

      vtkSmartPointer<vtkPolyData> polydata(reader->GetOutput());

      // Number of vertices
      unsigned int num_vertices = polydata->GetNumberOfPoints();

      // Create Mesh data
      std::shared_ptr<MeshData> result(new MeshData);
      result->vertices.resize(num_vertices);

      // Then read in first and subsequent lines
      for (unsigned int i(0); i < num_vertices; ++i)
      {
        double* point_coord = polydata->GetPoints()->GetPoint(i);
        result->vertices[i] = {point_coord[0], point_coord[1], point_coord[2]};

        log::Logger::Log<log::Trace, log::Singleton>("Vertex %i at %e, %e, %e",
                                                     i,
                                                     result->vertices[i].x(),
                                                     result->vertices[i].y(),
                                                     result->vertices[i].z());
      }

      // Read facet indices
      unsigned int num_facets = polydata->GetNumberOfCells();
      result->facets.resize(num_facets);

      for (unsigned int i(0); i < num_facets; ++i)
      {
        vtkCell* triangle = polydata->GetCell(i);
        assert(triangle->GetCellType() == VTK_TRIANGLE);
        result->facets[i][0] = triangle->GetPointId(0);
        result->facets[i][1] = triangle->GetPointId(1);
        result->facets[i][2] = triangle->GetPointId(2);
        log::Logger::Log<log::Trace, log::Singleton>("Facet %i with %i, %i, %i",
                                                     i,
                                                     result->facets[i][0],
                                                     result->facets[i][1],
                                                     result->facets[i][2]);
      }

      log::Logger::Log<log::Debug, log::Singleton>("Read %i vertices and %i triangular facets",
                                                   num_vertices,
                                                   num_facets);

      return std::make_tuple(result, polydata);
    }

    auto VTKMeshIO::read(Storage mode, std::string const &filename_or_data, bool fixFacetOrientation) const -> MeshPtr {
      MeshIO::MeshPtr meshData;
      VTKMeshIO::PolyDataPtr polyData;
      std::tie(meshData, polyData) = readUnoriented(mode, filename_or_data);

      if (fixFacetOrientation) {
	unsigned numSwapped = orientFacets(*meshData, *polyData);
	log::Logger::Log<log::Debug, log::Singleton>("Swapped %d facets", numSwapped);
      }

      return meshData;
    }

    std::string VTKMeshIO::write(Storage m, std::string const &filename,
                                 MeshData::Vertices const& vertices, MeshData::Facets const& facets,
                                 util::UnitConverter const* c, PointScalarData const& pt_scalar_fields) const {
        VtkErrorsThrow t;

        // Build the vtkPolyData
        auto pd = vtkSmartPointer<vtkPolyData>::New();

        // First, the points/vertices
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToDouble();
        points->SetNumberOfPoints(ssize(vertices));
        for (auto&& [i, v_lat]: util::enumerate(vertices)) {
            if (c) {
                auto const v = c->ConvertPositionToPhysicalUnits(v_lat);
                points->SetPoint(i, v.m_values.data());
            } else {
                points->SetPoint(i, v_lat.m_values.data());
            }
        }
        pd->SetPoints(points);

        // Second, the polys/facets/triangles
        auto tris = vtkSmartPointer<vtkCellArray>::New();
        tris->AllocateExact(ssize(facets), 3);
        for (auto&& tri: facets) {
            tris->InsertNextCell({tri[0], tri[1], tri[2]});
        }
        pd->SetPolys(tris);

        // Third, the scalar point data
        for (auto&& field: pt_scalar_fields) {
            auto&& name = field.first;
            auto&& data = field.second;
            auto da = vtkSmartPointer<vtkDoubleArray>::New();
            da->SetName(name.c_str());
            da->SetNumberOfComponents(1);
            // This let's VTK "borrow" the data inside our vector.  We
            // know we're not modifying the polydata, so the const cast is
            // OK.
            da->SetArray(const_cast<double*>(data.data()), ssize(data), 1);
            pd->GetPointData()->AddArray(da);
        }

        // Now, the writer
        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        // The vtkPolyData we just made
        writer->SetInputData(pd);

        // Based on the type of write, configure the writer, run it and
        // return any output string.
        switch (m) {
            case Storage::file:
                log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell to %s",
                                                             filename.c_str());
                writer->WriteToOutputStringOff();
                writer->SetFileName(filename.c_str());
                writer->Write();
                return {};

            case Storage::string:
                writer->WriteToOutputStringOn();
                writer->Write();
                return writer->GetOutputString();

        }
    }

}
