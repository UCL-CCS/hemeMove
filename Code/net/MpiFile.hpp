// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIFILE_HPP
#define HEMELB_NET_MPIFILE_HPP

#include "net/MpiFile.h"

namespace hemelb::net
{

    template<typename T, std::size_t N>
    void MpiFile::Read(std::span<T, N> buffer, MPI_Status* stat)
    {
        MpiCall{MPI_File_read}(*filePtr, buffer.data(), buffer.size(), MpiDataType<T>(), stat);
    }
    template<typename T, std::size_t N>
    void MpiFile::ReadAt(MPI_Offset offset, std::span<T, N> buffer, MPI_Status* stat)
    {
	MpiCall{MPI_File_read_at}(*filePtr, offset, buffer.data(), buffer.size(), MpiDataType<T>(), stat);
    }
    template<typename T, std::size_t N>
    void MpiFile::ReadAtAll(MPI_Offset offset, std::span<T, N> buffer, MPI_Status* stat)
    {
      MpiCall{MPI_File_read_at_all}(*filePtr, offset, buffer.data(), buffer.size(), MpiDataType<T>(), stat);

    }

    template<typename T, std::size_t N>
    void MpiFile::Write(std::span<T const, N> buffer, MPI_Status* stat)
    {
      MpiCall{MPI_File_write}(*filePtr, buffer.data(), buffer.size(), MpiDataType<T>(), stat);
    }
    template<typename T, std::size_t N>
    void MpiFile::WriteAt(MPI_Offset offset, std::span<T const, N> buffer, MPI_Status* stat)
    {
      MpiCall{MPI_File_write_at}(*filePtr, offset, buffer.data(), buffer.size(), MpiDataType<T>(), stat);
    }
}

#endif
