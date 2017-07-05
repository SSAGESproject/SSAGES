/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    file.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Block decompose and distribute file as string on MPI communicator.
 *
 * This is not a substitute for MPI_File functionality.
 * TODO:
 * - [ ] Implement proper MPI_File functions for parallel reading of files
 * - [ ] HDF5 ?
 */

#ifndef MXX_FILE_HPP
#define MXX_FILE_HPP

#include <mpi.h>

// C++ includes
#include <string>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <iostream>
#include <iomanip>

// mxx includes
#include "partition.hpp"

namespace mxx {

std::ifstream::pos_type get_filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

class rangebuf: public std::streambuf {
public:
    rangebuf(std::streampos start,
                    size_t size,
                    std::streambuf* sbuf):
        size_(size), sbuf_(sbuf), buf_(new char[64])
    {
        sbuf->pubseekpos(start, std::ios_base::in);
    }
    int underflow() {
        size_t r(this->sbuf_->sgetn(this->buf_,
            std::min<size_t>(sizeof(this->buf_), this->size_)));
        this->size_ -= r;
        this->setg(this->buf_, this->buf_, this->buf_ + r);
        return this->gptr() == this->egptr()
            ? traits_type::eof()
            : traits_type::to_int_type(*this->gptr());
    }

    ~rangebuf()
    {
        delete [] this->buf_;
    }
protected:
    size_t size_;
    std::streambuf* sbuf_;
    char* buf_;
};

std::string file_block_decompose(const char* filename, MPI_Comm comm = MPI_COMM_WORLD, std::size_t max_local_size = 0)
{
    // TODO: handle error if file doesn't exist

    // get size of input file
    std::size_t file_size = get_filesize(filename);

    // get communication parameters
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // restrict max local size (assuming that it is the same parameter on each
    // processor)
    if (max_local_size > 0 && file_size / p > max_local_size)
        file_size = p*max_local_size;
    partition::block_decomposition<std::size_t> part(file_size, p, rank);
    // block decompose
    std::size_t local_size = part.local_size();
    std::size_t offset = part.excl_prefix_size();

    // open file
    std::ifstream t(filename);
    // wrap in our custom range buffer (of type std::streambuf)
    rangebuf rb(offset, local_size, t.rdbuf());

    // read file (range) buffer into string stream
    std::stringstream ss;
    ss << &rb;

    std::string local_str(ss.str());

    return local_str;
}

template <typename _Iterator>
void write_files(const std::string& filename, _Iterator begin, _Iterator end, MPI_Comm comm = MPI_COMM_WORLD)
{
    // get MPI Communicator properties
    int rank, p;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // get max rank string length:
    std::stringstream sslen;
    sslen << p;
    int rank_slen = sslen.str().size();

    // concat rank at end of filename
    std::stringstream ss;
    ss << filename << "." << std::setfill('0') << std::setw(rank_slen) << p << "." << std::setfill('0') << std::setw(rank_slen) << rank;

    // open file with stream
    //std::cerr << "writing to file " << ss.str() << std::endl;
    std::ofstream outfs(ss.str());

    // write the content into the file, sep by newline
    while (begin != end)
    {
        outfs << *(begin++) << std::endl;
    }
    outfs.close();
}

} // namespace mxx

#endif // MXX_FILE_HPP
