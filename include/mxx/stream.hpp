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


#include "comm_fwd.hpp"
#include "collective.hpp"
#ifdef MXX_COLLECTIVE_DONE


#ifndef MXX_STREAM_HPP
#define MXX_STREAM_HPP

#include <iostream>
#include <sstream>

namespace mxx {


// stream class that buffers on each process till the explicit `sync_flush`
// is called
// then all data is collectively send to rank 0 and there added to the
// wrapped stream object
template <typename CharT, class Traits>
class sync_basic_ostream : public std::basic_ostream<CharT, Traits> {
protected:
    const mxx::comm& m_comm;
    std::basic_ostream<CharT, Traits>* m_stream;
    int m_root;
    typedef std::basic_stringbuf<CharT, Traits, std::allocator<CharT>> strbuf_t;
    std::unique_ptr<strbuf_t> m_buf;
public:
    typedef std::basic_ostream<CharT, Traits> base_stream;

    // for rank 0
    sync_basic_ostream(const mxx::comm& comm, int root, base_stream& stream) : m_comm(comm), m_stream(&stream), m_root(root), m_buf(new strbuf_t()) {
        MXX_ASSERT(0 <= root && root < comm.size());
        this->rdbuf(m_buf.get());
    }

    sync_basic_ostream(const mxx::comm& comm, int root) : m_comm(comm), m_stream(nullptr), m_root(root), m_buf(new strbuf_t()) {
        MXX_ASSERT(0 <= root && root < comm.size());
        MXX_ASSERT(root != comm.rank()); // the root node can't have a null stream
        this->rdbuf(m_buf.get());
    }

    sync_basic_ostream(sync_basic_ostream&& o)
        : m_comm(o.m_comm), m_stream(o.m_stream), m_root(o.m_root),
          m_buf(std::move(o.m_buf)) {
        o.setstate(std::ios_base::badbit);
        o.rdbuf();
        this->rdbuf(m_buf.get());
    }

    sync_basic_ostream(const sync_basic_ostream& o) = delete;

    void sync_flush() {
        // communicate all data to rank `root`
        std::basic_string<CharT, Traits> str = m_buf->str();
        std::vector<size_t> recv_counts = mxx::gather(str.size(), m_root, m_comm);
        std::vector<CharT> strings = mxx::gatherv(&str[0], str.size(), recv_counts, m_root, m_comm);
        // on `root`: output all strings
        if (m_comm.rank() == m_root) {
            typename std::vector<CharT>::iterator begin = strings.begin();
            for (int i = 0; i < m_comm.size(); ++i) {
                std::basic_string<CharT, Traits> recv_string(begin, begin + recv_counts[i]);
                *m_stream << recv_string;
                begin += recv_counts[i];
            }
        }
        // clear buffer content
        m_buf->str("");
    }

    // sync upon destruction
    virtual ~sync_basic_ostream() {
        sync_flush();
    }
};

using sync_ostream = sync_basic_ostream<char>;

inline sync_ostream sync_cout(const mxx::comm& comm, int root = 0) {
    return comm.rank() == root ? sync_ostream(comm, root, std::cout) : sync_ostream(comm, root);
}

inline sync_ostream sync_cerr(const mxx::comm& comm, int root = 0) {
    return comm.rank() == root ? sync_ostream(comm, root, std::cerr) : sync_ostream(comm, root);
}

} // namespace mxx


#define MXX_STREAM_DONE
#include "comm_def.hpp"

#endif // MXX_STREAM_HPP
#endif
