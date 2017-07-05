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
 * @file    env.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements a wrapper for a MPI environment.
 */

#ifndef MXX_ENV_HPP
#define MXX_ENV_HPP

#include <mpi.h>
#include <exception>
#include <stdexcept>
#include <string>

namespace mxx {

void my_mpi_errorhandler(MPI_Comm* comm, int* error_code, ...) {
    char buf[MPI_MAX_ERROR_STRING];
    int strlen;
    int error_class;

    MPI_Error_class(*error_code, &error_class);
    MPI_Error_string(error_class, buf, &strlen);
    std::string classstr(buf);
    //fprintf(stderr, "%3d: %s\n", my_rank, error_string);
    MPI_Error_string(*error_code, buf, &strlen);
    std::string errorstr(buf);
    //fprintf(stderr, "%3d: %s\n", my_rank, error_string);
    // throw exception, enables gdb stack trace analysis
    std::terminate();
    throw std::runtime_error("[MPI Error] class: " + classstr + "\n[MPI Error] " + errorstr);
}

class env {
public:

    env() {
        if (!env::initialized()) {
            MPI_Init(NULL, NULL);
        }
    }

    env(int& argc, char**& argv) {
        if (!env::initialized()) {
            MPI_Init(&argc, &argv);
        }
    }

    // TODO: add threading level constructors

    virtual ~env() {
        if (env::initialized() && !env::finalized()) {
            MPI_Finalize();
        }
    }

    /**
     * @brief Returns true if the MPI environment has been initilized with
     *        MPI_Init(_thread).
     */
    static inline bool initialized() {
        int init;
        MPI_Initialized(&init);
        return init != 0;
    }

    /**
     * @brief Returns true if the MPI environment has been finalized with
     *        `MPI_Finalize()`
     */
    static inline bool finalized() {
        int fin;
        MPI_Finalized(&fin);
        return fin != 0;
    }

    static void set_exception_on_error() {
        // set custom error handler (for debugging with working stack-trace on gdb)
        MPI_Errhandler errhandler;
#if MPI_VERSION >= 2
        MPI_Comm_create_errhandler(&my_mpi_errorhandler, &errhandler);
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
#else
        MPI_Errhandler_create(&my_mpi_errorhandler, &errhandler);
        MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);
#endif
    }
};


} // namespace mxx

#endif // MXX_ENV_HPP
