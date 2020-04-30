#pragma once

#include <mxx/env.hpp>
#include <mxx/comm.hpp>
#include <gtest/gtest.h>

#define private public
#define protected public
#include "Snapshot.h"

// Test accuracy up to 10^-10.
constexpr double eps = 1e-10;

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    return RUN_ALL_TESTS();
}
