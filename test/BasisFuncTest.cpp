// In case you want/have to test or access private members of
// your Method class, uncomment the following line:
//#define private public

#include "../src/Methods/BasisFunc.h"
#include "gtest/gtest.h"

using namespace SSAGES;

class BasisFuncTest : public ::testing::Test {
protected:
    virtual void SetUp()
    {
    }
    virtual void TearDown() 
    {
    }

    Basis* basistest;
};

TEST_F(BasisFuncTest, DummyTest)
{
    ASSERT_TRUE(true);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

