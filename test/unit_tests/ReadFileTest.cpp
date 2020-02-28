#include <iostream>
#include <fstream>
#include "gtest/gtest.h"
#include "../src/Utility/ReadFile.h"

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

class ReadFileTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {
        reader = new ReadFile();
        //Set up xyz files
        std::ofstream myfile;

        //Normal file with extra line (should pass)
        myfile.open (filexyz1);
        myfile << "3\n";
        myfile << "Comments here\n";
        myfile << "1 0 0 0\n";
        myfile << "2 1 4 5\n";
        myfile << "3 1 1 1";
        myfile.close();

        //Incorrect number of atoms
        myfile.open (filexyz2);
        myfile << "2\n";
        myfile << "Comments here\n";
        myfile << "1 0 0 0\n";
        myfile << "2 1 1 1\n";
        myfile << "3 1 1 1\n";
        myfile.close();

        //Incorrect number of atoms
        myfile.open (filexyz3);
        myfile << "4\n";
        myfile << "Comments here\n";
        myfile << "1 0 0 0\n";
        myfile << "2 1 1 1\n";
        myfile << "3 1 1 1\n";
        myfile.close();

        // Bad number of atoms
        myfile.open (filexyz4);
        myfile << "garbage\n";
        myfile << "Comments here\n";
        myfile << "1 0 0 0\n";
        myfile << "2 1 1 1\n";
        myfile << "3 1 1 1\n";
        myfile.close();

        //Bad file format
        myfile.open (filexyz5);
        myfile << "3\n";
        myfile << "Comments here\n";
        myfile << "1 0 0 0\n";
        myfile << "2 1 4 5\n";
        myfile << "3 bad 1 1\n";
        myfile.close();

    }

    virtual void TearDown() 
    {
        std::remove(filexyz1.c_str());
        std::remove(filexyz2.c_str());
        std::remove(filexyz3.c_str());
        std::remove(filexyz4.c_str());
        std::remove(filexyz5.c_str());

    }

    ReadFile *reader;
    std::string filexyz1 = "test1.xyz";
    std::string filexyz2 = "test2.xyz";
    std::string filexyz3 = "test3.xyz";
    std::string filexyz4 = "test4.xyz";
    std::string filexyz5 = "test5.xyz";

};

TEST_F(ReadFileTest, ReadXYZ)
{
    std::vector<std::array<double,4>> Values;

    Values = reader->ReadXYZ(filexyz1);
    EXPECT_EQ(Values[1][0], 2);
    EXPECT_EQ(Values[1][1], 1);
    EXPECT_EQ(Values[1][2], 4);
    EXPECT_EQ(Values[1][3], 5);

    EXPECT_THROW(reader->ReadXYZ(filexyz2), std::runtime_error);
    EXPECT_THROW(reader->ReadXYZ(filexyz3), std::runtime_error);
    EXPECT_THROW(reader->ReadXYZ(filexyz4), std::runtime_error);
    EXPECT_THROW(reader->ReadXYZ(filexyz5), std::runtime_error);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
