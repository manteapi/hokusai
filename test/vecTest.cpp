#include <iostream>
#include <hokusai/Vec.hpp>
#include <gtest/gtest.h>

using namespace hokusai;

TEST(Vec3Test, Constructor)
{
    Vec3<double> vd;
    EXPECT_EQ(0, vd[0]);
    EXPECT_EQ(0, vd[1]);
    EXPECT_EQ(0, vd[2]);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
