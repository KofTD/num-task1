#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "tsae.hpp"
class TDMA : public ::testing::Test {
 protected:
  using fp_type = float;

  std::vector<fp_type> a_ = {-1, -1};
  std::vector<fp_type> b_ = {-1, -1};
  std::vector<fp_type> c_ = {-4, -4};
  fp_type kappa1_ = 0.25;
  fp_type kappa2_ = 0.25;
  std::vector<fp_type> rhs_ = {0.5, -4, -6, 13.0 / 4};
};

TEST_F(TDMA, CanUseTDMA) {
  TridiagonalMatrix<fp_type> matrix(4);
  matrix.a = a_;
  matrix.b = b_;
  matrix.c = c_;
  matrix.kappa1 = kappa1_;
  matrix.kappa2 = kappa2_;
  TSAE<fp_type> system;
  system.matrix = matrix;
  system.rhs = rhs_;
  EXPECT_TRUE(canUseTdma(system));
}

TEST_F(TDMA, CalcTest) {
  TridiagonalMatrix<fp_type> matrix(4);
  matrix.a = a_;
  matrix.b = b_;
  matrix.c = c_;
  matrix.kappa1 = kappa1_;
  matrix.kappa2 = kappa2_;
  TSAE<fp_type> system;
  system.matrix = matrix;
  system.rhs = rhs_;
  auto result = tdma(system);
  std::vector<fp_type> standard = {1, 2, 3, 4};
  for (auto [lhs, rhs] : std::views::zip(result, standard)) {
    EXPECT_NEAR(lhs, rhs, std::numeric_limits<fp_type>::epsilon());
  }
}
