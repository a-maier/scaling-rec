#include "firefly/Reconstructor.hpp"

class BlackBoxUser: public firefly::BlackBoxBase<BlackBoxUser> {
public:
  BlackBoxUser(
    std::string const & file,
    std::vector<std::string> const & vars
  ): par_(file, vars) {}

  template<typename FFIntTemp>
  std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
    return par_.evaluate_pre(values);
  }

  void prime_changed() {
    par_.precompute_tokens();
  }

private:
  firefly::ShuntingYardParser par_;
};

using namespace firefly;

int main() {
  {
    std::cerr << "Reconstructing aajamp\n";
    BlackBoxUser bb{
      "../scaling-rec/data/aajamp",
      {"x12", "x23", "x34", "x45", "x51"}
    };
    firefly::Reconstructor<BlackBoxUser> rec{5, 1, bb};
    rec.reconstruct();
  }

  {
    std::cerr << "Reconstructing coeff_prop_4l\n";
    BlackBoxUser bb{
      "../scaling-rec/data/coeff_prop_4l",
      {"z", "d"}
    };
    firefly::Reconstructor<BlackBoxUser> rec{2, 1, bb};
    rec.reconstruct();
  }
}
