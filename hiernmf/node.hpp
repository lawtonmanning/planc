#include "distnmf/distr2.hpp"
#include <vector>

namespace planc {
template <class INPUTMATTYPE>
class Node {
  private:
    Node * lchild = NULL;
    Node * rchild = NULL;
    Node * parent = NULL;
    DistR2<INPUTMATTYPE> nmf = NULL;
    INPUTMATTYPE A;
    double error;
    bool activated = false;

    double compute_error() {
      return 0.0;
    }

  public:
    Node(const INPUTMATTYPE & A, std::vector<int> cols, Node * parent) {
      this->A = A.cols(cols);
      this->parent = parent;
      this->error = compute_error();
    }

    bool activate(DistR2<INPUTMATTYPE> nmf) {
      if (this->activated) {
        return true;
      }
      this->nmf = nmf;
      return true;
    }

};
}
