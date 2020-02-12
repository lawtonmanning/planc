#include "distnmf/distr2.hpp"
#include <vector>

namespace planc {
template <class INPUTMATTYPE>
class Node {
  private:
    Node * lchild = NULL;
    Node * rchild = NULL;
    Node * parent = NULL;
    INPUTMATTYPE A;
    INPUTMATTYPE W;
    double sigma;
    double score;
    bool activated = false;
    std::vector<int> cols;

    double compute_sigma() {
      return 0.0;
    }

    double compute_score() {
      if (this->lchild == NULL || this->rchild == NULL) {
        return -1.0;
      }

      return 0.0;

    }

  public:
    Node(const INPUTMATTYPE & A, std::vector<int> cols, Node * parent) {
      this->cols = cols;
      this->A = A.cols(cols);
      this->parent = parent;
      this->sigma = compute_sigma();
    }

    bool split() {
      return true;
    }

    bool accept() {
      return true;
    }


};
}
