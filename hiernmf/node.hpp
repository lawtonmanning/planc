#include <vector>
#include "distnmf/distr2.hpp"
#include "hiernmf/matutils.hpp"

namespace planc {
template <class INPUTMATTYPE>
class Node {
  protected:
    Node * lchild = NULL;
    Node * rchild = NULL;
    Node * parent = NULL;
    INPUTMATTYPE A0;
    INPUTMATTYPE A;
    INPUTMATTYPE W;
    double sigma;
    double score;
    bool activated = false;
    UVEC cols;

    void allocate_A() {
#ifdef BUILD_SPARSE
      int n_cols = this->cols.n_elem;

      arma::umat locs(2,n_cols);
      for (int i = 0; i < n_cols; i++) {
        locs(0,i) = i;
        locs(1,i) = this->cols(i);
      }
      VEC vals(n_cols,1);

      //SP_MAT S(locs,vals,this->A0.n_rows,n_cols);

      this->A = this->A0;
#else
      this->A = this->A0.cols(this->cols);
#endif
      print(this->A,"A");
    }

    void compute_sigma() {
      this->sigma = 0.0;
    }

    double compute_score() {
      if (this->lchild == NULL || this->rchild == NULL) {
        return -1.0;
      }

      return 0.0;

    }

  public:
    Node(INPUTMATTYPE & A, UVEC & cols, Node * parent) {
      this->cols = cols;
      this->A0 = A;
      this->parent = parent;
      this->allocate_A();
      this->compute_sigma();
    }

    bool split() {
      return true;
    }

    bool accept() {
      return true;
    }


};

template <class INPUTMATTYPE>
class RootNode : public Node<INPUTMATTYPE> {
  protected:
    void allocate_A() {
      this->A = this->A0;
    }

    void compute_sigma() {
      this->sigma = 0.0;
    }

  public:
    RootNode(INPUTMATTYPE & A, UVEC & cols) : Node<INPUTMATTYPE>(A, cols, NULL) {
    }
};
}
