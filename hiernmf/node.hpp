#include <vector>
#include "distnmf/distr2.hpp"
#include "hiernmf/matutils.hpp"

namespace planc {
  template <class INPUTMATTYPE>
    class Node {
      public:  
        Node * lchild = NULL;
        Node * rchild = NULL;
        Node * parent = NULL;
        INPUTMATTYPE A0;
        INPUTMATTYPE A;
        MAT W;
        MAT H;
        double sigma;
        double score;
        bool activated = false;
        UVEC cols;

        void allocate() {
#ifdef BUILD_SPARSE
          int n_cols = this->cols.n_elem;

          arma::umat locs(2,n_cols);
          for (int i = 0; i < n_cols; i++) {
            locs(1,i) = i;
            locs(0,i) = this->cols(i);
          }
          VEC vals(n_cols);
          vals.fill(1);

          SP_MAT S(locs,vals,this->A0.n_cols,n_cols);

          this->A = this->A0*S;
#else
          this->A = this->A0.cols(this->cols);
#endif
          this->W = arma::randu<MAT>(this->A.n_rows,2);
          this->H = arma::randu<MAT>(this->A.n_cols,2);
        }

        void compute_sigma() {
          this->sigma = powIter(this->A);
        }

        double compute_score() {
          if (this->lchild == NULL || this->rchild == NULL) {
            return -1.0;
          }

          return 0.0;

        }

        Node() {
        }

        Node(INPUTMATTYPE & A, UVEC & cols, Node * parent) {
          this->cols = cols;
          this->A0 = A;
          this->parent = parent;
          this->allocate();
          this->compute_sigma();
        }

        bool split() {
          print(this->H, "H");
          UVEC left = this->H.col(0) > this->H.col(1);
          print(left, "left");
          UVEC lcols = this->cols(find(left == 1));
          UVEC rcols = this->cols(find(left == 0));
          print(this->cols, "cols");
          print(lcols, "lcols");
          print(rcols, "rcols");

          return true;
        }

        bool accept() {
          return true;
        }


    };

  template <class INPUTMATTYPE>
    class RootNode : public Node<INPUTMATTYPE> {
      public:  
        void allocate() {
          this->A = this->A0;
          this->W = arma::randu<MAT>(this->A.n_rows,2);
          this->H = arma::randu<MAT>(this->A.n_cols,2);
        }

        void compute_sigma() {
          this->sigma = 0.0;
        }

        RootNode(INPUTMATTYPE & A, UVEC & cols) : Node<INPUTMATTYPE>() {
          this->cols = cols;
          this->A0 = A;
          this->parent = NULL;
          this->allocate();
          this->compute_sigma();
        }
    };
}
