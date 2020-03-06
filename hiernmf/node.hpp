#include <queue>
#include <vector>
#include "common/parsecommandline.hpp"
#include "distnmf/distr2.hpp"
#include "distnmf/mpicomm.hpp"
#include "hiernmf/matutils.hpp"

namespace planc {
  template <class INPUTMATTYPE>
    class Node {
      public:  
        Node * lchild = NULL;
        Node * rchild = NULL;
        bool lvalid,rvalid;
        Node * parent = NULL;
        INPUTMATTYPE A0;
        INPUTMATTYPE A;
        MAT W;
        MAT H;
        double sigma;
        double score;
        UVEC cols;
        bool accepted = false;
        int index;

        MPICommunicator * mpicomm;
        ParseCommandLine * pc;

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
        }

        void compute_sigma() {
          this->sigma = powIter(this->A);
        }

        void compute_score() {
          if (this->lvalid && this->rvalid) {
            this->score = (this->lchild->sigma+this->rchild->sigma) - this->sigma;
          }
          else if (this->lvalid) {
            this->score = this->lchild->sigma - this->sigma;
          }
          else if (this->rvalid) {
            this->score = this->rchild->sigma - this->sigma;
          }
          else {
            this->score = 0;
          }
        }

        Node() {
        }

        Node(INPUTMATTYPE & A, UVEC & cols, Node * parent) {
          this->cols = cols;
          this->A0 = A;
          this->parent = parent;
          this->mpicomm = parent->mpicomm;
          this->pc = parent->pc;
          this->allocate();
          this->compute_sigma();
        }

        void split() {
          MPI_Barrier(MPI_COMM_WORLD);
          char * title = (char *)malloc(20);
          sprintf(title, "cols(%d,%d):", this->index, mpicomm->rank());
          printf("A(%d): %dx%d\n",mpicomm->rank(),A.n_rows,A.n_cols);
          cols.t().print(title);
          free(title);
          MPI_Barrier(MPI_COMM_WORLD);
          MAT Wo = arma::randu<MAT>(itersplit(A.n_rows,pc->pc(),mpicomm->col_rank()),2);
          MAT Ho = arma::randu<MAT>(itersplit(A.n_cols,pc->pr(),mpicomm->row_rank()),2);

          DistR2<INPUTMATTYPE> nmf(A, Wo, Ho, *mpicomm, 1);;
          nmf.num_iterations(pc->iterations());
          nmf.compute_error(pc->compute_error());
          nmf.algorithm(R2);
          nmf.regW(pc->regW());
          nmf.regH(pc->regH());
          MPI_Barrier(MPI_COMM_WORLD);
          try {
            mpitic();
            nmf.computeNMF();
            double temp = mpitoc();

            if (this->mpicomm->rank() == 0) printf("NMF took %.3lf secs.\n", temp);
          } catch (std::exception &e) {
            printf("Failed rank %d: %s\n", this->mpicomm->rank(), e.what());
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
          this->W = nmf.getLeftLowRankFactor();
          this->H = nmf.getRightLowRankFactor();
          /*
          Ho = nmf.getRightLowRankFactor();
          this->H.zeros(2,A.n_cols);
          int sendcnt = Ho.n_rows*Ho.n_cols;
          int * recvcnts = (int *)malloc(this->mpicomm->size()*sizeof(int));
          int * displs = (int *)malloc(this->mpicomm->size()*sizeof(int));
          recvcnts[0] = itersplit(A.n_cols,this->mpicomm->size(),0)*2;
          displs[0] = 0;
          for (int i = 0; i < this->mpicomm->size(); i++) {
            recvcnts[i] = itersplit(A.n_cols,this->mpicomm->size(),i)*2;
            displs[i] = displs[i-1]+recvcnts[i-1];
          }
          MPI_Allgatherv(Ho.memptr(),recvcnts[this->mpicomm->rank()],MPI_DOUBLE,H.memptr(),recvcnts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
          this->H = this->H.t();
          */

          /*
          MPI_Barrier(MPI_COMM_WORLD);
          title = (char *)malloc(20);
          sprintf(title, "H(%d): %dx%d", mpicomm->rank(), H.n_rows, H.n_cols);
          H.print(title);
          free(title);
          MPI_Barrier(MPI_COMM_WORLD);
          */


          UVEC left = this->H.col(0) > this->H.col(1);

          UVEC lcols = this->cols(find(left == 1));
          UVEC rcols = this->cols(find(left == 0));

          
          this->lvalid = !lcols.is_empty();
          this->rvalid = !rcols.is_empty();
          
          if (lvalid) {
            this->lchild = new Node(this->A0, lcols, this);
            this->lchild->index = 2*this->index+1;
            printf("node(%d,%d) %f\n",lchild->index,mpicomm->rank(),lchild->sigma);
          }
          if (rvalid) {
            this->rchild = new Node(this->A0, rcols, this);
            this->rchild->index = 2*this->index+2;
            printf("node(%d,%d) %f\n",rchild->index,mpicomm->rank(),rchild->sigma);
          }

          this->compute_score();
        }

        void accept() {
          this->accepted = true;
          if (this->lvalid) {
            this->lchild->split();
          }
          if (this->rvalid) {
            this->rchild->split();
          }
        }
        
        template<class QUEUE>
        void enqueue(QUEUE & queue) {
          if (this->lvalid) {
            queue.push(this->lchild);
          }
          if (this->rvalid) {
            queue.push(this->rchild);
          }
        }

        bool operator > (const Node<INPUTMATTYPE> & rhs) const {
          return (this->score > rhs.score);
        }

        bool operator < (const Node<INPUTMATTYPE> & rhs) const {
          return (this->score < rhs.score);
        }
    };

  class ScoreCompare {
    public:
    template <typename T>
    bool operator()(T * a, T * b) {
      return a->score > b->score;
    }
  };


  template <class INPUTMATTYPE>
    class RootNode : public Node<INPUTMATTYPE> {
      public:  
        RootNode(INPUTMATTYPE & A, UVEC & cols, MPICommunicator * mpicomm, ParseCommandLine * pc) : Node<INPUTMATTYPE>() {
          this->cols = cols;
          this->A0 = A;
          this->parent = NULL;
          this->mpicomm = mpicomm;
          this->pc = pc;
          this->A = this->A0;
          this->sigma = 0.0;
          this->index = 0;
        }
    };
}
