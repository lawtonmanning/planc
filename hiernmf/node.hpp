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
          this->W.zeros(this->A.n_rows,2);
          this->H.zeros(this->A.n_cols,2);
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
          this->accepted = true;
          
          arma::arma_rng::set_seed_random();
          this->W = arma::randu<MAT>(itersplit(A.n_rows,pc->pc(),mpicomm->col_rank()),2);

          arma::arma_rng::set_seed_random();
          this->H = arma::randu<MAT>(itersplit(A.n_cols,pc->pr(),mpicomm->row_rank()),2);
          

          DistR2<INPUTMATTYPE> nmf(A, this->W, this->H, *mpicomm, 1);;
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
          
          
          UVEC lleft = this->H.col(0) > this->H.col(1);
          UVEC left(this->cols.n_elem,arma::fill::zeros);
          int * recvcnts = (int *)malloc(this->mpicomm->size()*sizeof(int));
          int * displs = (int *)malloc(this->mpicomm->size()*sizeof(int));
          recvcnts[0] = itersplit(A.n_cols,pc->pr(),0);
          displs[0] = 0;
          for (int i = 0; i < this->mpicomm->size(); i++) {
            recvcnts[i] = itersplit(A.n_cols,pc->pr(),i);
            displs[i] = displs[i-1]+recvcnts[i-1];
          }
          MPI_Allgatherv(lleft.memptr(),lleft.n_elem,MPI_UNSIGNED_LONG_LONG,left.memptr(),recvcnts,displs,MPI_UNSIGNED_LONG_LONG,MPI_COMM_WORLD);


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
          if (this->lvalid) {
            this->lchild->split();
          }
          if (this->rvalid) {
            this->rchild->split();
          }
        }
        
        template<class QUEUE>
        // TODO: change enqueue to enqueue_children
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
      return a->score < b->score;
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
          this->W.zeros(this->A.n_rows,2);
          this->H.zeros(this->A.n_cols,2);
          this->sigma = 0.0;
          this->index = 0;
        }
    };
}
