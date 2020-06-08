#include <queue>
#include <vector>
#include "common/parsecommandline.hpp"
#include "distnmf/distr2.hpp"
#include "distnmf/mpicomm.hpp"
#include "hiernmf/matutils.hpp"

namespace planc {
  struct NodeTimings {
    double NMF;
    double sigma;
    double top_words;
  };
  template <class INPUTMATTYPE>
    class Node {
      public:  
        Node * lchild = NULL;
        Node * rchild = NULL;
        bool lvalid,rvalid;
        Node * parent = NULL;
        INPUTMATTYPE A0;
        INPUTMATTYPE A;
        VEC W;
        double sigma;
        double score;
        UVEC cols;
        bool accepted = false;
        int index;

        int global_m, global_n;

        UVEC top_words;

        MPICommunicator * mpicomm;
        ParseCommandLine * pc;

        NodeTimings timings;

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
          this->sigma = powIter(this->A,this->pc->iterations(),this->pc->tolerance());
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

        void compute_top_words() {
#ifdef BUILD_SPARSE
          int k = this->pc->words();
#else
          int k = this->pc->globalm();
#endif
          if (k == 0) {
            return;
          }
          int p = this->mpicomm->size();
          VEC locWm = maxk(W, k);
          UVEC locWi = maxk_idx(W, k) + startidx(this->global_m, p, this->mpicomm->rank());

          int * kcounts = (int *)malloc(p*sizeof(int));
          int n = (int)(locWm.n_elem);
          MPI_Allgather(&n, 1, MPI_INT, kcounts, 1, MPI_INT, MPI_COMM_WORLD);
          int ktotal = 0;
          for (int i = 0; i < p; i++) {
            ktotal += kcounts[i];
          }

          int * kdispls = (int *)malloc(p*sizeof(int));
          kdispls[0] = 0;
          for (int i = 1; i < p; i++)
          {
            kdispls[i] = kdispls[i - 1] + kcounts[i - 1];
          }

          VEC gloWm = arma::zeros<VEC>(ktotal);
          MPI_Allgatherv(locWm.memptr(), locWm.n_elem, MPI_DOUBLE, gloWm.memptr(), kcounts, kdispls, MPI_DOUBLE, MPI_COMM_WORLD);


          UVEC gloWi = arma::zeros<UVEC>(ktotal);
          MPI_Allgatherv(locWi.memptr(), locWi.n_elem, MPI_UNSIGNED_LONG_LONG, gloWi.memptr(), kcounts, kdispls, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

          UVEC Wmi = maxk_idx(gloWm, k);
          UVEC Wi = gloWi.elem(Wmi);

          this->top_words = Wi;
        }

        Node() {
        }

        Node(INPUTMATTYPE & A, VEC W, UVEC & cols, Node * parent, int index) {
          this->cols = cols;
          this->A0 = A;
          this->global_m = parent->global_m;
          this->global_n = parent->global_n;
          this->W = W;
          this->parent = parent;
          this->index = index;
          this->mpicomm = parent->mpicomm;
          this->pc = parent->pc;
          this->allocate();
          mpitic();
          this->compute_sigma();
          this->timings.sigma = mpitoc();
          mpitic();
          this->compute_top_words();
          this->timings.top_words = mpitoc();
          this->lvalid = false;
          this->rvalid = false;
        }

        void split() {
          this->accepted = true;

          arma::arma_rng::set_seed_random();
          MAT W = arma::randu<MAT>(itersplit(A.n_rows,pc->pc(),mpicomm->col_rank()),2);
          MAT H = arma::randu<MAT>(itersplit(A.n_cols,pc->pr(),mpicomm->row_rank()),2);

          DistR2<INPUTMATTYPE> nmf(A, W, H, *mpicomm, 1);;
          nmf.num_iterations(pc->iterations());
          nmf.compute_error(pc->compute_error());
          nmf.tolerance(pc->tolerance());
          nmf.algorithm(R2);
          nmf.regW(pc->regW());
          nmf.regH(pc->regH());

          mpitic();
          nmf.computeNMF();
          this->timings.NMF = mpitoc();

          W = nmf.getLeftLowRankFactor();
          H = nmf.getRightLowRankFactor();

          UVEC lleft = H.col(0) > H.col(1);
          UVEC left(A.n_cols,arma::fill::zeros);
          int * recvcnts = (int *)malloc(this->mpicomm->size()*sizeof(int));
          int * displs = (int *)malloc(this->mpicomm->size()*sizeof(int));
          recvcnts[0] = itersplit(A.n_cols,pc->pr(),0);
          displs[0] = 0;
          for (int i = 1; i < this->mpicomm->size(); i++) {
            recvcnts[i] = itersplit(A.n_cols,pc->pr(),i);
            displs[i] = displs[i-1]+recvcnts[i-1];
          }
          MPI_Allgatherv(lleft.memptr(),lleft.n_elem,MPI_UNSIGNED_LONG_LONG,left.memptr(),recvcnts,displs,MPI_UNSIGNED_LONG_LONG,MPI_COMM_WORLD);


          UVEC lcols = this->cols(find(left == 1));
          UVEC rcols = this->cols(find(left == 0));

          this->lvalid = !lcols.is_empty();
          this->rvalid = !rcols.is_empty();
          
          if (this->lvalid) {
            this->lchild = new Node(this->A0, W.col(0), lcols, this, 2*this->index+1);
          }
          if (this->rvalid) {
            this->rchild = new Node(this->A0, W.col(1), rcols, this, 2*this->index+2);
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
        RootNode(INPUTMATTYPE & A, int global_m, int global_n, UVEC & cols, MPICommunicator * mpicomm, ParseCommandLine * pc) : Node<INPUTMATTYPE>() {
          this->cols = cols;
          this->A0 = A;
          this->global_m = global_m;
          this->global_n = global_n;
          this->parent = NULL;
          this->mpicomm = mpicomm;
          this->pc = pc;
          this->A = this->A0;
          this->sigma = 0.0;
          this->index = 0;
          this->lvalid = false;
          this->rvalid = false;
        }
    };
}
