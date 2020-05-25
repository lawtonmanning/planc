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
        VEC W;
        double sigma;
        double score;
        UVEC cols;
        bool accepted = false;
        int index;

        UVEC top_words;

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

        void compute_top_words() {
#ifdef BUILD_SPARSE
          int k = this->pc->words();
#else
          int k = this->pc->globalm();
#endif
          printf("words %d\n", k);
          if (k == 0) {
            return;
          }
          int p = this->mpicomm->size();
          VEC locWm = maxk(W, k);
          UVEC locWi = maxk_idx(W, k) + startidx(this->pc->globalm(), p, this->mpicomm->rank()) - 1;
          locWi.t().print("indices");
          printf("1\n");

          int * kcounts = (int *)malloc(p*sizeof(int));
          printf("1.1\n");
          int n = (int)(locWm.n_elem);
          MPI_Allgather(&n, 1, MPI_INT, kcounts, p, MPI_INT, MPI_COMM_WORLD);
          printf("1.2\n");
          int ktotal = 0;
          for (int i = 0; i < p; i++) {
            ktotal += kcounts[i];
          }
          printf("2\n");

          int * kdispls = (int *)malloc(p*sizeof(int));
          kdispls[0] = 0;
          for (int i = 1; i < p; i++)
          {
            kdispls[i] = kdispls[i - 1] + kcounts[i - 1];
          }
          printf("3\n");

          VEC gloWm = arma::zeros<VEC>(ktotal);
          MPI_Allgatherv(locWm.memptr(), locWm.n_elem, MPI_DOUBLE, gloWm.memptr(), kcounts, kdispls, MPI_DOUBLE, MPI_COMM_WORLD);

          UVEC gloWi = arma::zeros<UVEC>(ktotal);
          MPI_Allgatherv(locWi.memptr(), locWi.n_elem, MPI_UNSIGNED_LONG_LONG, gloWi.memptr(), kcounts, kdispls, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
          printf("4\n");
          //VEC Wm = maxk(gloWm,k);
          UVEC Wmi = maxk_idx(gloWm, k);
          UVEC Wi = gloWi.elem(Wmi);
          printf("5\n");

          this->top_words = Wi;
          printf("done computing words\n");
        }

        Node() {
        }

        Node(INPUTMATTYPE & A, VEC W, UVEC & cols, Node * parent, int index) {
          this->cols = cols;
          this->A0 = A;
          this->W = W;
          this->parent = parent;
          this->index = index;
          this->mpicomm = parent->mpicomm;
          this->pc = parent->pc;
          this->allocate();
          this->compute_sigma();
          this->compute_top_words();
          this->lvalid = false;
          this->rvalid = false;
        }

        void split() {
          this->accepted = true;

          MAT W = arma::randu<MAT>(itersplit(A.n_rows,pc->pc(),mpicomm->col_rank()),2);
          MAT H = arma::randu<MAT>(itersplit(A.n_cols,pc->pr(),mpicomm->row_rank()),2);
          /*
          if (mpicomm->rank() == 0) {
            W.eye();
            H.eye();
          }
          else {
            W.zeros();
            H.zeros();
          }*/

          DistR2<INPUTMATTYPE> nmf(A, W, H, *mpicomm, 1);;
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

          printf("node %d cols -- total:%d  left:%d  right:%d\n",this->index,this->cols.n_elem,arma::sum(lleft == 1),arma::sum(lleft == 0));

          this->lvalid = !lcols.is_empty();
          this->rvalid = !rcols.is_empty();
          
          if (this->lvalid) {
            this->lchild = new Node(this->A0, W.col(0), lcols, this, 2*this->index+1);
            printf("node(%d,%d) %f\n",lchild->index,mpicomm->rank(),lchild->sigma);
          }
          if (this->rvalid) {
            this->rchild = new Node(this->A0, W.col(1), rcols, this, 2*this->index+2);
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
          printf("size %d\n",queue.size());
          if (this->lvalid) {
            queue.push(this->lchild);
            printf("size %d\n",queue.size());
          }
          if (this->rvalid) {
            queue.push(this->rchild);
            printf("size %d\n",queue.size());
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
          this->sigma = 0.0;
          this->index = 0;
          this->lvalid = false;
          this->rvalid = false;
        }
    };
}
