#include <string>
#include <queue>
#include "common/distutils.hpp"
#include "common/parsecommandline.hpp"
#include "common/utils.hpp"
#include "distnmf/distr2.hpp"
#include "distnmf/distio.hpp"
#include "distnmf/mpicomm.hpp"
#include "hiernmf/node.hpp"

using namespace planc;

class HierNMFDriver {
  private:
    int m_argc;
    char **m_argv;
    int m_k;
    UWORD m_globalm, m_globaln;
    std::string m_Afile_name;
    std::string m_outputfile_name;
    int m_num_it;
    int m_pr;
    int m_pc = 1;
    FVEC m_regW;
    FVEC m_regH;
    double m_sparsity;
    iodistributions m_distio;
    uint m_compute_error;
    int m_num_k_blocks;
    static const int kprimeoffset = 17;
    normtype m_input_normalization;
    MPICommunicator * mpicomm;
    ParseCommandLine * pc;

#ifdef BUILD_SPARSE
    RootNode<SP_MAT> * root;
#else
    RootNode<MAT> * root;
#endif

    template <class INPUTMATTYPE>
      bool split(Node<INPUTMATTYPE> * node) {
        MAT W = arma::randu<MAT>(node->A.n_rows/this->m_pc,2);
        MAT H = arma::randu<MAT>(node->A.n_cols/this->m_pr,2);
        DistR2<INPUTMATTYPE> nmf(node->A, W, H, *this->mpicomm, 1);
        printf("Created NMF\n");
        nmf.num_iterations(this->m_num_it);
        nmf.compute_error(this->m_compute_error);
        nmf.regW(this->m_regW);
        nmf.regH(this->m_regH);
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
        node->W = nmf.getLeftLowRankFactor();
        MAT Ho = nmf.getRightLowRankFactor();
        H.zeros(2,node->A.n_cols);
        int sendcnt = Ho.n_rows*Ho.n_cols;
        int * recvcnts = (int *)malloc(this->mpicomm->size()*sizeof(int));
        int * displs = (int *)malloc(this->mpicomm->size()*sizeof(int));
        recvcnts[0] = itersplit(node->A.n_cols,this->mpicomm->size(),0)*2;
        displs[0] = 0;
        for (int i = 0; i < this->mpicomm->size(); i++) {
          recvcnts[i] = itersplit(node->A.n_cols,this->mpicomm->size(),i)*2;
          displs[i] = displs[i-1]+recvcnts[i-1];
        }
        MPI_Allgatherv(Ho.memptr(),recvcnts[this->mpicomm->rank()],MPI_DOUBLE,H.memptr(),recvcnts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
        node->H = H.t();
        return node->split();
      }


    void parseCommandLine() {
      pc = new ParseCommandLine(this->m_argc, this->m_argv);
      pc->parseplancopts();
      this->m_k = pc->lowrankk();
      this->m_Afile_name = pc->input_file_name();
      this->m_pr = pc->pr();
      this->m_pc = 1;
      this->m_sparsity = pc->sparsity();
      this->m_num_it = pc->iterations();
      this->m_distio = TWOD;
      this->m_regW = pc->regW();
      this->m_regH = pc->regH();
      this->m_num_k_blocks = 1;
      this->m_globalm = pc->globalm();
      this->m_globaln = pc->globaln();
      this->m_compute_error = pc->compute_error();
      this->m_input_normalization = pc->input_normalization();
      pc->printConfig();
    }

    void buildTree() {
      std::string rand_prefix("rand_");
      this->mpicomm = new MPICommunicator(this->m_argc, this->m_argv, this->m_pr, this->m_pc);

#ifdef BUILD_SPARSE    
      DistIO<SP_MAT> dio(*mpicomm, m_distio);
#else 
      DistIO<MAT> dio(*mpicomm, m_distio);
#endif

      dio.readInput(m_Afile_name, this->m_globalm, this->m_globaln, this->m_k,
          this->m_sparsity, this->m_pr, this->m_pc,
          this->m_input_normalization);

#ifdef BUILD_SPARSE
      SP_MAT A(dio.A());
#else 
      MAT A(dio.A());
#endif
      int rank = this->mpicomm->rank();
      arma::uvec cols(this->m_globaln);
      for (unsigned int i = 0; i < this->m_globaln; i++) {
        cols[i] = i;
      }

#ifdef BUILD_SPARSE
      this->root = new RootNode<SP_MAT>(A, cols, this->mpicomm, this->pc);
#else
      this->root = new RootNode<MAT>(A, cols, this->mpicomm, this->pc);
#endif

#ifdef BUILD_SPARSE
      std::priority_queue<Node<SP_MAT> *> leaves;
      Node<SP_MAT> * leaf;

      std::queue<Node<SP_MAT> *> nodes;
#else
      std::priority_queue<Node<MAT> *> leaves;
      Node<MAT> * leaf;

      std::queue<Node<MAT> *> nodes;
#endif
      nodes.push(this->root);
      this->root->split();
      this->root->accept();
      this->root->enqueue(leaves);
      this->root->enqueue(nodes);
      

      while (leaves.size() < this->m_k) {
        leaf = leaves.top();
        leaves.pop();
        leaf->accept();
        leaf->enqueue(leaves);
        leaf->enqueue(nodes);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (this->mpicomm->rank() == 0) {
        
      }
      MPI_Barrier(MPI_COMM_WORLD);



      delete this->mpicomm;
    }
  public:
    HierNMFDriver(int argc, char * argv[]) {
      this->m_argc = argc;
      this->m_argv = argv;
      this->parseCommandLine();
      this->buildTree();
    }
};

int main(int argc, char * argv[]) {
  HierNMFDriver(argc,argv);
}

