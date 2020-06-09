#include <string>
#include <queue>
#include <vector>
#include <fstream>
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


    void parseCommandLine() {
      pc = new ParseCommandLine(this->m_argc, this->m_argv);
      pc->parseplancopts();
      this->m_k = 2;
      this->m_Afile_name = pc->input_file_name();
      this->m_outputfile_name = pc->output_file_name();
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

      if (m_Afile_name.compare(0, rand_prefix.size(), rand_prefix) != 0) {
        int localm = A.n_rows;
        int localn = A.n_cols;
        int m,n;
        MPI_Allreduce(&localm,&m,1,MPI_INT,MPI_SUM,mpicomm->commSubs()[0]);
        MPI_Allreduce(&localn,&n,1,MPI_INT,MPI_SUM,mpicomm->commSubs()[1]);
        this->m_globalm = m;
        this->m_globaln = n;
      }

      int rank = this->mpicomm->rank();
      arma::uvec cols(this->m_globaln);
      for (unsigned int i = 0; i < this->m_globaln; i++) {
        cols[i] = i;
      }

#ifdef BUILD_SPARSE
      std::priority_queue<Node<SP_MAT> *, std::vector<Node<SP_MAT> *>, ScoreCompare> frontiers;
      Node<SP_MAT> * frontier;
      Node<SP_MAT> * node;

      std::queue<Node<SP_MAT> *> nodes;
#else
      std::priority_queue<Node<MAT> *, std::vector<Node<MAT> *>, ScoreCompare> frontiers;
      Node<MAT> * frontier;
      Node<MAT> * node;

      std::queue<Node<MAT> *> nodes;
#endif

      mpitic();

#ifdef BUILD_SPARSE
      this->root = new RootNode<SP_MAT>(A, this->m_globalm, this->m_globaln, cols, this->mpicomm, this->pc);
#else
      this->root = new RootNode<MAT>(A, this->m_globalm, this->m_globaln, cols, this->mpicomm, this->pc);
#endif

      nodes.push(root);

      this->root->split();
      this->root->accept();
      this->root->enqueue(frontiers);
      this->root->enqueue(nodes);


      for (int i = 2; i < pc->nodes(); i++) {
        if (frontiers.empty()) {
          break;
        }
        frontier = frontiers.top();
        frontiers.pop();
        frontier->accept();
        frontier->enqueue(frontiers);
        frontier->enqueue(nodes);
      }

      double temp = mpitoc();
      if (this->mpicomm->rank() == 0) {
        printf("H2NMF took %.3lf secs.\n", temp);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (this->mpicomm->rank() == 0) {
        std::ofstream ofs(this->m_outputfile_name, std::ofstream::out);
        
        printf("idx NMF sigma top_words\n");
        while (!nodes.empty()) {
          node = nodes.front();
          nodes.pop();
          
          printf("%d %f %f %f\n",node->index,node->timings.NMF,node->timings.sigma,node->timings.top_words);
          
          if (node->index == 0) {
            continue;
          }
          ofs << node->index << " ";
#ifdef BUILD_SPARSE
          ofs << node->top_words.t();
#else
          ofs << node->cols.t();
#endif
        }

        ofs.close();
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

