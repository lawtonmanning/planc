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

#ifdef BUILD_SPARSE
    RootNode<SP_MAT> * root;
#else
    RootNode<MAT> * root;
#endif

    void parseCommandLine() {
      ParseCommandLine pc(this->m_argc, this->m_argv);
      pc.parseplancopts();
      this->m_k = pc.lowrankk();
      this->m_Afile_name = pc.input_file_name();
      this->m_pr = pc.pr();
      this->m_pc = 1;
      this->m_sparsity = pc.sparsity();
      this->m_num_it = pc.iterations();
      this->m_distio = TWOD;
      this->m_regW = pc.regW();
      this->m_regH = pc.regH();
      this->m_num_k_blocks = 1;
      this->m_globalm = pc.globalm();
      this->m_globaln = pc.globaln();
      this->m_compute_error = pc.compute_error();
      this->m_distio = TWOD;
      this->m_input_normalization = pc.input_normalization();
      pc.printConfig();
    }

    void buildTree() {
      std::string rand_prefix("rand_");
      MPICommunicator mpicomm(this->m_argc, this->m_argv, this->m_pr, this->m_pc);

#ifdef BUILD_SPARSE    
      DistIO<SP_MAT> dio(mpicomm, m_distio);
#else 
      DistIO<MAT> dio(mpicomm, m_distio);
#endif

      if (m_Afile_name.compare(0, rand_prefix.size(), rand_prefix) == 0) {
        dio.readInput(m_Afile_name, this->m_globalm, this->m_globaln, this->m_k,
            this->m_sparsity, this->m_pr, this->m_pc,
            this->m_input_normalization);
      } else {
        dio.readInput(m_Afile_name);
      }

#ifdef BUILD_SPARSE
      SP_MAT A(dio.A());
#else 
      MAT A(dio.A());
#endif
      
      arma::uvec cols(this->m_globaln);
      for (unsigned int i = 0; i < this->m_globaln; i++) {
        cols[i] = i;
      }

      print(A,"A0");

#ifdef BUILD_SPARSE
      this->root = new RootNode<SP_MAT>(A,cols);
#else
      this->root = new RootNode<MAT>(A,cols);
#endif

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

