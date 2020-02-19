#include <string>
#include "common/distutils.hpp"
#include "common/parsecommandline.hpp"
#include "common/utils.hpp"
#include "distnmf/distr2.hpp"
#include "distnmf/mpicomm.hpp"

using namespace planc;

class HierNMFDriver {
  private:
    int m_argc;
    char **m_argv;
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

    void parseCommandLine() {
      ParseCommandLine pc(this->m_argc, this->m_argv);
      pc.parseplancopts();
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
      MPICommunicator mpicomm(this->m_argc, this->m_argv, this->m_pr, this->m_pc);

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

