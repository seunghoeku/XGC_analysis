void load_init(const std::string &filename);
void load_finalize();

adios2::StepStatus load_data(std::vector<long> &igid, std::vector<int> &iflag,
                             std::vector<float> &idw,
                             std::vector<float> &iphase,
                             std::vector<long> &egid, std::vector<int> &eflag,
                             std::vector<float> &edw,
                             std::vector<float> &ephase);