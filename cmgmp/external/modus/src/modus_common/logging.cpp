#include <modus/common/logging.hpp>

using namespace modus;


EnumerateCSModesInputPtr
modus::CreateInput(const Eigen::MatrixXd& N,
                   const Eigen::VectorXd& d,
                   double eps)
{
    EnumerateCSModesInputPtr input(new EnumerateCSModesInput);
    input->N = N;
    input->d = d;
    input->eps = eps;

    return input;
}

EnumerateSSModesInputPtr
modus::CreateInput(const Eigen::MatrixXd& N,
                   const Eigen::VectorXd& d,
                   const Eigen::MatrixXd& T,
                   const std::string& cs_mode,
                   double eps)
{
    EnumerateSSModesInputPtr input(new EnumerateSSModesInput);
    input->N = N;
    input->d = d;
    input->T = T;
    input->cs_mode = cs_mode;
    input->eps = eps;

    return input;
}

InputPtr modus::NextInputFromLog(std::istream& is) {
    // Read input stream line by line until we reach a JSON block. The block is
    // delimited by the input name.
    std::string line;
    bool found_json_block = false;
    int type;
    while (std::getline(is, line)) {
        // std::cout << line << std::endl;
        if (line == EnumerateCSModesInput::Name()) {
            found_json_block = true;
            type = Input::ENUMERATE_CS_MODES;
            break;
        }
        else if (line == EnumerateSSModesInput::Name()) {
            found_json_block = true;
            type = Input::ENUMERATE_SS_MODES;
            break;
        }
    }
    if (found_json_block) {
        std::string out;
        while (std::getline(is, line)) {
            out += line;
            if (line == "}") {
                break;
            }
            out += "\n";
        }

        // std::cout << out << std::endl;

        std::istringstream str_stream(out, std::ios::in);
        cereal::JSONInputArchive archive(str_stream);
        if (type == Input::ENUMERATE_CS_MODES) {
            EnumerateCSModesInputPtr input(new EnumerateCSModesInput);
            archive(CEREAL_NVP(input));
            return input;
        }
        else if (type == Input::ENUMERATE_SS_MODES) {
            EnumerateSSModesInputPtr input(new EnumerateSSModesInput);
            archive(CEREAL_NVP(input));
            return input;
        }
    }
    return nullptr;
}