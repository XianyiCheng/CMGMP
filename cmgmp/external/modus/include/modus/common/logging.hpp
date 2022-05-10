#pragma once
#include <memory>
#include <sstream>
// #include <glog/logging.h>
#include <modus/common/serialization.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/memory.hpp>


namespace modus
{

struct Input {
    enum {
        ENUMERATE_CS_MODES = 0,
        ENUMERATE_SS_MODES = 1,
    };
    int type;

    virtual ~Input() = default;

    template <typename Archive>
    void save(Archive& ar) const {
    }

    template <typename Archive>
    void load(Archive& ar) {
    }
};

using InputPtr = std::shared_ptr<Input>;

struct EnumerateCSModesInput : public Input {
    Eigen::MatrixXd N;
    Eigen::VectorXd d;
    double eps;
    // ModeEnumerationOptions options;

    EnumerateCSModesInput() {
        type = ENUMERATE_CS_MODES;
    }

    static std::string Name() {
        return "Enumerate CS Modes Input";
    }

    template <typename Archive>
    void save(Archive& ar) const {
        ar(CEREAL_NVP(N));
        ar(CEREAL_NVP(d));
        ar(CEREAL_NVP(eps));
        ar(CEREAL_NVP(type));
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(CEREAL_NVP(N));
        ar(CEREAL_NVP(d));
        ar(CEREAL_NVP(eps));
        ar(CEREAL_NVP(type));
    }
};

using EnumerateCSModesInputPtr = std::shared_ptr<EnumerateCSModesInput>;

struct EnumerateSSModesInput : public Input {
    Eigen::MatrixXd N;
    Eigen::VectorXd d;
    Eigen::MatrixXd T;
    std::string cs_mode;
    double eps;
    // ModeEnumerationOptions options;

    EnumerateSSModesInput() {
        type = ENUMERATE_SS_MODES;
    }

    static std::string Name() {
        return "Enumerate SS Modes Input";
    }

    template <typename Archive>
    void save(Archive& ar) const {
        ar(CEREAL_NVP(N));
        ar(CEREAL_NVP(d));
        ar(CEREAL_NVP(T));
        ar(CEREAL_NVP(cs_mode));
        ar(CEREAL_NVP(eps));
        ar(CEREAL_NVP(type));
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(CEREAL_NVP(N));
        ar(CEREAL_NVP(d));
        ar(CEREAL_NVP(T));
        ar(CEREAL_NVP(cs_mode));
        ar(CEREAL_NVP(eps));
        ar(CEREAL_NVP(type));
    }
};

using EnumerateSSModesInputPtr = std::shared_ptr<EnumerateSSModesInput>;

EnumerateCSModesInputPtr CreateInput(const Eigen::MatrixXd& N,
                                     const Eigen::VectorXd& d,
                                     double eps);
EnumerateSSModesInputPtr CreateInput(const Eigen::MatrixXd& N,
                                     const Eigen::VectorXd& d,
                                     const Eigen::MatrixXd& T,
                                     const std::string& cs_mode,
                                     double eps);

template <typename InputPtr>
// void WriteInputToLog(const InputPtr input) {
//     std::ostringstream strstream;
//     {
//         cereal::JSONOutputArchive archive(strstream);
//         archive(CEREAL_NVP(input));
//     }
//     LOG(INFO) << "\n" << input->Name() << "\n" << strstream.str();
//     LOG(INFO) << "EOI";
//     google::FlushLogFiles(google::INFO);
// }

InputPtr NextInputFromLog(std::istream& is);


}