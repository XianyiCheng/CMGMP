#pragma once
#include <vector>
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <modus/common/eigen.hpp>

namespace cereal
{
//   template <class Archive, class Derived> inline
//     typename std::enable_if<traits::is_output_serializable<BinaryData<typename Derived::Scalar>, Archive>::value, void>::type
//     save(Archive & ar, Eigen::PlainObjectBase<Derived> const & m){
//       typedef Eigen::PlainObjectBase<Derived> ArrT;
//       if(ArrT::RowsAtCompileTime==Eigen::Dynamic) ar(m.rows());
//       if(ArrT::ColsAtCompileTime==Eigen::Dynamic) ar(m.cols());
//       ar(binary_data(m.data(),m.size()*sizeof(typename Derived::Scalar)));
//     }

//   template <class Archive, class Derived> inline
//     typename std::enable_if<traits::is_input_serializable<BinaryData<typename Derived::Scalar>, Archive>::value, void>::type
//     load(Archive & ar, Eigen::PlainObjectBase<Derived> & m){
//       typedef Eigen::PlainObjectBase<Derived> ArrT;
//       Eigen::Index rows=ArrT::RowsAtCompileTime, cols=ArrT::ColsAtCompileTime;
//       if(rows==Eigen::Dynamic) ar(rows);
//       if(cols==Eigen::Dynamic) ar(cols);
//       m.resize(rows,cols);
//       ar(binary_data(m.data(),static_cast<std::size_t>(rows*cols*sizeof(typename Derived::Scalar))));
//     }

    template <class Archive, class Derived> inline
    void save(Archive& ar, const Eigen::PlainObjectBase<Derived>& m) {
        std::vector<std::vector<typename Derived::Scalar>> m_;
        m_.resize(m.rows());
        for (size_t i = 0; i < m.rows(); i++) {
            m_[i].resize(m.cols());
            for (size_t j = 0; j < m.cols(); j++) {
                m_[i][j] = m(i,j);
            }
        }
        ar(m_);
    }

    template <class Archive, class Derived> inline
    void load(Archive& ar, Eigen::PlainObjectBase<Derived>& m) {
        std::vector<std::vector<typename Derived::Scalar>> m_;
        ar(m_);
        m.resize(m_.size(), m_[0].size());
        for (size_t i = 0; i < m.rows(); i++) {
            m_[i].resize(m.cols());
            for (size_t j = 0; j < m.cols(); j++) {
                m(i,j) = m_[i][j];
            }
        }
    }
}

namespace modus
{

}