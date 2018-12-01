#ifndef DLT_H
#define DLT_H 1

#include <vector>
#include <Eigen/Geometry>
#include <Eigen/SVD>


Eigen::Matrix3f naiveProjective4pts(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts);

Eigen::Matrix3f projectiveDLT(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts);

Eigen::Matrix3f normalizedDLT(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts);


#endif /* ifndef DLT_H */
