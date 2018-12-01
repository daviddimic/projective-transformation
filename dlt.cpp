#include <iostream>
#include "dlt.hpp"


static Eigen::Matrix3f unitProjective4pts(const std::vector<Eigen::Vector3f>& inPts) {
	Eigen::Matrix3f A;
	for (unsigned i = 0; i < inPts.size() - 1; i++) {
		A.col(i) = inPts[i];
	}
	Eigen::Vector3f lambda = A.colPivHouseholderQr().solve(inPts.back());

	Eigen::Matrix3f P1;
	for (unsigned i = 0; i < inPts.size() - 1; i++) {
		P1.col(i) = inPts[i] * lambda[i];
	}

	return P1;
}

Eigen::Matrix3f naiveProjective4pts(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts) {

	if(inPts.size() != 4 || outPts.size() != 4) {
		std::cerr << "number of input or output points is not 4!" << std::endl;
		return {};
    }

    Eigen::Matrix3f P1 = unitProjective4pts(inPts);
	Eigen::Matrix3f P2 = unitProjective4pts(outPts);

	Eigen::Matrix3f P;
	Eigen::Matrix3f	P1t = P1.inverse();
	P = P2 * P1t;

	return P;
}



static Eigen::MatrixXf getMatrix29f(const Eigen::Vector3f& inPoint, const Eigen::Vector3f& outPoint) {
	Eigen::MatrixXf A(2,9);

	Eigen::MatrixXf a1(2, 3), a2(2, 3), a3(2, 3);
	a1.setZero();
	a1.row(1) =  outPoint[2] * inPoint;

	a2.setZero();
	a2.row(0) = -a1.row(1);

	a3.row(0) = outPoint[1] * inPoint;
	a3.row(1) = -outPoint[0] * inPoint;

	A.row(0) << a1.row(0), a2.row(0), a3.row(0);
	A.row(1) << a1.row(1), a2.row(1), a3.row(1);
	return A;
}

Eigen::Matrix3f projectiveDLT(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts) {

	Eigen::MatrixXf A(2 * inPts.size(), 9);

	for (unsigned i = 0; i < inPts.size(); i++) {
		auto a = getMatrix29f(inPts[i], outPts[i]);

		A.row(2*i) << a.row(0);
		A.row(2*i+1) << a.row(1);
	}

	//SVD decomposition
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXf V = svd.matrixV();

	//last col of matrix V is solution
	Eigen::MatrixXf lastCol = V.rightCols(1);

	Eigen::Matrix3f P;
	P.row(0) =  lastCol.block(0, 0, 3, 1).transpose();
	P.row(1) = lastCol.block(3, 0, 3, 1).transpose();
	P.row(2) = lastCol.block(6, 0, 3, 1).transpose();

	return P;
}


static Eigen::Vector2f centreOfPoints(const std::vector<Eigen::Vector2f>& inPts) {
	float sx = 0.0, sy = 0.0;
	for (unsigned i = 0; i < inPts.size(); i++) {
		sx += inPts[i](0);
		sy += inPts[i](1);
	}
	Eigen::Vector2f center;
	center << sx / inPts.size(),
			sy / inPts.size();

	return center;
}

static float euclidDist2D(const Eigen::Vector2f& center, const Eigen::Vector2f& point) {
	Eigen::Vector2f diff = point - center;
	return diff.norm();
}

float averageDist(const Eigen::Vector2f& center, const std::vector<Eigen::Vector2f>& inPts) {
	float sumDist = 0.0;
	for (unsigned i = 0; i < inPts.size(); i++) {
		sumDist += euclidDist2D(center, inPts[i]);
	}
	return sumDist / inPts.size();
}

static Eigen::Transform<float, 2, Eigen::Affine> normalize(const std::vector<Eigen::Vector2f>& inPts) {
	//Translation matrix G
	Eigen::Vector2f center = centreOfPoints(inPts);
	Eigen::Vector2f trans = {0 - center(0), 0 - center(1)};
	Eigen::Translation2f G = Eigen::Translation2f(trans);

	//Scale matrix S
    float scale = float(std::sqrt(2.0)) / averageDist(center, inPts);
	Eigen::UniformScaling<float> S(scale);

	//Translate and scale matrix T = SG
	Eigen::Transform<float, 2, Eigen::Affine> T = S*G;

	return T;
}


Eigen::Matrix3f normalizedDLT(const std::vector<Eigen::Vector3f>& inPts, const std::vector<Eigen::Vector3f>& outPts) {
	std::vector<Eigen::Vector2f> inPtsNorm2f,  outPtsNorm2f;

	//make 3rd cord z = 1
	for (unsigned i = 0; i < inPts.size(); i++) {
		Eigen::Vector2f point;

		point << inPts[i](0) / inPts[i](2),
				inPts[i](1) / inPts[i](2);
		inPtsNorm2f.push_back(point);

		point << outPts[i](0) / outPts[i](2),
				outPts[i](1) / outPts[i](2);
		outPtsNorm2f.push_back(point);
	}

	Eigen::Transform<float, 2, Eigen::Affine> T = normalize(inPtsNorm2f);
	Eigen::Transform<float, 2, Eigen::Affine> Tp = normalize(outPtsNorm2f);

	std::vector<Eigen::Vector3f> inPtsNorm3f, outPtsNorm3f;
	//Apply T on points
	//and add 3. cord z = 1
	for (unsigned i = 0; i < inPtsNorm2f.size(); i++) {
		Eigen::Vector2f point2f;
		Eigen::Vector3f point3f;

		point2f = T * inPtsNorm2f[i];
        point3f << point2f(0),
				point2f(1),
				1;
		inPtsNorm3f.push_back(point3f);

		point2f = Tp * outPtsNorm2f[i];
		point3f << point2f(0),
				point2f(1),
				1;
		outPtsNorm3f.push_back(point3f);
	}

	Eigen::Matrix3f Pp = projectiveDLT(inPtsNorm3f, outPtsNorm3f);

	Eigen::Matrix3f P = (Tp.inverse() * Pp * T).matrix();

	return P;
}
