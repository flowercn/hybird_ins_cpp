#ifndef KF_STATE_H
#define KF_STATE_H

#include <Eigen/Dense>
#include "glv.h"

using namespace Eigen;

class KFAlignVN {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // 0-2: phi(3), 3-5: dv(3), 6-8: eb(3), 9-11: db(3)
    static const int n = 12; 
    static const int m = 3;  // 观测: vn

    VectorXd xk;       // 12x1
    MatrixXd Pxk;      // 12x12
    MatrixXd Phikk_1;  // 状态转移矩阵 Fk
    MatrixXd Hk;       // 观测矩阵
    MatrixXd Qk;       // 过程噪声
    MatrixXd Rk;       // 测量噪声

    // 中间变量
    VectorXd yk; MatrixXd Kk;

    KFAlignVN();
    void Init(double nts, const GLV& glv, 
              double phi0, double web_psd, double wdb_psd, 
              double eb_sigma, double db_sigma, double wvn_err);

    void UpdatePhi(const Vector3d& dvn, const Matrix3d& Cnb, const Vector3d& wnie, double nts);
    void Update(const Vector3d& vn_meas);
    void Feedback(Quaterniond& qnb, Vector3d& vn, Vector3d& eb, Vector3d& db, Matrix3d& Cnb);
};

#endif // KF_STATE_H