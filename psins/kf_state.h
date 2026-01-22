#ifndef KF_STATE_H
#define KF_STATE_H

#include <Eigen/Dense>
#include "ins_state.h"
#include "glv.h"

using namespace Eigen;

class KFAlignVN {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 复刻 alignvn.m: 12维状态
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

    // 初始化 (对应 avnkfinit)
    void Init(double nts, const GLV& glv, 
              double phi0, double web_psd, double wdb_psd, 
              double eb_sigma, double db_sigma, double wvn_err);

    // 状态转移矩阵更新 (对应 loop 中的 kf.Phikk_1 更新)
    // 需要传入 dvn 和 Cnb 用于构建矩阵块
    void UpdatePhi(const Vector3d& dvn, const Matrix3d& Cnb, const Vector3d& wnie, double nts);

    // 测量更新 (对应 kfupdate)
    void Update(const Vector3d& vn_meas);

    // 反馈 (对应 alignvn.m 底部独特的 0.91/0.09 反馈逻辑)
    // 注意：alignvn.m 修改了输入的 qnb 和 vn
    void Feedback(Quaterniond& qnb, Vector3d& vn);
};

#endif // KF_STATE_H