#include "ins_math.h"
#include <cmath>
#include <iostream>

namespace INSMath {

    // 反对称矩阵 [0 -z y; z 0 -x; -y x 0]
    Matrix3d askew(const Vector3d& v) {
        Matrix3d m;
        m << 0.0, -v(2), v(1),
             v(2), 0.0, -v(0),
             -v(1), v(0), 0.0;
        return m;
    }

    // 旋转矢量 -> 旋转矩阵 (罗德里格斯公式)
    Matrix3d rv2m(const Vector3d& rv) {
        Matrix3d I = Matrix3d::Identity();
        double n2 = rv.squaredNorm();
        
        // 小角度近似 (防止除以0)
        if (n2 < 1.0e-8) {
            return I + askew(rv) + 0.5 * askew(rv) * askew(rv);
        } else {
            double n = std::sqrt(n2);
            double sn = std::sin(n);
            double cn = std::cos(n);
            Matrix3d sk = askew(rv);
            return I + (sn / n) * sk + ((1.0 - cn) / n2) * sk * sk;
        }
    }

    // 旋转矢量 -> 四元数
    Quaterniond rv2q(const Vector3d& rv) {
        double n2 = rv.squaredNorm();
        if (n2 < 1.0e-8) {
            // q = [1, 0.5*rv]
            Vector3d v = 0.5 * rv;
            return Quaterniond(1.0, v(0), v(1), v(2)).normalized();
        } else {
            double n = std::sqrt(n2);
            double half_n = 0.5 * n;
            double s = std::sin(half_n) / n;
            return Quaterniond(std::cos(half_n), s * rv(0), s * rv(1), s * rv(2));
        }
    }

    Quaterniond a2qua(const Vector3d& att) {
        Vector3d s, c;
        const Vector3d att2 = att * 0.5;

        s(0) = std::sin(att2(0)); c(0) = std::cos(att2(0)); // pitch
        s(1) = std::sin(att2(1)); c(1) = std::cos(att2(1)); // roll
        s(2) = std::sin(att2(2)); c(2) = std::cos(att2(2)); // yaw

        double w = c(0)*c(1)*c(2) - s(0)*s(1)*s(2);
        double x = s(0)*c(1)*c(2) - c(0)*s(1)*s(2);
        double y = c(0)*s(1)*c(2) + s(0)*c(1)*s(2);
        double z = c(0)*c(1)*s(2) + s(0)*s(1)*c(2);
        
        if (w < 0) {
            w = -w; x = -x; y = -y; z = -z;
        }

        return Quaterniond(w, x, y, z).normalized();
    }

    Matrix3d a2mat(const Vector3d& att) {
        double si = std::sin(att(0)), ci = std::cos(att(0)); // pitch
        double sj = std::sin(att(1)), cj = std::cos(att(1)); // roll
        double sk = std::sin(att(2)), ck = std::cos(att(2)); // yaw

        double sisj = si * sj, sicj = si * cj;
        Matrix3d Cnb;
        Cnb(0,0) = cj*ck - sisj*sk;  Cnb(0,1) = -ci*sk;  Cnb(0,2) = sj*ck + sicj*sk;
        Cnb(1,0) = cj*sk + sisj*ck;  Cnb(1,1) =  ci*ck;  Cnb(1,2) = sj*sk - sicj*ck;
        Cnb(2,0) = -ci*sj;           Cnb(2,1) =  si;     Cnb(2,2) = ci*cj;

        return Cnb;
    }

    Vector3d m2att(const Matrix3d& Cnb) {
        Vector3d att;
        double c21 = Cnb(2, 1);
        att(0) = std::asin(c21);

        if (std::abs(c21) < 0.9999999999) {
            att(1) = -std::atan2(Cnb(2, 0), Cnb(2, 2)); // roll
            att(2) = -std::atan2(Cnb(0, 1), Cnb(1, 1)); // yaw
        } else {
            att(1) = std::atan2(Cnb(0, 2), Cnb(0, 0));  // roll
            att(2) = 0.0;                               // yaw 置零
        }
        return att;
    }

    // --- 补充实现 (Eigen Wrapper) ---
    Matrix3d q2mat(const Quaterniond& qnb) {
        return qnb.toRotationMatrix();
    }

    Quaterniond m2qua(const Matrix3d& Cnb) {
        return Quaterniond(Cnb);
    }
    // -----------------------------

    Vector3d q2rv(const Quaterniond& q) {
        Quaterniond qt = q;
        if (qt.w() < 0) qt.coeffs() *= -1.0;
        double n2 = acos(qt.w());
        double k = (n2 > 1e-20) ? (2.0 * n2 / sin(n2)) : 2.0;
        return k * qt.vec();
    }

    Quaterniond qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in) {
        // Q_new = Q_n_correction * Q_old * Q_b_update
        Quaterniond q_ib = rv2q(rv_ib);
        Quaterniond q_in = rv2q(-rv_in); // 注意取反
        
        // Eigen乘法顺序: q_in * qnb0 * q_ib
        return (q_in * qnb0 * q_ib).normalized();
    }

    Vector3d aa2phi(const Vector3d& att1, const Vector3d& att0) {
        Matrix3d Cn1b = a2mat(att1);
        Matrix3d Cn0b = a2mat(att0);
        Matrix3d Cnn = Cn1b * Cn0b.transpose(); // Cn1b * Cbn0
        Vector3d phi;
        phi << (Cnn(1, 2) - Cnn(2, 1)), 
               (Cnn(2, 0) - Cnn(0, 2)), 
               (Cnn(0, 1) - Cnn(1, 0));
        return phi / 2.0;
    }

    Vector3d r2dms(double rad, const GLV& glv) {
        double s = (rad >= 0) ? 1.0 : -1.0;
        rad = std::abs(rad);
        double deg = rad / glv.deg;
        double d0 = std::floor(deg);
        double min = (deg - d0) * 60.0;
        double m0 = std::floor(min);
        double sec = (min - m0) * 60.0;
        return Vector3d(s * d0, m0, sec);
    }
}