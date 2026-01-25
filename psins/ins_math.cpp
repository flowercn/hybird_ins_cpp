#include "ins_math.h"

Matrix3d INSMath::askew(const Vector3d& v) {
    Matrix3d m;
    m <<   0.0,  -v(2),   v(1),
           v(2),   0.0,  -v(0),
          -v(1),   v(0),   0.0;
    return m;
}
 
Matrix3d INSMath::rv2m(const Vector3d& rv) {
    Matrix3d m;
    double n2 = rv.squaredNorm();
    double a, b;

    if (n2 < 1.0e-8) {
        a = 1.0 - n2 * (1.0 / 6.0 - n2 / 120.0);
        b = 0.5 - n2 * (1.0 / 24.0 - n2 / 720.0);
    } else {
        double n = std::sqrt(n2);
        a = std::sin(n) / n;
        b = (1.0 - std::cos(n)) / n2;
    }

    double arvx = a * rv(0), arvy = a * rv(1), arvz = a * rv(2);
    double bxx = b * rv(0) * rv(0), byy = b * rv(1) * rv(1), bzz = b * rv(2) * rv(2);
    double bxy = b * rv(0) * rv(1), bxz = b * rv(0) * rv(2), byz = b * rv(1) * rv(2);

    m(0, 0) = 1.0 - byy - bzz;  m(0, 1) = -arvz + bxy;      m(0, 2) = arvy + bxz;
    m(1, 0) = arvz + bxy;       m(1, 1) = 1.0 - bxx - bzz;  m(1, 2) = -arvx + byz;
    m(2, 0) = -arvy + bxz;      m(2, 1) = arvx + byz;       m(2, 2) = 1.0 - bxx - byy;

    return m;
}

Matrix3d INSMath::q2mat(const Quaterniond& qnb) {
    double n2 = qnb.squaredNorm();
    double s = (n2 > 0.0) ? (2.0 / n2) : 0.0;
    double qw = qnb.w(), qx = qnb.x(), qy = qnb.y(), qz = qnb.z();

    double q12 = s*qw*qx, q13 = s*qw*qy, q14 = s*qw*qz;
    double q22 = s*qx*qx, q23 = s*qx*qy, q24 = s*qx*qz;
    double q33 = s*qy*qy, q34 = s*qy*qz, q44 = s*qz*qz;

    Matrix3d Cnb;
    Cnb(0,0) = 1.0-q33-q44;    Cnb(0,1) = q23-q14;        Cnb(0,2) = q24+q13;
    Cnb(1,0) = q23+q14;        Cnb(1,1) = 1.0-q22-q44;    Cnb(1,2) = q34-q12;
    Cnb(2,0) = q24-q13;        Cnb(2,1) = q34+q12;        Cnb(2,2) = 1.0-q22-q33;

    return Cnb;
}

Quaterniond INSMath::a2qua(const Vector3d& att) {
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

Matrix3d INSMath::a2mat(const Vector3d& att) {
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

Vector3d INSMath::m2att(const Matrix3d& Cnb) {
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

Quaterniond INSMath::rv2q(const Vector3d& rv) {
    double n2 = rv.squaredNorm();
    double w, s;

    if (n2 < 1.0e-8) {
        w = 1.0 - n2 * (1.0 / 8.0 - n2 / 384.0);
        s = 0.5 - n2 * (1.0 / 48.0 - n2 / 3840.0);
    } else {
        double n_half = std::sqrt(n2) * 0.5;
        w = std::cos(n_half);
        s = std::sin(n_half) / (n_half * 2.0); 
    }

    Vector3d vec_part = s * rv;
    return Quaterniond(w, vec_part.x(), vec_part.y(), vec_part.z());
}

Quaterniond INSMath::qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in) {
    // Q_new = Q_n_correction * Q_old * Q_b_update
    Quaterniond q_ib = rv2q(rv_ib);
    Quaterniond q_in = rv2q(-rv_in); // 注意取反
    
    // Eigen乘法顺序: q_in * qnb0 * q_ib
    return (q_in * qnb0 * q_ib).normalized();
}

Vector3d INSMath::aa2phi(const Vector3d& att1, const Vector3d& att0) {
    Matrix3d Cn1b = a2mat(att1);
    Matrix3d Cn0b = a2mat(att0);
    Matrix3d Cnn = Cn1b * Cn0b.transpose(); // Cn1b * Cbn0
    Vector3d phi;
    phi << (Cnn(1, 2) - Cnn(2, 1)), 
           (Cnn(2, 0) - Cnn(0, 2)), 
           (Cnn(0, 1) - Cnn(1, 0));
    return phi / 2.0;
}

Vector3d INSMath::r2dms(double rad, const GLV& glv) {
    double s = (rad >= 0) ? 1.0 : -1.0;
    rad = std::abs(rad);
    double deg = rad / glv.deg;
    double d0 = std::floor(deg);
    double min = (deg - d0) * 60.0;
    double m0 = std::floor(min);
    double sec = (min - m0) * 60.0;
    return Vector3d(s * d0, m0, sec);
}

Vector3d INSMath::q2rv(const Quaterniond& q) {
    Quaterniond qt = q;
    if (qt.w() < 0) qt.coeffs() *= -1.0;
    double n2 = acos(qt.w());
    double k = (n2 > 1e-20) ? (2.0 * n2 / sin(n2)) : 2.0;
    return k * qt.vec();
}
