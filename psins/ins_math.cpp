#include "ins_math.h"

// ---------------------------------------------------------
// 辅助函数：严格对应 MATLAB askew.m
// ---------------------------------------------------------
Matrix3d INSMath::askew(const Vector3d& v) {
    Matrix3d m;
    // MATLAB: m = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
    // C++ Index: 0, 1, 2
    m <<  0.0,   -v(2),  v(1),
          v(2),   0.0,   -v(0),
         -v(1),   v(0),   0.0;
    return m;
}
 
// ---------------------------------------------------------
// 严格对应 MATLAB rv2m.m
// ---------------------------------------------------------
Matrix3d INSMath::rv2m(const Vector3d& rv) {
    Matrix3d m;
    double xx = rv(0) * rv(0);
    double yy = rv(1) * rv(1);
    double zz = rv(2) * rv(2);
    double n2 = xx + yy + zz;
    double a, b;

    // 对应 MATLAB: if n2<1.e-8 ...
    if (n2 < 1.0e-8) {
        // 泰勒级数近似
        // a = 1-n2*(1/6-n2/120);
        a = 1.0 - n2 * (1.0 / 6.0 - n2 / 120.0);
        // b = 0.5-n2*(1/24-n2/720);
        b = 0.5 - n2 * (1.0 / 24.0 - n2 / 720.0);
    } else {
        double n = std::sqrt(n2);
        a = std::sin(n) / n;
        b = (1.0 - std::cos(n)) / n2;
    }

    double arvx = a * rv(0);
    double arvy = a * rv(1);
    double arvz = a * rv(2);

    double bxx = b * xx;
    double bxy = b * rv(0) * rv(1);
    double bxz = b * rv(0) * rv(2);
    double byy = b * yy;
    double byz = b * rv(1) * rv(2);
    double bzz = b * zz;

    // MATLAB 代码直接给 m(1)~m(9) 赋值 (列优先)
    // m(1)=1-byy-bzz; m(4)=-arvz+bxy;  m(7)=arvy+bxz;
    // m(2)=arvz+bxy;  m(5)=1-bxx-bzz;  m(8)=-arvx+byz;
    // m(3)=-arvy+bxz; m(6)=arvx+byz;   m(9)=1-bxx-byy;
    
    // 转换到 C++ (row, col):
    // Col 1
    m(0, 0) = 1.0 - byy - bzz; 
    m(1, 0) = arvz + bxy;      
    m(2, 0) = -arvy + bxz;     

    // Col 2
    m(0, 1) = -arvz + bxy;
    m(1, 1) = 1.0 - bxx - bzz;
    m(2, 1) = arvx + byz;

    // Col 3
    m(0, 2) = arvy + bxz;
    m(1, 2) = -arvx + byz;
    m(2, 2) = 1.0 - bxx - byy;

    return m;
}

// ---------------------------------------------------------
// 严格对应 MATLAB q2mat.m
// ---------------------------------------------------------
Matrix3d INSMath::q2mat(const Quaterniond& qnb) {
    // 强制归一化 (对应 MATLAB q2mat 第一行虽然没写但隐含的要求，或者为了稳健性)
    // 但为了 100% 复刻，如果 MATLAB 没做归一化，我们通常也要小心。
    // PSINS 的 q2mat.m 只有: if length(qnb)==3...
    // 我们假设输入已经是四元数。
    
    // 映射: qnb(1)->w, qnb(2)->x, qnb(3)->y, qnb(4)->z
    double q11 = qnb.w() * qnb.w(); 
    double q12 = qnb.w() * qnb.x(); 
    double q13 = qnb.w() * qnb.y(); 
    double q14 = qnb.w() * qnb.z(); 
    
    double q22 = qnb.x() * qnb.x(); 
    double q23 = qnb.x() * qnb.y(); 
    double q24 = qnb.x() * qnb.z();     
    
    double q33 = qnb.y() * qnb.y(); 
    double q34 = qnb.y() * qnb.z();  
    
    double q44 = qnb.z() * qnb.z();

    Matrix3d Cnb;
    // MATLAB Code:
    // Cnb = [ q11+q22-q33-q44,  2*(q23-q14),     2*(q24+q13);
    //         2*(q23+q14),      q11-q22+q33-q44, 2*(q34-q12);
    //         2*(q24-q13),      2*(q34+q12),     q11-q22-q33+q44 ];
    
    Cnb(0,0) = q11+q22-q33-q44;  Cnb(0,1) = 2*(q23-q14);      Cnb(0,2) = 2*(q24+q13);
    Cnb(1,0) = 2*(q23+q14);      Cnb(1,1) = q11-q22+q33-q44;  Cnb(1,2) = 2*(q34-q12);
    Cnb(2,0) = 2*(q24-q13);      Cnb(2,1) = 2*(q34+q12);      Cnb(2,2) = q11-q22-q33+q44;

    return Cnb;
}

// ---------------------------------------------------------
// 严格对应 MATLAB a2qua.m
// ---------------------------------------------------------
Quaterniond INSMath::a2qua(const Vector3d& att) {
    // att = [pitch; roll; yaw]
    Vector3d att2 = att / 2.0;
    double sp = sin(att2(0)); double cp = cos(att2(0)); // pitch
    double sr = sin(att2(1)); double cr = cos(att2(1)); // roll
    double sy = sin(att2(2)); double cy = cos(att2(2)); // yaw
    
    // MATLAB Code:
    // qnb = [ cp*cr*cy - sp*sr*sy;  <- q0 (w)
    //         sp*cr*cy - cp*sr*sy;  <- q1 (x)
    //         cp*sr*cy + sp*cr*sy;  <- q2 (y)
    //         cp*cr*sy + sp*sr*cy ]; <- q3 (z)
    
    double w = cp*cr*cy - sp*sr*sy;
    double x = sp*cr*cy - cp*sr*sy;
    double y = cp*sr*cy + sp*cr*sy;
    double z = cp*cr*sy + sp*sr*cy;

    // MATLAB logic: if qnb(1)<0, qnb=-qnb;
    if (w < 0) {
        w = -w; x = -x; y = -y; z = -z;
    }
    
    // Eigen构造函数顺序: w, x, y, z
    return Quaterniond(w, x, y, z).normalized();
}

// ---------------------------------------------------------
// 严格对应 MATLAB a2mat.m (通过 a2qua->q2mat 实现，或手动)
// 这里我们手动实现以匹配你之前提供的代码
// ---------------------------------------------------------
Matrix3d INSMath::a2mat(const Vector3d& att) {
    double si = sin(att(0)), ci = cos(att(0)); // pitch (theta)
    double sj = sin(att(1)), cj = cos(att(1)); // roll (gamma)
    double sk = sin(att(2)), ck = cos(att(2)); // yaw (psi)

    Matrix3d Cnb;
    // MATLAB 3-1-2 Order from PSINS common usage (verify with your specific a2mat.m if available)
    // Assuming standard PSINS 3-1-2:
    // Cnb = [ cj*ck-si*sj*sk, -ci*sk,  sj*ck+si*cj*sk;
    //         cj*sk+si*sj*ck,  ci*ck,  sj*sk-si*cj*ck;
    //        -ci*sj,           si,     ci*cj ];
    
    Cnb(0,0) = cj*ck - si*sj*sk;  Cnb(0,1) = -ci*sk;   Cnb(0,2) = sj*ck + si*cj*sk;
    Cnb(1,0) = cj*sk + si*sj*ck;  Cnb(1,1) =  ci*ck;   Cnb(1,2) = sj*sk - si*cj*ck;
    Cnb(2,0) = -ci*sj;            Cnb(2,1) =  si;      Cnb(2,2) = ci*cj;

    return Cnb;
}

// ---------------------------------------------------------
// 严格对应 MATLAB m2att.m
// ---------------------------------------------------------
Vector3d INSMath::m2att(const Matrix3d& Cnb) {
    Vector3d att;
    // MATLAB: pitch = asin(Cnb(3,2)); -> C++ (2,1)
    att(0) = asin(Cnb(2, 1)); 

    if (abs(Cnb(2, 1)) <= 0.999999) {
        // MATLAB: roll = -atan2(Cnb(3,1), Cnb(3,3)); -> C++ (2,0), (2,2)
        att(1) = -atan2(Cnb(2, 0), Cnb(2, 2)); 
        // MATLAB: yaw  = -atan2(Cnb(1,2), Cnb(2,2)); -> C++ (0,1), (1,1)
        att(2) = -atan2(Cnb(0, 1), Cnb(1, 1)); 
    } else {
        // MATLAB: roll = atan2(Cnb(1,3), Cnb(1,1)); -> C++ (0,2), (0,0)
        att(1) = atan2(Cnb(0, 2), Cnb(0, 0));  
        att(2) = 0.0;
    }
    return att;
}

// ---------------------------------------------------------
// 严格对应 MATLAB rv2q.m
// ---------------------------------------------------------
Quaterniond INSMath::rv2q(const Vector3d& rv) {
    double n2 = rv.squaredNorm();
    double w, s;
    
    // PSINS 阈值 1.0e-8
    if (n2 < 1.0e-8) { 
        // 泰勒展开: 
        // w = 1 - n2*(1/8 - n2/384)
        w = 1.0 - n2 * (1.0 / 8.0 - n2 / 384.0);
        // s = 1/2 - n2*(1/48 - n2/3840)
        s = 0.5 - n2 * (1.0 / 48.0 - n2 / 3840.0);
    } else {
        double n = sqrt(n2);
        w = cos(n / 2.0);
        s = sin(n / 2.0) / n;
    }
    
    // 构造四元数: [w, s*rv]
    return Quaterniond(w, s * rv(0), s * rv(1), s * rv(2));
}

// ---------------------------------------------------------
// 对应 MATLAB qupdt2.m
// ---------------------------------------------------------
Quaterniond INSMath::qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in) {
    // MATLAB: qnb = qmul(rv2q(phi), qpb); 
    // 在 qupdt2 中: qnb = qmul(rv2q(-wnin*nts), qmul(qnb, rv2q(wm)));
    // 这里的顺序是: Q_new = Q_nav_correction * Q_old * Q_body_update
    
    Quaterniond q_ib = rv2q(rv_ib);
    Quaterniond q_in = rv2q(-rv_in); // 注意取反
    
    // Eigen乘法顺序: q_in * qnb0 * q_ib
    return (q_in * qnb0 * q_ib).normalized();
}

// ---------------------------------------------------------
// 对应 MATLAB aa2phi.m
// ---------------------------------------------------------
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