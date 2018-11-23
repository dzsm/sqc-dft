#include "solid_real_sphericalharmonics_maxl5.h"


inline void solid_real_spherical_harmonics_l0(double *srsh,double x,double y,double z)
{
    srsh[0]  = 0.2820947917738781;                                // 0  0
}

inline void solid_real_spherical_harmonics_l1(double *srsh,double x,double y,double z)
{

    srsh[0]  = 0.2820947917738781;                                // 0  0

    srsh[1]  = 0.4886025119029199*y;                              // 1 -1
    srsh[2]  = 0.4886025119029199*z;                              // 1  0
    srsh[3]  = 0.4886025119029199*x;                              // 1  1

}

inline void solid_real_spherical_harmonics_l2(double *srsh,double x,double y,double z)
{
    double o1,o2,o3,o4,o5,o6;

    o1  = z*z;
    o2 = 3.*o1;
    o3 = -1. + o2;
    o4  = -1.*y;
    o5 = o4 + x;
    o6 = x + y;

    srsh[0]  = 0.2820947917738781;                                // 0  0

    srsh[1]  = 0.4886025119029199*y;                              // 1 -1
    srsh[2]  = 0.4886025119029199*z;                              // 1  0
    srsh[3]  = 0.4886025119029199*x;                              // 1  1

    srsh[4]  = 1.092548430592079*x*y;                             // 2 -2
    srsh[5]  = 1.092548430592079*y*z;                             // 2 -1
    srsh[6]  = 0.31539156525252*o3;                               // 2  0
    srsh[7]  = 1.092548430592079*x*z;                             // 2  1
    srsh[8]  = 0.5462742152960395*o5*o6;                          // 2  2

}

inline void solid_real_spherical_harmonics_l3(double *srsh,double x,double y,double z)
{
    double o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14;

    o1  = z*z;
    o2 = 3.*o1;
    o3 = -1. + o2;
    o4  = -1.*y;
    o5 = o4 + x;
    o6 = x + y;
    o7 = x*x;
    o8 = -3.*o7;
    o9  = y*y;
    o10 = o8 + o9;
    o11 = 5.*o1;
    o12 = -1. + o11;
    o13 = -3.*o9;
    o14 = o13 + o7;

    srsh[0]  = 0.2820947917738781;                                // 0  0

    srsh[1]  = 0.4886025119029199*y;                              // 1 -1
    srsh[2]  = 0.4886025119029199*z;                              // 1  0
    srsh[3]  = 0.4886025119029199*x;                              // 1  1

    srsh[4]  = 1.092548430592079*x*y;                             // 2 -2
    srsh[5]  = 1.092548430592079*y*z;                             // 2 -1
    srsh[6]  = 0.31539156525252*o3;                               // 2  0
    srsh[7]  = 1.092548430592079*x*z;                             // 2  1
    srsh[8]  = 0.5462742152960395*o5*o6;                          // 2  2

    srsh[9]  = -0.5900435899266435*o10*y;                         // 3 -3
    srsh[10] = 2.890611442640554*x*y*z;                           // 3 -2
    srsh[11] = 0.4570457994644657*o12*y;                          // 3 -1
    srsh[12] = 0.3731763325901154*(-3. + o11)*z;                  // 3  0
    srsh[13] = 0.4570457994644657*o12*x;                          // 3  1
    srsh[14] = 1.445305721320277*o5*o6*z;                         // 3  2
    srsh[15] = 0.5900435899266435*o14*x;                          // 3  3

}

inline void solid_real_spherical_harmonics_l4(double *srsh,double x,double y,double z)
{
    double o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,
           o16,o17,o18,o19,o20,o21,o22,o23,o24,o25,o26,o27,o28;

    o1  = z*z;
    o2 = 3.*o1;
    o3 = -1. + o2;
    o4  = -1.*y;
    o5 = o4 + x;
    o6 = x + y;
    o7 = x*x;
    o8 = -3.*o7;
    o9  = y*y;
    o10 = o8 + o9;
    o11 = 5.*o1;
    o12 = -1. + o11;
    o13 = -3.*o9;
    o14 = o13 + o7;
    o15 = 7.*o1;
    o16 = -1. + o15;
    o17 = -3. + o15;
    o18 = o1*o1;
    o19 = 2.*o9;
    o20 = o1 + o19;
    o21 = -1. + o20;
    o22 = o9*o9;
    o23 = 8.*o22;
    o24 = -1. + o1;
    o25 = o24*o9;
    o26 = 8.*o25;
    o27 = o24*o24;
    o28 = o23 + o26 + o27;

    srsh[0]  = 0.2820947917738781;                                // 0  0
    srsh[1]  = 0.4886025119029199*y;                              // 1 -1
    srsh[2]  = 0.4886025119029199*z;                              // 1  0
    srsh[3]  = 0.4886025119029199*x;                              // 1  1
    srsh[4]  = 1.092548430592079*x*y;                             // 2 -2
    srsh[5]  = 1.092548430592079*y*z;                             // 2 -1
    srsh[6]  = 0.31539156525252*o3;                               // 2  0
    srsh[7]  = 1.092548430592079*x*z;                             // 2  1
    srsh[8]  = 0.5462742152960395*o5*o6;                          // 2  2
    srsh[9]  = -0.5900435899266435*o10*y;                         // 3 -3
    srsh[10] = 2.890611442640554*x*y*z;                           // 3 -2
    srsh[11] = 0.4570457994644657*o12*y;                          // 3 -1
    srsh[12] = 0.3731763325901154*(-3. + o11)*z;                  // 3  0
    srsh[13] = 0.4570457994644657*o12*x;                          // 3  1
    srsh[14] = 1.445305721320277*o5*o6*z;                         // 3  2
    srsh[15] = 0.5900435899266435*o14*x;                          // 3  3
    srsh[16] = 2.503342941796705*o5*o6*x*y;                       // 4 -4
    srsh[17] = -1.770130769779931*o10*y*z;                        // 4 -3
    srsh[18] = 0.94617469575756*o16*x*y;                          // 4 -2
    srsh[19] = 0.6690465435572892*o17*y*z;                        // 4 -1
    srsh[20] = 0.1057855469152043*(3. - 30.*o1 + 35.*o18);        // 4  0
    srsh[21] = 0.6690465435572892*o17*x*z;                        // 4  1
    srsh[22] = -0.47308734787878*o16*o21;                         // 4  2
    srsh[23] = 1.770130769779931*o14*x*z;                         // 4  3
    srsh[24] = 0.6258357354491761*o28;                            // 4  4

}

inline void solid_real_spherical_harmonics_l5(double *srsh,double x,double y,double z)
{
    double o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,
           o16,o17,o18,o19,o20,o21,o22,o23,o24,o25,o26,o27,o28,
           o29,o30,o31,o32,o33,o34,o35,o36;

    o1  = z*z;
    o2 = 3.*o1;
    o3 = -1. + o2;
    o4  = -1.*y;
    o5 = o4 + x;
    o6 = x + y;
    o7 = x*x;
    o8 = -3.*o7;
    o9  = y*y;
    o10 = o8 + o9;
    o11 = 5.*o1;
    o12 = -1. + o11;
    o13 = -3.*o9;
    o14 = o13 + o7;
    o15 = 7.*o1;
    o16 = -1. + o15;
    o17 = -3. + o15;
    o18 = o1*o1;
    o19 = 2.*o9;
    o20 = o1 + o19;
    o21 = -1. + o20;
    o22 = o9*o9;
    o23 = 8.*o22;
    o24 = -1. + o1;
    o25 = o24*o9;
    o26 = 8.*o25;
    o27 = o24*o24;
    o28 = o23 + o26 + o27;
    o29 = 16.*o22;
    o30 = 4.*o9;
    o31 = 9.*o1;
    o32 = -1. + o31;
    o33 = -14.*o1;
    o34 = 21.*o18;
    o35 = o33 + o34;
    o36 = 1. + o35;

    srsh[0]  = 0.2820947917738781;                                // 0  0
    srsh[1]  = 0.4886025119029199*y;                              // 1 -1
    srsh[2]  = 0.4886025119029199*z;                              // 1  0
    srsh[3]  = 0.4886025119029199*x;                              // 1  1
    srsh[4]  = 1.092548430592079*x*y;                             // 2 -2
    srsh[5]  = 1.092548430592079*y*z;                             // 2 -1
    srsh[6]  = 0.31539156525252*o3;                               // 2  0
    srsh[7]  = 1.092548430592079*x*z;                             // 2  1
    srsh[8]  = 0.5462742152960395*o5*o6;                          // 2  2
    srsh[9]  = -0.5900435899266435*o10*y;                         // 3 -3
    srsh[10] = 2.890611442640554*x*y*z;                           // 3 -2
    srsh[11] = 0.4570457994644657*o12*y;                          // 3 -1
    srsh[12] = 0.3731763325901154*(-3. + o11)*z;                  // 3  0
    srsh[13] = 0.4570457994644657*o12*x;                          // 3  1
    srsh[14] = 1.445305721320277*o5*o6*z;                         // 3  2
    srsh[15] = 0.5900435899266435*o14*x;                          // 3  3
    srsh[16] = 2.503342941796705*o5*o6*x*y;                       // 4 -4
    srsh[17] = -1.770130769779931*o10*y*z;                        // 4 -3
    srsh[18] = 0.94617469575756*o16*x*y;                          // 4 -2
    srsh[19] = 0.6690465435572892*o17*y*z;                        // 4 -1
    srsh[20] = 0.1057855469152043*(3. - 30.*o1 + 35.*o18);        // 4  0
    srsh[21] = 0.6690465435572892*o17*x*z;                        // 4  1
    srsh[22] = -0.47308734787878*o16*o21;                         // 4  2
    srsh[23] = 1.770130769779931*o14*x*z;                         // 4  3
    srsh[24] = 0.6258357354491761*o28;                            // 4  4
    srsh[25] = 0.6563820568401701*(20.*o25 + 5.*o27 + o29)*y;     // 5 -5
    srsh[26] = 8.302649259524165*o5*o6*x*y*z;                     // 5 -4
    srsh[27] = -0.4892382994352504*(-3. + o2 + o30)*o32*y;        // 5 -3
    srsh[28] = 4.793536784973324*o3*x*y*z;                        // 5 -2
    srsh[29] = 0.4529466511956969*o36*y;                          // 5 -1
    srsh[30] = 0.1169503224534236*(15. - 70.*o1 + 63.*o18)*z;     // 5  0
    srsh[31] = 0.4529466511956969*o36*x;                          // 5  1
    srsh[32] = -2.396768392486662*o21*o3*z;                       // 5  2
    srsh[33] = -0.4892382994352504*(-1. + o1 + o30)*o32*x;        // 5  3
    srsh[34] = 2.075662314881041*o28*z;                           // 5  4
    srsh[35] = 0.6563820568401701*(12.*o25 + o27 + o29)*x;        // 5  5

}

inline void solid_real_spherical_harmonics_l5_zero(double *srsh)
{
    srsh[0]  = 0.2820947917738781;
    srsh[1]  = 0.0;
    srsh[2]  = 0.0;
    srsh[3]  = 0.0;
    srsh[4]  = 0.0;
    srsh[5]  = 0.0;
    srsh[6]  = 0.0;
    srsh[7]  = 0.0;
    srsh[8]  = 0.0;
    srsh[9]  = 0.0;
    srsh[10] = 0.0;
    srsh[11] = 0.0;
    srsh[12] = 0.0;
    srsh[13] = 0.0;
    srsh[14] = 0.0;
    srsh[15] = 0.0;
    srsh[16] = 0.0;
    srsh[17] = 0.0;
    srsh[18] = 0.0;
    srsh[19] = 0.0;
    srsh[20] = 0.0;
    srsh[21] = 0.0;
    srsh[22] = 0.0;
    srsh[23] = 0.0;
    srsh[24] = 0.0;
    srsh[25] = 0.0;
    srsh[26] = 0.0;
    srsh[27] = 0.0;
    srsh[28] = 0.0;
    srsh[29] = 0.0;
    srsh[30] = 0.0;
    srsh[31] = 0.0;
    srsh[32] = 0.0;
    srsh[33] = 0.0;
    srsh[34] = 0.0;
    srsh[35] = 0.0;
}

// srsh[36] and l*(1 + l) + m; If |(x,y,z)|=1; it is regular

int solid_real_spherical_harmonics_maxl5::eval(double x,double y,double z,int maxl)
{

    if (x==0.0 && y==0.0 && z==0.0)
    {
        solid_real_spherical_harmonics_l5_zero(srsh);
        return 5;
    }

    switch(maxl)
    {
    case 0 :
        solid_real_spherical_harmonics_l0(srsh,x,y,z);
        break;
    case 1 :
        solid_real_spherical_harmonics_l1(srsh,x,y,z);
        break;
    case 2 :
        solid_real_spherical_harmonics_l2(srsh,x,y,z);
        break;
    case 3 :
        solid_real_spherical_harmonics_l3(srsh,x,y,z);
        break;
    case 4 :
        solid_real_spherical_harmonics_l4(srsh,x,y,z);
        break;
    case 5 :
        solid_real_spherical_harmonics_l5(srsh,x,y,z);
        break;
    default :
        throw __FUNCTION__;
    }
    return maxl;
}

int solid_real_spherical_harmonics_maxl5::eval(double x,double y,double z,double r,int maxl)
{
    //maxl = 5;
    if (r==0.0)
    {
        solid_real_spherical_harmonics_l5_zero(srsh);
        return 5;
    }

    x=x/r;
    y=y/r;
    z=z/r;

    switch(maxl)
    {
    case 0 :
        solid_real_spherical_harmonics_l0(srsh,x,y,z);
        break;
    case 1 :
        solid_real_spherical_harmonics_l1(srsh,x,y,z);
        break;
    case 2 :
        solid_real_spherical_harmonics_l2(srsh,x,y,z);
        break;
    case 3 :
        solid_real_spherical_harmonics_l3(srsh,x,y,z);
        break;
    case 4 :
        solid_real_spherical_harmonics_l4(srsh,x,y,z);
        break;
    case 5 :
        solid_real_spherical_harmonics_l5(srsh,x,y,z);
        break;
    default :
        throw __FUNCTION__;
    }
    return maxl;
}

