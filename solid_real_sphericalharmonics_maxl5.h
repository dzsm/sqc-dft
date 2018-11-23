#ifndef SOLID_REAL_SPHERICAL_HARMONICS_MAXL5
#define SOLID_REAL_SPHERICAL_HARMONICS_MAXL5

// srsh[36] and l*(1 + l) + m; If |(x,y,z)|=1; it is regular
struct solid_real_spherical_harmonics_maxl5
{
    double srsh[36];

    int eval(double x,double y,double z,int maxl);
    int eval(double x,double y,double z,double r,int maxl);

    inline double get (int l,int m) {return srsh[l*(1 + l) + m];}
};

#endif // SOLID_REAL_SPHERICAL_HARMONICS_MAXL5
