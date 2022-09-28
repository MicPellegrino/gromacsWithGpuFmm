#ifndef _BK_data_type_hpp
#define _BK_data_type_hpp


//memory intensive, for tree depths < 6

#define M2M_SOA_OPTIMIZATION

//gain ~10% performance but only before VOLTA
#define M2L_SOA_OPTIMIZATION

#define L2L_SOA_OPTIMIZATION

#ifndef GMX_FMM_DOUBLE
//#define GMX_FMM_DOUBLE
#endif

#ifndef GMX_FMM_DOUBLE

#define REAL float
#define INTTYPE int
#define REAL2 float2
#define REAL3 float3
#define REAL4 float4

#define make_real2(x,y) make_float2(x,y)
#define make_real3(x,y,z) make_float3(x,y,z)
#define make_real4(x,y,z,w) make_float4(x,y,z,w)

#else

#ifdef M2L_SOA_OPTIMIZATION
#undef M2L_SOA_OPTIMIZATION
#endif

#define REAL double
#define INTTYPE long long int
#define REAL2 double2
#define REAL3 double3
#define REAL4 double4

#define make_real2(x,y) make_double2(x,y)
#define make_real3(x,y,z) make_double3(x,y,z)
#define make_real4(x,y,z,w) make_double4(x,y,z,w)

#endif

#ifndef GMX_MATH_UNITS_H

#define ANGSTROM         (1e-10)                           /* Old...	*/
#define KILO             (1e3)                             /* Thousand	*/
#define NANO             (1e-9)                            /* A Number	*/
#define PICO             (1e-12)                           /* A Number	*/
#define A2NM             (ANGSTROM/NANO)                   /* NANO	        */
#define NM2A             (NANO/ANGSTROM)                   /* 10.0		*/
#define RAD2DEG          (180.0/M_PI)                      /* Conversion	*/
#define DEG2RAD          (M_PI/180.0)                      /* id		*/
#define CAL2JOULE        (4.184)                           /* id		*/
#define E_CHARGE         (1.602176565e-19)                 /* Coulomb, NIST 2010 CODATA */

#define AMU              (1.660538921e-27)                 /* kg, NIST 2010 CODATA  */
#define BOLTZMANN        (1.3806488e-23)                   /* (J/K, NIST 2010 CODATA */
#define AVOGADRO         (6.02214129e23)                   /* no unit, NIST 2010 CODATA */
#define RGAS             (BOLTZMANN*AVOGADRO)              /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)                       /* (kJ/(mol K)) */
#define FARADAY          (E_CHARGE*AVOGADRO)               /* (C/mol)      */
#define ELECTRONVOLT     (E_CHARGE*AVOGADRO/KILO)          /* (kJ/mol)   */
#define PLANCK1          (6.62606957e-34)                  /* J s, NIST 2010 CODATA */
#define PLANCK           (PLANCK1*AVOGADRO/(PICO*KILO))    /* (kJ/mol) ps */

#define EPSILON0_SI      (8.854187817e-12)                 /* F/m,  NIST 2010 CODATA */
/* Epsilon in our MD units: (e^2 / Na (kJ nm)) == (e^2 mol/(kJ nm)) */
#define EPSILON0         ((EPSILON0_SI*NANO*KILO)/(E_CHARGE*E_CHARGE*AVOGADRO))

#define SPEED_OF_LIGHT   (2.99792458E05)                   /* nm/ps, NIST 2010 CODATA */
#define ATOMICMASS_keV   (931494.061)                      /* Atomic mass in keV, NIST 2010 CODATA   */
#define ELECTRONMASS_keV (510.998928)                      /* Electron mas in keV, NIST 2010 CODATA  */

#define RYDBERG          (1.0973731568539e-02)             /* nm^-1, NIST 2010 CODATA */

#define ONE_4PI_EPS0     (1.0/(4.0*M_PI*EPSILON0))
#define FACEL            (10.0*ONE_4PI_EPS0)

#endif

#endif
