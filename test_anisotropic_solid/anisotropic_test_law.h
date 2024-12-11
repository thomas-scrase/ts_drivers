#ifndef ANISOTROPIC_TEST_LAW_HEADER
#define ANISOTROPIC_TEST_LAW_HEADER

#include "anisotropic_constitutive.h"

namespace oomph
{

class AnisotropicGeneralisedMooneyRivlin : public AnisotropicStrainEnergyFunction
{
public:
  /// Constructor takes the pointers to the constitutive parameters:
  /// Poisson's ratio, the Mooney-Rivlin parameter. Young's modulus is set
  /// to 1, implying that it has been used to scale the stresses
  AnisotropicGeneralisedMooneyRivlin(double* nu_pt, double* c1_pt)
    : AnisotropicStrainEnergyFunction(),
      Nu_pt(nu_pt),
      C1_pt(c1_pt),
      E_pt(new double(1.0)),
      Must_delete_e(true)
  {
   N_Principal_Vectors_Of_Anisotropy = 1;
   N_Additional_Strain_Invariants = 1;
  }

  /// Constructor takes the pointers to the constitutive parameters:
  /// Poisson's ratio, the Mooney-Rivlin parameter and Young's modulus
  AnisotropicGeneralisedMooneyRivlin(double* nu_pt, double* c1_pt, double* e_pt)
    : AnisotropicStrainEnergyFunction(),
      Nu_pt(nu_pt),
      C1_pt(c1_pt),
      E_pt(e_pt),
      Must_delete_e(false)
  {
   N_Principal_Vectors_Of_Anisotropy = 1;
   N_Additional_Strain_Invariants = 1;
  }


  /// Virtual destructor
  virtual ~AnisotropicGeneralisedMooneyRivlin()
  {
    if (Must_delete_e) delete E_pt;
  }

  void I(const DenseMatrix<double>& g,
        const DenseMatrix<double>& g_up,
        const DenseMatrix<double>& G,
        const DenseMatrix<double>& G_up,
        const double& detg,
        const double& getG,
        const Vector<Vector<double>>& a,
        Vector<double>& I,
        Vector<DenseMatrix<double>>& dIdG)
 {
  const unsigned dim = g.ncol();
  // The only additional strain invariant this model uses is a_iG_{ij}a_j
  I[0] = 0.0;
  for(unsigned i=0; i<dim; i++)
  {
   for(unsigned j=0; j<dim; j++)
   {
    I[0] += a[0][i] * G(i,j) * a[0][j];

    dIdG[0](i,j) = a[0][i] * a[0][j];
   }
  }
 }

  /// Return the strain energy in terms of strain tensor
  double W(const DenseMatrix<double>& gamma, const Vector<Vector<double>> a)
  {
    return AnisotropicStrainEnergyFunction::W(gamma, a);
  }


  /// Return the strain energy in terms of the strain invariants
  double W(const Vector<double>& I)
  {
    double G = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));
    return 0.5 * ((*C1_pt) * (I[0] - 3.0) + (G - (*C1_pt)) * (I[1] - 3.0) +
                  ((*C1_pt) - 2.0 * G) * (I[2] - 1.0) +
                  (1.0 - (*Nu_pt)) * G * (I[2] - 1.0) * (I[2] - 1.0) /
                    (2.0 * (1.0 - 2.0 * (*Nu_pt))))
           + pow(I[3] - 1.0, 2.0);
  }


  /// Return the derivatives of the strain energy function with
  /// respect to the strain invariants
  void derivatives(Vector<double>& I, Vector<double>& dWdI)
  {
    double G = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));
    dWdI[0] = 0.5 * (*C1_pt);
    dWdI[1] = 0.5 * (G - (*C1_pt));
    dWdI[2] = 0.5 * ((*C1_pt) - 2.0 * G +
                     2.0 * (1.0 - (*Nu_pt)) * G * (I[2] - 1.0) /
                       (2.0 * (1.0 - 2.0 * (*Nu_pt))));
    dWdI[3] = 2.0 * (I[3] - 1.0);
  }


  /// Pure virtual function in which the user must declare if the
  /// constitutive equation requires an incompressible formulation
  /// in which the volume constraint is enforced explicitly.
  /// Used as a sanity check in PARANOID mode. False.
  bool requires_incompressibility_constraint()
  {
    return false;
  }

private:
  /// Poisson's ratio
  double* Nu_pt;

  /// Mooney-Rivlin parameter
  double* C1_pt;

  /// Young's modulus
  double* E_pt;

  /// Boolean flag to indicate if storage for elastic modulus
  /// must be deleted in destructor
  bool Must_delete_e;
};



// A nonsense strain energy function to illustrate how to use anisotropic strain energy functions
class AnisotropicTestStrainEnergy : public AnisotropicStrainEnergyFunction
{
public:
 AnisotropicTestStrainEnergy(double* c1_pt, double* c2_pt, double* c3_pt, double* c4_pt) : AnisotropicStrainEnergyFunction(),
                                                                            C1_pt(c1_pt),
                                                                            C2_pt(c2_pt),
                                                                            C3_pt(c3_pt),
                                                                            C4_pt(c4_pt)
 {
  N_Principal_Vectors_Of_Anisotropy = 1;
  // One additional strain invariant arises from the fibre orientation
  N_Additional_Strain_Invariants = 1;
 }

 // Tell the constitutive law how to compute the additional strain invariants and the derivatives of the strain invariants with respect to the deformed metric tensor
 void I(const DenseMatrix<double>& g,
        const DenseMatrix<double>& g_up,
        const DenseMatrix<double>& G,
        const DenseMatrix<double>& G_up,
        const double& detg,
        const double& getG,
        const Vector<Vector<double>>& a,
        Vector<double>& I,
        Vector<DenseMatrix<double>>& dIdG)
 {
  const unsigned dim = g.ncol();
  // The only additional strain invariant this model uses is a_iG_{ij}a_j
  I[0] = 0.0;
  for(unsigned i=0; i<dim; i++)
  {
   for(unsigned j=0; j<dim; j++)
   {
    I[0] += a[0][i] * G(i,j) * a[0][j];

    dIdG[0](i,j) = a[0][i] * a[0][j];
   }
  }
 }

 // Compute the strain energy
 double W(const Vector<double>& I)
 {
  // return (*C1_pt)*(I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0) + (*C3_pt)*pow(I[2]-1.0, 2.0)
  //        + (*C4_pt)*pow(I[3]-1.0, 2.0);

  double Q = (*C3_pt)*I[3]*I[3];
  return (*C1_pt)/2.0*(exp(Q - 1.0)) + (*C2_pt)/2.0*(I[2]-1.0)*log(I[2]);
 }

 /// Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 void derivatives(Vector<double>& I, Vector<double>& dWdI)
 {
   // dWdI[0] = (*C1_pt);
   // dWdI[1] = (*C2_pt);
   // dWdI[2] = (*C3_pt)*2.0*(I[2]-1.0);
   // dWdI[3] = 0.0;//(*C4_pt)*2.0*(I[3]-1.0);

   double Q = (*C3_pt)*I[3]*I[3];
   dWdI[0] = 0.0;
   dWdI[1] = 0.0;
   dWdI[2] = (*C2_pt)/2.0*(log(I[2]) + 1.0 - 1.0/I[2]);
   dWdI[3] = (*C1_pt)*(*C3_pt)*exp(Q-1.0)*I[3];
 }

 bool requires_incompressibility_constraint()
 {
  return false;
 }

private:
 double* C1_pt;
 double* C2_pt;
 double* C3_pt;
 double* C4_pt;
};


// A nonsense strain energy function to illustrate how to use anisotropic strain energy functions with active stress
class AnisotropicActiveStressTestStrainEnergy : public AnisotropicStrainEnergyFunction
{
public:
 AnisotropicActiveStressTestStrainEnergy(double* c1_pt, double* c2_pt, double* c3_pt) : AnisotropicStrainEnergyFunction(),
                                                                                        C1_pt(c1_pt),
                                                                                        C2_pt(c2_pt),
                                                                                        C3_pt(c3_pt)
 {
  // A second principal vector of anisotropy is used to store the active stress component
  N_Principal_Vectors_Of_Anisotropy = 2;
  // A second additional strain invariant is included as a dummy invariant to store the active stress along the fibres
  N_Additional_Strain_Invariants = 2;                                                      
 }

 // Tell the constitutive law how to compute the additional strain invariants
 void I(const DenseMatrix<double>& g,
        const DenseMatrix<double>& g_up,
        const DenseMatrix<double>& G,
        const DenseMatrix<double>& G_up,
        const double& detg,
        const double& getG,
        const Vector<Vector<double>>& a,
        Vector<double>& I,
        Vector<DenseMatrix<double>>& dIdG)
 {
  const unsigned dim = g.ncol();
  // The only additional strain invariant this model uses is a_iG_{ij}a_j
  I[0] = 0.0;
  // We use a dummy strain invariant to encode the active stress along the fibres - we get the stress along the fibres from the first element in the dummy additional principal vector of anisotropy
  I[1] = a[1][0];
  for(unsigned i=0; i<dim; i++)
  {
   for(unsigned j=0; j<dim; j++)
   {
    I[0] += a[0][i] * G(i,j) * a[0][j];
    dIdG[0](i,j) = a[0][i] * a[0][j];

    // The second "strain invariant" is actually just used to store the active stress along the fibre direction
    // so it has zero derivatives with respect to the deformed metric tensor
    dIdG[1](i,j) = 0.0;
   }
  }
 }

 // Compute the strain energy - adds an additional term which represents the active stress in the direction of the fibres
 double W(const Vector<double>& I)
 {
  return (*C1_pt)*(I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0)
         + (*C3_pt + I[4]) * I[3];
 }

 /// Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 /// - adds an additional term which represents the active stress in the direction of the fibres
 void derivatives(Vector<double>& I, Vector<double>& dWdI)
 {
   dWdI[0] = (*C1_pt);
   dWdI[1] = (*C2_pt);
   dWdI[2] = 0.0;
   dWdI[3] = (*C3_pt) + I[4];
   dWdI[4] = 0.0; // This is a dummy strain invariant so it should not be included in the fill-in of the 2nd Piola-Kirchhoff stress tensor
 }

 bool requires_incompressibility_constraint()
 {
  return true;
 }

private:
 double* C1_pt;
 double* C2_pt;
 double* C3_pt;
};

};

#endif