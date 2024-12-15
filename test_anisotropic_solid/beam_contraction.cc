#include "generic.h"

#include "anisotropic_solid.h"
#include "anisotropic_test_law.h"

#include "constitutive.h"

#include "meshes/simple_cubic_mesh.h"

using namespace std;
using namespace oomph;

// Define the solid mesh
template<class ELEMENT>
class ElasticSimpleCubicMesh : public virtual SimpleCubicMesh<ELEMENT>,
                               public virtual SolidMesh
{
public:
 ElasticSimpleCubicMesh(const unsigned nx, const unsigned ny, const unsigned nz,
                        TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper) :
                        SimpleCubicMesh<ELEMENT>(nx, ny, nz, 4.0, 1.0, 1.0, time_stepper_pt),
                        SolidMesh()
 {
  set_lagrangian_nodal_coordinates();
 }
};


namespace Parameters
{
 // Define the principal vector of anisotropy.
 // In the x-y plane and defined by some angle from the positive x direction.
 // Additional dummy principal vector of anisotropy to contain the active stress in the fibre direction
 double Angle = 0.0;
 double Active_stress = 0.0;
 void principal_vectors_of_anisotropy(const Vector<double>& s,
                                      const Vector<double>& xi,
                                      Vector<Vector<double>>& a)
 {
  // There are 2 principal vectors of anisotropic
  a.resize(2, Vector<double>(3,0.0));
  // The fibres are aligned with the x axis
  a[0][0] = cos(Angle);
  a[0][1] = sin(Angle);
  a[0][2] = 0.0;

  // The active stress has only one component - that in the direction of the fibres
  a[1][0] = Active_stress;
  a[1][1] = 0.0;
  a[1][2] = 0.0;
 }

 double Lambda_sq = 1.0;

 double C1 = 0.25;
 double C2 = 1.0;
 double C3 = 1.0; // Scaling of the stress generated due to deformation in the fibre direction
 double C4 = 1.0; // Scaling of the active stress component
 // The pointers to the strain energy function and the constitutive law
 AnisotropicStrainEnergyFunction* anisotropic_strain_energy_function_pt;
 AnisotropicConstitutiveLaw* anisotropic_constutitive_law_pt;
}

template<class ELEMENT>
class CubicAnisotropicSolidProblem : public Problem
{
public:
 CubicAnisotropicSolidProblem(const unsigned& nx, const unsigned& ny, const unsigned& nz);
 ~CubicAnisotropicSolidProblem() {}
 void actions_before_newton_solve(){}
 void actions_after_newton_solve(){}
 void doc_solution(const unsigned& nplot);
};

template<class ELEMENT>
CubicAnisotropicSolidProblem<ELEMENT>::CubicAnisotropicSolidProblem(const unsigned& nx, const unsigned& ny, const unsigned& nz)
{
 Problem::max_residuals() = 100.0;
 Problem::max_newton_iterations() = 50;

 Problem::mesh_pt() = new ElasticSimpleCubicMesh<ELEMENT>(nx, ny, nz);

 Parameters::anisotropic_strain_energy_function_pt = new ActiveAnisotropicGeneralisedMooneyRivlin(&Parameters::C1, &Parameters::C2, &Parameters::C3, &Parameters::C4);
 Parameters::anisotropic_constutitive_law_pt = new AnisotropicStrainEnergyFunctionConstitutiveLaw(Parameters::anisotropic_strain_energy_function_pt);

 const unsigned n_element = Problem::mesh_pt()->nelement();
 for(unsigned i=0; i<n_element; i++)
 {
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Problem::mesh_pt()->element_pt(i));
  elem_pt->anisotropic_constitutive_law_pt() = Parameters::anisotropic_constutitive_law_pt;
  // elem_pt->body_force_fct_pt() = &Parameters::body_force;
  elem_pt->principal_vectors_of_anisotropy_fct_pt() = &Parameters::principal_vectors_of_anisotropy;
  elem_pt->lambda_sq_pt() = &Parameters::Lambda_sq;
 }

 // Apply boundary conditions - pin all of them on one side
 const unsigned dirichlet_boundary = 4;
 const unsigned n_bound_node = Problem::mesh_pt()->nboundary_node(dirichlet_boundary);
 for(unsigned i=0; i<n_bound_node;i++)
 {
  SolidNode* node_pt = dynamic_cast<SolidNode*>(Problem::mesh_pt()->boundary_node_pt(dirichlet_boundary,i));

  for(unsigned i=0; i<3; i++)
  {
   node_pt->variable_position_pt()->pin(i);
  }
 }

 cout << "There are " << assign_eqn_numbers() << " degrees of freedom" << std::endl;
}

template<class ELEMENT>
void CubicAnisotropicSolidProblem<ELEMENT>::doc_solution(const unsigned& label)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solution with specified number of plot points per element
 sprintf(filename,"RESLT_BEAM_CONTRACTION/soln%i.dat", label);
 some_file.open(filename);
 Problem::mesh_pt()->output(some_file, npts);
 some_file.close();
}

int main()
{
 unsigned nx = 5;
 unsigned ny = 5;
 unsigned nz = 5;
 
 CubicAnisotropicSolidProblem<AnisotropicQPVDElement<3, 3>> problem(nx, ny, nz);

 // Check whether the problem can be solved
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0)
 {
  cout << "passed: Problem can be solved." << std::endl;
 }
 else
 {
  throw OomphLibError("failed!",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }

 //Output solution
 unsigned out_num = 0;
 problem.doc_solution(out_num++);

 // Increase then decrease the active stress component
 unsigned N_Active = 40;
 double Active_Magnitude = 0.5;
 for(unsigned i=1; i<N_Active; i++)
 {
  Parameters::Active_stress = Active_Magnitude * sin((double)i/(double)N_Active * 2.0 * MathematicalConstants::Pi);
  // Solve the problem with this Sign
  problem.newton_solve();
  //Output solution
  problem.doc_solution(out_num++);
 }
}