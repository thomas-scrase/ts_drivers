#include "generic.h"

#include "ts_src.h"

#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

namespace ExactSoln
{
 // The exact solution
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = 0.0;
 }

 void source_function(const Vector<double>& x, double& source)
 {
  source = 0.0;
 }

 const double Omega = 1.0;
 const double Gamma = 0.0;

 const double Domain_Length = 1.0;
} // end of namespace


template<class ELEMENT>
class OneDMonodomainProblem : public Problem
{
public:
 OneDMonodomainProblem(const unsigned& n_element);
 ~OneDMonodomainProblem(){}
 void actions_before_newton_solve();
 void actions_after_newton_solve(){}
 void doc_solution(const unsigned& nplot);
};

template<class ELEMENT>
OneDMonodomainProblem<ELEMENT>::OneDMonodomainProblem(const unsigned& n_element)
{
 // Build the 1D mesh
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element, ExactSoln::Domain_Length);

 //Pin the values on the boundaries, at the ends of the mesh
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 mesh_pt()->boundary_node_pt(1,0)->pin(0);

 for(unsigned i=0;i<n_element;i++)
 {
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
  elem_pt->source_fct_pt() = &ExactSoln::source_function;
 }

 cout << "There are " << assign_eqn_numbers() << " degrees of freedom" << std::endl;
};

template<class ELEMENT>
void OneDMonodomainProblem<ELEMENT>::actions_before_newton_solve()
{
 Node* left_node_pt=mesh_pt()->boundary_node_pt(0,0);

 Vector<double> x(1);
 x[0] = left_node_pt->x(0);

 Vector<double> u(1);
 ExactSoln::get_exact_u(x,u);

 left_node_pt->set_value(0,u[0]);

 Node* right_node_pt=mesh_pt()->boundary_node_pt(1,0);
 
 x[0] = right_node_pt->x(0);

 ExactSoln::get_exact_u(x,u);

 right_node_pt->set_value(0,u[0]);
}

template<class ELEMENT>
void OneDMonodomainProblem<ELEMENT>::doc_solution(const unsigned& label)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solution with specified number of plot points per element
 sprintf(filename,"soln%i.dat",label);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 sprintf(filename,"exact_soln%i.dat",label);
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,20*npts,ExactSoln:get_exact_u);
 some_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 sprintf(filename,"error%i.dat",label);
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,ExactSoln::get_exact_u,
                          error,norm);
 some_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl;
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;
}


int main()
{
 unsigned n_element=40;
 OneDMonodomainProblem<QMonodomainElement<1,4>> problem(n_element);

 
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

 // Solve the problem with this Sign
 problem.newton_solve();

 //Output solution for this case (label output files with "0")
 problem.doc_solution(0);
}
