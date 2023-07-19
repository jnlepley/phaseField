// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){

    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index
    
    // The initial condition is a set of overlapping circles/spheres defined
    // by a hyperbolic tangent function. The center of each circle/sphere is
    // given by "center" and its radius is given by "radius".
    
    /*
     List of scalar variables
     Index: 0; var. n0
     Index: 1; var. n1
     Index: 2; var. n2
     Index: 3; var. n3
     Index: 4; var. n4
     Index: 5; var. n5
     Index: 6; var. psi
     Index: 7; var. mu0
     Index: 8; var. mu1
     Index: 9; var. mu2
     Index: 10; var. mu3
     Index: 11; var. mu4
     Index: 12; var. mu5
     Index: 13; var. mupsi
     Index: 14; var. cM
     Index: 15; var. cP
     Index: 16; var. Phi
    */
    
    unsigned int N=6;
    
//    double posx=p[0];
//    double posy=p[1];
//    double cx=0.5*userInputs.domain_size[0];
//    double cy=userInputs.domain_size[1];
//    double rad=std::sqrt((posx-cx)*(posx-cx) + (posy-cy)*(posy-cy));
//    double psipro = 0.5*(1.0 + std::tanh((0.025*userInputs.domain_size[0] - rad)/deltaV));
    
//    if (index == 2*N){
//        scalar_IC = psipro;
    /*} else*/
    
    if (index == 2*N+3){
        scalar_IC = 1000.0;
    } else {
        scalar_IC = 0.0;
    }
    // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).
	/*
    if (index == 2){
        if (direction == 1){
            double x=p[0];
            double y=p[1];
            scalar_BC=std::sin(y/7.0);
        }
    }
	*/
    // -------------------------------------------------------------------------
}

