// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
    
    // Variable index: 0
    set_variable_name                (0,"psi");
    set_variable_type                (0,SCALAR);
    set_variable_equation_type        (0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n0, n1, n2, n3, n4, n5, psi, grad(psi), cM, Phi");
    set_dependencies_gradient_term_RHS(0, "n0, n1, n2, n3, n4, n5, psi, cM, Phi, grad(mupsi)");
    
    //Number of order parameters
    int N=6;
    
    // Variable indices: [1,N]
    for (int var_index=0; var_index<N; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));
        
        set_variable_name                (var_index+1,var_name);
        set_variable_type                (var_index+1,SCALAR);
        set_variable_equation_type       (var_index+1,EXPLICIT_TIME_DEPENDENT);
        
        set_dependencies_value_term_RHS(var_index+1, "n0, n1, n2, n3, n4, n5, grad(psi), cM, Phi");
        set_dependencies_gradient_term_RHS(var_index+1, "n0, n1, n2, n3, n4, n5, psi, cM, Phi, grad(mu0), grad(mu1), grad(mu2), grad(mu3), grad(mu4), grad(mu5)");
        
    }

    // Variable indices: [N+1, N+N-1+1] or [N+1, 2N]
    for (int var_index=0; var_index<N; var_index++){
        std::string var_name2 = "mu";
        var_name2.append(std::to_string(var_index));
        
        set_variable_name                (N+var_index+1,var_name2);
        set_variable_type                (N+var_index+1,SCALAR);
        set_variable_equation_type       (N+var_index+1,EXPLICIT_TIME_DEPENDENT);
        
        set_dependencies_value_term_RHS(N+var_index+1, "n0, n1, n2, n3, n4, n5, psi");
        set_dependencies_gradient_term_RHS(N+var_index+1, "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5)");
    }

    // Variable index: 2N+1
    set_variable_name                (2*N+1,"mupsi");
    set_variable_type                (2*N+1,SCALAR);
    set_variable_equation_type        (2*N+1,EXPLICIT_TIME_DEPENDENT);
    
    set_dependencies_value_term_RHS(2*N+1, "n0, n1, n2, n3, n4, n5, psi");
    set_dependencies_gradient_term_RHS(2*N+1, "grad(psi)");
    
    // Variable index: 2N+2
    set_variable_name				(2*N+2,"cM");
    set_variable_type				(2*N+2,SCALAR);
    set_variable_equation_type		(2*N+2,EXPLICIT_TIME_DEPENDENT);
    
    set_dependencies_value_term_RHS(2*N+2, "n0, n1, n2, n3, n4, n5, cM, grad(cM), psi, grad(psi), Phi, grad(Phi)");
    set_dependencies_gradient_term_RHS(2*N+2, "cM, grad(cM), grad(Phi)");
    
    // Variable index: 2N+3
    set_variable_name				(2*N+3,"cP");
    set_variable_type				(2*N+3,SCALAR);
    set_variable_equation_type		(2*N+3,EXPLICIT_TIME_DEPENDENT);
    
    set_dependencies_value_term_RHS(2*N+3, "cP, grad(cP), psi, grad(psi), grad(Phi)");
    set_dependencies_gradient_term_RHS(2*N+3, "cP, grad(cP), grad(Phi)");
    
    
    // Variable index: 2N+4
    set_variable_name				(2*N+4,"Phi");
    set_variable_type				(2*N+4,SCALAR);
    set_variable_equation_type		(2*N+4,TIME_INDEPENDENT);
    
    set_dependencies_value_term_RHS(2*N+4, "n0, n1, n2, n3, n4, n5, cM, grad(psi), Phi");
    set_dependencies_gradient_term_RHS(2*N+4, "psi, grad(Phi), cM, cP, grad(cM), grad(cP)");
    set_dependencies_value_term_LHS(2*N+4, "n0, n1, n2, n3, n4, n5, grad(psi), cM, cP, Phi, change(Phi)");
    set_dependencies_gradient_term_LHS(2*N+4, "psi, grad(change(Phi))");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    //Number of order parameters
    unsigned int N=6;
    
    //Individual corrosion currents
    std::vector<scalarvalueType> icorrVi;
    icorrVi.resize(N);
    //Defining corrosion currents for each order parameter
    icorrVi[0] = constV(icorrV_0);
    icorrVi[1] = constV(icorrV_1);
    icorrVi[2] = constV(icorrV_2);
    icorrVi[3] = constV(icorrV_3);
    icorrVi[4] = constV(icorrV_4);
    icorrVi[5] = constV(icorrV_5);
    
    //Initializing fields for order parameters, chemical potentials and their derivatives
    std::vector<scalarvalueType> n;
   	n.resize(N);
    std::vector<scalargradType> nx;
    nx.resize(N);
    //std::vector<scalarvalueType> mu;
    //mu.resize(N);
    std::vector<scalargradType> mux;
    mux.resize(N);
    
    //Timestep
    scalarvalueType delt = constV(userInputs.dtValue);
    
    // The domain parameter and its derivatives
    scalarvalueType psi = variable_list.get_scalar_value(0);
    scalargradType psix = variable_list.get_scalar_gradient(0);
    
    for (unsigned int i=0; i<N; i++){
        n[i] = variable_list.get_scalar_value(i+1);
        nx[i] = variable_list.get_scalar_gradient(i+1);
        //mu[i] = variable_list.get_scalar_value(N+i);
        mux[i] = variable_list.get_scalar_gradient(N+i+1);
    }
    
    
    //The domain parameter chemical potential and it's derivatives
    //scalarvalueType mupsi = variable_list.get_scalar_value(2*N+1);
    scalargradType mupsix = variable_list.get_scalar_gradient(2*N+1);
    
    // The concentration of metal ion and its derivatives
    scalarvalueType cM = variable_list.get_scalar_value(2*N+2);
    scalargradType cMx = variable_list.get_scalar_gradient(2*N+2);
    
    // The concentration of supporting cation and its derivatives
    scalarvalueType cP = variable_list.get_scalar_value(2*N+3);
    scalargradType cPx = variable_list.get_scalar_gradient(2*N+3);
    
    // The electrolite potential and its derivatives
    scalarvalueType Phi = variable_list.get_scalar_value(2*N+4);
    scalargradType Phix = variable_list.get_scalar_gradient(2*N+4);
    
 	//Capping order parameter and domain parameter
    //Capping n to lower threshold bound and upper bound of 1
    
    //Initializing Capped fields
    std::vector<scalarvalueType> ncp;
    ncp.resize(N);
    for (unsigned int i=0; i<N; i++){
        ncp[i] = constV(0.0);
    }
    scalarvalueType psicp = constV(0.0);

    for (unsigned int j=0; j<psi.n_array_elements;j++){
        for (unsigned int i=0; i<N; i++){
            ncp[i][j]=n[i][j];
            if (n[i][j] < lthresh)
                ncp[i][j] = lthresh;
            if (n[i][j] > 1.0)
                ncp[i][j] = 1.0;
        }
            psicp[j]=psi[j];
            if (psi[j] < lthresh)
                psicp[j] = lthresh;
            if (psi[j] > 1.0)
                psicp[j] = 1.0;
    }
    n = ncp;
    psi = psicp;
    //End of capping
    
    //Calculating bulk part of free energy
    //Initializing capped fields
    std::vector<scalarvalueType> fnV;
    fnV.resize(N);
    
    //Calculating the i-th order parameter contribution to fnV
    scalarvalueType part_nsq_sum = constV(0.0);
    for (unsigned int i=0; i<N; i++){
        part_nsq_sum = constV(0.0);
        	for (unsigned int j=0; j<N; j++){
                if (j != i){
                    part_nsq_sum = part_nsq_sum + n[j]*n[j];
                }
            }
        fnV[i] = constV(W)*n[i]*(n[i]*n[i] - constV(1.0) + constV(2.0*gamma)*(part_nsq_sum + psi*psi));
    }
    
    //Calculating the fpsiV
    scalarvalueType fpsiV = constV(0.0);
    scalarvalueType nsq_sum = constV(0.0);
    for (unsigned int i=0; i<N; i++){
        nsq_sum = nsq_sum + n[i]*n[i];
    }
    fpsiV = constV(W)*psi*(psi*psi - constV(1.0) + constV(2.0*gamma)*nsq_sum);
    
    //Magnitude of the gradient of the domain parameter
    scalarvalueType magpsix = constV(0.0);
    for (int i=0; i<dim; i++){
        magpsix = magpsix + psix[i]*psix[i];
    }
    magpsix = std::sqrt(magpsix);
    
    //Miscellaneous dot products
    scalarvalueType psixcMx = constV(0.0);
    scalarvalueType psixcPx = constV(0.0);
    scalarvalueType psixcMPhix = constV(0.0);
    scalarvalueType psixcPPhix = constV(0.0);
    
    for (int i=0; i<dim; i++){
        psixcMx = psixcMx + psix[i]*cMx[i];
        psixcPx = psixcPx + psix[i]*cPx[i];
        psixcMPhix = psixcMPhix + psix[i]*cM*Phix[i];
        psixcPPhix = psixcPPhix + psix[i]*cP*Phix[i];
    }
    
    //Inverse of psi
    scalarvalueType invpsi=constV(1.0)/psi;
    
    //Overpotential
    scalarvalueType eta=constV(VsV-EcorrV)-Phi;

    //Calculation of reaction current
    //Calculating weight factors (xi_i) of order parameters
    scalarvalueType sum_op = constV(0.0);
    for (unsigned int i=0; i<N; i++){
        sum_op = sum_op + n[i];
    }
    
    //Capping the smallest value of sum_op to avoid division by zero
    scalarvalueType sum_op_cp= constV(0.0);
    for (unsigned int j=0; j<psi.n_array_elements;j++){
        sum_op_cp[j] = sum_op[j];
        if (sum_op[j] < lthresh)
            sum_op_cp[j] = lthresh;
    }
    sum_op = sum_op_cp;
    
    std::vector<scalarvalueType> xi;
    xi.resize(N);
    for (unsigned int i=0; i<N; i++){
        xi[i] =n[i]/sum_op;
    }
                               
    //Reaction Current
    scalarvalueType irxni = constV(0.0);
    scalarvalueType irxn = constV(0.0);
    
    scalarvalueType irxn_gb = constV(0.0);    //reaction current at grain boundaries (used only if n[i]*n[i]*n[j]*n[j] != 0)
    scalarvalueType c_gb = constV(alpha);
    //Total reaction current
    for (unsigned int i=0; i<N; i++){
        irxni=icorrVi[i]*(constV(1.0)-cM/constV(cMsatV))*std::exp(constV(zMV*(1.0-betaV)*FarC/(RV*TV))*eta);
        for (unsigned int j = i+1; j<N; j++) {
            irxn_gb += xi[i]*xi[i]*xi[j]*xi[j];
        }
        irxn = irxn + irxni*xi[i] + irxn_gb*c_gb;
    }
    
    //Velocity Field
    scalarvalueType v = -constV(VMV/(zMV*FarC))*irxn;
    
    //Interface mobility
    scalarvalueType MnV = 2.0*(VMV*irxn/(zMV*FarC))*std::sqrt(2.0*epssqV);
    
    // The residuals
    std::vector<scalarvalueType> rnV;
    rnV.resize(N);
    std::vector<scalargradType> rnxV;
    rnxV.resize(N);
    std::vector<scalarvalueType> rmuV;
    rmuV.resize(N);
    std::vector<scalargradType> rmuxV;
    rmuxV.resize(N);
    
    for (unsigned int i=0; i<N; i++){
        rnV[i] = n[i] + v*delt*magpsix;
        rnxV[i] = -MnV*psi*delt*mux[i];
        
        rmuV[i] = fnV[i];
        rmuxV[i] = constV(epssqV)*nx[i];
    }
    
    scalarvalueType rpsiV = psi - v*delt*magpsix;
    scalargradType rpsixV = -MnV*psi*delt*mupsix;
    
    scalarvalueType rmupsiV = fpsiV;
    scalargradType rmupsixV = constV(epssqV)*psix;
    
    scalarvalueType rcMV = cM + delt*(constV(DMV)*invpsi*psixcMx + constV(DMV*zMV*FarC/(RV*TV))*psixcMPhix*invpsi + constV(1.0/(zMV*FarC))*invpsi*magpsix*irxn);
    scalargradType rcMxV = -delt*(DMV*cMx + constV(DMV*zMV*FarC/(RV*TV))*cM*Phix);
    
    scalarvalueType rcPV = cP + delt*(constV(DPV)*invpsi*psixcPx + constV(DPV*zPV*FarC/(RV*TV))*psixcPPhix*invpsi);
    scalargradType rcPxV = -delt*(DPV*cPx + constV(DPV*zPV*FarC/(RV*TV))*cP*Phix);
    
    // Residuals for the equation to evolve the domain parameter
    variable_list.set_scalar_value_term_RHS(0,rpsiV);
    variable_list.set_scalar_gradient_term_RHS(0,rpsixV);
    
    for (unsigned int i=0; i<N; i++){
        // Residuals for the equation to evolve the order parameter
        variable_list.set_scalar_value_term_RHS(i+1,rnV[i]);
        variable_list.set_scalar_gradient_term_RHS(i+1,rnxV[i]);
        
        // Residuals for the equation to evolve the chemical potential
        variable_list.set_scalar_value_term_RHS(N+i+1,rmuV[i]);
        variable_list.set_scalar_gradient_term_RHS(N+i+1,rmuxV[i]);
    }
    
    // Residuals for the equation to evolve the chemical potential of the domain parameter
    variable_list.set_scalar_value_term_RHS(2*N+1,rmupsiV);
    variable_list.set_scalar_gradient_term_RHS(2*N+1,rmupsixV);
    
    // Residuals for the equation to evolve the concentration of metal ion
    variable_list.set_scalar_value_term_RHS(2*N+2,rcMV);
    variable_list.set_scalar_gradient_term_RHS(2*N+2,rcMxV);
    
    // Residuals for the equation to evolve the concentration of supporting cation
    variable_list.set_scalar_value_term_RHS(2*N+3,rcPV);
    variable_list.set_scalar_gradient_term_RHS(2*N+3,rcPxV);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    //Number of order parameters
    unsigned int N=6;
    
    //Individual corrosion currents
    std::vector<scalarvalueType> icorrVi;
    icorrVi.resize(N);
    //Defining corrosion currents for each order parameter
    icorrVi[0] = constV(icorrV_0);
    icorrVi[1] = constV(icorrV_1);
    icorrVi[2] = constV(icorrV_2);
    icorrVi[3] = constV(icorrV_3);
    icorrVi[4] = constV(icorrV_4);
    icorrVi[5] = constV(icorrV_5);
    
    
 // --- Getting the values and derivatives of the model variables ---
    
    //Initializing fields for order parameters, chemical potentials and their derivatives
    
    // The domain parameter and its derivatives
    scalarvalueType psi = variable_list.get_scalar_value(0);
    scalargradType psix = variable_list.get_scalar_gradient(0);
    
    std::vector<scalarvalueType> n;
    n.resize(N);
    
    for (unsigned int i=0; i<N; i++){
        n[i] = variable_list.get_scalar_value(i+1);
    }
    
    // The concentration of metal ion and its derivatives
    scalarvalueType cM = variable_list.get_scalar_value(2*N+2);
    scalargradType cMx = variable_list.get_scalar_gradient(2*N+2);
    
    // The concentration of supporting cation and its derivatives
    scalarvalueType cP = variable_list.get_scalar_value(2*N+3);
    scalargradType cPx = variable_list.get_scalar_gradient(2*N+3);
    
    // The electrolite potential and its derivatives
    scalarvalueType Phi = variable_list.get_scalar_value(2*N+4);
    scalargradType Phix = variable_list.get_scalar_gradient(2*N+4);
    
    //Capping order parameter and domain parameter
    //Capping n to lower threshold bound and upper bound of 1

    //Initializing Capped fields
    std::vector<scalarvalueType> ncp;
    ncp.resize(N);
    for (unsigned int i=0; i<N; i++){
        ncp[i] = constV(0.0);
    }
    scalarvalueType psicp = constV(0.0);
    
    for (unsigned int j=0; j<psi.n_array_elements;j++){
        for (unsigned int i=0; i<N; i++){
            ncp[i][j]=n[i][j];
            if (n[i][j] < lthresh)
                ncp[i][j] = lthresh;
            if (n[i][j] > 1.0)
                ncp[i][j] = 1.0;
        }
        psicp[j]=psi[j];
        if (psi[j] < lthresh)
            psicp[j] = lthresh;
        if (psi[j] > 1.0)
            psicp[j] = 1.0;
    }
    n = ncp;
    psi = psicp;
    //End of capping
    
    //Overpotential
    scalarvalueType eta = constV(VsV-EcorrV) - Phi;
    
    //Calculation of reaction current
    //Calculating weight factors (xi_i) of order parameters
    scalarvalueType sum_op = constV(0.0);
    for (unsigned int i=0; i<N; i++){
        sum_op = sum_op + n[i];
    }
    
    //Capping the smallest value of sum_op to avoid division by zero
    scalarvalueType sum_op_cp= constV(0.0);
    for (unsigned int j=0; j<psi.n_array_elements;j++){
        sum_op_cp[j] = sum_op[j];
        if (sum_op[j] < lthresh)
            sum_op_cp[j] = lthresh;
    }
    sum_op = sum_op_cp;
    
    std::vector<scalarvalueType> xi;
    xi.resize(N);
    for (unsigned int i=0; i<N; i++){
        xi[i] =n[i]/sum_op;
    }
    
    //Reaction Current
    scalarvalueType irxni = constV(0.0);
    scalarvalueType irxn = constV(0.0);
    
    scalarvalueType irxn_gb = constV(0.0);    //reaction current at grain boundaries (used only if n[i]*n[i]*n[j]*n[j] != 0)
    scalarvalueType c_gb = constV(alpha);
    //Total reaction current
    for (unsigned int i=0; i<N; i++){
        irxni=icorrVi[i]*(constV(1.0)-cM/constV(cMsatV))*std::exp(constV(zMV*(1.0-betaV)*FarC/(RV*TV))*eta);
        for (unsigned int j = i+1; j<N; j++) {
            irxn_gb += xi[i]*xi[i]*xi[j]*xi[j];
        }
        irxn = irxn + irxni*xi[i] + irxn_gb*c_gb;
    }
   
    //Magnitude of the gradient of the domain parameter
    scalarvalueType magpsix = constV(0.0);
    for (int i=0; i<dim; i++){
        magpsix = magpsix + psix[i]*psix[i];
    }
    magpsix = std::sqrt(magpsix);
    
    scalarvalueType kappa = constV(FarC*FarC/(RV*TV))*(constV(zMV*(zMV*DMV-znV*DnV))*cM + constV(zPV*(zPV*DPV-znV*DnV))*cP);
    scalarvalueType Msum = constV(FarC*zMV*(DnV-DMV));
    scalarvalueType Psum = constV(FarC*zPV*(DnV-DPV));
    
    scalarvalueType rPhi = -magpsix*irxn;
    scalargradType rPhix = psi*kappa*Phix -  Msum*psi*cMx - Psum*psi*cPx;
    
    // Residuals
    variable_list.set_scalar_value_term_RHS(2*N+4,rPhi);
    variable_list.set_scalar_gradient_term_RHS(2*N+4,rPhix);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    //Number of order parameters
    unsigned int N=6;
    // --- Getting the values and derivatives of the model variables ---

    //Individual corrosion currents
    std::vector<scalarvalueType> icorrVi;
    icorrVi.resize(N);
    //Defining corrosion currents for each order parameter
    icorrVi[0] = constV(icorrV_0);
    icorrVi[1] = constV(icorrV_1);
    icorrVi[2] = constV(icorrV_2);
    icorrVi[3] = constV(icorrV_3);
    icorrVi[4] = constV(icorrV_4);
    icorrVi[5] = constV(icorrV_5);
    
    //Initializing fields for order parameters, chemical potentials and their derivatives
    
    // The domain parameter and its derivatives
    scalarvalueType psi = variable_list.get_scalar_value(0);
    scalargradType psix = variable_list.get_scalar_gradient(0);
    
    std::vector<scalarvalueType> n;
    n.resize(N);
    
    for (unsigned int i=0; i<N; i++){
        n[i] = variable_list.get_scalar_value(i+1);
    }
    
    // The concentration of metal ion and its derivatives
    scalarvalueType cM = variable_list.get_scalar_value(2*N+2);
    
    // The concentration of supporting cation and its derivatives
    scalarvalueType cP = variable_list.get_scalar_value(2*N+3);
    
    // The electrolite potential and its derivatives
    scalarvalueType Phi = variable_list.get_scalar_value(2*N+4);
    
    // The change in potential in the electrode and its derivatives
    scalarvalueType DPhi = variable_list.get_change_in_scalar_value(2*N+4);
    scalargradType DPhix = variable_list.get_change_in_scalar_gradient(2*N+4);

    //Capping order parameter and domain parameter
    //Capping n to lower threshold bound and upper bound of 1
    
    //Initializing Capped fields
    std::vector<scalarvalueType> ncp;
    ncp.resize(N);
    for (unsigned int i=0; i<N; i++){
        ncp[i] = constV(0.0);
    }
    scalarvalueType psicp = constV(0.0);
    
    for (unsigned int j=0; j<psi.n_array_elements;j++){
        for (unsigned int i=0; i<N; i++){
            ncp[i][j]=n[i][j];
            if (n[i][j] < lthresh)
                ncp[i][j] = lthresh;
            if (n[i][j] > 1.0)
                ncp[i][j] = 1.0;
        }
        psicp[j]=psi[j];
        if (psi[j] < lthresh)
            psicp[j] = lthresh;
        if (psi[j] > 1.0)
            psicp[j] = 1.0;
    }
    n = ncp;
    psi = psicp;
    //End of capping

    scalarvalueType kappa = constV(FarC*FarC/(RV*TV))*(constV(zMV*(zMV*DMV-znV*DnV))*cM + constV(zPV*(zPV*DPV-znV*DnV))*cP);

    //Magnitude of the gradient of the domain parameter
    scalarvalueType magpsix = constV(0.0);
    for (int i=0; i<dim; i++){
        magpsix = magpsix + psix[i]*psix[i];
    }
    magpsix = std::sqrt(magpsix);
    
    //Overpotential
    scalarvalueType eta=constV(VsV-EcorrV)-Phi;
    
    //Calculation of the derivative of the reaction current with respect to the potential
    //Calculating weight factors (xi_i) of order parameters
    scalarvalueType sum_op = constV(0.0);
    for (unsigned int i=0; i<N; i++){
        sum_op = sum_op + n[i];
    }
    
    //Capping the smallest value of sum_op to avoid division by zero
    scalarvalueType sum_op_cp= constV(0.0);
    for (unsigned int j=0; j<psi.n_array_elements;j++){
        sum_op_cp[j] = sum_op[j];
        if (sum_op[j] < lthresh)
            sum_op_cp[j] = lthresh;
    }
    sum_op = sum_op_cp;
    
    std::vector<scalarvalueType> xi;
    xi.resize(N);
    for (unsigned int i=0; i<N; i++){
        xi[i] =n[i]/sum_op;
    }

    //Reaction Current derivarive
    scalarvalueType irxnip = constV(0.0);
    scalarvalueType irxnp = constV(0.0);
    //Total reaction current
    for (unsigned int i=0; i<N; i++){
        irxnip = -icorrVi[i]*(constV(1.0)-cM/constV(cMsatV))*constV(zMV*(1.0-betaV)*FarC/(RV*TV))*std::exp(constV(zMV*(1.0-betaV)*FarC/(RV*TV))*eta);
        irxnp = irxnp + irxnip*xi[i];
    }

    scalarvalueType rDPhi = magpsix*irxnp*DPhi;
    scalargradType rDPhix = -psi*kappa*DPhix;
    
    // Residuals
    variable_list.set_scalar_value_term_LHS(2*N+4,rDPhi);
    variable_list.set_scalar_gradient_term_LHS(2*N+4,rDPhix);
}


/*
// =================================================================================
// thresholdField: a function particular to this app
// =================================================================================
//Method that caps the value of the order parameter and the domain parameter
template <int dim,int degree>
void customPDE<dim,degree>::capFields(dealii::VectorizedArray<double> & ncp, dealii::VectorizedArray<double> & psicp,
                                      dealii::VectorizedArray<double> n, dealii::VectorizedArray<double> psi) const {
    //Capping n to lower threshold bound and upper bound of 1
    for (int j=0; j<ncp.n_array_elements;j++){
        ncp[j]=n[j];
        if (n[j] < lthresh)
            ncp[j] = lthresh;
        if (n[j] > 1.0)
            ncp[j] = 1.0;
    }
    //Capping psi to lower threshold bound and upper bound of 1
    for (int j=0; j<ncp.n_array_elements;j++){
        psicp[j]=psi[j];
        if (psi[j] < lthresh)
            psicp[j] = lthresh;
        if (psi[j] > 1.0)
            psicp[j] = 1.0;
    }
}
*/
