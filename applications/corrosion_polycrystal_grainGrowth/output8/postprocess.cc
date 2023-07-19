// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
    
    // Variable 0
    set_variable_name                (0,"irxn");
    set_variable_type                (0,SCALAR);
    
    set_dependencies_value_term_RHS(0, "n0, n1, n2, n3, n4, n5, psi, cM, Phi");
    set_dependencies_gradient_term_RHS(0, "");
    
    set_output_integral             (0,false);
}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
 
    //Number of order parameters
    unsigned int N=6;
    
    //Individual corrosion currents
    std::vector<scalarvalueType> icorrVi;
    icorrVi.resize(N);
    icorrVi[0] = constV(icorrV_0);
    icorrVi[1] = constV(icorrV_1);
    icorrVi[2] = constV(icorrV_2);
    icorrVi[3] = constV(icorrV_3);
    icorrVi[4] = constV(icorrV_4);
    icorrVi[5] = constV(icorrV_5);
    
    scalarvalueType psi = variable_list.get_scalar_value(0);
    
    std::vector<scalarvalueType> n;
    n.resize(N);
    
    for (unsigned int i=0; i<N; i++){
        n[i] = variable_list.get_scalar_value(i+1);
    }
    
    //cM
    scalarvalueType cM = variable_list.get_scalar_value(2*N+2);
    
    //Phi
    scalarvalueType Phi = variable_list.get_scalar_value(2*N+4);

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
// Residuals for the equation to evolve the order parameter
pp_variable_list.set_scalar_value_term_RHS(0, irxn);
}
