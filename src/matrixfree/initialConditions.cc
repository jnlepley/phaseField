//methods to apply initial conditions

#include "../../include/matrixFreePDE.h"
#include "../../include/initialConditions.h"
#include "../../include/IntegrationTools/PField.hh"
#include "../../include/OrderParameterRemapper.h"


template <int dim>
class InitialConditionPField : public Function<dim>
{
public:
  unsigned int index;
  Vector<double> values;
  typedef PRISMS::PField<double*, double, dim> ScalarField;
  ScalarField &inputField;

  InitialConditionPField (const unsigned int _index, ScalarField &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}

  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
	  double scalar_IC;

	  double coord[dim];
	  for (unsigned int i = 0; i < dim; i++){
		  coord[i] = p(i);
	  }

	  scalar_IC = inputField(coord);

	  return scalar_IC;
  }
};


//methods to apply initial conditions
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::applyInitialConditions(){

    if (userInputs.load_grain_structure){
        // Create the dummy field
        vectorType grain_index_field;

        // Clear the order parameter fields
        unsigned int op_list_index = 0;
        for (unsigned int var_index=0; var_index < userInputs.number_of_variables; var_index++){
            if (op_list_index < userInputs.variables_for_remapping.size()){
                if (var_index == userInputs.variables_for_remapping.at(op_list_index)){
                    *solutionSet[var_index] = 0.0;
                    op_list_index++;
                }
            }
        }

        simplified_grain_representations.clear();

        // Get the index of one of the scalar fields
        unsigned int scalar_field_index = 0;
        for (unsigned int var=0; var<userInputs.number_of_variables; var++){
            if (userInputs.var_type.at(var) == SCALAR){
                scalar_field_index = var;
                break;
            }
        }

        matrixFreeObject.initialize_dof_vector (grain_index_field, scalar_field_index);

        // Declare the PField types and containers
        typedef PRISMS::PField<double*, double, dim> ScalarField;
        typedef PRISMS::Body<double*, dim> Body;
        Body body;

        // Create the filename of the the file to be loaded

        // Load the data from the file using a PField
        std::string filename = userInputs.grain_structure_filename;
        filename += ".vtk";

        body.read_vtk(filename);
        ScalarField &id_field = body.find_scalar_field(userInputs.grain_structure_variable_name);

        pcout << "Applying PField initial condition...\n";

        VectorTools::interpolate (*dofHandlersSet[scalar_field_index], InitialConditionPField<dim>(0,id_field), grain_index_field);

        grain_index_field.update_ghost_values();

        // Get the max and min grain ids
        unsigned int max_id = (unsigned int)grain_index_field.linfty_norm();
        unsigned int min_id = 0;

        // Now locate all of the grains and create simplified representations of them
        QGaussLobatto<dim> quadrature2 (degree+1);

        //////////////////////////
        // *FESet.at(scalar_field_index) is initialized here, and passed in again (duplicated) on
        // line 110 (103). (its just a pointer so it shouldn't mater)
        // I know FESet is a vector of pointer-to-FESystem objects,
        // but I'm not sure how that vector is structured. What is .at(scalar_field_index) "pointing"
        // to?
        //////////////////////////

        FloodFiller<dim, degree> flood_filler(*FESet.at(scalar_field_index), quadrature2);

        pcout << "Locating the grains...\n";

        // Make a new vector of GrainSet objects
        std::vector<GrainSet<dim>> grain_sets;
        for (unsigned int id=min_id; id<max_id+1; id++){
            pcout << "Locating grain " << id << "...\n";

            // For each feature ID, make a new vector of GrainSets (that all will have the same grain ID)
            std::vector<GrainSet<dim>> grain_sets_single_id;

            // Calculate the grain sets that share this ID
            flood_filler.calcGrainSets(*FESet.at(scalar_field_index), *dofHandlersSet_nonconst.at(scalar_field_index), &grain_index_field, (double)id - userInputs.order_parameter_threshold, (double)id + userInputs.order_parameter_threshold, 0, grain_sets_single_id);

            for (unsigned int g=0; g<grain_sets_single_id.size(); g++){
                grain_sets_single_id.at(g).setGrainIndex(id);
            }

            grain_sets.insert(grain_sets.end(), grain_sets_single_id.begin(), grain_sets_single_id.end());
        }

        pcout << "Generating simplified representations of the grains...\n";
        for (unsigned int g=0; g<grain_sets.size(); g++){
            SimplifiedGrainRepresentation<dim> simplified_grain_representation(grain_sets.at(g));

            if (dim == 2){
                pcout << "Grain: " << simplified_grain_representation.getGrainId() << " " << simplified_grain_representation.getOrderParameterId() << " Center: " << simplified_grain_representation.getCenter()(0) << " " << simplified_grain_representation.getCenter()(1) << "  Radius: " << simplified_grain_representation.getRadius() << std::endl;
            }
            else {
                pcout << "Grain: " << simplified_grain_representation.getGrainId() << " " << simplified_grain_representation.getOrderParameterId() << " Center: " << simplified_grain_representation.getCenter()(0) << " " << simplified_grain_representation.getCenter()(1) << " " << simplified_grain_representation.getCenter()(2) << "  Radius: " << simplified_grain_representation.getRadius() << std::endl;
            }

            simplified_grain_representations.push_back(simplified_grain_representation);
        }

        // Delete grains with very small radii that correspond to a single interfacial element
        for (unsigned int g=0; g<simplified_grain_representations.size(); g++){
            if (simplified_grain_representations.at(g).getRadius() < userInputs.min_radius_for_loading_grains){
                simplified_grain_representations.erase(simplified_grain_representations.begin()+g);
                g--;
            }
        }

        pcout << "Reassigning the grains to new order parameters...\n";
        SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
        simplified_grain_manipulator.reassignGrains(simplified_grain_representations, userInputs.buffer_between_grains, userInputs.variables_for_remapping);

        pcout << "After reassignment: " << std::endl;
        for (unsigned int g=0; g<simplified_grain_representations.size(); g++){
            if (dim == 2){
                pcout << "Grain: " << simplified_grain_representations.at(g).getGrainId() << " " << simplified_grain_representations.at(g).getOrderParameterId() << " Center: " << simplified_grain_representations.at(g).getCenter()(0) << " " << simplified_grain_representations.at(g).getCenter()(1) << std::endl;
            }
            else {
                pcout << "Grain: " << simplified_grain_representations.at(g).getGrainId() << " " << simplified_grain_representations.at(g).getOrderParameterId() << " Center: " << simplified_grain_representations.at(g).getCenter()(0) << " " << simplified_grain_representations.at(g).getCenter()(1) << " " << simplified_grain_representations.at(g).getCenter()(2) << std::endl;
            }

        }

        pcout << "Placing the grains in their new order parameters...\n";
        OrderParameterRemapper<dim> order_parameter_remapper;
        order_parameter_remapper.remap_from_index_field(simplified_grain_representations, &grain_index_field, solutionSet, *dofHandlersSet_nonconst.at(scalar_field_index), FESet.at(scalar_field_index)->dofs_per_cell, userInputs.buffer_between_grains);

        // Smooth the order parameters
        double dt_for_smoothing = dealii::GridTools::minimal_cell_diameter(triangulation)/1000.0;

        op_list_index = 0;
        for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
            if (op_list_index < userInputs.variables_for_remapping.size()){
                if ( fieldIndex == userInputs.variables_for_remapping.at(op_list_index)){

                    for (unsigned int cycle=0; cycle<userInputs.num_grain_smoothing_cycles; cycle++){
                        computeLaplaceRHS(fieldIndex);

                        unsigned int invM_size = invM.local_size();
                        for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
                            solutionSet[fieldIndex]->local_element(dof)=solutionSet[fieldIndex]->local_element(dof)-
                            invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof)*dt_for_smoothing;
                        }

                        solutionSet[fieldIndex]->update_ghost_values();
                    }

                    op_list_index++;
                }
            }
        }

    }

    unsigned int op_list_index = 0;
    for (unsigned int var_index=0; var_index < userInputs.number_of_variables; var_index++){

        bool is_remapped_op = false;
        if (op_list_index < userInputs.variables_for_remapping.size()){
            if (var_index == userInputs.variables_for_remapping.at(op_list_index)){
                is_remapped_op = true;
                op_list_index++;
            }
        }


        if (!is_remapped_op || !userInputs.load_grain_structure){

            if (userInputs.load_ICs[var_index] == false){
                pcout << "Applying non-PField initial condition...\n";

                if (userInputs.var_type[var_index] == SCALAR){
                    VectorTools::interpolate (*dofHandlersSet[var_index], InitialCondition<dim,degree>(var_index,userInputs,this), *solutionSet[var_index]);
                }
                else {
                    VectorTools::interpolate (*dofHandlersSet[var_index], InitialConditionVector<dim,degree>(var_index,userInputs,this), *solutionSet[var_index]);
                }
            }

            else {
                // Declare the PField types and containers
                typedef PRISMS::PField<double*, double, dim> ScalarField;
                typedef PRISMS::Body<double*, dim> Body;
                Body body;

                // Create the filename of the the file to be loaded
                std::string filename;
                if (userInputs.load_parallel_file[var_index] == false){
                    filename = userInputs.load_file_name[var_index] + ".vtk";
                }
                else {
                    int proc_num = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
                    std::ostringstream conversion;
                    conversion << proc_num;
                    filename = userInputs.load_file_name[var_index] + "." + conversion.str() + ".vtk";
                }

                // Load the data from the file using a PField
                body.read_vtk(filename);
                ScalarField &conc = body.find_scalar_field(userInputs.load_field_name[var_index]);

                if (userInputs.var_type[var_index] == SCALAR){
                    pcout << "Applying PField initial condition...\n";
                    VectorTools::interpolate (*dofHandlersSet[var_index], InitialConditionPField<dim>(var_index,conc), *solutionSet[var_index]);
                }
                else {
                    std::cout << "PRISMS-PF Error: Cannot load vector fields. Loading initial conditions from file is currently limited to scalar fields" << std::endl;
                }

            }

            pcout << "Application of initial conditions for field number " << var_index << " complete \n";
        }
    }
}




// =================================================================================

// I don't think vector fields are implemented in PFields yet
//template <int dim>
//class InitialConditionPFieldVec : public Function<dim>
//{
//public:
//  unsigned int index;
//  Vector<double> values;
//  typedef PRISMS::PField<double*, double, 2> ScalarField2D;
//  ScalarField2D &inputField;
//
//  InitialConditionPFieldVec (const unsigned int _index, ScalarField2D &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}
//
//  void vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
//  {
//	  double coord[dim];
//	  for (unsigned int i = 0; i < dim; i++){
//		  coord[i] = p(i);
//	  }
//
//	  vector_IC = inputField(coord);
//  }
//};

#include "../../include/matrixFreePDE_template_instantiations.h"
