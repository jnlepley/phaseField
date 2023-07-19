#include "../../include/FloodFiller.h"
#include "../../include/matrixFreePDE.h"

#include <numeric>
#include <iostream>
#include <vector>
#include <queue>
#include <map>

template <int dim, int degree>
void FloodFiller<dim, degree>::calcGrainSets(dealii::FESystem<dim> & fe, dealii::DoFHandler<dim> &dof_handler, vectorType* solution_field, double threshold_lower, double threshold_upper, unsigned int order_parameter_index, std::vector<GrainSet<dim>> & grain_sets){

    // unsigned int grain_index = 0;

    // I would use this, although userInputs is a protected member variable (probably for good reason)
    //unsigned int level = MatrixFreePDE::userInputs.refine_factor;

    unsigned int level = 7;

    // Loop through the whole mesh and set the user flags to false (so everything is considered unmarked)
    typename dealii::DoFHandler<dim>::cell_iterator di = dof_handler.begin(level);
    // Iterate through every 'di' cell in the mesh
    while (di != dof_handler.end(level)) {
        if (!di->has_children()){
            // Clear each cell of the 'marked' marking
            di->clear_user_flag();
        }
        ++di;
    }

    // Make a new datatype called 'grain_set', a way of representing a single grain
    GrainSet<dim> grain_set;

    // Add that set to a list of grains, for use in initialConditions.cc
    grain_sets.push_back(grain_set);

    // Set that grain's order parameter index to 0 (in inital conditions.cc)
    // Not sure why it's zero in initialConditions.cc
    grain_sets.back().setOrderParameterIndex(order_parameter_index);

    //////////////////////////
    // I feel like the below while loop is making many unnecisary calculations
    // Idealy, I'd call recusiveFloodFill once per grain, and not once per 'di'
    // cell_iterator. If there are 8 grains in the scene and ~10,000 di's per grain, 
    // well, the computation will be unnecisarily slow
    // The issue is that we have no good way of picking out specific grains
    // and begining the flood fill there
    //////////////////////////
    // The flood fill loop
    di = dof_handler.begin(level);
    unsigned int numberOfCellsIterared = 0;
    while (di != dof_handler.end(level))
    {
        // If the cell doesn't have children,
        if (!di->has_children()){

            bool grain_assigned = false;
            // Flood fill in that grain
            queueFloodFill<typename dealii::DoFHandler<dim>::cell_iterator>(di, dof_handler.end(level), solution_field, threshold_lower, threshold_upper, grain_sets, grain_assigned);

            // Flood filler will indicate if the grain is assigned or not
            if (grain_assigned){
                // Get the grain set initialized for the next grain to be found
                // grain_index++;
                // Add it to the list of grains (grain_sets)
                GrainSet<dim> new_grain_set;
                new_grain_set.setOrderParameterIndex(order_parameter_index);
                grain_sets.push_back(new_grain_set);
            }
        } 

        ++di;
        ++numberOfCellsIterared;
    }

    // std::cout << "  Total cells: " << numberOfCellsIterared << "\n";

    // Generate global list of the grains, merging grains split between multiple processors
    if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1) {

        // Send the grain set info to all processors so everyone has the full list
        createGlobalGrainSetList(grain_sets);

        // Merge grains that are split across processors
        mergeSplitGrains(grain_sets);
	}
}

template <int dim, int degree>
template <typename T>
void FloodFiller<dim, degree>::queueFloodFill(T di, T di_end, vectorType* solution_field, double threshold_lower, double threshold_upper, std::vector<GrainSet<dim>> & grain_sets, bool & grain_assigned){

    // If the current di's ele_val doesn't fall in the thresholds (not equal to target), then leave 
    if (!checkCell(di, di_end, solution_field, threshold_lower, threshold_upper)){
        return;
    }

    // Print the cell's begining location to console
    std::cout << "  Grain found at: " << di->vertex(0) << "\n";

    // Make a queue used for the flood fill
    std::queue<T> floodQueue;
    floodQueue.emplace(di);

    // Mark the current (first cell in this grain_id)
    di->set_user_flag();
    std::vector<dealii::Point<dim>> vertex_list;
    for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
        vertex_list.push_back(di->vertex(v));
    }
    grain_sets.back().addVertexList(vertex_list);
    vertex_list.clear();

    // Keep track of the number of cells iterated
    unsigned int numberOfCellsInFill = 1;

    // queue loop:
    while (floodQueue.size() > 0) {

        // Dequeue the current cell
        T currentCell = floodQueue.front();
        floodQueue.pop();

        // Imediatly find currentCell's possible neighbors (4 max if 2D, and 6 max if 3D)
        std::vector<T> possible_neighbors;
        for (unsigned int n=0; n<2*dim; n++){ 
            // Only include cells that should be marked
            if(checkCell(currentCell->neighbor(n), di_end, solution_field, threshold_lower, threshold_upper)) {
                possible_neighbors.push_back(currentCell->neighbor(n));
            }
        }

        // Marking the enqueued cells, as we enqueue more
        for(unsigned int i=0; i<possible_neighbors.size(); i++) {
            // Marking the cell
            possible_neighbors[i]->set_user_flag();
            std::vector<dealii::Point<dim>> vertex_list;
            for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
                vertex_list.push_back(possible_neighbors[i]->vertex(v));
            }
            grain_sets.back().addVertexList(vertex_list);
            // Enqueuing it 
            floodQueue.emplace(possible_neighbors[i]);
            
            ++numberOfCellsInFill;
        }

    }

    std::cout << "  Grain filled: found " << numberOfCellsInFill << " cells in grain.\n";

    grain_assigned = true;

}

// checkCell returns true IF:
// 1) The cell's 'ele_val' is in the thresholds AND
// 2) di isn't already marked AND
// 3) it isn't past the last cell iterator AND
// 4) the di is owned by this processor
template <int dim, int degree>
template <typename T>
bool FloodFiller<dim, degree>::checkCell(T di, T di_end, vectorType* solution_field, double threshold_lower, double threshold_upper){

    // If the cell in question isn't the end, isn't set, and is locally owned...
    if (di != di_end && !di->user_flag_set() && !di->has_children() && di->is_locally_owned()){

        // A vector that contains the values (of feature IDs)
        // of every quadrature point. Unsure how they are structured.
        std::vector<double> var_values(num_quad_points);

        // Get the most common value for the element
        // I don't know exacly what reinit does, but I know it is necisary for 
        // geting the values at out 'di'
        fe_values.reinit(di);

        // Return the values of a finite element function restricted to the current cell
        // that was last selected by the reinit function
        // solution_field is A vector of values that describes (globally) the finite element function 
        // that this function should evaluate at the quadrature points of the current cell. 
        // var_values[q] will contain the value of the field described by fe_function at the qth quadrature point.
        fe_values.get_function_values(*solution_field, var_values);

        // A loop that gives the most common feature ID on all quadrature points
        // Values on quadrature points is stored in var_values[q_point]
        // the double is the key type, the vals stored on quadreature points
        // the int is the value type, and tells how many times the associated double has been seen
        // Note, indexing into maps is like this:  map[key] ---> returns associated value, and is modifiable
        std::map<double, int> quadratureValues;

        // Variables for keeping track of the most common value seen at all quadrature points
        // and for keeping track *how many times* that most common value *has* been seen
        int maxNumberSeen = 0;
        double mostCommonQPointValue = -1;
        for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
            // Add the number of times that var_values[q_point] has been seen
            ++quadratureValues[var_values[q_point]];
            if (quadratureValues[var_values[q_point]] > maxNumberSeen) {
                maxNumberSeen = quadratureValues[var_values[q_point]];
                mostCommonQPointValue = var_values[q_point];
            }
        }
        //std::cout << mostCommonQPointValue << "\n";

        if (threshold_lower < mostCommonQPointValue && mostCommonQPointValue < threshold_upper){
            return true;
        }
    }
    return false;

}

// =================================================================================
// All-to-all communication of the grain sets
// =================================================================================
template <int dim, int degree>
void FloodFiller<dim, degree>::createGlobalGrainSetList (std::vector<GrainSet<dim>> & grain_sets) const
{
    int numProcs=dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	//int thisProc=dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    unsigned int num_grains_local = grain_sets.size();

    // Convert the grain_set object into a group of vectors
    std::vector<unsigned int> order_parameters;
    std::vector<unsigned int> num_elements;
    std::vector<double> vertices;

    for (unsigned int g=0; g<grain_sets.size(); g++){
        order_parameters.push_back(grain_sets.at(g).getOrderParameterIndex());

        std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[g].getVertexList();
        num_elements.push_back(vertex_list.size());

        for (unsigned int c=0; c<num_elements[g]; c++){
            for (unsigned int v=0; v<dealii::Utilities::fixed_power<dim>(2.0); v++){
                for (unsigned int d=0; d<dim; d++){
                    vertices.push_back(vertex_list[c][v][d]);
                }
            }
        }
    }

    unsigned int num_vertices = 0;
    for (unsigned int g=0; g<grain_sets.size(); g++){
        num_vertices += num_elements[g] * dealii::Utilities::fixed_power<dim>(2) * dim;
    }

    // Communicate how many grains each core has
    std::vector<int> num_grains_per_core(numProcs,0);

    MPI_Allgather(&num_grains_local, 1, MPI_INT, &num_grains_per_core[0], 1, MPI_INT, MPI_COMM_WORLD);

    int num_grains_global = std::accumulate(num_grains_per_core.begin(), num_grains_per_core.end(), 0);

    // Communicate the order_parameters
    std::vector<int> offset(numProcs,0);
    for (int n=1; n<numProcs; n++){
        offset[n] = offset[n-1] + num_grains_per_core[n-1];
    }

    std::vector<unsigned int> order_parameters_global(num_grains_global,0);

    MPI_Allgatherv(&order_parameters[0], num_grains_local, MPI_UNSIGNED, &order_parameters_global[0], &num_grains_per_core[0], &offset[0], MPI_UNSIGNED, MPI_COMM_WORLD);

    // Communicate the number of elements
    std::vector<unsigned int> num_elements_global(num_grains_global,0);

    MPI_Allgatherv(&num_elements[0], num_grains_local, MPI_UNSIGNED, &num_elements_global[0], &num_grains_per_core[0], &offset[0], MPI_UNSIGNED, MPI_COMM_WORLD);

    // Communicate the vertices
    unsigned int total_elements = std::accumulate(num_elements_global.begin(), num_elements_global.end(), 0);
    int num_vertices_global = (unsigned int)total_elements * dealii::Utilities::fixed_power<dim>(2) * (unsigned int)dim;
    std::vector<double> vertices_global(num_vertices_global,0);

    std::vector<int> num_vertices_per_core;

    unsigned int g = 0;
    for (int i=0; i<numProcs; i++){
        int num_vert_single_core = 0;
        for (int j=0; j<num_grains_per_core.at(i); j++){
            num_vert_single_core += num_elements_global.at(g) * dealii::Utilities::fixed_power<dim>(2) * (unsigned int)dim;
            g++;
        }
        num_vertices_per_core.push_back(num_vert_single_core);
    }

    offset.at(0) = 0;
    for (int n=1; n<numProcs; n++){
        offset[n] = offset[n-1] + num_vertices_per_core[n-1];
    }

    MPI_Allgatherv(&vertices[0], num_vertices, MPI_DOUBLE, &vertices_global[0], &num_vertices_per_core[0], &offset[0], MPI_DOUBLE, MPI_COMM_WORLD);

    // Put the GrainSet objects back together
    grain_sets.clear();

    for (int g=0; g<num_grains_global; g++){
        GrainSet<dim> new_grain_set;
        for (unsigned int c=0; c<num_elements_global.at(g); c++){
            std::vector<dealii::Point<dim>> verts;
            for (unsigned int v=0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){
                double coords[dim];
                for (unsigned int d=0; d<dim; d++){
                    coords[d] = vertices_global.front();
                    vertices_global.erase(vertices_global.begin());
                }
                dealii::Tensor<1,dim> tensor_coords(coords);
                dealii::Point<dim> vert(tensor_coords);
                verts.push_back(vert);
            }
            new_grain_set.addVertexList(verts);
        }
        new_grain_set.setOrderParameterIndex(order_parameters_global.at(g));
        grain_sets.push_back(new_grain_set);
    }


}

// =================================================================================
// Check to see if any grains on different processors share vertices
// =================================================================================

template <int dim, int degree>
void FloodFiller<dim, degree>::mergeSplitGrains (std::vector<GrainSet<dim>> & grain_sets) const
{

    // Loop though each vertex in the base grain "g"
    for (unsigned int g=0; g<grain_sets.size(); g++){

        std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[g].getVertexList();

        // Now cycle through the other grains to find overlapping elements
        for (unsigned int g_other=g+1; g_other<grain_sets.size(); g_other++){
            bool matching_vert = false;

            std::vector<std::vector<dealii::Point<dim>>> vertex_list_other = grain_sets[g_other].getVertexList();

            for (unsigned int c = 0; c < vertex_list.size(); c++){
                for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){

                    for (unsigned int c_other = 0; c_other < vertex_list_other.size(); c_other++){
                        for (unsigned int v_other = 0; v_other < dealii::Utilities::fixed_power<dim>(2.0); v_other++){
                            // Check if the vertices match
                            if (vertex_list[c][v] == vertex_list_other[c_other][v_other]){
                                matching_vert = true;
                                break;
                            }
                            if (matching_vert){
                                break;
                            }
                        }
                        if (matching_vert){
                            break;
                        }
                    }
                    if (matching_vert){
                        break;
                    }
                }
                if (matching_vert){
                    break;
                }
            }

            if (matching_vert){

                for (unsigned int c_base = 0; c_base < vertex_list.size(); c_base++){
                    grain_sets[g_other].addVertexList(vertex_list.at(c_base));
                }
                grain_sets.erase(grain_sets.begin()+g);
                g--;
                break;
            }
        }
    }
}

// Template instantiations
template class FloodFiller<2,1>;
template class FloodFiller<2,2>;
template class FloodFiller<2,3>;
template class FloodFiller<2,4>;
template class FloodFiller<3,1>;
template class FloodFiller<3,2>;
template class FloodFiller<3,3>;
template class FloodFiller<3,4>;
