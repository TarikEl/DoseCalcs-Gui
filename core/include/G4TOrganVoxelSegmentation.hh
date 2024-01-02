/*
 * G4TOrganVoxelSegmentation.hh
 *
 *  Created on: Feb 9, 2019
 *      Author: tarik
 */

#ifndef INCLUDE_G4TORGANVOXELSEGMENTATION_HH_
#define INCLUDE_G4TORGANVOXELSEGMENTATION_HH_


#include "globals.hh"
#include <vector>



class G4TOrganVoxelSegmentation {

public:
	// constructor and destructor
	G4TOrganVoxelSegmentation();
	~G4TOrganVoxelSegmentation();

public:

	// defined as a constants values
	G4int number_of_z_real_values; // the number of z real of a phantom have to be the maximum possible value
	G4int number_of_y_real_values; // the number of y real of a phantom have to be the maximum possible value
	G4int number_of_x_real_values; // the number of x real of a phantom have to be the maximum possible value
	G4int length_z_real;
	G4int length_y_real;
	G4int length_x_real;
	G4double real_thickness_z;
	G4double real_thickness_y;
	G4double real_thickness_x;


	// sended by detector construction class
	G4int number_of_z_simulated_values; // is small than z real
	G4int number_of_y_simulated_values; // is small than y real
	G4int number_of_x_simulated_values; // is small than x real
	G4int length_z;
	G4int length_y;
	G4int length_x;
	G4double thickness_z;
	G4double thickness_y;
	G4double thickness_x;

	void get_number_of_x_y_z_simulated(G4double nz,G4double ny, G4double nx ){
		number_of_z_simulated_values = nz;
		number_of_y_simulated_values = ny;
		number_of_x_simulated_values = nx;
	}

	void get_lenght_of_x_y_z_Dim_simulated(G4double lenZ,G4double lenY, G4double lenX ){
		length_z = lenZ;
		length_y = lenY;
		length_x = lenX;
	}

	void get_thickness_of_x_y_z_voxel_simulated(G4double thiZ,G4double thiY, G4double thiX ){
		thickness_z = thiZ;
		thickness_y = thiY;
		thickness_x = thiX;
	}


	// // sended by detector class
	std::vector<G4int> voxels_materials_indices;
	std::vector<G4int> voxels_density_values;

	// filled by function fill_vector_of_organ_indices()
	std::vector<G4int> voxels_organs_indices;

	void get_material_voxels_indices(std::vector<G4int> vector_mat_ind ){
		voxels_materials_indices = vector_mat_ind;
	}

	void get_density_voxels_values(std::vector<G4int> vector_den_val ){
		voxels_density_values = vector_den_val;
	}

	std::vector<G4int> get_organ_voxels_indices(){
		return voxels_organs_indices;
	}


	G4int real_z_value;
	G4int real_x_value;
	G4int real_y_value;
	G4int simulated_z_value;
	G4int simulated_x_value;
	G4int simulated_y_value;


	std::vector< G4int > voxels_x;
	std::vector< std::vector<G4int> > voxels_y;
	std::vector< std::vector< std::vector<G4int> > > voxels_z;

	// this method used for passing values from detectorConstruction to this object
	void get_phantom_voxels_data();

	// for the methods of segmentation, the z used is the real, then to segement the simulated phantom, we have to convert the simulated z values to the real z values for using methods of segmentation
	// the method of segmentation are based on z real
	G4int convert_from_z_simulated_to_z_real(G4int );
	G4int convert_from_y_simulated_to_y_real(G4int );
	G4int convert_from_x_simulated_to_x_real(G4int );

	void fill_vector_of_organ_indices();
	G4int get_organ_indice(G4int,G4int,G4int); // it's use as argument real zi, yi and xi.
	void convert_z_x_y_vectors_index_to_vector();

	G4int emptyIndice = 0;
	G4int HeartIndice = 1;
	G4int BrainIndice = 2;



private:


};






#endif /* INCLUDE_G4TORGANVOXELSEGMENTATION_HH_ */
