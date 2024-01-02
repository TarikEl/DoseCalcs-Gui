/*
 * G4TOrganVoxelSegmentation.cc
 *
 *  Created on: Feb 9, 2019
 *      Author: tarik
 */


#include "G4TOrganVoxelSegmentation.hh"
#include "globals.hh"

extern  size_t* MateIDs; // index of material of each voxel unsigned int* fMateIDs; // index of material of each voxel

G4TOrganVoxelSegmentation::G4TOrganVoxelSegmentation():
number_of_z_real_values(400) ,number_of_z_simulated_values(100), number_of_x_real_values(50),
length_z_real(180), length_y_real(40), length_x_real(20),
real_thickness_z(0.1), real_thickness_y(0.1), real_thickness_x(0.1),
real_z_value(0),real_y_value(0),real_x_value(0),
number_of_x_simulated_values(0), number_of_y_real_values(0) ,number_of_y_simulated_values(0),
length_z(170),length_y(30),length_x(30),
simulated_z_value(0),simulated_x_value(0),simulated_y_value(0),
thickness_z(0.1), thickness_y(0.1), thickness_x(0.1)
{

}

G4TOrganVoxelSegmentation::~G4TOrganVoxelSegmentation(){

}

//called from
// it call
void G4TOrganVoxelSegmentation::get_phantom_voxels_data(){


}

// it call
G4int G4TOrganVoxelSegmentation::convert_from_z_simulated_to_z_real(G4int zi){

	G4int real_zi;
	length_z_real = real_thickness_z * number_of_z_real_values;
	real_zi = (zi/length_z)*length_z_real;
	// some code to convert simulated_zi to real_zi
	return real_zi;
}

//called from G4TOrganVoxelSegmentation::fill_vector_of_organ_indices()
// it call
G4int G4TOrganVoxelSegmentation::convert_from_y_simulated_to_y_real(G4int yi){
	G4int real_yi;
	// some code to convert real_zi to simulated_zi
	return real_yi;
}

//called from G4TOrganVoxelSegmentation::fill_vector_of_organ_indices()
// it call
G4int G4TOrganVoxelSegmentation::convert_from_x_simulated_to_x_real(G4int xi){
	G4int real_xi;
	// some code to convert real_zi to simulated_zi
	return real_xi;
}


//called from DICOMDetectorConsruction OrgSeg->fill_vector_of_organ_indices();
// it call convert_from_z_simulated_to_z_real(zzi);  convert_from_y_simulated_to_y_real(yyi); convert_from_x_simulated_to_x_real(xxi); get_organ_indice( realZ , realY , realX );
void G4TOrganVoxelSegmentation::fill_vector_of_organ_indices(){

	G4int aiia = 0;

	for(G4int zzi = 0; zzi< number_of_z_simulated_values; zzi++){

		for(G4int yyi = 0; yyi< number_of_y_simulated_values; yyi++){

			for(G4int xxi = 0; xxi< number_of_x_simulated_values; xxi++){

				G4int realZ = convert_from_z_simulated_to_z_real(zzi);
				G4int realY = convert_from_y_simulated_to_y_real(yyi);
				G4int realX = convert_from_x_simulated_to_x_real(xxi);

				G4int organ_indice = get_organ_indice( realZ , realY , realX );

				voxels_organs_indices[aiia] = organ_indice ;
				aiia++ ;
			}
		}
	}
}


//called from
// it call
void G4TOrganVoxelSegmentation::convert_z_x_y_vectors_index_to_vector(){

	G4int aiia = 0;
	for(G4int zzi = 0; zzi< number_of_z_simulated_values; zzi++){

		for(G4int yyi = 0; yyi< number_of_y_simulated_values; yyi++){

			for(G4int xxi = 0; xxi< number_of_x_simulated_values; xxi++){

				//voxels_organs_indices.push_back(voxels_z[zzi][yyi][xxi]);
				//aiia++;
			}
		}
	}
}



G4int organIndice = 0 ; // organ indice initialization

// it's use as argument real zi, yi and xi.
// called from G4TOrganVoxelSegmentation::fill_vector_of_organ_indices()
// it call
G4int G4TOrganVoxelSegmentation::get_organ_indice(G4int zed, G4int why,G4int exe){

	if (zed == 12) { // for heart
		// some code if we want to determine the voxels_heart in z and y and x axe using a various methods with more precision befor we decide that is heart
		if(why == 13 ){
			if(exe == 1){

				organIndice = HeartIndice;
				//voxels_organs_indices.push_back(HeartIndice);
				//voxels_z[zed][why][exe] = HeartIndice;
			}
		}
	}
	if(zed == 4){ // for brain
		// some code if we want to determine the voxels_brain in z and y and x axe using a various methods with more precision befor we decide that is brain

		if(why == 13 ){
			if(exe == 1){
				organIndice = BrainIndice;
			}
		}
	}else { // in case in empty from materiatls or just air or...
		organIndice = emptyIndice ;
	}

	return organIndice;

}

// -- first compare total body mass and volume and size to the one of ICRP 110
// -- we get the ICRP correspondant phantom size, then for each organ we delimit is
// with resulted size, and we get all the delimited voxels density, to test if it's belonged to the interval
// of organ intervall then its saved if not the voxels is removed.

// when removing or adding a voxel to an organ, we have to check the new Organ mass, Volume Limits x y z



