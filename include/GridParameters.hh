#ifndef grid_h
#define grid_h
#include <G4ThreeVector.hh>
#include "G4SystemOfUnits.hh"
namespace Grid 
{
	
	//grid_data_generator
	G4double CAD_GRID_RELATIVE_OFFSETS[32][3]={ //measured from cad file, relative to the first detector
		{0., 0., 0.,},
		{  0.,   0., -23.},
		{ -0.,   0., -50.},
		{  0.,   0., -73.},
		{  0., -25.,  23.},
		{  0., -25.,   0.},
		{  0., -25., -23.},
		{  0., -25., -50.},
		{  0., -25., -73.},
		{ -0., -25., -96.},
		{ -0., -50.,  37.},
		{  0., -50.,  14.},
		{  0., -50.,  -9.},
		{  0., -50., -64.},
		{  0., -50., -87.},
		{   0.,  -50., -110.},
		{  0., -75.,  37.},
		{  0., -75.,  14.},
		{  0., -75.,  -9.},
		{  0., -75., -64.},
		{  0., -75., -87.},
		{   0.,  -75., -110.},
		{   0., -100.,   23.},
		{   0., -100.,    0.},
		{   0., -100.,  -23.},
		{  -0., -100.,  -50.},
		{  -0., -100.,  -73.},
		{   0., -100.,  -96.},
		{   0., -125.,    0.},
		{   0., -125.,  -23.},
		{   0., -125.,  -50.},
		{   0., -125.,  -73.},	
	};
	//above are nominal grid parameters

	G4double CAD_GEO_CENTER[3]={21.0624, 165.599, 163.2};
	//measured from cad
	//calculation see the notebook in /home/xiaohl/FHNW/STIX/CAD/steps/material_calculation.inpy


	G4ThreeVector getDetectorCenterCoordsCAD(int i){
		G4ThreeVector v(CAD_GEO_CENTER[0]*mm+CAD_GRID_RELATIVE_OFFSETS[i][0]*mm, 
				CAD_GEO_CENTER[1]*mm+CAD_GRID_RELATIVE_OFFSETS[i][1]*mm, 
				CAD_GEO_CENTER[2]*mm+CAD_GRID_RELATIVE_OFFSETS[i][2]*mm
				);
		return v;
	};
	G4ThreeVector getGridCenterCAD(int i, int layer=0){
		//layer=0 rear, 1 for front
		G4double offset=-1;
		G4double x0[2]={-34.2226+offset , -579.523+ offset}; //measured from cad drawing

		G4double grid0Center[3]={0, 165.599, 164.0};
		G4ThreeVector v(x0[layer]*mm+CAD_GRID_RELATIVE_OFFSETS[i][0]*mm, 
				grid0Center[1]*mm+CAD_GRID_RELATIVE_OFFSETS[i][1]*mm, 
				grid0Center[2]*mm+CAD_GRID_RELATIVE_OFFSETS[i][2]*mm
				);
		return v;
	}


};

#endif
