#ifndef SELECTIONSHWW_H
#define SELECTIONSHWW_H

#include "Enums.h"
#include "TypeDefs.h"

class SmurfTree;

bool hww_topveto(SmurfTree *tree);
bool hww_dy_selection(SmurfTree *tree);
bool hww_sfmet_selection(SmurfTree *tree);
bool hww_sfdymva_selection(SmurfTree *tree, const float& dymva);
bool hww_vbf_selection(SmurfTree *tree);
bool hww_pass_wwBaseline(SmurfTree *tree, const Option option);
bool hww_pass_wwSelection(SmurfTree *tree, const Option option);
bool hww_pass_wwSSSelection(SmurfTree *tree, const Option option);
bool hww_pass_wgammaSelection(SmurfTree *tree, const Option option);
bool hww_pass_wgammaSSSelection(SmurfTree *tree, const Option option);
bool hww_pass_wwPassFailSelection(SmurfTree *tree, const Option option);
bool hww_pass_wwSSPassFailSelection(SmurfTree *tree, const Option option);
void evaluate_hww_cuts(const float &analysis, float &lep1ptCut, float &lep2ptCut, 
		       float &mllCut, float &dPhiCut, float &mtLowCut, float &mtHighCut, float &mllLooseCut) ;
bool hww_pass_cutSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin );
bool hww_pass_mvaSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin );
bool hww_pass_2DSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin, const Option option );

// for assigning events to correct sources
bool hww_assign_this_event(SmurfTree *tree, DataType dataType);

#endif

