/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//Includes the routine that reads the total-field/scattered-field (TFSF) definitions

#include "headers.h"

#include "read_tfsf.h"

//definition of PWDataType needed
#include "Cpw.h"

//definition of FBDataType needed
#include "Cfb.h"

//definition of HGBDataType needed
#include "Chgb.h"

//definition of GSMBDataType needed
#include "Cgsmb.h"

//definition of KBDataType needed
#include "Ckohler.h"

//definition of Ctfsf needed
#include "Ctfsf.h"

//definition of Cwfs needed
#include "waveforms/Cwfs.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;

extern int rank;

extern void add_slash_to_path(string& path);
extern int create_path(const string& path);


void read_tfsf(Ctfsf &TFSF, const Config& fdtdconfig, const Config& validsettings, const Cwfs& Waveforms)
{
	//Read TFSF settings
	string tfsf_setting_path = "TFSF";
	if (fdtdconfig.exists(tfsf_setting_path))
	{
		Setting& TFSFsettings = fdtdconfig.lookup(tfsf_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(TFSFsettings,validsettings);

		//Read plane-wave settings
		string pw_setting_path = "PlaneWaves";
		if (TFSFsettings.exists(pw_setting_path))
		{
			Setting& PWCollec = read_list_from_group(TFSFsettings,pw_setting_path);

			int num_of_pws = PWCollec.getLength();
			for (int pwindex=0; pwindex<num_of_pws; pwindex++)
			{
				Setting& PWsettings = PWCollec[pwindex];	//go to the pwindex'th plane-wave setting
				//check group for invalid settings
				CheckAngoraGroupSetting(PWsettings,validsettings);

				if (SettingEnabledForGrid(PWsettings))		//apply only if enabled for this grid
				{
					PWDataType PWData;	//call the constructor for each new plane wave

					read_value_from_group<double>(PWsettings,"THETA",PWData.THETA);
					PWData.THETA *= M_PI/180;	//convert to radians

					read_value_from_group<double>(PWsettings,"PHI",PWData.PHI);
					PWData.PHI *= M_PI/180;	//convert to radians

					read_value_from_group<double>(PWsettings,"PSI",PWData.PSI);
					PWData.PSI *= M_PI/180;	//convert to radians

					read_optional_value_from_group<int>(PWsettings,"PWMarginBackX",PWData.PWMarginBackX);
					read_optional_value_from_group<int>(PWsettings,"PWMarginFrontX",PWData.PWMarginFrontX);
					read_optional_value_from_group<int>(PWsettings,"PWMarginLeftY",PWData.PWMarginLeftY);
					read_optional_value_from_group<int>(PWsettings,"PWMarginRightY",PWData.PWMarginRightY);
					read_optional_value_from_group<int>(PWsettings,"PWMarginLowerZ",PWData.PWMarginLowerZ);
					read_optional_value_from_group<int>(PWsettings,"PWMarginUpperZ",PWData.PWMarginUpperZ);

					if (read_optional_value_from_group<double>(PWsettings,"PWOriginX",PWData.PWOriginX))
					{
						//if given in the config file, it is with respect to the grid origin
						PWData.PWOriginX += OriginX;
					}

					if (read_optional_value_from_group<double>(PWsettings,"PWOriginY",PWData.PWOriginY))
					{
						//if given in the config file, it is with respect to the grid origin
						PWData.PWOriginY += OriginY;
					}

					if (read_optional_value_from_group<double>(PWsettings,"PWOriginZ",PWData.PWOriginZ))
					{
						//if given in the config file, it is with respect to the grid origin
						PWData.PWOriginZ += OriginZ;
					}

					read_optional_value_from_group<double>(PWsettings,"L_req",PWData.L_req);
					read_optional_value_from_group<bool>(PWsettings,"DisplayWarnings",PWData.DisplayWarnings);
					read_optional_value_from_group<bool>(PWsettings,"DisplayCellsPerLambda",PWData.DisplayCellsPerLambda);

					//read the extra amplitude factor
					read_optional_value_from_group<double>(PWsettings,"E0",PWData.E0);

					//read the index of the plane-wave waveform in the Waveforms object
					string waveform_tag;
					read_value_from_group<string>(PWsettings,"waveform_tag",waveform_tag);
					const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
					//If not in check mode, add the plane wave to the TFSF object
					if (!check_mode)
					{
						//the pointer to the right waveform object
						PWData.waveform = waveformptr;

						//Add the plane wave to the TFSF object
						TFSF.AddPlaneWave(PWData);
					}
				}
			}
		}

		//Read focused-beam settings
		string fp_setting_path = "FocusedBeams";
		if (TFSFsettings.exists(fp_setting_path))
		{
			Setting& FBCollec = read_list_from_group(TFSFsettings,fp_setting_path);

			int num_of_fps = FBCollec.getLength();
			for (int fpindex=0; fpindex<num_of_fps; fpindex++)
			{
				Setting& FBsettings = FBCollec[fpindex];	//go to the fpindex'th focused-beam setting
				//check group for invalid settings
				CheckAngoraGroupSetting(FBsettings,validsettings);

				if (SettingEnabledForGrid(FBsettings))		//apply only if enabled for this grid
				{
					FBDataType FBData;	//call the constructor for each new focused beam

					if(!read_optional_value_from_group<string>(FBsettings,"angular_discretization",FBData.angular_discretization))
					{
						FBData.angular_discretization = "cartesian";
					}
					else if ((FBData.angular_discretization!="cartesian")&&(FBData.angular_discretization!="radial"))
					{
						throw AngoraInvalidSettingValueException(FBsettings["angular_discretization"],"should be either \"cartesian\" or \"radial\"");
					}

					read_value_from_group<double>(FBsettings,"theta_max",FBData.theta_max);
					FBData.theta_max *= M_PI/180;	//convert to radians

					if (FBData.angular_discretization=="radial")
					{
						read_value_from_group<int>(FBsettings,"n_1",FBData.n_1);
						read_value_from_group<int>(FBsettings,"n_2",FBData.n_2);
					}
					else
					{
						//then these are optional (left at their default sentinel values)
						read_optional_value_from_group<int>(FBsettings,"n_1",FBData.n_1);
						read_optional_value_from_group<int>(FBsettings,"n_2",FBData.n_2);
					}

					read_value_from_group<double>(FBsettings,"f",FBData.f);

					read_value_from_group<double>(FBsettings,"pw_pol",FBData.pw_pol);
					FBData.pw_pol *= M_PI/180;	//convert to radians

					read_optional_value_from_group<int>(FBsettings,"FBMarginBackX",FBData.FBMarginBackX);
					read_optional_value_from_group<int>(FBsettings,"FBMarginFrontX",FBData.FBMarginFrontX);
					read_optional_value_from_group<int>(FBsettings,"FBMarginLeftY",FBData.FBMarginLeftY);
					read_optional_value_from_group<int>(FBsettings,"FBMarginRightY",FBData.FBMarginRightY);
					read_optional_value_from_group<int>(FBsettings,"FBMarginLowerZ",FBData.FBMarginLowerZ);
					read_optional_value_from_group<int>(FBsettings,"FBMarginUpperZ",FBData.FBMarginUpperZ);

					if (read_optional_value_from_group<double>(FBsettings,"FBOriginX",FBData.FBOriginX))
					{
						//if given in the config file, it is with respect to the grid origin
						FBData.FBOriginX += OriginX;
					}

					if (read_optional_value_from_group<double>(FBsettings,"FBOriginY",FBData.FBOriginY))
					{
						//if given in the config file, it is with respect to the grid origin
						FBData.FBOriginY += OriginY;
					}

					if (read_optional_value_from_group<double>(FBsettings,"FBOriginZ",FBData.FBOriginZ))
					{
						//if given in the config file, it is with respect to the grid origin
						FBData.FBOriginZ += OriginZ;
					}

					read_optional_value_from_group<double>(FBsettings,"L_req",FBData.L_req);
					read_optional_value_from_group<bool>(FBsettings,"DisplayWarnings",FBData.DisplayWarnings);
					read_optional_value_from_group<bool>(FBsettings,"DisplayCellsPerLambda",FBData.DisplayCellsPerLambda);

					//read the extra amplitude factor
					read_optional_value_from_group<double>(FBsettings,"E0",FBData.E0);

					//read the index of the focused-beam waveform in the Waveforms object
					string waveform_tag;
					read_value_from_group<string>(FBsettings,"waveform_tag",waveform_tag);
					const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
					//If not in check mode, add the focused beam to the TF/SF object
					if (!check_mode)
					{
						//the pointer to the right waveform object
						FBData.waveform = waveformptr;

						//Add a focused beam to the TF/SF object
						TFSF.AddFocusedBeam(FBData);
					}
				}
			}
		}

		//Read Hermite-Gaussian-beam settings
		string hgp_setting_path = "HermiteGaussianBeams";
		if (TFSFsettings.exists(hgp_setting_path))
		{
			Setting& HGBCollec = read_list_from_group(TFSFsettings,hgp_setting_path);

			int num_of_hgps = HGBCollec.getLength();
			for (int hgpindex=0; hgpindex<num_of_hgps; hgpindex++)
			{
				Setting& HGBsettings = HGBCollec[hgpindex];	//go to the hgpindex'th Hermite-Gaussian-beam setting
				//check group for invalid settings
				CheckAngoraGroupSetting(HGBsettings,validsettings);

				if (SettingEnabledForGrid(HGBsettings))		//apply only if enabled for this grid
				{
					HGBDataType HGBData;	//call the constructor for each new Hermite-Gaussian beam
					//Add Hermite-Gaussian beam(s) to the TFSF object
					read_value_from_group<string>(HGBsettings,"direction",HGBData.direction);

					read_value_from_group<double>(HGBsettings,"beam_halfwidth",HGBData.beam_halfwidth);

					read_value_from_group<int>(HGBsettings,"x_order",HGBData.x_order);
					read_value_from_group<int>(HGBsettings,"y_order",HGBData.y_order);

					read_value_from_group<double>(HGBsettings,"polarization",HGBData.polarization);
					HGBData.polarization *= M_PI/180;	//convert to radians

					read_optional_value_from_group<int>(HGBsettings,"HGBMarginBackX",HGBData.HGBMarginBackX);
					read_optional_value_from_group<int>(HGBsettings,"HGBMarginFrontX",HGBData.HGBMarginFrontX);
					read_optional_value_from_group<int>(HGBsettings,"HGBMarginLeftY",HGBData.HGBMarginLeftY);
					read_optional_value_from_group<int>(HGBsettings,"HGBMarginRightY",HGBData.HGBMarginRightY);
					read_optional_value_from_group<int>(HGBsettings,"HGBMarginLowerZ",HGBData.HGBMarginLowerZ);
					read_optional_value_from_group<int>(HGBsettings,"HGBMarginUpperZ",HGBData.HGBMarginUpperZ);

					if (read_optional_value_from_group<double>(HGBsettings,"HGBOriginX",HGBData.HGBOriginX))
					{
						//if given in the config file, it is with respect to the grid origin
						HGBData.HGBOriginX += OriginX;
					}

					if (read_optional_value_from_group<double>(HGBsettings,"HGBOriginY",HGBData.HGBOriginY))
					{
						//if given in the config file, it is with respect to the grid origin
						HGBData.HGBOriginY += OriginY;
					}

					if (read_optional_value_from_group<double>(HGBsettings,"HGBOriginZ",HGBData.HGBOriginZ))
					{
						//if given in the config file, it is with respect to the grid origin
						HGBData.HGBOriginZ += OriginZ;
					}

					read_optional_value_from_group<double>(HGBsettings,"L_req",HGBData.L_req);
					read_optional_value_from_group<bool>(HGBsettings,"DisplayWarnings",HGBData.DisplayWarnings);
					read_optional_value_from_group<bool>(HGBsettings,"DisplayCellsPerLambda",HGBData.DisplayCellsPerLambda);

					//read the extra amplitude factor
					read_optional_value_from_group<double>(HGBsettings,"E0",HGBData.E0);

					//read the index of the Hermite-Gaussian-beam waveform in the Waveforms object
					string waveform_tag;
					read_value_from_group<string>(HGBsettings,"waveform_tag",waveform_tag);
					const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
					//If not in check mode, add the Hermite-Gaussian beam to the TF/SF object
					if (!check_mode)
					{
						//the pointer to the right waveform object
						HGBData.waveform = waveformptr;

						//Add a Hermite-Gaussian beam to the TF/SF object
						TFSF.AddHermiteGaussianBeam(HGBData);
					}
				}
			}
		}

		//Read Gaussian-Schell-model-beam settings
		string gsmp_setting_path = "GaussianSchellModelBeams";
		if (TFSFsettings.exists(gsmp_setting_path))
		{
			Setting& GSMBCollec = read_list_from_group(TFSFsettings,gsmp_setting_path);

			int num_of_gsmps = GSMBCollec.getLength();
			for (int gsmpindex=0; gsmpindex<num_of_gsmps; gsmpindex++)
			{
				Setting& GSMBsettings = GSMBCollec[gsmpindex];	//go to the gsmpindex'th Gaussian-Schell-model-beam setting
				//check group for invalid settings
				CheckAngoraGroupSetting(GSMBsettings,validsettings);

				if (SettingEnabledForGrid(GSMBsettings))		//apply only if enabled for this grid
				{
					GSMBDataType GSMBData;	//call the constructor for each new Gaussian Schell-model beam
					//Add Gaussian Schell-model beam(s) to the TFSF object
					read_value_from_group<string>(GSMBsettings,"direction",GSMBData.direction);

					read_value_from_group<int>(GSMBsettings,"n_mode_x",GSMBData.n_mode_x);
					read_value_from_group<int>(GSMBsettings,"n_mode_y",GSMBData.n_mode_y);

					read_value_from_group<double>(GSMBsettings,"beam_halfwidth",GSMBData.beam_halfwidth);

					read_value_from_group<double>(GSMBsettings,"correlation_length",GSMBData.correlation_length);

					read_value_from_group<double>(GSMBsettings,"polarization",GSMBData.polarization);
					GSMBData.polarization *= M_PI/180;	//convert to radians

					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginBackX",GSMBData.GSMBMarginBackX);
					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginFrontX",GSMBData.GSMBMarginFrontX);
					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginLeftY",GSMBData.GSMBMarginLeftY);
					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginRightY",GSMBData.GSMBMarginRightY);
					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginLowerZ",GSMBData.GSMBMarginLowerZ);
					read_optional_value_from_group<int>(GSMBsettings,"GSMBMarginUpperZ",GSMBData.GSMBMarginUpperZ);

					if (read_optional_value_from_group<double>(GSMBsettings,"GSMBOriginX",GSMBData.GSMBOriginX))
					{
						//if given in the config file, it is with respect to the grid origin
						GSMBData.GSMBOriginX += OriginX;
					}

					if (read_optional_value_from_group<double>(GSMBsettings,"GSMBOriginY",GSMBData.GSMBOriginY))
					{
						//if given in the config file, it is with respect to the grid origin
						GSMBData.GSMBOriginY += OriginY;
					}

					if (read_optional_value_from_group<double>(GSMBsettings,"GSMBOriginZ",GSMBData.GSMBOriginZ))
					{
						//if given in the config file, it is with respect to the grid origin
						GSMBData.GSMBOriginZ += OriginZ;
					}

					read_optional_value_from_group<double>(GSMBsettings,"L_req",GSMBData.L_req);
					read_optional_value_from_group<bool>(GSMBsettings,"DisplayWarnings",GSMBData.DisplayWarnings);
					read_optional_value_from_group<bool>(GSMBsettings,"DisplayCellsPerLambda",GSMBData.DisplayCellsPerLambda);

					//read the extra amplitude factor
					read_optional_value_from_group<double>(GSMBsettings,"E0",GSMBData.E0);

					//read the index of the Gaussian-Schell-model-beam waveform in the Waveforms object
					string waveform_tag;
					read_value_from_group<string>(GSMBsettings,"waveform_tag",waveform_tag);
					const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
					//If not in check mode, add the Gaussian Schell-model beam to the TF/SF object
					if (!check_mode)
					{
						//the pointer to the right waveform object
						GSMBData.waveform = waveformptr;

						//Add a Gaussian Schell-model beam to the TF/SF object
						TFSF.AddGaussianSchellModelBeam(GSMBData);
					}
				}
			}
		}

		//Read Kohler-beam settings
		string ki_setting_path = "KohlerBeams";
		if (TFSFsettings.exists(ki_setting_path))
		{
			Setting& KBCollec = read_list_from_group(TFSFsettings,ki_setting_path);

			int num_of_kis = KBCollec.getLength();
			for (int kiindex=0; kiindex<num_of_kis; kiindex++)
			{
				Setting& KBsettings = KBCollec[kiindex];	//go to the kiindex'th Kohler-beam setting
				//check group for invalid settings
				CheckAngoraGroupSetting(KBsettings,validsettings);

				if (SettingEnabledForGrid(KBsettings))	//apply only if enabled for this grid
				{
					KBDataType KBData;	//create the Kohler-beam data structure

//					int random_seed;	//common random seed among nodes
//					read_optional_value_from_group<int>(KBsettings,"random_seed",KBData.random_seed);
//
//					read_value_from_group<double>(KBsettings,"theta_max",KBData.theta_max);
//					KBData.theta_max *= M_PI/180;	//convert to radians
//
//					read_value_from_group<int>(KBsettings,"n_theta",KBData.n_theta);
//					read_value_from_group<int>(KBsettings,"n_phi",KBData.n_phi);
//
//					read_value_from_group<double>(KBsettings,"f",KBData.f);

					if(!read_optional_value_from_group<string>(KBsettings,"simulation_type",KBData.simulation_type))
					{
						KBData.simulation_type = "deterministic";
					}
					else if ((KBData.simulation_type!="deterministic")&&(KBData.simulation_type!="stochastic"))
					{
						throw AngoraInvalidSettingValueException(KBsettings["simulation_type"],"should be either \"deterministic\" or \"stochastic\"");
					}

					if(!read_optional_value_from_group<string>(KBsettings,"angular_discretization",KBData.angular_discretization))
					{
						KBData.angular_discretization = "cartesian";
					}
					else if ((KBData.angular_discretization!="cartesian")&&(KBData.angular_discretization!="radial"))
					{
						throw AngoraInvalidSettingValueException(KBsettings["angular_discretization"],"should be either \"cartesian\" or \"radial\"");
					}

					read_value_from_group<double>(KBsettings,"theta_max",KBData.theta_max);
					KBData.theta_max *= M_PI/180;	//convert to radians

					if (KBData.angular_discretization=="radial")
					{
						read_value_from_group<int>(KBsettings,"n_1",KBData.n_1);
						read_value_from_group<int>(KBsettings,"n_2",KBData.n_2);
					}
					else
					{
						//then these are optional (left at their default sentinel values)
						read_optional_value_from_group<int>(KBsettings,"n_1",KBData.n_1);
						read_optional_value_from_group<int>(KBsettings,"n_2",KBData.n_2);
					}

//					read_value_from_group<double>(KBsettings,"f",KBData.f);

					read_value_from_group<double>(KBsettings,"pw_pol",KBData.pw_pol);
					KBData.pw_pol *= M_PI/180;	//convert to radians

					read_optional_value_from_group<int>(KBsettings,"KBMarginBackX",KBData.KBMarginBackX);
					read_optional_value_from_group<int>(KBsettings,"KBMarginFrontX",KBData.KBMarginFrontX);
					read_optional_value_from_group<int>(KBsettings,"KBMarginLeftY",KBData.KBMarginLeftY);
					read_optional_value_from_group<int>(KBsettings,"KBMarginRightY",KBData.KBMarginRightY);
					read_optional_value_from_group<int>(KBsettings,"KBMarginLowerZ",KBData.KBMarginLowerZ);
					read_optional_value_from_group<int>(KBsettings,"KBMarginUpperZ",KBData.KBMarginUpperZ);

					if (read_optional_value_from_group<double>(KBsettings,"KBOriginX",KBData.KBOriginX))
					{
						//if given in the config file, it is with respect to the grid origin
						KBData.KBOriginX += OriginX;
					}

					if (read_optional_value_from_group<double>(KBsettings,"KBOriginY",KBData.KBOriginY))
					{
						//if given in the config file, it is with respect to the grid origin
						KBData.KBOriginY += OriginY;
					}

					if (read_optional_value_from_group<double>(KBsettings,"KBOriginZ",KBData.KBOriginZ))
					{
						//if given in the config file, it is with respect to the grid origin
						KBData.KBOriginZ += OriginZ;
					}

					read_optional_value_from_group<double>(KBsettings,"L_req",KBData.L_req);
					read_optional_value_from_group<bool>(KBsettings,"DisplayWarnings",KBData.DisplayWarnings);
					read_optional_value_from_group<bool>(KBsettings,"DisplayCellsPerLambda",KBData.DisplayCellsPerLambda);

					//read the extra amplitude factor
					read_optional_value_from_group<double>(KBsettings,"E0",KBData.E0);

					//read the index of the Kohler-beam waveform in the Waveforms object
					string waveform_tag;
					read_value_from_group<string>(KBsettings,"waveform_tag",waveform_tag);
					const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];

					//If not in check mode, add the Kohler beam to the TF/SF object
					if (!check_mode)
					{
						//the pointer to the right waveform object
						KBData.waveform = waveformptr;

						//Add a Kohler beam to the TF/SF object
						TFSF.AddKohlerBeam(KBData);
					}
				}
			}
		}
	}
}
