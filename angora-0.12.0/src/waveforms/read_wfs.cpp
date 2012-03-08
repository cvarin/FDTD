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

//Includes the routine that reads the time-waveform definitions

#include "headers.h"

#include "read_wfs.h"

//definition of Cwfs needed
#include "Cwfs.h"

extern string config_filename;


void read_wfs(Cwfs &Waveforms, const Config& fdtdconfig, const Config& validsettings)
{
	//Read waveform settings
	string waveform_setting_path = "Waveforms";
	if (fdtdconfig.exists(waveform_setting_path))
	{
		Setting& Waveformsettings = fdtdconfig.lookup(waveform_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(Waveformsettings,validsettings);

		//Read gaussian-waveform settings
		string gaussianwaveform_setting_path = "GaussianWaveforms";
		if (Waveformsettings.exists(gaussianwaveform_setting_path))
		{
			Setting& GaussianWaveformCollec = read_list_from_group(Waveformsettings,gaussianwaveform_setting_path);

			int num_of_gaussianwaveforms = GaussianWaveformCollec.getLength();
			for (int gaussianwaveformindex=0; gaussianwaveformindex<num_of_gaussianwaveforms; gaussianwaveformindex++)
			{
				Setting& Gaussianwaveformsettings = GaussianWaveformCollec[gaussianwaveformindex];	//go to the gaussianwaveformindex'th gaussian-waveform setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Gaussianwaveformsettings,validsettings);

				if (SettingEnabledForGrid(Gaussianwaveformsettings))		//apply only if enabled for this grid
				{
					//read tag string
					string NewWaveformTag;
					read_value_from_group<string>(Gaussianwaveformsettings,"tag",NewWaveformTag);
					double tau;
					read_value_from_group<double>(Gaussianwaveformsettings,"tau",tau);
					double delay;
					try{read_value_from_group<double>(Gaussianwaveformsettings,"delay",delay);}
					catch (AngoraSettingNotFoundException&)
					{
						delay = 0;
					}
					double amplitude;
					try{read_value_from_group<double>(Gaussianwaveformsettings,"amplitude",amplitude);}
					catch (AngoraSettingNotFoundException&)
					{
						amplitude = 1;
					}

					//Add the waveform to the waveform collector object (even in check mode, since other objects refer to these waveforms)
					Waveforms.AddGaussianWaveform(amplitude,tau,delay,NewWaveformTag);
				}
			}
		}

		//Read differentiated-Gaussian-waveform settings
		string diffgaussianwaveform_setting_path = "DifferentiatedGaussianWaveforms";
		if (Waveformsettings.exists(diffgaussianwaveform_setting_path))
		{
			Setting& DiffgaussianWaveformCollec = read_list_from_group(Waveformsettings,diffgaussianwaveform_setting_path);

			int num_of_diffgaussianwaveforms = DiffgaussianWaveformCollec.getLength();
			for (int diffgaussianwaveformindex=0; diffgaussianwaveformindex<num_of_diffgaussianwaveforms; diffgaussianwaveformindex++)
			{
				Setting& Diffgaussianwaveformsettings = DiffgaussianWaveformCollec[diffgaussianwaveformindex];	//go to the diffgaussianwaveformindex'th gaussian-waveform setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Diffgaussianwaveformsettings,validsettings);

				if (SettingEnabledForGrid(Diffgaussianwaveformsettings))		//apply only if enabled for this grid
				{
					//read tag string
					string NewWaveformTag;
					read_value_from_group<string>(Diffgaussianwaveformsettings,"tag",NewWaveformTag);
					double tau;
					read_value_from_group<double>(Diffgaussianwaveformsettings,"tau",tau);
					double delay;
					try{read_value_from_group<double>(Diffgaussianwaveformsettings,"delay",delay);}
					catch (AngoraSettingNotFoundException&)
					{
						delay = 0;
					}
					double amplitude;
					try{read_value_from_group<double>(Diffgaussianwaveformsettings,"amplitude",amplitude);}
					catch (AngoraSettingNotFoundException&)
					{
						amplitude = 1;
					}
					int ndiff;
					read_value_from_group<int>(Diffgaussianwaveformsettings,"n",ndiff);

					//Add the waveform to the waveform collector object (even in check mode, since other objects refer to these waveforms)
					Waveforms.AddDiffGaussianWaveform(amplitude,tau,delay,ndiff,NewWaveformTag);
				}
			}
		}

		//Read modulated-Gaussian-waveform settings
		string modulatedgaussianwaveform_setting_path = "ModulatedGaussianWaveforms";
		if (Waveformsettings.exists(modulatedgaussianwaveform_setting_path))
		{
			Setting& ModulatedgaussianWaveformCollec = read_list_from_group(Waveformsettings,modulatedgaussianwaveform_setting_path);

			int num_of_modulatedgaussianwaveforms = ModulatedgaussianWaveformCollec.getLength();
			for (int modulatedgaussianwaveformindex=0; modulatedgaussianwaveformindex<num_of_modulatedgaussianwaveforms; modulatedgaussianwaveformindex++)
			{
				Setting& Modulatedgaussianwaveformsettings = ModulatedgaussianWaveformCollec[modulatedgaussianwaveformindex];	//go to the modulatedgaussianwaveformindex'th gaussian-waveform setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Modulatedgaussianwaveformsettings,validsettings);

				if (SettingEnabledForGrid(Modulatedgaussianwaveformsettings))		//apply only if enabled for this grid
				{
					//read tag string
					string NewWaveformTag;
					read_value_from_group<string>(Modulatedgaussianwaveformsettings,"tag",NewWaveformTag);

					double tau;
					read_value_from_group<double>(Modulatedgaussianwaveformsettings,"tau",tau);
					double delay;
					try{read_value_from_group<double>(Modulatedgaussianwaveformsettings,"delay",delay);}
					catch (AngoraSettingNotFoundException&)
					{
						delay = 0;
					}
					double amplitude;
					try{read_value_from_group<double>(Modulatedgaussianwaveformsettings,"amplitude",amplitude);}
					catch (AngoraSettingNotFoundException&)
					{
						amplitude = 1;
					}
					double f_0;
					read_value_from_group<double>(Modulatedgaussianwaveformsettings,"f_0",f_0);
					double phase;
					try{read_value_from_group<double>(Modulatedgaussianwaveformsettings,"phase",phase);} //read in degrees!
					catch (AngoraSettingNotFoundException&)
					{
						phase = 0;
					}
					//convert from degrees to radians (since Cwf classes accept radians)
					phase *= (M_PI/180.0);

					string modulation_type;
					read_value_from_group<string>(Modulatedgaussianwaveformsettings,"modulation_type",modulation_type);
					if ((modulation_type!="sine")&&(modulation_type!="cosine"))
					{
						throw AngoraInvalidSettingValueException(Modulatedgaussianwaveformsettings, "should be \"sine\" or \"cosine\"");
					}

					//determine if the modulated pulse is also differentiated
					bool differentiated;
					try{read_value_from_group<bool>(Modulatedgaussianwaveformsettings,"differentiated",differentiated);}
					catch (AngoraSettingNotFoundException&)
					{
						differentiated=false;
					}

					//add the waveform to the collection
					if (!differentiated) //if not differentiated
					{
						if (modulation_type=="sine")
						{
							Waveforms.AddSineModulatedGaussianWaveform(amplitude,tau,f_0,phase,delay,NewWaveformTag);
						}
						else if (modulation_type=="cosine")
						{
							Waveforms.AddCosineModulatedGaussianWaveform(amplitude,tau,f_0,phase,delay,NewWaveformTag);
						}
					}
					else
					{
						if (modulation_type=="sine")
						{
							Waveforms.AddDiffSineModulatedGaussianWaveform(amplitude,tau,f_0,phase,delay,NewWaveformTag);
						}
						else if (modulation_type=="cosine")
						{
							Waveforms.AddDiffCosineModulatedGaussianWaveform(amplitude,tau,f_0,phase,delay,NewWaveformTag);
						}
					}
				}
			}
		}
	}
}
