% Copyright (c) 2011,2012 Ilker R Capoglu 
% Returns true if the dataset dataset_name exists in the HDF5 file filename.

function dataset_exists = hdf5_exists(filename,dataset_name)

% filename: The full path of the HDF5 file (string)
% dataset_name: The full path of the array dataset in the HDF5 file (string)

dataset_exists = true;

fid=H5F.open(filename);

try
    arrayid=H5D.open(fid,['/' dataset_name]);
catch err
   errmsg = err.message;
   if(strfind(errmsg, 'not found'))
      dataset_exists = false;
   end
end

H5F.close(fid);