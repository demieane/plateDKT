%==========================================================================
%   READ ME
%
%   Last review: 20/10/2023
%
%
%   What is the geometry of the wing? 
%           -> Create DeLaunay triangulation mesh (PDE Toolkit)
%   Boundary conditions
%   Material properties
%==========================================================================

% PHD VERIFICATION SECTION
% (UPDATED) mainDKT_eig.m 

% (UPDATED) mainDKT_static_dynamic.m (Static load Navier solution at the end) 

% TROUBLESHOOT

%The libstdc++ library shipped with Ubuntu 18.04 may be older and may not contain the required GLIBCXX version that MATLAB requires.
%To resolve this issue you may upgrade the libstdc++ libraries by running the following command as they appear:
%sudo add-apt-repository ppa:ubuntu-toolchain-r/test
%sudo apt-get update
%sudo apt-get install gcc-4.9
%sudo apt-get upgrade libstdc++6
%OR 
%Use MATLAB's libraries by setting the environment variable:
%export LD_LIBRARY_PATH=/usr/local/MATLAB/R2021a/sys/os/glnxa64
%Then attempt to open MATLAB again. 

%The errors are caused by a mismatch between the libstdc++ library shipped with Ubuntu 22.04 and MATLAB R2022b.
%As a workaround, you can direct MATLAB to use the system's libstdc++ library rather than the MATLAB shipped version. This can be achieved by setting the following environment variable prior to starting MATLAB:
%export LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
%You can then force load the system libstdc++ with the expected GLIBCXX version.

%export LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6
%alias... matlab
%In the first case the environment is set only for the MATLAB process, 
%in the second case it is set for the current shell session.