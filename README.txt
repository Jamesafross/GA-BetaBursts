on manjaro needed to install exec stack

yay -S execstack

then

execstack -c ~/.MathWorks/ServiceHost/-mw_shared_installs/*/bin/glnxa64/libmwfoundation_crash_handling.so
execstack -c ~/.MathWorks/ServiceHost/-mw_shared_installs/*/bin/glnxa64/mathworksservicehost/rcf/matlabconnector/serviceprocess/rcf/service/libmwmshrcfservice.so
cd /usr/local/MATLAB/R2024b
sudo find . -name "*.so" -exec execstack -c {} +


check in CMakeLIsts.txt: 

set(MATLAB_ROOT /usr/local/MATLAB/R2024b)

corresponds with system matlab version.  