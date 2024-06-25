function [Input]=Generate_Input

addpath(genpath('/home/ran/C-projects/YAML-MATLAB'))
yaml_file='par.yaml';
y_str = ReadYaml(yaml_file);

time_sim  = y_str.time_sim;  dt=y_str.dt;
steps_tot = round(time_sim / dt);

Tstart   = y_str.Tstart;
Tend     = y_str.Tend;
AmpInp   = y_str.AmpInp;

Input    = zeros(steps_tot,1);
Input(round(Tstart/dt):round(Tend/dt))=AmpInp;




dlmwrite('Input.txt',Input,'delimiter','\t');
