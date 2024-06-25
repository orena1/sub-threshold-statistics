function Create_Connectivity_Matrix_From_Data(PyrCell,FSCell)
% Create_Connectivity_Matrix_From_Data(’L5MouseR1’, ‘BasektR1’)


yaml_file='par.yaml';
y_str = ReadYaml(yaml_file);

ne=y_str.ne;ni=y_str.ni;ke=y_str.ke;ki=y_str.ki;
gee=y_str.gee;gei=y_str.gei;gie=y_str.gie;gii=y_str.gii;tee=y_str.tee;tei=y_str.tei;tie=y_str.tie;tii=y_str.tii;

g_S=y_str.g_S;
Cvk=y_str.Cvk;


JEE=gee/sqrt(ke)/tee;JEI=gei/sqrt(ki)/tei;JIE=gie/sqrt(ke)/tie;JII=gii/sqrt(ki)/tii;


%%%%% Create Balance matrix J%%%%%%%%%%%%%%
[J_EE]=make_sparse_normalInDeg(ne,ne,ke,JEE,Cvk);

[J_IE]=make_sparse_normalInDeg(ni,ne,ke,JIE,Cvk);

[J_EI]=make_sparse_normalInDeg(ne,ni,ki,JEI,Cvk);

[J_II]=make_sparse_normalInDeg(ni,ni,ki,JII,Cvk);

J_B=[[J_EE,J_EI];[J_IE, J_II]];
  
%%%%% Create Specific matrix J%%%%%%%%%%%%%%

[J_EE_S]=zeros(ne,ne);
[J_IE_S]=zeros(ni,ne);

[J_EI_S_LR]=make_sparse_mat(ne/2,ni/2,sqrt(ki),g_S*JEI);
[J_EI_S_RL]=make_sparse_mat(ne/2,ni/2,sqrt(ki),g_S*JEI);
[J_II_S_LR]=make_sparse_mat(ni/2,ni/2,sqrt(ki),g_S*JII);
[J_II_S_RL]=make_sparse_mat(ni/2,ni/2,sqrt(ki),g_S*JII);


J_EI_S=[[zeros(ne/2,ni/2),J_EI_S_LR];[J_EI_S_RL, zeros(ne/2,ni/2)]];
J_II_S=[[zeros(ni/2,ni/2),J_II_S_LR];[J_II_S_RL, zeros(ni/2,ni/2)]];

J_S=[[J_EE_S,J_EI_S];[J_IE_S, J_II_S]];

J=J_B+J_S;
J=J-diag(diag(J));
zIndeg=J(1,:)>0;
zIndeg(1:10)

zOutdeg=J(:,1)>0;
zOutdeg(1:10)

%%%%%%%%%%%%%%%%% Save Out Deg %%%%%%%%%%%%%%%%%%
[ind_i,ind_j,val] = find(J);
ind_i=ind_i-1; %C++ format 0 to n-1
ind_j=ind_j-1; %C++ format 0 to n-1
J_transpose=[ind_j,ind_i,val];
dlmwrite('connectivityOutDeg.txt',J_transpose,'delimiter','\t');

%%%%%%%%%%%%%%%%% Save In Deg  %%%%%%%%%%%%%%%%%%
% [s_ind_i,sort_ind]=sort(ind_i);
% s_ind_j=ind_j(sort_ind); %C++ format 0 to n-1
% s_val=val(sort_ind); %C++ format 0 to n-1
[ind_i,ind_j,val] = find(J');
ind_i=ind_i-1; %C++ format 0 to n-1
ind_j=ind_j-1; %C++ format 0 to n-1

% Generate rho,alpha for E neurons
Total_E_neruons_in_list=find(ind_j==ne,1)-1;
% Layer 5 Rat
if strcmp(PyrCell,'L5RatR1')
z= readtable('Rm4000_Rin23_norm.csv');
end


% Layer 4 Spiny
if strcmp(PyrCell,'SpinyR1')
z= readtable('Rm1000_Rin45_norm.csv');
end


if strcmp(FSCell,'BasketR1')
z= readtable('/home/ran/Dropbox (HHMI)/VoltageWholeCellProject-2018/MultiCompProject/rhos/NMO_130658/Rm3000_Rin36_norm.csv');
end
% Generate rho,alpha for I neurons

[rhoE,alphaE]=ChooseAlphaRhoFromDataEI(z,Total_E_neruons_in_list);
disp(['Exc alpha: ' num2str(mean(alphaE)) ' rho: ' num2str(mean(rhoE))] )

Total_I_neruons_in_list=length(val)-Total_E_neruons_in_list;
[rhoI,alphaI]=ChooseAlphaRhoFromDataEI(z,Total_I_neruons_in_list);
disp(['Inh alpha: ' num2str(mean(alphaI)) ' rho: ' num2str(mean(rhoI))] )

rho=[rhoE;rhoI];
alpha=[alphaE;alphaI];
figure(50)
subplot(2,1,1)
hist(rhoE,100)
subplot(2,1,2)
hist(rhoI,100)

% rho=ones(size(val));
% alpha=ones(size(val));
% J=[s_ind_i,s_ind_j,s_val,rho,alpha];
J_direct=[ind_j,ind_i,val,rho,alpha];

dlmwrite('connectivityInDeg.txt',J_direct,'delimiter','\t');


    function [J]=make_sparse_mat(N1,N2,K,J_0)
%         plotpplo
        % Create Balance matrix J      
        J=rand(N1,N2);
        J(J>=K/N2)=0;
        J=ceil(J);
        J=J_0.*J;
        J        = sparse(J);
            
    end



    function [J]=make_sparse_normalInDeg(N1,N2,K,J_0,Cvk)
        J=zeros(N1,N2);
        if Cvk==100
            Cvk=sqrt(K)/K
        end
        tmpvec=(1+Cvk*randn(N1,1));
        kindeg=round((    tmpvec.* (tmpvec>0   ) ).*K); % truncate by 0
        kindeg(kindeg>N2)=N2; %Truncate by maximal indegree
        figure;hist(kindeg,100)
        % Create Balance matrix J     
        for nn1=1:N1
            J(nn1,randperm(N2,kindeg(nn1)))=J_0;
        end
        J=sparse(J);
            
    end


end