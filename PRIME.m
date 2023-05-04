
function
%--------------------------------------------------------------------------
%   A program for a process-based 3 Parameter Continuous Mixing Model.
%   The model requires measured sap flux density, soil water and tree sap
%   isotopic ompositions, and soil and leaf water potentials.
%   Written by Bing Cheng Si and Eric J Neil, 
%   University of Saskatchewan, October 7, 2022
%--------------------------------------------------------------------------

clear all;
global BESTX BESTF ICALL PX PF;
global threshold;

%You can select from the following models.
AllmodelTypeArray=...
    ["Kumaraswamy" , "TrCauchy" , "Beta" , "TrGamma" , "TrWeibull",...
     "TrLogistic" , "Tr2Normal" , "Tr2SNormal" , "Lognormal" , "Simplex" ,...
     "Lindley", "MixLindley" , "TrExponent" , "3MixUniform" , "UnitLindley" , ...
     "MixULindley" , "3Tr2SNormal" , "3Bimodel" , "3Lognormal" , "3UnitWeibull" ,...
     "LogitNorm" , "FoldedNorm" , "Bimodel" , "Wigner" , "LogLogistic",...
     "UnitWeibull" , "UnitGamma1" , "UnitGamma2" , "UnitFoldN", "UnitLogistic" ,...
     "UnitLindley" , "Pert1" , "Pert2" , "UnitInvGaus" , "ToppLeone" ,... 
     "JohnsonSb" , "UnitGompertz"];       
       
%Select which model to use. Select 1 for running all models
modelTypeArray = AllmodelTypeArray(1,:);

%You can select one from the following algorithm for optimization.
AllAlgorithmArray=...
   ["lsqnonlin" , "fminunc" , "fmincon" , "patternsearch" ,"fminsearch",...
    "multis_lsq" , "multis_unc" , "multis_con" , "global" , "simulanneal",...
    "particleswarm", "ga", "purecmaes" , "DE"  , "SCE_UA" ,...       
    "aeDE" , "SPUCI", "surrogateopt" , "paretosearch" , "bayesopt"];

Algorithm = AllAlgorithmArray(1,:);            % Section of Algorithms 

%--------------------------------------------------------------------------
%             Read all the data for the analysis
%--------------------------------------------------------------------------
fileLoc='C:\Users\...;
fileName='PRIME input data - July 31.xlsx';
Delta=readmatrix([fileLoc,fileName],'Sheet','Isotopes','Range','B4:F11');
SoilPsi=readmatrix([fileLoc,fileName],'Sheet','Soil Water Potential','Range','B4:C11');
SapFlux=readmatrix([fileLoc,fileName],'Sheet','Tree Sap Flux','Range','B4:C15');
numVars=3; varNames={'Time','Psi'}; varTypes={'double','double'}; 
dataRange='B4:C11'; dataSheet='Leaf_Potential';
opts = spreadsheetImportOptions('NumVariables',numVars,'VariableNames',varNames, ...
    'VariableTypes',varTypes,'Sheet',dataSheet,'DataRange',dataRange);
LeafPsi=readmatrix([fileLoc,fileName],opts);

%--------------------------------------------------------------------------
LeafPsi(:,1)=LeafPsi(:,1)*24;   % convert dunit day to hours.
SapFlux(:,1)=SapFlux(:,1)*24;   % convert dunit day to hours.

mm=1500;                        %Number of points for summations or integration
Zr=250;                         %Zr=max(Delta(:,1)); % Maximum root depth
%weighted heavier to 18O and Sap flux density realtive 2H in optimization
Scaling=[1 20 1];
leafStime=min(LeafPsi(:,1));    %Time started Leaf potential measurement
leafEndtime=max(LeafPsi(:,1));  %Time ended Leaf potential measurement
psidepth=zeros(mm,1);           %Vector for storing mm depth intervals
leafPsitime=zeros(mm,1);        %Vector for storing time of leaf potentials
for i=1:mm
    psidepth(i)=(i-1)*Zr/(mm-1);
    leafPsitime(i)=leafStime+(i-1)*(leafEndtime-leafStime)/(mm-1);
end

%--------------------------------------------------------------------------
% The following are for interpolating leaf water potential, soil potential,
% and soil isotopes values for integration
%--------------------------------------------------------------------------
Interpolation="matlab";
InterpMethod=["nearestinterp" "pchip" "smoothingspline" "cubicinterp"];
fitType=InterpMethod(1);
switch Interpolation
    case 'matlab'      %make sure  leaf potential and sap flux
        % measurements have similar the start and end time
        % also make sure soil isotope and soil
        % potenital measurements start at soil surface
        % and end at the Zr depth
        C2H_M=fit(Delta(:,1),(Delta(:,2)+1000)/6.420135,fitType);
        C18O_M=fit(Delta(:,1),(Delta(:,3)+1000)*2.0052,fitType);
        soilPsi=fit(SoilPsi(:,1),SoilPsi(:,2),fitType);
        leafPsi= fit(LeafPsi(:,1),LeafPsi(:,2),fitType);

    case 'segMent2Continu'
        C2H_M=@(x)segMent2Continu(x,Delta(:,1),(Delta(:,2)+1000)/6.420135);
        C18O_M=@(x)segMent2Continu(x,Delta(:,1),(Delta(:,3)+1000)*2.0052);
        soilPsi=@(x)segMent2Continu(x,SoilPsi(:,1),SoilPsi(:,2));
        leafPsi= @(t)segMent2Continu(t,LeafPsi(:,1),LeafPsi(:,2));
end

%--------------------------------------------------------------------------
%    Store interpolated values in a vector for repeated use
%           in SSSE_rev3--the objective function
%--------------------------------------------------------------------------
mData.soilPsiv = arrayfun(@(z) (soilPsi(z)-(5.5+z)*0.00981),psidepth);  %Total soil potential
mData.leafPsiv=arrayfun(@(t) leafPsi(t),leafPsitime);
mData.leafPsiAtSapTv = arrayfun(@(t) leafPsi(t),SapFlux(:,1));
mData.C2H_Mv=arrayfun(@(z) C2H_M(z),psidepth);
mData.C18O_Mv=arrayfun(@(z) C18O_M(z),psidepth);
mData.psidepth=psidepth;
mData.SapFlux=SapFlux;
mData.Delta=Delta;

%--------------------------------------------------------------------------
%           Options settings to control
%        different optimization algorithms
%--------------------------------------------------------------------------
opts_lsq =  optimoptions('lsqnonlin','UseParallel',false,'Display','off'...
    ,'Algorithm', 'levenberg-marquardt','MaxFunctionEvaluations',1000, ...
    'MaxFunctionEvaluations',1000,'OutputFcn', @myoutput);
opts_fmin = optimoptions('fmincon','UseParallel',false,'Display','off', ...
    'Algorithm','interior-point','FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',1000,'HessianFcn', ...
    [],'FunctionTolerance',0.0001,'OutputFcn', @myoutput);
opts_sqp = optimoptions('fmincon','UseParallel',false,'Display', ...
    'off','Algorithm','sqp','FiniteDifferenceType','central','OutputFcn', ...
    @myoutput);
opts_unc = optimoptions('fminunc','UseParallel',false,'Display', ...
    'off','Algorithm','Quasi-Newton','OutputFcn', @myoutput);

%--------------------------------------------------------------------------
%                     Start the optimization
%--------------------------------------------------------------------------
NofModels=length(modelTypeArray);
NofAlgorithms=length(Algorithm);
ThetaSSUnc=cell(NofModels,NofAlgorithms);
ThetaSSCon=cell(NofModels,NofAlgorithms);
ThetaErrSSUnc=cell(NofModels,NofAlgorithms);
ThetaErrSSCon=cell(NofModels,NofAlgorithms);
C2H18O_unc=cell(NofModels,NofAlgorithms);
C2H18O_con=cell(NofModels,NofAlgorithms);
Qttunc_con=cell(NofModels,NofAlgorithms);
rootwater_con=cell(NofModels,NofAlgorithms);
Theta_All=cell(NofAlgorithms,1);
%stats=cell(NofModels,1);
depth=zeros(251,1);
for i=1:251
    depth(i)=(i-1)*Zr/250;
end

%--------------------------------------------------------------------------
%         Iterate through the selected Algorithms and models
%--------------------------------------------------------------------------
for j=1:NofAlgorithms
    for k=1:NofModels
         clear Theta0 ub lb A b;
%--------------------------------------------------------------------------
%                      Assign initial values
%                         for optimization
%--------------------------------------------------------------------------
        ModelType = modelTypeArray(k);
        switch ModelType
            case { 'Lindley','UnitLindley', 'TrExponent','Pert1'}
                Theta0 =   [  0,  0.001;
                             10,    0.2;
                            0.4,   0.02;
                              0,    0.0;           %Aineq  =%bineq
                              1,    1.0 ];
            case '3MixUniform'
                Theta0 =  [   0,     0,      0 ;     % Lower Bound< theta < upperB
                              1,     1,     0.2;     %upper Bound
                            0.8,   0.0,   0.015;     % theta0=   % Initial guess values
                              1,     1,      0 ;     %Aineq  < bineq
                              1,    40,     90 ];
            case {'TrWeibull','Kumaraswamy','Beta','LogitNorm', ...
                    'UnitWeibull' ,'UnitInvGaus','ToppLeone'}
                Theta0 =   [  0,  0.01,     0.0;
                             Zr,  Zr/2,     0.1;
                              2,     7,    0.02;
                              0,     0,     0  ;          %Aineq  < bineq
                              1,     1,     1.0];
            case {'Tr2Normal', 'TrGamma'   ,'TrCauchy'  , 'TrLogistic' ,...
                    'Pert2'   , 'FoldNorm'  ,'Tr2SNormal', 'Bimodel'    ,...
                 'FoldedNorm', 'Lognormal' ,'Wigner'    , 'LogLogistic',...
                 'UnitFoldN' , 'UnitGamma1', 'UnitGamma2','UnitLogistic',...
                 'UnitGompertz'}
                Theta0 =  [0.01,   0.1,   0.001;
                             Zr,  Zr/2,     0.2;
                             10,     2,    0.04;
                              0,     0,     0.0;           %Aineq  =%bineq
                              1,     1,      1.0 ];
            case {'MixLindley','MixULindley','3Tr2SNormal', '3Bimodel',...
                    '3Lognormal','3UnitWeibull'}
                Theta0 = [0.01, 0.01,  0.0, 0.001;
                            Zr, Zr/2,   10,   0.2;
                            2,     4,  0.1,  0.02;
                            0,     0,  0.0,   0.0;           %Aineq  =%bineq
                            1,     1,  1.0,   1.0 ];
        
            case { 'JohnsonSb'}
                Theta0 = [-500,  0.01,     0.001;
                            Zr,  Zr/2,       0.2;
                             7,     4,      0.02;
                             0,     0,       0.0;           %Aineq  =%bineq
                             1,     1,       1.0 ];
            case {'Simplex'}
                Theta0 =   [ 0,   0.0,     0.001;
                             1,   800,       0.2;
                           0.1,   0.8,      0.02;
                             0,     0,       0.0;           %Aineq  =%bineq
                             1,     1,       1.0 ];
            
        end
        %------------------------------------------------------------------
        lb=Theta0(1,:);ub=Theta0(2,:);xx0=Theta0(3,:);
        
        if strcmp(ModelType,'3MixUniform')
            A=Theta0(4,:); b=Theta0(5,1); Seg_D=Theta0(5,2:3);
        else
            A=[]; b=[]; Seg_D=[];
        end
        fprintf('The %s PDF model with  %s algorithm\n', ...
            ModelType,Algorithm(j));
        if(~isequal(Algorithm(j),'SCE_UA') ...
                && ~isequal(Algorithm(j),'SPUCI'))
            for jjjj=1:length(xx0)
                fprintf('      theta(%s)',num2str(jjjj));
            end
            fprintf('      SSE\n');
        end
%------------------------------------------------------------------
       
%--------------------------------------------------------------------------
%      Calculate the difference between measured
%             and predicted for lsqnonlin
%--------------------------------------------------------------------------
Diff=@(theta)Cook_SSE(theta,mData,ModelType,Scaling,Seg_D);

%--------------------------------------------------------------------------
%    SSE for fmincon,fminunc,generic algorithm and simulated annealing

SSE = @(theta) sum(Diff(theta).^2);

%--------------------------------------------------------------------------
nvar=length(xx0); 

x=optimvar("x",1,nvar,'LowerBound',lb,'UpperBound',ub); x0.x=xx0;
    funSSE=fcn2optimexpr(@(x) SSE(x),x,...
            'ReuseEvaluation',true,'OutputSize',[1 1] );
 probSSE=optimproblem('Objective',funSSE);
      funDiff=fcn2optimexpr(@(x) Diff(x),x,Analysis=="on");
% % % %       funDiff=fcn2optimexpr(@(x) Diff(x),x);
 probDiff=optimproblem('Objective',sum(funDiff.^2));
%-------------------------------------------------------------------------- 
 y=optimvar("y",1,nvar);  y0.y=xx0; %unconstrained problem
        funUnc=fcn2optimexpr(@(y) SSE(y),y,...
            'ReuseEvaluation',true,'OutputSize',[1 1]);
 probSSEUnc=optimproblem('Objective',funUnc);
 
%--------------------------------------------------------------------------
%    rng default % For reproducibility
%    ms = ...
%     MultiStart('UseParallel',false,OutputFcn=@multimyoutput,Display='off');
%--------------------------------------------------------------------------       
        switch Algorithm(j)
            case 'lsqnonlin'
                [sollsq,fvallsq,eflaglsq,outputlsq] = ...
                  solve(probDiff,x0,'Solver','lsqnonlin','Options',opts_lsq);
                Theta=sollsq.x;
            case 'fmincon'
                [solcon,fvalcon,eflagcon,outputcon] = ...
                    solve(probSSE,x0,'Solver','fmincon','Options',opts_fmin);
                Theta=solcon.x;
            case 'fminunc'
                [solunc,fvalunc,eflagunc,outputunc] =...
                    solve(probSSEUnc,y0,"Solver","fminunc","Options",opts_unc);
                Theta=solunc.y;
            case 'fminsearch'
                Theta = fminsearch(@(Theta) SSE(Theta),xx0, ...
                    optimset('Display','off','OutputFcn', @myoutput));
            case 'patternsearch'  % can use latin hypercube(searchlht) or searchga,
                opts = optimoptions('patternsearch','UseParallel',false,'Display', ...
                    'off','MaxFunctionEvaluations',300,'MaxIterations',100, ...
                    'SearchFcn',{@GPSPositiveBasis2N},'OutputFcn',@pamyoutput);
                [solpat,fvalpat,eflagpat,outputpat] = solve(probSSE,x0,...
                    "Solver","patternsearch","Options",opts); 
                Theta=solpat.x;
            case 'paretosearch'   % outputfcn is the same as patternsearch
                 opts = optimoptions(@paretosearch,'UseParallel',false,...
                     'Display','off','OutputFcn',@pamyoutput);
                [solpto,fvalpto,eflagpto,outputpto] = ...
                    solve(probSSE,x0,"Solver","paretosearch","Options",opts);
                Theta=solpto.x;
            case 'surrogateopt'
                opts = optimoptions("surrogateopt",'UseParallel',false,...
                    "MaxFunctionEvaluations",1000,'OutputFcn',@myoutput);
                [solsur,fvalsur,eflagsur,outputsur] = ...
                    solve(probSSE,x0,"Solver","surrogateopt");%,"Options",opts);
                Theta=solsur.x;
            case 'particleswarm'
                opts = optimoptions("particleswarm",'UseParallel',...
                    false,'MaxIterations',1000,'OutputFcn', @psmyoutput);
                [solpar,fvalpar,eflagpar,outputpar] =...
                    solve(probSSE,x0,"Solver","particleswarm","Options",opts)
                Theta=solpar.x;
            
            case 'global'      % Requires the global optimization tools box
                gs = GlobalSearch('OutputFcn', @multimyoutput,'Display',['' ...
                    'off'],'MaxTime',100);
                [solglo,fvalglo,eflagglo,outputglo] =...
                    solve(probSSE,x0,gs,Options=opts_fmin);
                Theta = solglo.x;
            case 'simulanneal'  % Requires the global optimization tools box
                options_Snea = optimoptions(@simulannealbnd,'HybridFcn', ...
                    {@fmincon,opts_fmin},'OutputFcn', @samyoutput);
                [solsim,fvalsim,eflagsim,outputsim] =...
                    solve(probSSE,x0,"Solver","simulannealbnd",...
                    "Options",options_Snea); 
                Theta = solsim.x;
            case 'ga'
                opts_ga = optimoptions(@ga,'HybridFcn',{@fmincon, ...
                    opts_fmin});   %,'OutputFcn', @gamyoutput);
                [solga,fvalga,eflagga,outputga] = ...
                    solve(probSSE,x0,"Solver","ga","Options",opts_ga); 
                Theta = solga.x;
            case 'gamultiobj'  % Requires the global optimization tools box
                [solgam,fvalgam,eflaggam,outputgam] =...
                    solve(probSSE,x0,"Solver","gamultiobj");
                Theta = solgam.x;
            case 'multis_lsq'
                [solmull,fvalmull,eflagmull,outputmull] = ...
                solve(probDiff,x0,ms,MinNumStartPoints=50, ...
                Options=opts_lsq);
                Theta = solmull.x;
            case 'multis_unc'
                [solmulu,fvalmulu,eflagmulu,outputmulu] = ...
                    solve(probSSEUnc,y0,ms,'Solver','fminunc', ...
                    'MinNumStartPoints',50,'Options',opts_unc);
                Theta = solmulu.y;
            case 'multis_con'
                [solmulc,fvalmulc,eflagmulc,outputmulc] =...
                    solve(probSSE,x0,ms,'Solver','fmincon', ...
                    'MinNumStartPoints',50,'Options',opts_fmin);
                Theta = solmulc.x;
            case 'bayesopt'
                fun=@(in) SSE_Bayesopt(in);
                X1=optimizableVariable('x',[lb(1) ub(1)]);
                X2=optimizableVariable('y',[lb(2) ub(2)]);
                X3=optimizableVariable('z',[lb(3) ub(3)]);
                vars=[X1,X2,X3];
                
                results= bayesopt(fun,vars,'UseParallel',false,'IsObjectiveDeterministic',true, ...
                    'AcquisitionFunctionName','expected-improvement-plus', ...
                    'MaxObjectiveEvaluations',100,'Verbose',0,'OutputFcn',@bayesoutputfun);
                Theta=table2array(results.XAtMinObjective);
                Fval=results.MinObjective;
                %ThetaB=results.XAtMinEstimatedObjective;
                
                %-------------------------------------------------------------------------
                %                The following Algorithms are not time consuming
                %                And additional gain by running these algorithm
                %                is not much!
                %-------------------------------------------------------------------------
            case 'SPUCI'            %constrained minimization using
                %the Shuffled Complex Evolution with
                % principal components analysis
                ngs=3; maxn=1000; kstop=10; pcento=0.1; peps=0.001; %iseed=10;
                iniflg=0; DSP='No';
                [Theta,~]= SPUCI(@(Theta) SSE(Theta), ...
                    xx0,lb,ub,maxn,kstop,pcento,peps,ngs,iniflg,DSP);
            case 'aeDE'                %constrained minimization using
                %Diferential Evolution algorithm
                format long
                threshold = 1e-4;               % Threshold for mutation scheme
                nvars = length(xx0);                     % Number of variables
                Objf  = @(Theta) SSE(Theta);% Objective function
                Const = @ (Theta)CookNonlinB(Theta);    % Constraint function
                Opts1.Popsize  = 100;       % Population size
                Opts1.tol      = 1e-5;      % Tolerance of stopping condition
                Opts1.Totalgen = 1000;      % Maximum iteration
                Opts1.eps2in   = 1.5;       % Parameter of penalty function
                Opts1.eps2max  = 6;         % Parameter of penalty function
                Opts1.Display  = 'No';      % Yes or No to display results
                [Theta,~,~] = aeDE(nvars,Objf,Const,lb,ub,Opts1);
            case 'purecmaes'  %unconstrained minimization using
                % Covariance Matrix Adaptation Evolution Strategy
                strfitnessfct = @(Theta) SSE(Theta);  % name of objective/fitness function
                Theta_P=purecmaes(strfitnessfct,xx0');
                Theta=Theta_P';
            case 'SCE_UA'
                ngs=3; maxn=10000; kstop=10; pcento=0.1; peps=0.001; DSP='No';
                iniflg=0;
                [Theta,~]= sceua(@(Theta) SSE(Theta), ...
                    xx0,lb,ub,maxn,kstop,pcento,peps,ngs,iniflg,DSP);
            case 'DE'
                format long
                threshold = 1e-4;               % Threshold for mutation scheme
                nvars = length(xx0);                  % Number of variables
                Objf  = @(Theta) SSE(Theta);          % Objective function
                Const = @ (Theta)CookNonlinB(Theta);  % Constraint function
                Options1.Popsize  = 100;          % Population size
                Options1.tol      = 1e-5;        % Tolerance of stopping condition
                Options1.Totalgen = 1000;        % Maximum iteration
                Options1.eps2in   = 1.5;         % Parameter of penalty function
                Options1.eps2max  = 6;           % Parameter of penalty function
                Options1.Display  = 'No';       % Yes or No to display 
                [Theta,~,~]  = DE2(nvars,Objf,Const,lb,ub,Options1);
        end                  % end of switch
        
        %--------------------------------------------------------------------------
        %               Calculate the Hessian for the parameter
        %                          uncertainty
        %------------------------------------------------------------------
        fprintf('                 fminunc \n');
        [Theta_unc,fval_unc,~,~,~,hessian_unc] = fminunc(@(Theta) ...
            SSE(Theta),Theta,opts_unc);
        fprintf('                 fmincon \n');
        [Theta_con,fval_con,~,~,~,~,hessian_con]= ...
            fmincon(@(Theta)SSE(Theta),Theta, ...
            A,b,[],[],lb,ub,[], opts_sqp);
        
        %------------------------------------------------------------------
        %      Store the output parameters and the SSE for each iteratation
        %------------------------------------------------------------------
        ErrUnc=sqrt(diag(inv(hessian_unc))).'; % the pseudo inverse of a matrix
        ErrCon=sqrt(diag(inv(hessian_con))).';
        
        if nvar<5 
            for iii=1:5-nvar
                Ins_zeros(iii)=0;
            end
        end
        ThetaSSUnc{k,j}=[Theta_unc(1:nvar-1),Ins_zeros,Theta_unc(nvar),fval_unc];
        ThetaSSCon{k,j}=[Theta_con(1:nvar-1),Ins_zeros,Theta_unc(nvar),fval_con];
              
        %------------------------------------------------------------------
        %             Calculate and store the parameter
        %               uncertainty from the hessian
        %------------------------------------------------------------------
        ThetaErrSSUnc{k,j}=[ErrUnc(1:nvar-1),Ins_zeros,ErrUnc(nvar),fval_unc];
        ThetaErrSSCon{k,j}=[ErrCon(1:nvar-1),Ins_zeros,ErrCon(nvar),fval_con];
         clear Ins_zeros; % because this vector change sizes in different iterations
         
        %------------------------------------------------------------------
        %            Output and store predicted 18O, 2H
        %            sap flux, and root water uptake PDF
        %------------------------------------------------------------------
        [~,C2H_unc,C18O_unc,Qtt_unc]=Cook_SSE(Theta_unc,mData,ModelType,Scaling,Seg_D);
        [~,C2H_con,C18O_con,Qtt_con]=Cook_SSE(Theta_con,mData,ModelType,Scaling,Seg_D);
        %     [~,C2H_unc,C18O_unc,Qtt_unc]=Cook_SSE(Theta_unc,ModelType,A);
        %     [~,C2H_con,C18O_con,Qtt_con]=Cook_SSE(Theta_con,ModelType,A);
        C2H18O_unc{k,j}=[C2H_unc,C18O_unc];
        C2H18O_con{k,j}=[C2H_con,C18O_con];
        Qttunc_con{k,j}=[SapFlux(:,1),Qtt_unc,Qtt_con];
        rootwater_con{k,j}=arrayfun(@(z) RootUptakePDF(Theta_con,z,Zr, ...
            ModelType,Seg_D),depth);
    end                       % end of loop k  of models

    %----------------------------------------------------------------------
    %                     Plotting Results
    %----------------------------------------------------------------------
    fig=cell(NofAlgorithms,1);
    fig{j}=figure ('Name',Algorithm(j));
    ColorChoice=[0.3467,    0.5360,    0.6907;...
        0.9153,    0.2816,    0.2878;...
        0.4416,    0.7490,    0.4322;...
        1.0000,    0.5984,    0.2000;...
        0.6769,    0.4447,    0.7114;...
        0     ,    0     ,    0     ];
    LinSty={'-','--',':','-.'};
    numLineStyle=3;
    tiledlayout(1,6);
    
    %----------------------------------------------------------------------
    ax1=nexttile;         %%%% Soil matric potential
    plot(SoilPsi(:,2),SoilPsi(:,1),'k:d',arrayfun(@(z) soilPsi(z), ...
        psidepth),psidepth,'k-','MarkerSize',6,'MarkerFaceColor', ...
        'k','LineWidth',2);
    set(gca, 'YDir','reverse')
    title(ax1,'Soil Potential')
    ylabel(ax1,'Depth (cm)')
    xlabel(ax1,'Matric Potential (MPa)')
    axis(ax1,[-inf 0 0 250])
    grid on
    legend('Measured','Interpolated')
    
    %----------------------------------------------------------------------
    nexttile;             %%%% measured plant leaf potential
    LeafPotentialTimeS=timeseries(LeafPsi(:,2),LeafPsi(:,1), ...
        "Name",'Measured');
    LeafPotentialTimeS.Name='Leaf \psi (MPa)';
    LeafPotentialTimeS.TimeInfo.Units='Hours';
    plot(LeafPotentialTimeS,'k:s','MarkerSize',6,'MarkerFaceColor', ...
        'k','LineWidth',2);
    hold on
    LeafPotentialTimeS1=timeseries(mData.leafPsiv, leafPsitime,"Name", ...
        'interpolated');
    LeafPotentialTimeS1.Name='Leaf \psi (MPa)';
    LeafPotentialTimeS1.TimeInfo.Units='Hours';
    plot(LeafPotentialTimeS1,'k-','MarkerSize',6,'MarkerFaceColor', ...
        'k','LineWidth',2);
    grid on
    legend('Measured','Interpolated')

    %----------------------------------------------------------------------
    nexttile;             %%%% Measured and predicted sap flux density
    legendinfo=cell(length(modelTypeArray)+2,1);

    for k=1:length(modelTypeArray)
        %modelType=modelTypeArray(k);
        %plot water flux from Fminunc
        Q=Qttunc_con{k,j};
        SapFluxTimeS=timeseries(Q(:,2),Q(:,1),"Name",modelTypeArray(k));
        SapFluxTimeS.Name='Sap Flux Density(cm h^-^1)';
        SapFluxTimeS.TimeInfo.Units='Hours';
        plot(SapFluxTimeS,'color', ...
            ColorChoice(1+fix(k/(numLineStyle+0.1)),:),'linestyle', ...
            LinSty(1+mod(k-1,numLineStyle)),'LineWidth',2);
        legendinfo{k}=modelTypeArray(k);
        hold on
    end
    SapFluxTimeS=timeseries(SapFlux(:,2),SapFlux(:,1),"Name",'Measured');
    SapFluxTimeS.Name='Sap Flux Density(cm h^-^1)';
    SapFluxTimeS.TimeInfo.Units='Hours';
    legendinfo{length(modelTypeArray)+1}="Measured";
    p1=plot(SapFluxTimeS,'k:v','MarkerSize',6,'MarkerFaceColor', ...
        'k','LineWidth',2);
    grid on
    legend(p1, legendinfo(end-1,:),'location','southeast')
    hold off
    
    %----------------------------------------------------------------------
    ax4=nexttile;        %%%% Soil water 2H composition as a function depth
                            % and the corresponding plant water 2H compositions
    hold on
    for k=1:length(modelTypeArray)
        %modelType=modelTypeArray(k);
        C=C2H18O_con{k,j};
        plot(C(1)*ones((length(depth)),1),depth,'color',ColorChoice(1+ ...
            fix(k/(numLineStyle+0.1)),:) ...
            ,'linestyle',LinSty(1+mod(k-1,numLineStyle)),'LineWidth',2);
        set(gca,'YDir','reverse')
        legendinfo{k} = modelTypeArray(k);
    end
    p1=plot(ax4,Delta(:,2),Delta(:,1),'k:*','LineWidth',2);
    legendinfo{length(modelTypeArray)+1}='Soil Water \delta^2H';
    p2=plot(ax4,Delta(1,4)*ones(length(depth),1),depth, ...
        'g:pentagram','MarkerSize',6, ...
        'MarkerIndices',1:30:length(depth), 'LineWidth',2);
    legendinfo{length(modelTypeArray)+2}='Measured Sap \delta^2H ';
    legend([p1 p2],legendinfo(end-1:end,:),'location','southwest')
    grid on
    title(ax4,'\delta^2H')
    ylabel(ax4,'Depth (cm)')
    xlabel(ax4,['\delta^2H(',char(8240),')'])
    hold off
    
    %----------------------------------------------------------------------
    ax5=nexttile;       %%%% Soil water 18O composition as a function depth
                           % and the corresponding plant water 18O compositions
    hold on
    for k=1:length(modelTypeArray)
        % modelType=modelTypeArray(k);
        C=C2H18O_con{k,j}; %plot water flux from Fminunc
        plot(C(2)*ones(length(depth),1),depth, 'color', ...
            ColorChoice(1+fix(k/(numLineStyle+0.1)),:), ...
            'linestyle',LinSty(1+mod(k-1,numLineStyle)),'LineWidth',2);
        set(gca,'YDir','reverse')
        legendinfo{k}=modelTypeArray(k);
    end
    p3 = plot(ax5,Delta(:,3),Delta(:,1),'k:*','LineWidth', 2);
    legendinfo{length(modelTypeArray)+1}='Soil water \delta^1^8O';
    p4=plot(ax5,Delta(1,5)*ones(length(depth),1),depth, ...
        'g:pentagram','MarkerSize',6, ...
        'MarkerIndices',1:30:length(depth), 'LineWidth',2);
    legendinfo{length(modelTypeArray)+2}='Measured Sap \delta^1^8O';
    legend([p3 p4],legendinfo(end-1:end,:),'location','southwest')
    grid on
    title(ax5,'\delta^1^8O')
    ylabel(ax5,'Depth (cm)')
    xlabel(ax5,['\delta^1^8O (', char(8240),')'])
    hold off
    
    %----------------------------------------------------------------------
    ax6=nexttile;        %%%% Root water uptake probability desnity functions
    hold on
    for k=1:length(modelTypeArray)
        %modelType=modelTypeArray(k);
        %plot water flux from Fminunc
        plot(rootwater_con{k,j},depth,'color', ...
            ColorChoice(1+fix(k/(numLineStyle+0.1)),:), ...
            'linestyle',LinSty(1+mod(k-1,numLineStyle)),'LineWidth',2);
        set(gca,'YDir','reverse');
        xlim([0 0.2])
        ylim([0 Zr])
        legendinfo{k}=modelTypeArray(k);
    end
    title(ax6,'RootWaterUptake Pattern')
    xlabel(ax6,'Root Water uptake Probability Density (cm^-^1)')
    ylabel(ax6,'Depth(cm)')
    legend (legendinfo(1:end-2,:), 'location','southeast')
    hold off
    
    %----------------------------------------------------------------------
    %                     Write the result to files
    %----------------------------------------------------------------------
    NoModels=length(modelTypeArray);
    algorName = cell(NoModels+1,1);
    algorName{1} ="Algorithm";
    for kk=1:NoModels
        algorName{kk+1}=Algorithm(j);
    end
    modelName= ["Model Name";modelTypeArray'];        %model name array
    
    %nvar=length(xx0);
        
    for kkk=1:5              %Assuming 5 parameters
        paramName1{kkk}=['Theta(',num2str(kkk),')']; %"Beta", "Ks","SSE"];         %Parameter name array
        paramName2{kkk}=['P_std(',num2str(kkk),')'];% "Beta-std", "Ks-std","SSE"]; %Parameter STD array
    end
    paramName1{6}='SSE__Con';
    paramName2{6}='SSE__Unc';
    %If you see complex numbers, the excel output results for SSE may not be reliable because
    %MatLab stores complex numbers in a different way than excel

        Theta_All{j}=[algorName(2:end,1),modelName(2:end,1),cell2mat(ThetaSSCon(:,j)),...
            cell2mat(ThetaErrSSCon(:,j)),cell2mat(ThetaSSUnc(:,j)),...
            cell2mat(ThetaErrSSUnc(:,j))];
        
        exportgraphics(fig{j},[fileName(1:end-5),'_figure5.gif'],'Append',true);
        
        
        writematrix(Theta_con, [fileName(1:end-5),'_parm.xlsx']);

        writecell(Qttunc_con, [fileName(1:end-5),'_PredictedFlux.xlsx']);
 
        writecell(rootwater_con, [fileName(1:end-5),'_rootwater_con.xlsx']);
        writecell(rootwater_unc, [fileName(1:end-5),'_rootwater_uncon.xlsx']);
 
        writecell(C2H18O_con, [fileName(1:end-5),'_PlantIsotopes.xlsx']);
 
        writematrix(modelName, [fileName(1:end-5),'_ModelName.xlsx']);
        
        

    %----------------------------------------------------------------------
end                %% end of loop j on Algorithm
%--------------------------------------------------------------------------
Headers={'Algorithm','Model Name',paramName1{1:end},paramName2{1:end},...
    paramName1{1:end},paramName2{1:end}};
writematrix([Headers;vertcat(Theta_All{:})], [fileName(1:end-5),'_parm6.xlsx']);
disp([Headers;vertcat(Theta_All{:})]);

%--------------------------------------------------------------------------
% %comparing different Matlab Algorithms
%--------------------------------------------------------------------------
%  sols = [sollsq.x;solcon.x; solpto.x;...
%     solsur.x; solunc.y; solpat.x; solpar.x; solglo.x;...
%     solsim.x; solga.x;solgam.x; solmull.x;...
%     solmulu.y ;solmulc.x];
% fvals = [fvallsq.x;fvalcon.x; fvalpto.x;...
%     fvalsur.x; fvalunc.y; fvalpat.x; fvalpar.x; fvalglo.x;...
%     fvalsim.x; fvalga.x;fvalgam.x; fvalmull.x;...
%     fvalmulu.y ;fvalmulc.x];
% fevals = [outputflsq.funcCount;  outputpcon.funccount;...
%     outputgpto.funccount; outputsur.funccount;outputunc.funccount;...
%     outputpat.funccount;  outputpar.funccount;outputglo.funccount; ...
%     outputsim.funccount;  outputga.funccount; outputgam.funccount; ...
%     outputmull.funccount; outputmulu.funccount; outputmulc.funccount];
% minfval=min(fvals);
% [row,col]=find(fvals,minfval,1);
% bestTheta=sols(row,col);
% if j = 1:N
%     stats(j).bestTheta=bestTheta;
%     stats{j}.minfvals=minfval;
%     stats{j}.modelType=ModelType;
% end
% for j=1:NofModels-1
%     for jj=j+1:NofModels
%         Theta1=stats(j).bestTheta; Theta2=stats(jj).bestTheta;
% SSE_comp=@(p) Cook_SSE_comp(p,Theta1,Theta2,stats(j).ModelType,stats(jj).ModelType,mData, ...
%          Scaling,Seg_depth);
%  p=0.5;
%  [p_unc,fval_unc,~,~,~,hessian_unc] = fminunc(@SSE_comp,p,opts_unc);
%  parm(j).mType1=stats(j).ModelType;     parm(j).mType2=stats(jj).ModelType;
%  parm(j).Theta1=Theta1;     parm(j).Theta1=Theta2;      parm(j).p=p_unc;
%  parm(j).psigma=1/hessian_unc(1,1);     parm(j).fval=fval_unc;
% stats.Properties.RowNames = ["lsqnonlin" "fmincon" "paretpsearch"  ...
%     "surrogateopt" "fminunc" "patternsearch" "particleswarm" "global" ...
%     "simulannealbnd" "ga" "gamultiobj" "multistart_lsq" "multistart_unc"...
%     "multistart_con"   ];
% stats{j}.Properties.VariableNames = ["Solution" "Objective" "# Fevals"];
% disp(stats{j})

% %--------------------------------------------------------------------------
% %     Objective function is presented as a nested function
% %--------------------------------------------------------------------------
function [Diff, C2H, C18O, sapFluxPv]= Cook_SSE(theta,mData,type,Scaling,Seg_depth)
   psidepth=mData.psidepth;   
   mm=size(psidepth,1);
   Zr=max(psidepth);
   RootWUv=arrayfun(@(z) RootUptakePDF(theta,z,Zr,type, ...
       Seg_depth),mData.psidepth);
   SPsi_RWv= mData.soilPsiv.*RootWUv;
   IntSPsiRw_dz = sum(SPsi_RWv)*Zr/mm; 
   IntRwCH_dz = sum(RootWUv.*mData.C2H_Mv)*Zr/mm;
   IntRwCO_dz = sum(RootWUv.*mData.C18O_Mv)*Zr/mm;
   IntSPsiRwCH_dz = sum(SPsi_RWv.*mData.C2H_Mv)*Zr/mm;
   IntSPsiRwCO_dz = sum(SPsi_RWv.*mData.C18O_Mv)*Zr/mm;
      
   %The following do not follow the formula presented in Cook, because
   %theta(3) and summation increment in time cancelled out when
   %normalization is done. So Qtt is the totle flux/theta(3)/delt_t;
   %Qt2H=At2H/theta(3)/delt_t
   Qtt   = sum(-mData.leafPsiv)/mm+IntSPsiRw_dz;
   Qt2H  = sum(-mData.leafPsiv).*IntRwCH_dz/mm+IntSPsiRwCH_dz;
   Qt18O = sum(-mData.leafPsiv).*IntRwCO_dz/mm+IntSPsiRwCO_dz; 
   
   C2H = Qt2H/Qtt*6.420135-1000;   % normalization and convert to del notion
   C18O= Qt18O/Qtt/2.0052-1000;    % normalization and convert to del notion
   sapFluxPv =  theta(end)*(-mData.leafPsiAtSapTv+IntSPsiRw_dz); %Predicted sap flux
   Diff=[(C2H-mData.Delta(1,4))*Scaling(1);(C18O-mData.Delta(1,5))*Scaling(2);   ...
       (sapFluxPv-mData.SapFlux(:,2))*Scaling(3)];  % magnifying the sap flux by scaling

 end              % end of the Objective function


 function fval= SSE_Bayesopt(in)
        theta(1)=in.x;
        theta(2)=in.y;
        theta(3)=in.z;
        fval=SSE(theta);
 end

%--------------------------------------------------------------------------
%   output function handler for particleswarm algorithm 
%--------------------------------------------------------------------------
    function stop = psmyoutput(optimValues,state)
        stop = false;
        if ~isequal(state,'iter')
            bestx=optimValues.bestx;
            bestf = optimValues.bestfval;
            for jjj=1:length(bestx)
                fprintf('     %8.4f',bestx(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
% output function for lsqnonlin,surrogateopt,fmincon,fminunc algorithm 
%--------------------------------------------------------------------------
    function stop = myoutput(x,optimValues,state)
        stop = false;
        if ~isequal(state,'iter')
            bestf = SSE(x);
            for jjj=1:length(x)
                fprintf('     %8.4f',x(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
%                 Output function for Surrogateopt,fmincon,fminunc 
%--------------------------------------------------------------------------
    function stop = fminmyoutput(x,optimValues,state)
        stop = false;
        if ~isequal(state,'iter')
            bestf=optimValues.fval;
            for jjj=1:length(x)
                fprintf('     %8.4f',x(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
%      output function for Multistart anf Global algorithm 
%--------------------------------------------------------------------------
    function stop = multimyoutput(optimValues,state)
        stop = false;
        if ~isequal(state,'iter')
            bestx = optimValues.bestx;
            bestf = optimValues.bestfval;
            for jjj=1:length(bestx)
                fprintf('     %8.4f',bestx(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
%   output function handler for Patternsearch and Paretosearch Algorithm
%--------------------------------------------------------------------------

    function [stop,options,optchanged] = pamyoutput(optimValues,options,flag)
        optchanged = false;
        stop = false;
        if ~isequal(flag,'iter')
            bestx = optimValues.x;
            bestf = optimValues.fval;
            for jjj=1:length(bestx)
                fprintf('     %8.4f',bestx(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
%     output function handler for Simulated Annealing Algorithm
%--------------------------------------------------------------------------

    function [state,options,optchanged] = samyoutput(options,optimvalues,flag)
        optchanged = false;
        state = false;
        if ~isequal(flag,'iter')
            bestx = optimvalues.bestx;
            bestf = optimvalues.bestfval;
            for jjj=1:length(bestx)
                fprintf('     %8.4f',bestx(jjj));
            end
            fprintf('     %8.5f\n',bestf);
        end
    end

%--------------------------------------------------------------------------
%             output function handler for Genetic Algorithm
%--------------------------------------------------------------------------
function [state, options,optchanged] = gamyoutput(options,state,flag)
   optchanged = false;
   if ~isequal(flag,'iter')
       ibest = state.Best(end);
       ibest = find(state.Score == ibest,1,'last');
       bestx = state.Population(ibest,:);
       bestf = SSE(bestx);
       for jjj=1:length(bestx)
            fprintf('     %8.4f',bestx(jjj));
       end
            fprintf('     %8.5f\n',bestf);
  end
end

function stop = bayesoutputfun(results,state)
   stop = false;
     if ~isequal(state,'iter')
        bestx = table2array(results.XAtMinObjective);
        bestf=results.MinObjective;
        for jjj=1:length(bestx)
            fprintf('     %8.4f',bestx(jjj));
        end
            fprintf('     %8.5f\n',bestf);
    end
end

%--------------------------------------------------------------------------
%                    Nonlinear constraint function handler
%                          for Weibull function
%--------------------------------------------------------------------------

    function [c,ceq] = nonlcon(x,modelType)
        if strcmp(modelType,'Weibull')
            c(1)=x(1)*gamma(1+1/x(2))-250;
            c(2)=10-x(1)*x(1)*(gamma(1+2/x(2))-(gamma(1+1/x(2)))^2);
            ceq=[];
        else
            c=[];
            ceq=[];
        end
    end

%--------------------------------------------------------------------------
    function ZZ = segMent2Continu(x,z,y) %x is a value and LeafPsi is a 2 b yn matrix
        nn =length(z);
        if x <( z(1)+z(2))/2
            ZZ=y(1);
            return
        end

        for ii = 2:(nn-1)
            if ((z(ii-1)+z(ii))/2 <= x)  && (x < (z(ii)+z(ii+1))/2)
                ZZ=y(ii);
                return
            end
        end

        if x >=(z(nn-1)+z(nn))/2
            ZZ=y(nn);
            return
        end
        ZZ=9999.99;
    end   %end of the function

end       %end of main program

%--------------------------------------------------------------------------
%                     End of the program
%--------------------------------------------------------------------------

