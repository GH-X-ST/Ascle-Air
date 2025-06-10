function [emax,W1,W2,Cinf] = simWCsyn_HinfLoopShaping(P,LoopShape_bounds,Weight_bounds,CondNo_bounds,WeightStruc,omega,Cinit)

% This function simultaneously synthesises by one algorithm loop-shaping
% weights 'W1/W2' and a robustly stabilising controller 'Cinf' which
% achieves a pre-specified performance level and maximises the robust
% stability margin 'emax' to normalised coprime factor uncertainty. Here,
% the notation is such that 'W2' pre-multiplies 'P' and 'W1'
% post-multiplies it to form the shaped plant.
% 
% USAGE  :  [emax,W1,W2,Cinf] = simWCsyn_HinfLoopShaping(P,LoopShape_bounds, ...
%                               Weight_bounds,CondNo_bounds,WeightStruc,omega,Cinit)
% 
% INPUTS :  P                 = SYSTEM matrix representing nominal plant.
%           LoopShape_bounds  = a 2x1 SYSTEM matrix such that the magnitude
%                               of the first (second) transfer function
%                               gives a upper (lower) bound for the desired
%                               loop-shape. 
% 
%           THE FOLLOWING INPUT ARGUMENTS ARE OPTIONAL
% 
%           Weight_bounds     = a 2x2 SYSTEM matrix such that the magnitude
%                               of each individual transfer function is an
%                               upper/lower bound for the singular values
%                               of W1/W2 as follows:
%                                                 W1    W2
%                                          upper [1,1   1,2]
%                                          lower [2,1   2,2].
%                               (DEFAULT value is '[1e10 1; 1e-10 1]').
%           CondNo_bounds     = a 1x2 SYSTEM matrix such that the magnitude
%                               of the first (second) transfer function gives
%                               an upper bound for the condition number of 
%                               W1 (W2). (DEFAULT value is '[20 20]').
%           WeightStruc       = 'DiagDiag' for diagonal loop-shaping weights W1/W2.
%                               'DiagFull' for diagonal W1 and non-diagonal W2.
%                               'FullDiag' for non-diagonal W1 and diagonal W2.
%                               'FullFull' for non-diagonal W1 and non-diagonal W2.
%                               (DEFAULT value is 'DiagDiag').
%           omega             = a vector containing the frequencies at
%                               which the internal LMI constraints will be
%                               checked. (DEFAULT value is 'logspace(-4,4,100)').
%           Cinit             = this is an internally stabilising controller 
%                               used to initialise the algorithm (for the
%                               unshaped plant).
% 
% OUTPUTS:  emax              = best attained robust stability margin by
%                               this algorithm.
%           W1                = loop-shaping weight pre-multiplying the
%                               nominal plant P.
%           W2                = loop-shaping weight post-multiplying the
%                               nominal plant P.
%           Cinf              = robust stabilising controller for shaped
%                               plant 'Ps = W2*P*W1' which achieves robust
%                               stability margin 'emax'. To compute the
%                               loop-shaping controller 'C' for nominal
%                               plant 'P', simply let 'C = W1*Cinf*W2'.
%
% Copyright: Alexander Lanzon - 17 July 2000.


% Checking number of inputs and outputs
if (nargin < 2 ) | (nargin > 7)
  disp('ERROR: Wrong number of inputs');
  return;
end
 
if (nargout == 0) | (nargout > 4)
  disp('ERROR: Wrong number of outputs');
  return;
end

% Checking required inputs
[mattype,m,n,s] = minfo(P);
if ~(strcmp(mattype,'syst'))
  disp('ERROR: Nominal Plant "P" not a SYSTEM matrix');
  return;
end

[mattype,rowd,cold,num] = minfo(LoopShape_bounds);
if ~(strcmp(mattype,'syst')) | (rowd ~= 2) | (cold ~= 1)
  disp('ERROR: Variable "LoopShape_bounds" not properly defined');
  return;
end

% Checking and setting optional inputs
if (nargin >= 3)
  [mattype,rowd,cold,num] = minfo(Weight_bounds);
  if ~(strcmp(mattype,'syst') | strcmp(mattype,'cons')) | (rowd ~= 2) | (cold ~= 2)
	disp('ERROR: Variable "Weight_bounds" not properly defined');
	return;
  end
else
  Weight_bounds = [1e10 1; 1e-10 1];
end

if (nargin >= 4)
  [mattype,rowd,cold,num] = minfo(CondNo_bounds);
  if ~(strcmp(mattype,'syst') | strcmp(mattype,'cons')) | (rowd ~= 1) | (cold ~= 2) 
	disp('ERROR: Variable "CondNo_bounds" not properly defined');
	return;
  end
else
  CondNo_bounds = [20 20];
end

if (nargin >= 5)
  if ~((strcmp(WeightStruc,'DiagDiag')) | (strcmp(WeightStruc,'DiagFull')) ...
        | (strcmp(WeightStruc,'FullDiag')) | (strcmp(WeightStruc,'FullFull')))
	disp('ERROR: Variable "WeightStruc" not properly defined');
	return;
  end
else
  WeightStruc = 'DiagDiag';
end

if (nargin >= 6)
  [mattype,rowd,cold,num] = minfo(omega);
  if ~(strcmp(mattype,'cons')) | (rowd ~= 1)
	disp('ERROR: Variable "omega" not properly defined');
	return;
  end
else
  omega = logspace(-4,4,100);
end

if (nargin >= 7)
  [mattype,rowd,cold,num] = minfo(Cinit);
  if ~(strcmp(mattype,'syst')) | (rowd ~= n) | (cold ~= m)
	disp('ERROR: Variable "Cinit" not properly defined');
	return;
  end
else
  Cinit = ncfsyn(P,1.2);
end
C = Cinit;

% Checking if Weight_bounds are [1;1]
W1_bounds_FLAG = 'false';
dummy = linfnorm(msub(sel(Weight_bounds,':',1),[1;1]));
if (norm(dummy(1:2)) == 0)
  W1_bounds_FLAG = 'true';
end

W2_bounds_FLAG = 'false';
dummy = linfnorm(msub(sel(Weight_bounds,':',2),[1;1]));
if (norm(dummy(1:2)) == 0)
  W2_bounds_FLAG = 'true';
end

if (strcmp(W1_bounds_FLAG,'true') & strcmp(W2_bounds_FLAG,'true'))
  disp('ERROR: It is not possible to have Weight_bounds = [1 1; 1 1]');
  return;
end

% Determining axis limits for loop-shape
loopshape_upper_axis_limit = var2con(vnorm(frsp(sel(LoopShape_bounds,2,1),omega(1))),omega(1))*1e3;
loopshape_lower_axis_limit = var2con(vnorm(frsp(sel(LoopShape_bounds,1,1),omega(end))),omega(end))/1e3;


% Plotting graphs of original data given
figure(1); hold off;
           HANDLEC = vplot('liv,m',frsp(1e-10,1));
           xlabel('Frequency (radians/sec)');
	       ylabel('Pointwise Robust Stability Margin');
		   axis([omega(1) omega(end) 0 1]); hold on; zoom on;

figure(2); vplot('liv,lm',vnorm(frsp(sel(LoopShape_bounds,1,1),omega)),'--', ...
				          vnorm(frsp(sel(LoopShape_bounds,2,1),omega)),'--', ...
				          frsp(1,omega),':', ...
				          vsvd(frsp(P,omega)),'-');
           xlabel('Frequency (radians/sec)');
		   ylabel('Singular Values of Nominal Plant & Loop-Shape Boundaries');
		   axis([omega(1) omega(end) loopshape_lower_axis_limit loopshape_upper_axis_limit]); 
		   zoom on;

figure(3); if strcmp(W1_bounds_FLAG,'true')
             subplot(1,1,1);
             vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,2),omega)),'--', ...
				            vnorm(frsp(sel(Weight_bounds,2,2),omega)),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for \sigma_i (W_2)'); zoom on;
		   elseif strcmp(W2_bounds_FLAG,'true')
             subplot(1,1,1);
			 vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,1),omega)),'--', ...
				            vnorm(frsp(sel(Weight_bounds,2,1),omega)),'--');
             xlabel('Frequency (radians/sec)');
			 ylabel('Boundaries for \sigma_i (W_1)'); zoom on;
		   else
			 subplot(2,1,1);
			 vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,1),omega)),'--', ...
				            vnorm(frsp(sel(Weight_bounds,2,1),omega)),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for \sigma_i (W_1)'); zoom on;

             subplot(2,1,2);
             vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,2),omega)),'--', ...
				            vnorm(frsp(sel(Weight_bounds,2,2),omega)),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for \sigma_i (W_2)'); zoom on;
		   end
 
figure(4); if strcmp(W1_bounds_FLAG,'true')
             subplot(1,1,1);
             vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,2),omega)),'--', ...
				           frsp(1,omega),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for k(W_2)'); zoom on;
		   elseif strcmp(W2_bounds_FLAG,'true')
             subplot(1,1,1);
             vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,1),omega)),'--', ...
				           frsp(1,omega),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for k(W_1)'); zoom on;
		   else
             subplot(2,1,1);
             vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,1),omega)),'--', ...
				           frsp(1,omega),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for k(W_1)'); zoom on;

             subplot(2,1,2);
             vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,2),omega)),'--', ...
				           frsp(1,omega),'--');
             xlabel('Frequency (radians/sec)');
		     ylabel('Boundaries for k(W_2)'); zoom on;
		   end

figure(5); vplot('liv,lm',frsp(1,omega),':', ...
				          vsvd(frsp(Cinit,omega)),'-');
           xlabel('Frequency (radians/sec)');
		   ylabel('Singular Values of Initial Controller'); zoom on; pause(.5);

		   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF ACTUAL PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set Uhat and Vhat for diagonal/non-diagonal weights
[Uhat,Vhat] = unitary_svd_approx(P,omega,WeightStruc);


% Begin iterations
ii = 0; flag = 'NOTEXIT'; emaxIN = -1; decvarsIN = 99999;
while ~(strcmp(flag,'exit'))
  ii = ii + 1;

  disp(' '); disp(' '); disp(['W-iteration: ' num2str(ii)]); disp(' ');
  [W1_tmp,W2_tmp,rho(ii,:),decvarsOUT] = W_iter(P,C,Uhat,Vhat,LoopShape_bounds, ...
												Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN);
  Ps = mmult(W2_tmp,P,W1_tmp);
  
  % Plotting graphs after W-iteration
  figure(2); vplot('liv,lm',vnorm(frsp(sel(LoopShape_bounds,1,1),omega)),'--', ...
				            vnorm(frsp(sel(LoopShape_bounds,2,1),omega)),'--', ...
				            frsp(1,omega),':', ...
				            vsvd(frsp(Ps,omega)),'-');
             xlabel('Frequency (radians/sec)');
		     ylabel('Singular Values of Shaped Plant');
			 axis([omega(1) omega(end) loopshape_lower_axis_limit loopshape_upper_axis_limit]);

  figure(3); if strcmp(W1_bounds_FLAG,'true')
	           subplot(1,1,1);
			   vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,2),omega)),'--', ...
				              vnorm(frsp(sel(Weight_bounds,2,2),omega)),'--', ...
				              vsvd(frsp(W2_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Singular Values of W_2');
			 elseif strcmp(W2_bounds_FLAG,'true')
	           subplot(1,1,1);
               vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,1),omega)),'--', ...
				              vnorm(frsp(sel(Weight_bounds,2,1),omega)),'--', ...
				              vsvd(frsp(W1_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Singular Values of W_1');
		     else
	           subplot(2,1,1);
               vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,1),omega)),'--', ...
				              vnorm(frsp(sel(Weight_bounds,2,1),omega)),'--', ...
				              vsvd(frsp(W1_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Singular Values of W_1');
			   
			   subplot(2,1,2);
			   vplot('liv,lm',vnorm(frsp(sel(Weight_bounds,1,2),omega)),'--', ...
				              vnorm(frsp(sel(Weight_bounds,2,2),omega)),'--', ...
				              vsvd(frsp(W2_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Singular Values of W_2');
			 end

  figure(4); if strcmp(W1_bounds_FLAG,'true')
	           subplot(1,1,1);
			   vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,2),omega)),'--', ...
				             frsp(1,omega),'--', ...
				             vcond(frsp(W2_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Condition Number of W_2'); 
             elseif strcmp(W2_bounds_FLAG,'true')
	           subplot(1,1,1);
               vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,1),omega)),'--', ...
				             frsp(1,omega),'--', ...
				             vcond(frsp(W1_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Condition Number of W_1');
		     else
	           subplot(2,1,1);
               vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,1),omega)),'--', ...
				             frsp(1,omega),'--', ...
				             vcond(frsp(W1_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Condition Number of W_1');
			   
			   subplot(2,1,2);
			   vplot('liv,m',vnorm(frsp(sel(CondNo_bounds,1,2),omega)),'--', ...
				             frsp(1,omega),'--', ...
				             vcond(frsp(W2_tmp,omega)),'-');
               xlabel('Frequency (radians/sec)');
		       ylabel('Condition Number of W_2'); 
			 end
			 pause(.5);

			 
  disp(' '); disp(' '); disp(['Controller Synthesis ' num2str(ii)]);
  [Cinf_tmp,emaxC_tmp(ii)] = ncfsyn(Ps,1);
  C = mmult(W1_tmp,Cinf_tmp,W2_tmp);
  clp = madd(mmult(daug(-eye(m),eye(n)),minv(sbs(abv(eye(m),Cinf_tmp),abv(Ps,eye(n))))),daug(eye(m),zeros(n)));
  
  % Plotting graphs after Controller Synthesis
  figure(1); set(HANDLEC,'visible','off');
             HANDLEC = vplot('liv,m',vinv(vnorm(frsp(clp,omega))),'-');
             xlabel('Frequency (radians/sec)');
	         ylabel('Pointwise Robust Stability Margin');
		     axis([omega(1) omega(end) 0 1]);

  figure(5); vplot('liv,lm',frsp(1,omega),':', ...
			                vsvd(frsp(Cinf_tmp,omega)),'-');
             xlabel('Frequency (radians/sec)');
		     ylabel('Singular Values of Controller C_\infty'); pause(.5);

			 
  disp(' '); disp(' '); disp(['Iteration ' num2str(ii) ' successfully completed!!!']);
  
  
  % Setting input variables for next iteration's "gevp" LMI optimisation
  decvarsIN = decvarsOUT;
  emaxIN    = emaxC_tmp(ii);
  
  % Asking whether to do another iteration or not!
  flag = input('Type "exit" if no more iterations required: ','s');
end


disp(' '); disp(' '); disp(['Successfully converged after ' num2str(ii) ' iterations']);

% Plotting graph for robust stability margin as iterations proceed
figure(6); plot(1:ii,emaxC_tmp,'*');
           xlabel('Iteration Number');
		   ylabel('Robust Stability Margin after Controller Synthesis');
		   axis([0 (ii+1) 0 1]); 
		   pause(.5); grid on; zoom on;
		   
		   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ACTUAL PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setting outputs
if (nargout >= 1)
  emax = emaxC_tmp(end);
end

if (nargout >= 2)
  W1 = W1_tmp;
end

if (nargout >= 3)
  W2 = W2_tmp;
end

if (nargout >= 4)
  Cinf = Cinf_tmp;
end

figure(1); hold off;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Uhat,Vhat] = unitary_svd_approx(P,omega,WeightStruc)

% ***** This is an INTERNAL function for use with "simWCsyn_HinfLoopShaping" *****
% 
% This function takes a continuous frequency-by-frequency Singular Value
% Decomposition of plant "P" and fits transfer function matrices "Uhat"
% and/or "Vhat" to the matrices of singular vectors.
% 
% Alexander Lanzon - 6 June 2000.


[mattype,m,n,num] = minfo(P);

if strcmp(WeightStruc,'DiagDiag')
  Uhat = eye(m);
  Vhat = eye(n);
  return;
end


% Take a continuous SVD
[U_g,S_g,V_g] = vsvdcont(frsp(P,omega));


% Approximate the continuous frequency matrices of singular vectors with
% system matrices in RL-infinity
disp(' '); disp(' '); disp('Computing rational approximations to matrices of singular vectors ...');

if strcmp(WeightStruc,'DiagFull')
  U_sys = rational_approx(U_g);
  V_sys = eye(n);
end

if strcmp(WeightStruc,'FullDiag')
  U_sys = eye(m);
  V_sys = rational_approx(V_g);
end

if strcmp(WeightStruc,'FullFull')
  U_sys = rational_approx(U_g);
  V_sys = rational_approx(V_g);
end

% If "U_sys" or "V_sys" is tall, then construct an all-pass dilation to it
% as given by ZDG Lemma 13.31 (valid also for RL-infinity transfer functions)
if (m > n)
  [A,B,C,D] = unpck(U_sys);
  X = lyap(A',C'*C);
  U_sys = pck(A,[B -pinv(X)*C'*null(D')],C,[D null(D')]);
elseif (n > m)
  [A,B,C,D] = unpck(V_sys);
  X = lyap(A',C'*C);
  V_sys = pck(A,[B -pinv(X)*C'*null(D')],C,[D null(D')]);
end

% Setting output variable
Uhat = U_sys;
Vhat = V_sys;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [U_sys,delta] = rational_approx(U_g)

% ***** This is an INTERNAL function for use with "unitary_svd_approx" *****
% 
% This function finds a rational approximation "U_sys" to 
% the VARYING matrix "U_g".
% 
% Alexander Lanzon - 6 June 2000.


[mattype,m,n,num] = minfo(U_g);
omega = getiv(U_g).';

% Element-by-element fit with a specified order.
disp(' ');

figure;
accum_vert = [];
for ii = 1:m
  accum_horiz = [];
  for jj = 1:n
	flag = 'y';
	while ~(strcmp(flag,'n'))
	  vplot('bode',sel(U_g,ii,jj),'y.');
	  subplot(2,1,1);	title('Bode Diagram for matrix of Singular Vectors'); 
	  zoom on; subplot(2,1,1); hold on; subplot(2,1,2); hold on;
	  
	  order = input(['Please input order of (' num2str(ii) ',' num2str(jj) ') approximation: ']);
	  approx_sys = fitsys(sel(U_g,ii,jj),order);
	  
	  vplot('bode',frsp(approx_sys,omega),'r.'); 
	  subplot(2,1,1);	title('Bode Diagram for matrix of Singular Vectors'); 
	  zoom on; subplot(2,1,1); hold off; subplot(2,1,2); hold off;
	  
	  flag = input('Press ENTER to re-approximate, "n" to go to next approximation: ','s');
	end
	accum_horiz = sbs(accum_horiz,approx_sys);
	disp(' ');
  end
  accum_vert = abv(accum_vert,accum_horiz);
end
U_fit = minrealsys(accum_vert,1e-7);
close;

% Model reducing the graph symbol of "U".
[nlcf_sys,sig] = sncfbal(U_fit);
red_nlcf_sys = hankmr(nlcf_sys,sig,sum(sig>.01),'d');
U_red = cf2sys(red_nlcf_sys);

% Setting outputs
U_sys = syscl(U_red);
delta = pkvnorm(msub(U_g,frsp(U_sys,omega)));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [W1,W2,rho,decvarsOUT] = W_iter(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN)

% ***** This is an INTERNAL function for use with "simWCsyn_HinfLoopShaping" *****
% 
% This function determines whether the plant "P" is TALL or FAT and prepares
% the input data for the function "W_iter_tall_P". The definition of all
% input/output arguments of this function is as defined in "simWCsyn_HinfLoopShaping".
% 
% Alexander Lanzon - 17 July 2000.


[mattype,m,n,s] = minfo(P);
if (m >= n)
  [W1_tmp,W2_tmp,rho_tmp,decvarsOUT_tmp] = ...
       W_iter_tall_P(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN,1);
else
  [dummy1,dummy2,rho_tmp,decvarsOUT_tmp] = ...
       W_iter_tall_P(transp(P),transp(C),cjt(transp(Vhat)),cjt(transp(Uhat)),LoopShape_bounds, ...
					 mmult(Weight_bounds,[0 1;1 0]),mmult(CondNo_bounds,[0 1;1 0]),omega,emaxIN,decvarsIN,0);
  W1_tmp = transp(dummy2);
  W2_tmp = transp(dummy1);
end

W1         = W1_tmp;
W2         = W2_tmp;
rho        = rho_tmp;
decvarsOUT = decvarsOUT_tmp;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [W1,W2,rho,decvarsOUT] = W_iter_tall_P(P,C,Uhat,Vhat,LoopShape_bounds,Weight_bounds,CondNo_bounds,omega,emaxIN,decvarsIN,tall)

% ***** This is an INTERNAL function for use with "W_iter" *****
% 
% This function performs the W-iteration for a TALL or SQUARE plant by
% defining the LMI optimisation problem, solving it and then constructing
% the weights through spectral factorisation. The definition of all
% input/output arguments of this function is as defined in "W_iter".
% 
% Alexander Lanzon - 17 July 2000.


global HANDLEW;
[mattype,m,n,s] = minfo(P);

% Checking if Weight_bounds are [1;1]
W1_bounds_FLAG = 'false';
dummy = linfnorm(msub(sel(Weight_bounds,':',1),[1;1]));
if (norm(dummy(1:2)) == 0)
  W1_bounds_FLAG = 'true';
end

W2_bounds_FLAG = 'false';
dummy = linfnorm(msub(sel(Weight_bounds,':',2),[1;1]));
if (norm(dummy(1:2)) == 0)
  W2_bounds_FLAG = 'true';
end


% Defining some commonly appearing constants
NonDiagCorrector     = daug(cjt(Uhat),minv(Vhat));
LoopShapeUpperConstr = mmult(NonDiagCorrector,abv(P,ndaug(sel(LoopShape_bounds,1,1),n)));
LoopShapeLowerConstr = mmult(NonDiagCorrector,abv(P,ndaug(sel(LoopShape_bounds,2,1),n)));
CostFunctLeftConstr  = mmult(NonDiagCorrector,sbs(zeros(m+n,m),abv(P,eye(n))));
CostFunctRightConstr = mmult(NonDiagCorrector,sbs(abv(eye(m),C),abv(P,eye(n))));

NonDiagCorrector_g     = frsp(NonDiagCorrector,omega);
LoopShapeUpperConstr_g = frsp(LoopShapeUpperConstr,omega);
LoopShapeLowerConstr_g = frsp(LoopShapeLowerConstr,omega);
CostFunctLeftConstr_g  = frsp(CostFunctLeftConstr,omega);
CostFunctRightConstr_g = frsp(CostFunctRightConstr,omega);

if strcmp(W1_bounds_FLAG,'true')
  mod_w1_upper_g = frsp(1.01,omega);
  mod_w1_lower_g = frsp(1/1.01,omega);
else
  mod_w1_upper_g = vnorm(frsp(sel(Weight_bounds,1,1),omega));
  mod_w1_lower_g = vnorm(frsp(sel(Weight_bounds,2,1),omega));
end
if strcmp(W2_bounds_FLAG,'true')
  mod_w2_upper_g = frsp(1.01,omega);
  mod_w2_lower_g = frsp(1/1.01,omega);
else
  mod_w2_upper_g = vnorm(frsp(sel(Weight_bounds,1,2),omega));
  mod_w2_lower_g = vnorm(frsp(sel(Weight_bounds,2,2),omega));
end

mod_k1_g       = vnorm(frsp(sel(CondNo_bounds,1,1),omega));
mod_k2_g       = vnorm(frsp(sel(CondNo_bounds,1,2),omega));


% Solve each optimisation problem at each individual frequency
jj = 0;
for ww = omega
  jj = jj + 1;
  
  % Define "big" constants for complex LMIs
  dummy = var2con(NonDiagCorrector_g,ww);
  bigNonDiagCorrector = [real(dummy) imag(dummy); -imag(dummy) real(dummy)];
  
  dummy = var2con(LoopShapeUpperConstr_g,ww);
  bigLoopShapeUpperConstr = [real(dummy) imag(dummy); -imag(dummy) real(dummy)];
  
  dummy = var2con(LoopShapeLowerConstr_g,ww);
  bigLoopShapeLowerConstr = [real(dummy) imag(dummy); -imag(dummy) real(dummy)];
  
  dummy = var2con(CostFunctLeftConstr_g,ww);
  bigCostFunctLeftConstr = [real(dummy) imag(dummy); -imag(dummy) real(dummy)];
  
  dummy = var2con(CostFunctRightConstr_g,ww);
  bigCostFunctRightConstr = [real(dummy) imag(dummy); -imag(dummy) real(dummy)];
  
  
  % Defining LMI optimisation problem
  setlmis([]);
  DD       = lmivar(3,diag([[1:(m+n)] [1:(m+n)]]));
  DDneg    = lmivar(3,diag([[[1:m] -[(m+1):(m+n)]] [[1:m] -[(m+1):(m+n)]]]));
  MUlower  = lmivar(3,daug((m+n+1)*eye(m),(m+n+2)*eye(n),(m+n+1)*eye(m),(m+n+2)*eye(n)));
  MUupper  = lmivar(3,daug((m+n+3)*eye(m),(m+n+4)*eye(n),(m+n+3)*eye(m),(m+n+4)*eye(n)));
  mu2lower = lmivar(3,(m+n+1));
  mu1lower = lmivar(3,(m+n+2));
  mu2upper = lmivar(3,(m+n+3));
  mu1upper = lmivar(3,(m+n+4));
  
  lmiterm([  1 1 1        0],1/(var2con(mod_w1_upper_g,ww)^2));
  lmiterm([ -1 1 1 mu1lower],1,1);
  lmiterm([  2 1 1 mu1upper],1,1);
  lmiterm([ -2 1 1        0],1/(var2con(mod_w1_lower_g,ww)^2));
  lmiterm([  3 1 1 mu1upper],1,1);
  lmiterm([ -3 1 1 mu1lower],var2con(mod_k1_g,ww)',var2con(mod_k1_g,ww));
  
  lmiterm([  4 1 1        0],(var2con(mod_w2_lower_g,ww)^2));
  lmiterm([ -4 1 1 mu2lower],1,1);
  lmiterm([  5 1 1 mu2upper],1,1);
  lmiterm([ -5 1 1        0],(var2con(mod_w2_upper_g,ww)^2));
  lmiterm([  6 1 1 mu2upper],1,1);
  lmiterm([ -6 1 1 mu2lower],var2con(mod_k2_g,ww)',var2con(mod_k2_g,ww));
    
  lmiterm([  7 1 1 MUlower],1,1);
  lmiterm([ -7 1 1      DD],bigNonDiagCorrector',bigNonDiagCorrector);
  lmiterm([  8 1 1      DD],bigNonDiagCorrector',bigNonDiagCorrector);
  lmiterm([ -8 1 1 MUupper],1,1);
  
  lmiterm([  9 1 1 DDneg],bigLoopShapeUpperConstr',bigLoopShapeUpperConstr);
  lmiterm([-10 1 1 DDneg],bigLoopShapeLowerConstr',bigLoopShapeLowerConstr);
  
  lmiterm([ 11 1 1  DD],bigCostFunctLeftConstr',bigCostFunctLeftConstr);
  lmiterm([-11 1 1  DD],bigCostFunctRightConstr',bigCostFunctRightConstr);

  lmisys = getlmis;
  
  
  % Solve LMI minimisation problem
  if emaxIN < 0.01
	[gammasqr(jj),decvarsOUT_tmp(:,jj)] = gevp(lmisys,1,[1e-3 200 0 20 1]);
  else
	[gammasqr(jj),decvarsOUT_tmp(:,jj)] = gevp(lmisys,1,[1e-3 200 0 20 1],1.2/(emaxIN*emaxIN),decvarsIN(:,jj),1.05);
  end
  
  
  % Seeing the progress of the iterations
  if (mod(jj,10) == 0)
	fprintf(1,'%d\n',jj)
  else
	fprintf(1,'%d.',jj)
  end
end

rho_tmp = 1./sqrt(gammasqr);  
modD2   = sqrt(decvarsOUT_tmp(1:m,:));
modD1  = 1./sqrt(decvarsOUT_tmp((m+1):(m+n),:));

figure(1); if (emaxIN ~= -1)
             set(HANDLEW,'visible','off');
		   end
           HANDLEW = vplot('liv,m',vpck(rho_tmp.',omega.'),':');
           xlabel('Frequency (radians/sec)');
	       ylabel('Pointwise Robust Stability Margin');
		   axis([omega(1) omega(end) 0 1]);


disp(' '); disp(' '); disp('Fitting transfer functions to magnitudes ...');
% Fitting transfer functions to magnitudes

% [dummy1,dummy2,dummy3,wgc1] = margin(pck2ss(minrealsys(sel(LoopShape_bounds,2,1))));
% [dummy1,dummy2,dummy3,wgc2] = margin(pck2ss(minrealsys(sel(LoopShape_bounds,1,1))));
% 
% lowpass_g  = vnorm(frsp(nd2sys([1/wgc1 1], ...
% 							     nconv([1/sqrt(wgc1*wgc2) 1],2)),omega));
% highpass_g = vnorm(frsp(nd2sys(conv([1/wgc1 0],[1/(1*wgc2) 1]), ...
% 							     nconv([1/sqrt(wgc1*wgc2) 1],2)),omega));
%
% THE ABOVE HAS BEEN COMMENTED OUT AS WE WISH A GOOD FIT WHERE THE ROBUST
% STABILITY MARGIN IS BAD!!!

lowpass_g  = vpck(gammasqr.',omega.');
highpass_g = vpck(gammasqr.',omega.');


if (tall == 1)
  alpha1 = '1';
  alpha2 = '2';
else
  alpha1 = '2';
  alpha2 = '1';
end

if strcmp(W1_bounds_FLAG,'true')
  D1 = eye(n);
else
  accum = [];
  for kk = 1:n
	modD1kk_g = vpck(modD1(kk,:).',omega.');
	disp(' '); disp(['Fitting D_' alpha1 ' (' num2str(kk) ',' num2str(kk) ')']); 
	figure(7); 
	sys = fitmag(modD1kk_g,mmult(lowpass_g,vinv(modD1kk_g)));
	accum = daug(accum,sys);
  end
  D1 = syscl(accum);
end

if strcmp(W2_bounds_FLAG,'true')
  D2 = eye(m);
else
  accum = [];
  for kk = 1:m
	modD2kk_g = vpck(modD2(kk,:).',omega.');
	disp(' '); disp(['Fitting D_' alpha2 ' (' num2str(kk) ',' num2str(kk) ')']); 
	figure(7); 
	sys = fitmag(modD2kk_g,mmult(highpass_g,vinv(modD2kk_g)));
	accum = daug(accum,sys);
  end
  D2 = syscl(accum);
end
close;

% Constructing Weights through Spectral Factorisation
dummy = linfnorm(msub(Uhat,eye(m)));
if (norm(dummy(1:2)) == 0)
  W2_tmp = D2;
else
  halfGsys     = minrealsys(mmult(Uhat,cjt(D2)),5e-8);
  Gsys         = mmult(halfGsys,cjt(halfGsys));
  [Gst,Gun]    = sdecomp(Gsys);
  [A1,B1,C1,D] = unpck(Gst);
  Ham          = [A1-B1*inv(D)*C1 -B1*inv(D)*B1'; C1'*inv(D)*C1 -(A1-B1*inv(D)*C1)'];
  [X1,X2]      = ric_schr(Ham);
  W2_tmp       = pck(A1,B1,sqrtm(inv(D))*(C1+B1'*X2*inv(X1)),sqrtm(D));
end

dummy = linfnorm(msub(Vhat,eye(n)));
if (norm(dummy(1:2)) == 0)
  W1_tmp = D1;
else
  halfGsys     = minrealsys(mmult(Vhat,D1),5e-8);
  Gsys         = mmult(halfGsys,cjt(halfGsys));
  [Gst,Gun]    = sdecomp(transp(Gsys));
  [A1,B1,C1,D] = unpck(Gst);
  Ham          = [A1-B1*inv(D)*C1 -B1*inv(D)*B1'; C1'*inv(D)*C1 -(A1-B1*inv(D)*C1)'];
  [X1,X2]      = ric_schr(Ham);
  W1_tmp       = transp(pck(A1,B1,sqrtm(inv(D))*(C1+B1'*X2*inv(X1)),sqrtm(D)));
end

W1         = W1_tmp;
W2         = W2_tmp;
rho        = rho_tmp;
decvarsOUT = decvarsOUT_tmp;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bigM] = big(M)

% This function computes   bigM := [real(M) imag(M); -imag(M) real(M)];
%
% USAGE  : [bigM] = big(M)
% INPUT  : M      = a CONSTANT complex/real matrix
% OUTPUT : bigM   = a CONSTANT real matrix
%
% Alexander Lanzon - 10 April 2001

if ((nargin ~= 1) | (nargout > 1))
  error('USAGE: [bigM] = big(M)');
end

bigM = [real(M) imag(M); -imag(M) real(M)];

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = nconv(IN,n)

% This function convolves input polynomial row vector 'IN' with itself
% 'n' times and outputs the result in variable 'OUT'.
% 
% USAGE  : [OUT] = nconv(IN,n)
% INPUT  : IN    = an input polynomial row vector
%          n     = how many times this vector will be convolved with itself
% OUTPUT : OUT   = resulting output polynomial row vector
%
% Alexander Lanzon - 15 May 2000

if (nargin > 2)
  disp('USAGE: [OUT] = nconv(IN,n)');
  return;
end

accum = IN;
for kk = 2:n
  accum = conv(accum,IN);
end
OUT = accum;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matout] = ndaug(matin,n)

% Stacks along the main diagonal "n" copies of "matin".
%
% USAGE  : [matout] = ndaug(matin,n)
% INPUT  : matin    = a SYSTEM/VARYING/CONSTANT matrix
%          n        = how many copies of "matin" will be diagonally augmented
% OUTPUT : matout   = resulting output SYSTEM/VARYING/CONSTANT matrix
%
% Alexander Lanzon - 5 May 2000

if (nargin == 0) | (nargin > 2)
  disp('USAGE: [matout] = ndaug(matin,n)');
  return;
end

accum = matin;
for kk = 2:n
  accum = daug(accum,matin);
end
matout = accum;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SSsys] = pck2ss(sys)

% Converts a SYSTEM matrix created by the 'pck' command to a STATE-SPACE
% object, usually created by the 'ss' command.
%
% USAGE : [SSsys] = pck2ss(sys)
% INPUT : sys     = SYSTEM matrix
% OUTPUT: SSsys   = STATE-SPACE object
%
% Alexander - 26/01/99

[a,b,c,d] = unpck(sys);
SSsys     = ss(a,b,c,d);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys] = ss2pck(SSsys)

% Converts a STATE-SPACE object created by the 'ss' command to a SYSTEM
% matrix, created by the 'pck' command.
%
% USAGE : [sys] = ss2pck(SSsys)
% INPUT : SSsys   = STATE-SPACE object
% OUTPUT: sys     = SYSTEM matrix
%
% Alexander - 26/01/99

[a,b,c,d] = ssdata(SSsys);
sys     = pck(a,b,c,d);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sysout = syscl(sys)

%  SYSCL scales a state space realization in order to
%  reduce errors in later calculations, using diagonal
%  similarity transformations, and deleting any rows
%  and columns which are zero (apart from the diagonal of A).
%  
%  *** Modified by "Alexander Lanzon - 7 June 2000" so that it 
%      does not display n1, m, p, n1.
%
%  Copyright (c) 1991-98 by MUSYN Inc. and The MathWorks, Inc.
%  $Revision: 1.3$

[a,b,c,d]=unpck(sys);
[n,n1]=size(a);
if n1~=n
  disp('A is not square')
  return
end
[n1,m]=size(b);
if n1~=n
  disp('dimensions do not match')
  return
end
[p,n1]=size(c);
if n1~=n
  disp('dimensions do not match')
  return
end
ib=16;bb=256;il=1;ih=n;

% remove and save the diagonal of a.

adiag=diag(a); a=a-diag(adiag,0);

% find zero rows in a & b and zero columns in a & c

perm=ones(1,n);for i=1:n,perm(i)=i;end;
j=ih;
while j>0,
  nozra=any((a(perm,perm))');
  nozb=any([b(perm,:),zeros(n,1)]')|nozra;
  nozca=any(a(perm,perm));
  nozc=any([c(:,perm);zeros(1,n)])|nozca;
  if ~nozb(j)|~nozc(j),
    for i=j+1:n,perm(i-1)=perm(i);end;
    n=n-1;ih=ih-1;j=ih+1;
    perm=perm(1:n);
  elseif j==1,break;
  end;
  j=j-1;
end;
j=ih;
while j>0,
  if ~nozra(j),t=perm(j);perm(j)=perm(ih);perm(ih)=t;
    t=nozra(j);nozra(j)=nozra(ih);nozra(ih)=t;
    ih=ih-1;j=ih+1;
  elseif j==1,break;
  end;
  j=j-1;
end;
j=il;
while j<=ih,
  if ~nozca(j),t=perm(j);perm(j)=perm(il);perm(il)=t;
    t=nozca(j);nozca(j)=nozca(il);nozca(il)=t;
    il=il+1;j=il-1;
  elseif j==ih,break;
  end;
  j=j+1;
end;
perm=perm(1:n);a=a(perm,perm);b=b(perm,:);c=c(:,perm);
adiag=adiag(perm);
%n,a,b,c,il,ih
% scaling

fail=1;
while fail,fail=0;
  for i=1:n,
    cosum=sum([abs(a(:,i));abs(c(:,i))]);
    rosum=sum([abs(a(i,:)), abs(b(i,:))]);
    f=1;g=rosum/ib;s=cosum+rosum;
    while cosum<g,f=f*ib;cosum=cosum*bb;end;
    g=rosum*ib;
    while cosum>=g,f=f/ib;cosum=cosum/bb;end;
    if(cosum+rosum)/f<0.95*s,
      g=1/f;fail=1;
      a(i,:)=a(i,:)*g;a(:,i)=a(:,i)*f;
      b(i,:)=b(i,:)*g;c(:,i)=c(:,i)*f;
    end;
  end;
end;
a=a+diag(adiag,0);
sysout=pck(a,b,c,d);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U_g,S_g,V_g] = vsvdcont(P_g)

% This function takes a CONTINUOUS frequency-by-frequency Singular Value
% Decomposition of the VARYING matrix "P_g" and outputs the results in the
% VARYING matrices "U_g", "S_g" and "V_g". It is very similar to MATLAB's
% VSVD command, but has some very important differences:
% (a) it works only on VARYING matrices,
% (b) it does not order the singular values in decreasing order,
% (c) if "P_g" have "m" outputs and "n" inputs and if "q = min{m,n}", 
%     then "U_g" has "m" outputs and "q" inputs, "S_g" has "q" outputs
%     and "q" inputs, and "V_g" has "n" outputs and "q" inputs. 
% This function also corrects the singular vectors in "U_g" and "V_g"
% with diagonal all-pass factors up to which they are unique at
% frequencies of non-repeated singular values. Hence the output
% variables "U_g", "S_g" and "V_g" are almost everywhere (ie: at all
% frequencies except those where repeated singular values occur) continuous
% functions of frequency.
% 
% USAGE  : [U_g,S_g,V_g] = vsvdcont(P_g)
% INPUT  : P_g = input  (m x n) VARYING matrix for Singular Value Decomposition.
% OUTPUT : U_g = output (m x q) VARYING matrix containing front Singular Vectors.
%          S_g = output (q x q) VARYING matrix containing Singular Values.
%          V_g = output (n x q) VARYING matrix containing back Singular Vectors.
%
% Alexander Lanzon - 25 May 2000


if ((nargin ~= 1) | (nargout ~= 3))
  disp('USAGE: [U_g,S_g,V_g] = vsvdcont(P_g)');
  return;
end


[mattype,m,n,num] = minfo(P_g);
q = min([m n]);
omega = getiv(P_g).';

flag = 'NOTok';
while ~(strcmp(flag,'ok'))
  [Utmp_g,Stmp_g,Vtmp_g] = vsvd(P_g);
  
  Utmp_g = sel(Utmp_g,1:m,1:q);
  Stmp_g = sel(Stmp_g,1:q,1:q);
  Vtmp_g = sel(Vtmp_g,1:n,1:q);
  
  % Plotting singular values of plant
  figure; vplot('liv,lm',vdiag(Stmp_g),'.');
          xlabel('Frequency (radians/sec)');
		  ylabel('Singular Values of Nominal Plant');
		  grid on; zoom on;
		  
  % Entering the frequency points at which singular values are repeated
  disp(' '); disp(' '); 
  disp('Enter a row vector containing the frequencies at which the');
  disp('following singular values intersect. If there are NO intersections,');
  disp('simply hit ENTER.'); 
  disp(' ');
  intersection_freqs = cell(q-1,1);
  for ii = 1:(q-1)
	intersection_freqs{ii} = input(['Singular Values ' num2str(ii) ' and ' num2str(ii+1) ': ']);
  end
  
  % Correcting ordering of singular values and singular vectors
  freq_points  = []; 
  what_to_swap = [];
  for ii = 1:(q-1)
	freq_points  = [freq_points  intersection_freqs{ii}];
	what_to_swap = [what_to_swap ii*ones(1,length(intersection_freqs{ii}))];
  end
  [freq_points,pos] = sort(freq_points);
  what_to_swap = what_to_swap(pos);
  
  Tr = eye(q); accum = Tr; 
  freq_points = [freq_points Inf];
  jj = 1; 
  for ii = 2:length(omega)
	if (omega(ii) < freq_points(jj))
	  accum = [accum; Tr];
	else
	  Tr = daug(eye(what_to_swap(jj)-1),[0 1;1 0],eye(q-what_to_swap(jj)-1))*Tr;
	  accum = [accum; Tr];
	  jj = jj + 1;
	end
  end
  Tr_g = vpck(accum,omega.');
  
  Utmp_g = mmult(Utmp_g,Tr_g);                % Note that Tr' = inv(Tr) at each frequency
  Stmp_g = mmult(vcjt(Tr_g),Stmp_g,Tr_g);
  Vtmp_g = mmult(Vtmp_g,Tr_g);
  
  % Correcting singular vectors with a diagonal all-pass factor
  U_old = var2con(Utmp_g,omega(1));
  V_old = var2con(Vtmp_g,omega(1));
  
  U_accum = U_old;
  V_accum = V_old;
  
  for ww = omega(2:end)
	U_new = var2con(Utmp_g,ww);
	V_new = var2con(Vtmp_g,ww);
	
	if (q == m)
	  DiagAllPass = diag(exp(-i*angle(diag(U_old'*U_new))));
	elseif (q == n)
	  DiagAllPass = diag(exp(-i*angle(diag(V_old'*V_new))));
	end
	
	U_accum = [U_accum; U_new*DiagAllPass];
	V_accum = [V_accum; V_new*DiagAllPass];
	
	U_old = U_new*DiagAllPass;
	V_old = V_new*DiagAllPass;
  end
  Utmp_g = vpck(U_accum,omega.');
  Vtmp_g = vpck(V_accum,omega.');
  
  % Plotting the resulting continuous singular values
  vplot('liv,lm',vdiag(Stmp_g),'.');
  xlabel('Frequency (radians/sec)');
  ylabel('Continuous Singular Values of Nominal Plant');
  grid on; zoom on;
  
  % Checking whether to repeat the whole procedure again
  disp(' '); 
  flag = input('Type "ok" if result is fine, anything else to repeat whole procedure: ','s');
  close;
end

% Setting output variables
U_g = Utmp_g;
S_g = Stmp_g;
V_g = Vtmp_g;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [msys] = minrealsys(sys,tol,dispflag)

% Minimal realisation for a SYSTEM matrix.
%
% This function removes all uncontrollable/unobservable states from 
% SYSTEM matrix "sys" to give a minimal SYSTEM realisation "msys".
%
% USAGE : [msys]   = minrealsys(sys,tol,dispflag)
% INPUT : sys      = non-minimal SYSTEM matrix.
%         tol      = tolerance for trucation of uncontrollable 
%                    and unobservable modes (DEFAULT = sqrt(eps)).
%         dispflag = 1 displays the number of states removed
%                    0 does not display the number of removed states.
%                    (DEFAULT = 1).
% OUTPUT: msys     = minimal SYSTEM matrix.
%
% Alexander Lanzon - 7 June 2000

if (nargin == 0) | (nargin > 3)
  disp('USAGE : [msys] = minrealsys(sys,tol,dispflag)');
  return;
end

if (nargin <= 1)
  tol = sqrt(eps);
end

if (nargin <= 2)
  dispflag = 1;
end

[a,b,c,d]     = unpck(syscl(sys));
SSsys         = ss(a,b,c,d);
mSSsys        = minreal(SSsys,tol,dispflag);
[am,bm,cm,dm] = ssdata(mSSsys);
msys          = syscl(pck(am,bm,cm,dm));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
