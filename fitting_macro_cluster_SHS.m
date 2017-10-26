###########################################################
%MIT License

%Copyright (c) 2017 Tomasz M. Stawski

%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
###########################################################


clear all; clf;
###########################################################
%dependencies: NAN, optim
###########################################################
pkg load nan;
pkg load optim;
graphics_toolkit gnuplot; %This is configuration-specific. 
%Run "available_graphics_toolkits" to choose the most appropirate one. 

%--------------------------------------------------------%
%input and output definitions
%--------------------------------------------------------%
%distribution histogram from MCSAS fit
distfile = 'distribution_data/histogram_linear.dat';
%the non-aggregated silica scattering pattern 
reffile = 'distribution_data/silica_initial.dat'; 
%set input directory
imdir = 'ResultBkgSub_selected/'; 
%output directory
fitting_results_dir='fitting_results/test/';
mkdir (fitting_results_dir);

%--------------------------------------------------------%
%constants
%--------------------------------------------------------%
%background level from the MCSAS fit
bkg = 0.269/100; 
%contrast from silica in water in [m^-2]
SLD = (1.8831e-5-9.469e-6)*1e16;
%interaction potential in kT
U = 2.5; 

###########################################################
%define model functions
###########################################################
%--------------------------------------------------------%
%Sphere form factor
%--------------------------------------------------------%
function PofQ = Sphere(x,R)
    PofQ = (3.*(sin(x*R')-(x*R').*cos(x*R'))./(x*R').^3).^2;
end

%--------------------------------------------------------%
%Volume in cm^3 when R in nm
%--------------------------------------------------------%
function VofR = V(R)
   VofR = ((4/3)*pi()*(R.*1e-7).^3);
end

%--------------------------------------------------------%
%Debye-Bueche structure factor - simplified
%--------------------------------------------------------%
function S_q = S_DBsimp(x,par_A,D)
			S_q=par_A./(x.^D);
end

%--------------------------------------------------------%
%Gaussian distribution
%--------------------------------------------------------%
function DD = GaussianDist(R,R_mean,sigma)
		DD = ((1/sqrt(2*pi.*sigma.*sigma)).*exp(-(R-R_mean).*(R-R_mean)./ ...
		(2*sigma.*sigma)));
end

%--------------------------------------------------------%
%sticky HS structure factor
%--------------------------------------------------------%
function S_q=S_SHS(x,R,v,U)
	par_width=0.1;  
	if (v==0)
		S_q=1;
	else
		kappa = 2*x.*R;
		%width of apotential well in nm
		%minimum delta is 0.15 nm
		if (par_width.*2.*R > 0.15) 
			delta = par_width.*2.*R; 
		else
			delta = 0.15;
		endif	
		%display(delta);fflush(stdout);
		ni = delta./(delta+2*R);
		Vtau = (1/12)*(ni.^(-1)).*exp(-U);
		tau = ones(rows(x),1)*Vtau;
		%display(tau);fflush(stdout);
		Veta = (v.*((2.*R+delta)./(2.*R)).^3); 
		eta = ones(rows(x),1)*Veta;  		
		%display(eta);fflush(stdout);
		epsilon = tau + (eta./(1-eta));	
		gamma = v.*((1+eta/2)./(3*(1-eta).^2));
		lambda = (6*eta.^(-1)).*(epsilon-sqrt((epsilon.^2)-gamma));
		miu = lambda.*eta.*(1-eta);
		beta = -(3*eta.*(2+eta).^2-2*miu.*(1+7*eta+eta.^2)+(miu.^2).*(2+eta)) ...
		./(2*(1-eta).^4);
		alpha = ((1+2*eta-miu).^2)./((1-eta).^4);
		A = 2*(eta.*lambda./kappa).*sin(kappa);
		B = 2*((eta.*lambda./kappa).^2).*(1-cos(kappa));
		C = alpha.*(kappa.^3).*(sin(kappa)-kappa.*cos(kappa)); 
		D = beta.*(kappa.^2).*(2*kappa.*sin(kappa)-((kappa.^2)-2).*cos(kappa)-2);
		E = 0.5*eta.*alpha.*((4*kappa.^3-24*kappa).*sin(kappa)-((kappa.^4) ...
		-12*(kappa.^2)+24).*cos(kappa)+24);
		F = 24*eta./(kappa.^6);
		
		C_q = A - B - (C+D+E).*F;
		S_q = 1./(1-C_q);
	endif		
end

%--------------------------------------------------------%
%polydisperse SHS structure factor with Gaussian distribution
%--------------------------------------------------------%
function S_q = S_SHS_poly(x,dr,delta_RHS_avg,sigma,v,steps,U)
		if (v==0)
			S_q = 1;
		else		
			%range of delta_RHS
			delta_RHS_min = 0;
			delta_RHS_max = delta_RHS_avg+8*sigma;
			d_delta_RHS= (delta_RHS_max-delta_RHS_min)./(steps-1);
			delta_RHS = [delta_RHS_min:d_delta_RHS:delta_RHS_max]';
			
			S_q= nansum(S_SHS(x,dr'+delta_RHS',v,U).*GaussianDist(delta_RHS',...
			delta_RHS_avg,sigma),2)./nansum(GaussianDist(delta_RHS',...
			delta_RHS_avg,sigma),2);
			
		endif
end

%--------------------------------------------------------%
%complete fitting function allowing for the polydispersity of delta_RHS 
%using Gaussian distribution
%--------------------------------------------------------%
function result = FF_SHS_poly(x,dr,dist,normal,bkg,SLD,v,par_A,D,steps,fit_par,U)	
	integrand = 0;
	for i = 1:columns(dr)	
		%display(i);fflush(stdout);
		integrand = integrand + Sphere(x,dr(i)).*V(dr(i)).*dist(i)...
		.*(S_SHS_poly(x,dr(i),fit_par(1),fit_par(2),v,steps,U)+S_DBsimp(x,par_A,D));
	endfor
	result = (bkg + (SLD.^2).*integrand./normal);
end

%--------------------------------------------------------%

###########################################################
%Form factor calcutions based on the size distribution
###########################################################
%--------------------------------------------------------%
%read in the distribution data
%--------------------------------------------------------%
A=real(dlmread(distfile));
%select relevant data from the original MCSAS output file: X, Y, dY
matrixDist = [A(2:end,1).*1e9,A(2:end,3),A(2:end,4)]; 

%--------------------------------------------------------%
%interpolate the MC distibution to given number of steps
%--------------------------------------------------------%
number_step = 250;
dr_step = (matrixDist(rows(matrixDist),1)-matrixDist(1,1))/(number_step-1);
dr = [matrixDist(1,1):dr_step:matrixDist(rows(matrixDist),1)];
dist = (interp1(matrixDist(:,1),matrixDist(:,2),dr,"pchip"))';
%re-normalization of the distribution when changing the number of 
%interpolation steps
normal = number_step/(rows(matrixDist(:,1))-1);


%--------------------------------------------------------%
%Read in bkg-corrected data files from the input directory
%--------------------------------------------------------%
%count the files in the input directory
Files=dir([imdir,'*.*']); %count the files in the input directory
%select data frames to processes
first_frame = 1; %select data frames to processes e.g. 490 to 3500
last_frame=length(Files);

%Time the processing
total_time=0;

try
%follow up on a broken fit
		p_out_poly=dlmread([fitting_results_dir,...
		'fitting_parameters_SHS_poly_frames.out']);
		p_out_poly_error=dlmread([fitting_results_dir,...
		'fitting_parameters_SHS_poly_errors_frames.out']);
		p_out_mono=dlmread([fitting_results_dir,...
		'fitting_parameters_SHS_mono_frames.out']);
		p_out_mono_error=dlmread([fitting_results_dir,...
		'fitting_parameters_SHS_mono_errors_frames.out']);
catch
%or create a new empty results matrix of a right size
    p_out_poly = zeros(5,last_frame); 
    p_out_poly_error = zeros(5,last_frame);
    p_out_mono = zeros(4,last_frame);
    p_out_mono_error = zeros(4,last_frame);
 end_try_catch 

###########################################################
%--------------------------------------------------------%
%fitting including the structure factor using a
%Local monodisperse approximation - LMA
%J. Appl. Cryst. (1994). 27, 595-608

%In parts 1 and 2 we fit the model using a non-smeared 
%(monodisperse) S_SHS(q) model. These fit yield 4 parameters: 
%v, R_eHS, A, and D which are used to obtain a final smeared
%(polydisperse) fit in part 3, hence yielding v, <R_eHS>, sigma, A and D.

%--------------------------------------------------------%
%reduce the input data to make fits faster

%loop through the all data in the directory
for k=first_frame:last_frame  
	tic; %start measuring time
	try
	
	%read in the data file for the k-th frame	
		data=real(dlmread([imdir,Files(k).name])); 
		display(["Processing frame ",int2str(k-first_frame+1)," out of ",...
		int2str(last_frame-first_frame+1)]);fflush(stdout);
	%vector of indepoendet variables
		XX = real(data(:,1));
	%vector of observed values
		YY = real(data(:,2));
		dYY = real(data(:,3));
	%reduced data set
	%interpolate the corresponding intensity and uncertainties
		number_steps = 500;
		XX_reduced = [XX(1):(XX(end)-XX(1))/(number_steps-1):XX(end)]';
		YY_reduced = interp1(XX(1:end),YY(1:end),XX_reduced,"linear");
		dYY_reduced = interp1(XX(1:end),dYY(1:end),XX_reduced,"linear");
	
	%reduce fitting range to 0.3 - 2.55 nm^-1  and the number of points to 
	%a chosen value interpolate the corresponding intensity and uncertainties
	%used for fitting of the local maximum 
		number_steps = 100;
		XX_reduced_narrow = [XX(100):(XX(1000)-XX(100))/(number_steps-1):XX(1000)]';
		YY_reduced_narrow = interp1(XX(100:1000),YY(100:1000),...
		XX_reduced_narrow,"linear");
		dYY_reduced_narrow = interp1(XX(100:1000),dYY(100:1000),...
		XX_reduced_narrow,"linear");

%******************************FITTING     STEP 1      ************************
	
	%--------------------------------------------------------%
	%fitting function - step 1: fitting the low-q part - SHS part is fixed and not 
	%smeared (Eq. 5), aggregate - part is fitted
	%--------------------------------------------------------%
	%fitting function
	v=0.0; %v-parameter in Eq. 5
	delta_RHS=0; %aReHS-parameter in Eq. 5
	FF_cluster = @(x,p) (bkg +  (SLD.^2).*nansum((Sphere(x,dr').*V(dr).*dist'...
	./normal.*(S_SHS(x,dr+delta_RHS,v,U)+S_DBsimp(x,p(1),p(2)))),2));
	
	%initial fitting paramaters
	% p(1) par_A: A in Eq. 11
	% p(2) D: D in Eq. 11
	par_A=0.1;
	D = 1.8;
	%parameter matrix
	pin = [par_A;D]; 
	%scalar tolerance
	stol = 0.0001;
	%max number of iterations
	niter = 50;
	%weights
	wt = dYY_reduced; %this is not optimal actually, but it works
	%bidirectional numeric gradient
	dp = 0.01*ones(size(pin));
	%function for gradient
	dFdp = "dfdp";
	%bounds
	bound_lower = [0;1.00001];
	bound_upper = [1;4];
	bounds = [bound_lower,bound_upper];
	options.bounds = bounds;
	%fitting
	[f_step1,p_step1,cvg_step1,iter_step1,corp_step1,covp_step1,covr_step1,...
	stdresid_step1,Z_step1,r2_step1] = leasqr(XX_reduced,YY_reduced,pin, ...
	FF_cluster, stol,niter,wt,dp,dFdp, options);
	
%******************************FITTING     STEP 2      ************************
	%--------------------------------------------------------%
	%fitting function - step 2: fitting the high-q part, SHS part is fitted and 
	%not smeared (Eq. 5), agg - part is fixed
	%--------------------------------------------------------%
	%fitting function
	FF_SHS = @(x,p) log(bkg +  (SLD.^2).*nansum((Sphere(x,dr').*V(dr).*dist'...
	./normal.*(S_SHS(x,dr+p(1),p(2),U)+S_DBsimp(x,p_step1(1),p_step1(2)))),2));
	% p(1) delta_RHS; aReHS in Eq. 5
	% p(2) v; v in Eq. 5
	delta_RHS=0.5;
	v=0.1;
	%parameter matrix
	pin = [delta_RHS;v];
	%par_A and D are taken from the fit in step 1.
	
	%scalar tolerance
	stol = 0.0001;
	%max number of iterations
	niter = 50;
	%weights
	wt = dYY_reduced_narrow;
	%bidirectional numeric gradient
	dp = 0.01*ones(size(pin));
	%function for gradient
	dFdp = "dfdp";
	%bounds
	bound_lower = [0;0.0];
	bound_upper = [5;0.65];
	bounds = [bound_lower,bound_upper];
	options.bounds = bounds;
	%fitting
	[f_step2,p_step2,cvg_step2,iter_step2,corp_step2,covp_step2,covr_step2,...
	stdresid_step2,Z_step2,r2_step2] = leasqr(XX_reduced_narrow,...
	log(YY_reduced_narrow),pin, FF_SHS, stol,niter,wt,dp,dFdp, options);
	
%--------------------------------------------------------%
%final monodisperse results - combined 1 and 2
%partial
%--------------------------------------------------------%
	p_monodisperse=[p_step2;p_step1];
	p_out_mono(:,k) = p_monodisperse;
	dlmwrite([fitting_results_dir,'fitting_parameters_SHS_mono_frames.out'],...
	p_out_mono);
	d_delta_RHS_mono = sqrt(diag(covp_step2))(1,:);
	d_v_mono = sqrt(diag(covp_step2))(2,:);
	d_par_A_mono=sqrt(diag(covp_step1))(1,:);
	d_D_mono= sqrt(diag(covp_step1))(2,:);
	p_monodisperse_error = [d_delta_RHS_mono;d_v_mono;d_par_A_mono;d_D_mono];
	p_out_mono_error(:,k) = p_monodisperse_error;
	dlmwrite([fitting_results_dir,'fitting_parameters_SHS_mono_errors_frames.out']...
	, p_out_mono_error);
	
	
	%******************************FITTING     STEP 3      **********************
	%--------------------------------------------------------%
	%polydisperse parameters - step 3
	%--------------------------------------------------------%
	%initial guess based on the non-smeared fit in step 2
	delta_RHS = p_monodisperse(1); 
	%initial guess
	sigma = 0.05*delta_RHS; 

	v = p_monodisperse(2);
	par_A = p_monodisperse(3);
	D = p_monodisperse(4);
	int_acc = 120;
	%fitting function
	FF_polydisperse = @(x,p) (log(FF_SHS_poly(x,dr,dist,normal,bkg,...
	SLD,v,par_A,D,int_acc,p,U)));
	%parameter matrix
	pin = [delta_RHS;sigma]; 
	%scalar tolerance
	stol = 0.0001;
	%max number of iterations
	niter = 50;
	%weights
	wt = dYY_reduced_narrow;
	%bidirectional numeric gradient
	dp = 0.01*ones(size(pin));
	%function for gradient
	dFdp = "dfdp";
	%bounds
	bound_lower = [delta_RHS*0.95;0];
	bound_upper = [delta_RHS*1.05;3];
	bounds = [bound_lower,bound_upper];
	options.bounds = bounds;
	%fitting
	[f_step3,p_step3,cvg_step3,iter_step3,corp_step3,covp_step3,covr_step3,...
	stdresid_step3,Z_step3,r2_step3] = leasqr(XX_reduced_narrow,...
	log(YY_reduced_narrow),pin, FF_polydisperse, stol,niter,wt,dp,dFdp, options);

%--------------------------------------------------------%
%final polydisperse parameters combined 1 to 3
%--------------------------------------------------------%
	p_polydisperse = [p_step3(1);p_step3(2);p_monodisperse(2);p_monodisperse(3);...
	p_monodisperse(4)];
	p_out_poly(:,k) = p_polydisperse;
	dlmwrite([fitting_results_dir,'fitting_parameters_SHS_poly_frames.out'], ...
	p_out_poly);
	
	%fittting uncertainties
	d_delta_RHS = sqrt(diag(covp_step3))(1,:);
	d_sigma = sqrt(diag(covp_step3))(2,:);
	d_v = sqrt(diag(covp_step2))(2,:);
	d_par_A=sqrt(diag(covp_step1))(1,:);
	d_D= sqrt(diag(covp_step1))(2,:);
	p_polydisperse_error = [d_delta_RHS;d_sigma;d_v;d_par_A;d_D];
	p_out_poly_error(:,k) = p_polydisperse_error;
	dlmwrite([fitting_results_dir,'fitting_parameters_SHS_poly_errors_frames.out']...
	, p_out_poly_error);
	
	figure(1);clf;
	errorbar(XX,YY,dYY);hold on;
	plot(XX,FF_SHS_poly(XX,dr,dist,normal,bkg,SLD,p_polydisperse(3),...
	p_polydisperse(4),p_polydisperse(5),40,[p_polydisperse(1);...
	p_polydisperse(2)],U),'r','LineWidth',2);
	title (["Frame: ",int2str(k)," Filename: ",Files(k).name]);
	xlabel ("q [nm^{-1}]");
	ylabel ("I(q) [cm^{-1}]");
	set(gca,'YScale','log','XScale','log');
	grid on;
	print([fitting_results_dir,'fit_frame_SHS_poly-',int2str(k),'.pdf']);
	catch
end_try_catch

	total_time = total_time+toc;
	
	%on-screen information
	display(["Completed: ",int2str((k-first_frame+1)...
	./(last_frame-first_frame+1)*100),"%"]);fflush(stdout);
	display(["Estimated remaining time: ",int2str(total_time./(k-first_frame+1)...
	.*(last_frame-k)./60)," minutes."]);fflush(stdout);
	display(["Elapsed time: ",int2str(total_time./60)," minutes."]);
	fflush(stdout);
	%summary
	figure(2);
	errorbar([first_frame:k],p_out_poly(1,first_frame:k),...
	p_out_poly_error(1,first_frame:k));
	title ("<aR_{eHS}>");
	xlabel ("frame number");
	ylabel ("<aR_{eHS}> [nm]");
  axis([0 1 0 3],"autox"); 
	grid on;
	print([fitting_results_dir,'aR_eHS_poly.pdf']);

	figure(3);
	errorbar([first_frame:k],p_out_poly(2,first_frame:k),...
	p_out_poly_error(2,first_frame:k));
	title ("sigma");
	xlabel ("frame number");
	ylabel ("sigma [nm]");
  axis([0 1 0 2],"autox"); 
	grid on;
	print([fitting_results_dir,'sigma_poly.pdf']);

	figure(4);
	errorbar([first_frame:k],p_out_poly(3,first_frame:k),...
	p_out_poly_error(3,first_frame:k));
	title ("local volume fraction");
	xlabel ("frame number");
	ylabel ("v");
  axis([0 1 0 0.7],"autox");
	grid on;
	print([fitting_results_dir,'v.pdf']);

	figure(5);
	errorbar([first_frame:k],p_out_poly(4,first_frame:k),...
	p_out_poly_error(4,first_frame:k));
	title ("cluster size");
	xlabel ("frame number");
	ylabel ("A");
	grid on;
	print([fitting_results_dir,'par_A.pdf']);

	figure(6);
  errorbar([first_frame:k],p_out_poly(5,first_frame:k),...
	p_out_poly_error(5,first_frame:k));
	title ("fractal dimension");
	xlabel ("frame number");
	ylabel ("D");
  axis([0 1 1 4],"autox");
	grid on;
	print([fitting_results_dir,'D.pdf']);

	endfor

display(total_time);








