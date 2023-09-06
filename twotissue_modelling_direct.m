function [V_T,k_1,k_2,k_3,k_4,Vb,bd] = twotissue_modelling(plasma_name,blood_name,tac_file)


%read in plasma data

%subject name placeholder

subjnamepl_T = split(tac_file,'.');

subjnamepl=subjnamepl_T{1};

plasma=readtable(plasma_name);

plasmtimes=table2array(plasma(:,1));
plasmvals=table2array(plasma(:,2));


%%%unweighted
plasmtimes=[0;plasmtimes];
plasmvals=[0;plasmvals];

%to deal with NaNs from the interpolation process

plasmtimesn=[-10^10;plasmtimes];
plasmvalsn=[0;plasmvals];
begin=plasmtimesn(find(plasmvalsn==max(plasmvalsn(:))));

%% define function to fit plasma data with...




f = @(B,x) min(( x<begin).*interp1(plasmtimesn,plasmvalsn,(x < begin).*x) + (x >=begin).*(B(1).*exp(-(x-begin).*log(2)./B(2))+B(3).*exp(-(x-begin).*log(2)./B(4))+ B(5).*exp(-(x-begin).*log(2)./B(6))),1000);

%% perform fitting
B0=[45,1.5,30,25,15,40];
W=[0.0001;1./(.05.*plasmvals(2:end))];

W(isinf(W))=0.0001;

opts1=optimset('display','off','MaxFunEvals',1000,'MaxIter',1000,'FinDiffType','central');
Bf=lsqnonlin(@(Beta) W.*(f(Beta,plasmtimes)-plasmvals),B0,[0 0 0 0 0 0],[100000.0 120 100000.0 20000 100000.0 20000],opts1);






%define initial parameters
%%% old fitting
%B0=[38,1.5,25,25,12.7,40];

%%perform initial fitting of plasma data...

%mdl=fitnlm(plasmtimes,plasmvals,f,B0);
%B1=mdl.Coefficients.Estimate;






%%% re-fit with weights using initial 'guess' from unweighted fitting


%mdl=fitnlm(plasmtimes,plasmvals,f,B1,'Weights',[0.0001;1./(.05.^2.*plasmvals(2:end).^2)]);


%Bf=mdl.Coefficients.Estimate;

%define plasma function



X=linspace(0,90,10000);
Y=f(Bf,X);

%%write plasma function to folder, also creating PNG of fit/data for QC purposes...

activity_model=array2table([X'.*60,Y'],'VariableNames',{'time_seconds','plasma'});

writetable(activity_model,strcat(subjnamepl,'_matlabParent_Plasma_Activity_Model.txt'),'Delimiter','\t');




fi=figure('visible','off')

fplot(@(x) f(Bf,x),[0,90])
hold on
scatter(plasmtimes, plasmvals, 'x');
saveas(fi,strcat(subjnamepl,'_plasmafit.png'),'png');
close(fi)


%read in blood data, define blood function as simple linear interpolation, and get blood values for the frame times through interpolation.

blood=readtable(blood_name);
bloodtimes=table2array(blood(:,1))*60;
bloodvals=table2array(blood(:,2));
bloodvals=[0; 0; bloodvals; 0];
bloodtimes=[ -10^6; 0; bloodtimes; 10^6];
bloodfunc=@(x) interp1(bloodtimes,bloodvals,x);



%% read in TAc data. Using the first column of TAC (assumed to be wholebrain!) create weights using Poisson distribution SD/mean. 

tac=readtable(tac_file);
tactimes=(table2array(tac(:,1)) + table2array(tac(:,2)))/2;
tacdur=(table2array(tac(:,2)) - table2array(tac(:,1)));
tac=tac(:,[3:end])
tac=[array2table(tactimes),tac];
brain=table2array(tac(:,2));
WB_unc = brain.*exp(-1*tactimes*0.693147/(60*20.334));
kNEC=tacdur.*WB_unc;
NSD=1./sqrt(kNEC);
NSD(1)=NSD(2);
NSDsq=(NSD.^2);
blood_measures = interp1(bloodtimes,bloodvals,table2array(tac(:,1)));




%%% simplifying TAC for use with fitting...

tac_array = table2array(tac(:,[2:end]))';
t1=size(tac_array);
t=t1(1);
tac_minusblood=tac_array;
times=table2array(tac(:,1));

%%preallocating

V_T=zeros(1,t);
k_1=zeros(1,t);
k_2=zeros(1,t);
k_3=zeros(1,t);
k_4=zeros(1,t);
vb=zeros(1,t);

%%writing nsd for manual checking

NSD_model=array2table([times,NSD],'VariableNames',{'time_seconds','weights'});

writetable(NSD_model,strcat(subjnamepl,'NSD.txt'),'Delimiter','\t');

%doing whole brain with blood delay....comments here.

%initial parameters. Same as PNEURO but /60 since we're using seconds
%seconds A0
%A0=[.001 .01 .001 .0001,0.05,0];


A0=[.1 .1 .01 .01 0.05 0];
i=1;


%defining function handles for plasma and blood to simplify things a bit...

bloodfunc_m=@(x) interp1(bloodtimes,bloodvals,x*60);
plasmafunc_m= @(x) f(Bf,x);




%% defining 2tcm fit function.


%old seconds function
%%componeintblood_delay = @(A,x) (1-A(5)).*(x>0).*max(nanmin(integral(@(s) (x.*plasmafunc(((x.*s)-A(6)))).*(A(1).*exp(-A(2).*(x-x.*s)) + A(3).*exp(-A(4).*(x-x.*s))),0,1,'ArrayValued',true,'AbsTol',10^(-2),'RelTol',10^-(2)),10000),-10000)+(A(5)*bloodfunc(x-A(6)));



%now in minutes!
componeintblood_delay_m = @(A,x) (1-A(5)).*(x>0).*max(nanmin(integral(@(s) (x.*plasmafunc_m(((x.*s)+A(6)))).*((A(1).*( ((A(2) + A(3) + A(4) +  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2) -A(3) - A(4))./(((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)).*exp(-((A(2) + A(3) + A(4) +  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2).*(x-x.*s)) + (A(1).*( ((A(2) + A(3) + A(4) -  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2) - A(3) - A(4)) ./ (-1* ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)).*exp(-((A(2) + A(3) + A(4) -  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2).*(x-x.*s))),0,1,'ArrayValued',true,'AbsTol',10^(-2),'RelTol',10^-(2)),10000),-10000)+(A(5)*bloodfunc_m(x+A(6)));

% performing nonlinear fit
%opts1=optimset('algorithm','levenberg-marquardt','display','off')
opts1=optimset('display','off','MaxFunEvals',10000,'MaxIter',10000,'FinDiffType','central');
tic
%old seconds version...
%A=lsqnonlin(@(Beta) (1./NSD).*(componeintblood_delay(Beta,times)-tac_minusblood(i,:)'),A0,[0 0 0 0 0 -10],[.1 .1 .1 .1 .1, 10],opts1);

[A,~,residual,~,~,~,jac]=lsqnonlin(@(Beta) (1./sqrt(NSD)).*(componeintblood_delay_m(Beta,times./60)-tac_minusblood(i,:)'),A0,[0 0 0 0 0 -.2],[8 8 8 8 .5 .2],opts1);



%you can compute CI's to error check things if you wish. This is sort of the same as computing standard error like PMOD reports, but arguably more useful....these *may not be accurate* when dealing with bounded coefficients as we are...
% CI{1} = nlparci(A,residual,'jacobian',jac);
toc
%setting blood delay and initial parameters...

bd=A(6);
A0=[A(1) A(2) A(3) A(4) A(5)];



%redefining function for all regions, including blood delay fit from whole brain!

componeint_m = @(A,x) (1-A(5)).*(x>0).*max(nanmin(integral(@(s) (x.*plasmafunc_m(((x.*s)+bd))).*((A(1).*( ((A(2) + A(3) + A(4) +  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2) -A(3) - A(4))./(((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)).*exp(-((A(2) + A(3) + A(4) +  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2).*(x-x.*s)) + (A(1).*( ((A(2) + A(3) + A(4) -  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2) - A(3) - A(4)) ./ (-1* ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)).*exp(-((A(2) + A(3) + A(4) -  ((A(2)+A(3)+A(4)).^2 -4.*A(2).*A(4)).^0.5)/2).*(x-x.*s))),0,1,'ArrayValued',true,'AbsTol',10^(-2),'RelTol',10^-(2)),10000),-10000)+(A(5)*bloodfunc_m(x+bd));

%%% redoing wholebrain fit WITHOUT BD for consistency -- and doing it seperately to generate a figure. 


[A,~,residual,~,~,~,jac]=lsqnonlin(@(Beta) (1./NSD).*(componeint_m(Beta,times./60)-tac_minusblood(i,:)'),A0,[0 0 0 0 0],[1 1 1 1 .15],opts1);

%you can compute CI's to error check things if you wish. This is sort of the same as computing standard error like PMOD reports, but arguably more useful....these *may not be accurate* when dealing with bounded coefficients as we are...
% CI{i} = nlparci(A,residual,'jacobian',jac);

k_1(i) = A(1);
k_2(i) = A(2);
k_3(i) = A(3);
k_4(i)= A(4);
Vb(i)=A(5);

V_T(i) = (k_1(i)/k_2(i))*(1 + (k_3(i)/k_4(i)));





%%% generating WB fit plot for QC purposes.

fi=figure('visible','off')

fplot(@(x) componeint_m(A,x),[0,90])
hold on
scatter(times./60, tac_minusblood(i,:)', 'x');
saveas(fi,strcat(subjnamepl,'_WBfit.png'),'png');
close(fi)









%resetting initial guess -- now we don't need a parameter for blood delay; using parameters from initial!

A0=[A(1) A(2) A(3) A(4) A(5)];

for i=2:t
warning('off','all')
tic


%A=nlinfit(times,tac_minusblood(i,:)',componeint,A0,'Weights',1./(NSD.^2),'Options',options);

[A,~,residual,~,~,~,jac]=lsqnonlin(@(Beta) (1./NSD).*(componeint_m(Beta,times./60)-tac_minusblood(i,:)'),A0,[0 0 0 0 0],[1 1 1 1 .15],opts1);

%you can compute CI's to error check things if you wish. This is sort of the same as computing standard error like PMOD reports, but arguably more useful....these *may not be accurate* when dealing with bounded coefficients as we are...
% CI{i} = nlparci(A,residual,'jacobian',jac);

k_1(i) = A(1);
k_2(i) = A(2);
k_3(i) = A(3);
k_4(i)= A(4);
Vb(i)=A(5);

V_T(i) = (k_1(i)/k_2(i))*(1 + (k_3(i)/k_4(i)));



if i<6
fi=figure('visible','off');

fplot(@(x) componeint_m(A,x),[0,90]);
hold on
scatter(times./60, tac_minusblood(i,:)', 'x');
saveas(fi,strcat(subjnamepl,'_',int2str(i),'_WBfit.png'),'png');
close(fi)
end


toc
end

























