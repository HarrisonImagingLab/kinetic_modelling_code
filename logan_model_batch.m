

function [V_T,Tstar,STE] = logan_model_batch(plasma_name,tac_file,varargin);


%function [V_T, Tstar, varargout] = logan_model_single_tac(times,vals,plasmatimes,plasmavals,varargin)


%%%%%%% This should be 
%% logan_model_batch(plasma_name,tac_file,maxerr (default .2, optional), Tstar (optional, default is to derive from TAC));
%% IF you input your own T*, note that the maxerr *won't do anything*. That's fine, but for input parsing you need to put in a filler.
%% plasma should be the FITTED, RESAMPLED plasmas. this function does not have plasma fitting. That can be accomplished by 2TCM functions or the independent plasma fitting function. 

plasma=readtable(plasma_name);
plasmatimes=table2array(plasma(:,1));
plasmavals=table2array(plasma(:,2));

tac=readtable(tac_file);
times=(table2array(tac(:,1)) + table2array(tac(:,2)))/2;
tac=tac(:,[3:end]);
tac=[array2table(times),tac];
tac_array = table2array(tac(:,[2:end]))';

t1=size(tac_array);
t=t1(1);


maxerr=0.2;
%settingdefaultvalues


V_T=zeros(1,t);
Tstar=zeros(1,t);

for i=1:t


vals=tac_array(i,:)';


if length(varargin)==1
	maxerr=varargin{1};
	Tstar(i)=getTstar(maxerr, plasmatimes,plasmavals,times,vals);
elseif length(varargin)==2
	maxerr=varargin{1};
	Tstar(i)=varargin{2};
else
Tstar(i)=getTstar(maxerr, plasmatimes,plasmavals,times,vals);
end

plasmaintegrated=cumtrapz(plasmatimes,plasmavals);

indexs=find(times >= .99*Tstar(i));

nonindexs=find(times <.99*Tstar(i));

skipped_indexs = length(times) - length(indexs);
times_to_use=times(indexs);
clear plasma_time_indexs;
for j=1:length(times_to_use)
[v,plasma_time_indexs(j)]=min(abs(plasmatimes-times_to_use(j)));
end

%times_to_not_use=times(nonindexs);
%for j=1:length(times_to_not_use)
%[v,plasma_time_not_use_indexs(j)]=min(abs(plasmatimes-times_to_not_use(j)));
%end



Y=cumtrapz(times,vals)./vals;
numerator=plasmaintegrated(plasma_time_indexs);
X=numerator(:)./vals(indexs,:);

%unusedX= plasmaintegrated(plasma_time_not_use_indexs)./vals(nonindexs);

% this would have worked, but we need standard error too...

p=polyfit(X,Y(indexs),1);


V_T(i)=p(1);




newT=table(X,Y(indexs));
lme=fitlme(newT,'Var2~X','FitMethod','REML');

SE = 100*lme.Coefficients.SE(2)/lme.Coefficients.Estimate(2);

STE(i)=SE;



end




%G=figure
%hold on
%scatter(X,Y(indexs),'o','g');
%scatter(unusedX,Y(nonindexs),'x','r');
%refline(p(1),p(2));

