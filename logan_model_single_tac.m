function [V_T, Tstar] = logan_model_single_tac(times,vals,plasmatimes,plasmavals,varargin)
%%%%%%% This should be 
%% logan_model_single_tac_nofig(times_of_tac,values_of_tac,plasmatimes,plasmavals,maxerr (default .2, optional), Tstar (optional, default is to derive from TAC));
%% IF you input your own T*, note that the maxerr *won't do anything*. That's fine, but for input parsing you need to put in a filler.



maxerr=0.2;
%settingdefaultvalues



if length(varargin)==1
	maxerr=varargin{1};
	Tstar=getTstar(maxerr, plasmatimes,plasmavals,times,vals);
elseif length(varargin)==2
	maxerr=varargin{1};
	Tstar=varargin{2};
else
Tstar=getTstar(maxerr, plasmatimes,plasmavals,times,vals);
end


plasmaintegrated=cumtrapz(plasmatimes,plasmavals);

indexs=find(times >= .99*Tstar);

nonindexs=find(times <.99*Tstar);

skipped_indexs = length(times) - length(indexs);
times_to_use=times(indexs);
for j=1:length(times_to_use)
[v,plasma_time_indexs(j)]=min(abs(plasmatimes-times_to_use(j)));
end





times_to_not_use=times(nonindexs);
for j=1:length(times_to_not_use)
[v,plasma_time_not_use_indexs(j)]=min(abs(plasmatimes-times_to_not_use(j)));
end



Y=cumtrapz(times,vals)./vals;
numerator=plasmaintegrated(plasma_time_indexs);
X=numerator(:)./vals(indexs,:);



p=polyfit(X,Y(indexs),1);


V_T=p(1);


G=figure
hold on
scatter(X,Y(indexs),'o','g');


if Tstar > min(times)
unusedX= plasmaintegrated(plasma_time_not_use_indexs)./vals(nonindexs);
scatter(unusedX,Y(nonindexs),'x','r');
end
refline(p(1),p(2));

