function [Tstar,DV]=getTstar(maxerr, plasma_times,plasma_vals,tac_times,tac_vals)

plasma_time_indexs=zeros(1,length(tac_times));

for i=1:length(tac_times)
[~,plasma_time_indexs(i)]=min(abs(plasma_times-tac_times(i)));
end


%%%generating left and right hand quantities for logan plot.
plasmaintegrated=cumtrapz(plasma_times,plasma_vals);
%%% sampling only at relevant time points
X=plasmaintegrated(plasma_time_indexs)./tac_vals;

%%%generating Y values
Y=cumtrapz(tac_times,tac_vals)./tac_vals;

%%computing T*

Tstar=inf;
for i=1:length(tac_vals)-1

[p,~]=polyfit(X([i:1:end]),Y([i:1:end]),1);
DV=p(1);
resids = abs(p(1)*X([i:1:end]) + p(2) - Y([i:1:end]));
residspercent= resids./abs(p(1)*X([i:1:end]) + p(2));

if max(residspercent(:)) < maxerr
Tstar=tac_times(i);
break
end

end


end
