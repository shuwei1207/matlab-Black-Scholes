
function Binomial_Value = Binomial_BS(S0,K,r,q,sigma,T,OptionType,ExerciseType,NT);
                          
% Parameters
%T=2;
%S0=50;
%K=52;
%r=0.05;
%q=0;
%sigma=0.3;
%OptionType='p'; %'p' for put option; otherwise call option
%ExerciseType='e'; %'a' for American option; otherwise European option

%Total number of time intervals %NT=252*T % if daily
%NT=252*T 
dt=T/NT;
u=exp(sigma*sqrt(dt));
d=1/u;
a=exp((r-q)*dt);
p=(a-d)/(u-d);

f=zeros(NT+1, NT+1);

% Option Prices at maturity
for j = 0:NT;
    if(OptionType=='p')
        % Put
        f(NT+1,j+1)=max(K-S0*(u^j)*(d^(NT-j)),0);
    else
        %Call
        f(NT+1,j+1)=max(S0*(u^j)*(d^(NT-j))-K,0);
    end;
end;
    

% Backward Induction
for i = (NT-1):-1:0; 
    for j = 0:i;
        if(OptionType=='p')
            if (ExerciseType=='a')
                %put american
                EV = max(K-S0*(u^j)*(d^i-j),0);
                f(i+1,j+1) = max(EV,exp(-r*dt)*(p*f(i+2,j+2)+(1-p)*f(i+2,j+1)));
            else
                %put european
                EV = 0;
                f(i+1,j+1) = max(EV,exp(-r*dt)*(p*f(i+2,j+2)+(1-p)*f(i+2,j+1)));
            end;
        else
            if (ExerciseType=='a')
                %call american
                EV = max(S0*(u^j)*(d^i-j)-K,0);
                f(i+1,j+1) = max(EV,exp(-r*dt)*(p*f(i+2,j+2)+(1-p)*f(i+2,j+1)));
            else
                %call european
                EV = 0;
                f(i+1,j+1) = max(EV,exp(-r*dt)*(p*f(i+2,j+2)+(1-p)*f(i+2,j+1)));
            end;
        end;
    end;
end;

Binomial_Value=f(1,1)
