function [H_perm,JSD,Comp_JS]=cecp_param(data,dim)
% Input: data: one dimensional vector N x 1.
%        dim : integer type, Embedding dimension (type=numeric) 
%               Commonly used value of dim ranges from 3 to 7
%
% Outpur: H_perm: Normalised Permutation Entropy based on ordinal
%                 patterns.
%         JSD: Jansen Shannon Divergence
%         Comp_JS=Complexity based on Jansen Shannon Divergence.
% (!=factorial).
% dim=3; 
on=[0:1:dim-1]';
% All possible permutation patterns
possible_patterns=perms(on);
result=0;
result(1:length(possible_patterns),1)=0;

for i=1:(length(data)-(dim-1))
    temp=data(i:(i+(dim-1)),1);
    tempseq=[0:1:dim-1]';
    tempdata=[temp tempseq];
    tempdata=(sortrows(tempdata));
    
    for j=1:length(possible_patterns)
        if((possible_patterns(j,:))==transpose(tempdata(:,2)))
            result(j,1)=result(j,1)+1;
    end
    end
end
op=result; % Ordinal Patterns

%% Computation of Normalised Permutation Entropy
entropy_max=log(length(op));

for k=1:length(op)
    p(k,1)=op(k,1)/sum(op);
end
H_perm=-sum(p(p>0).*log(p(p>0)))/entropy_max; % Normalised Permutation Entropy

%% Computation of Complexity based on Jansen Shannon Divergence.

constant1= (0.5+((1 - 0.5)/length(op)))*log(0.5+((1 - 0.5)/length(op)));
constant2=((1 - 0.5)/length(op))*log((1 - 0.5)/length(op))*(length(op) - 1) ;
constant3= 0.5*log(length(op));
Q_o= -1/(constant1+constant2+constant3);

temp_op_prob=op/sum(op);
temp_op_prob2=(0.5*temp_op_prob)+(0.5*(1/length(op)));

H_temp_op_prob2= -sum(temp_op_prob2(temp_op_prob2>0).*log(temp_op_prob2(temp_op_prob2>0)));
H_temp_op_prob= -sum(temp_op_prob(temp_op_prob>0).*log(temp_op_prob(temp_op_prob>0)));

JSD=H_temp_op_prob2-(0.5*H_temp_op_prob)-(0.5*log(length(op))); % Jansen Shannon Divergence

Comp_JS=Q_o*JSD*H_perm; % Complexity based on Jansen Shannon Divergence

% Complexity Entropy Causality coordinates computation
cecp=[Comp_JS H_perm];


end

