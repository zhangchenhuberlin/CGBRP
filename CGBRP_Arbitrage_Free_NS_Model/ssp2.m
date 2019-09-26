
% This procedure is the standard Linear Gaussian State Space Model
% Using Kalman Filter
% 
function [XNS Xsam Ptt] = ssp2(yy,Zt,mu,F,Q,R,smo)

if (nargin <= 5); smo = 0; end

bigt  = size(yy,1);
K  = size(Zt,2);
yn = size(yy,2);

DD = eye(K);
%% Start Recursion (ref: P378 Hamilton)
XNS = zeros(K,bigt);
Ptt = zeros(K^2,bigt);
Ptl = zeros(K^2,bigt);
Xini = zeros(K,1);   % use 0 as the initial news
Xini(1,1)=0;
Xt_t = Xini;
vecP = (eye(K^2) - kron(F,F))\eye(K^2) * Q(:);  % vec(P_1|0)
%vecP=ones(K*K,1)*10^7;
Pt_t = reshape(vecP,K,K);    % P_1|0

for t = 1:bigt;
    %% Predicting
    Xt_lag = F * Xt_t;
    Pt_lag = F*Pt_t*F' + DD*Q*DD';
    itat_lag = yy(t,:)' - Zt * Xt_lag - mu;
    ft_lag = Zt*Pt_lag*Zt' + R * eye(yn);
        
    %% updating:
    Kalg = Pt_lag* Zt'*(ft_lag\eye(yn));  %% Kalman gain
    Xt_t = Xt_lag + Kalg * itat_lag;
    Pt_t = Pt_lag - Kalg * Zt* Pt_lag;
    XNS(:,t) = Xt_t;
    Pt_t = (Pt_t + Pt_t')/2;
    Ptt(:,t) = Pt_t(:);
    Ptl(:,t) = Pt_lag(:);
end;
   %Ptt
if smo == 1
    XSM = XNS;
    PtT = Ptt;
         for sm = 1 : (bigt - 1)
             Pt1_t = reshape(Ptl(:,bigt+1 - sm),K,K);
             Pt0_t = reshape(Ptt(:,bigt   - sm),K,K);
             Pt1_T = reshape(PtT(:,bigt+1 - sm),K,K);
             XSM(:, bigt - sm) = XNS(:,bigt-sm) + Pt0_t * F' * inv(Pt1_t)* (XSM(:, bigt + 1 - sm) - F * XNS(:,bigt-sm)); 
             PT = Pt0_t + Pt0_t * F' * inv(Pt1_t)* (Pt1_T - Pt1_t) * inv(Pt1_t)' * F * Pt0_t';
%             PT = chol(PT)' * chol(PT);
             PT = (PT + PT') / 2;
             PtT(:, bigt - sm) = PT(:);
         end
    XNS = real(XSM);
    Ptt = real(PtT);
end
    Xsam = XNS;
for j = 1:bigt
    Xsam(:,j) = normrnd(XNS(:,j),real(sqrt(diag(reshape(Ptt(:,j),K,K)))));
    %Xsam(:,j) = mvnrnd(XNS(:,j),reshape(Ptt(:,j),K,K));
end
    