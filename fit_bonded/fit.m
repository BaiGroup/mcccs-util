%% import data s.t. x=data(:,1) and y=data(:,4)

%% Torsion
v0=[-213.7684 -513.1350  -32.6299 -119.4727  879.0060]; weightT=-1e4;

func=@(v)objFunc(v,data,weightT);
opts=optimset('FinDiffType','central','FunValCheck','on','Display','off','TolFun',1e-14,'TolX',1e-6,'Algorithm','interior-point');

%% Local Minimum %%
opts=optimset(opts,'Display','iter');%,'Algorithm',{'levenberg-marquardt',.005});
[v,fv,ef,out,grad,hess]=fminunc(func,v0,opts)
% [v,fv,ef,out]=fminsearch(func,v0,opts)

%% Convert to OPLS convention and plot
[y0,i0]=min(torsion(v,data(:,1)));
[y0,x0]=minTorsion(v,data(i0,1));
x0/pi*180
v(end)=v(end)-y0
figure; plot(data(:,1),data(:,4),'b<'); hold on; plot(data(:,1),torsion(v,data(:,1)),'c'); xlabel('\phi'); ylabel('E [K]');% plot(data(:,1),data(:,2),'rd');
if numel(v)==4
    opls=[v(2)+v(4)-v(1)-v(3),v(1),-v(2),v(3)]
elseif numel(v)==5
    opls=[v(2)+v(4)+v(5)-v(1)-v(3),v(1),-v(2),v(3),-v(4)]
end
