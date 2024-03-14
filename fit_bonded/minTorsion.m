function [y,x]=minTorsion(v,x0)

  func=@(x)torsion(v,x);
  % y0=func(data(1,1));
  opts=optimset('FinDiffType','central','FunValCheck','on','Display','off','TolFun',1e-14,'TolX',1e-6,'Algorithm','interior-point');
  [x,fv,ef,out,lambda,grad,hessian]=fmincon(func,x0,[],[],[],[],[x0-pi],[x0+pi],[],opts);
  y=func(x);
end
