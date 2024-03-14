function f=objFunc(v,data,annealT)

  y=torsion(v,data(:,1));
  [y0,i0]=min(y);
  y=y-minTorsion(v,data(i0,1));
  if annealT<0
    f=sum(abs(y-data(:,4)).*exp(-data(:,4)/annealT));
  else
    f=sum((y-data(:,4)).^2);
  end
end
