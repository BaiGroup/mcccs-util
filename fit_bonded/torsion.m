function f=torsion(v,x)

  x=x/180*pi;
  if numel(v)==4
    f=v(1)*cos(x)+v(2)*cos(2*x)+v(3)*cos(3*x);
  elseif numel(v)==5
    f=v(1)*cos(x)+v(2)*cos(2*x)+v(3)*cos(3*x)+v(4)*cos(4*x);
  end
  f=f+v(end);
end
