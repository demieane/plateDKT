x1=xc_fem_data(1,end)
x2=xc_fem_data(1,1)
y1 = min(min(th_fem_data))
y2 = max(max(th_fem_data))

slope = (y1-y2)/(x1-x2);
ct = y1-slope*x1

tx_linear = slope *xc_fem_data  + ct;%0.015;
tx_modified = th_fem_data.*0 + tx_linear;


