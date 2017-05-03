% William Page (587000) - rads all the daq files in a single 2d variable

function [u1,v1,T1,P1] = read_datadaqs(path,ndaq)

% Intialise matrix of daqs
data = zeros(ndaq,2);
for i=1:ndaq
    data(i,:) = load([i,'.daq']);
end