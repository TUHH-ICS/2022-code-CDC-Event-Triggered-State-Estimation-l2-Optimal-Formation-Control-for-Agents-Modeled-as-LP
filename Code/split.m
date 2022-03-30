function [out0,out1] = split(in1,in2,vt)
%SPLIT Summary of this function goes here
%   Detailed explanation goes here
out0 = (vt(2)*in1-vt(1)*in2)/(vt(2)-vt(1));
out1 = (in2-in1)/(vt(2)-vt(1));
end

