function [xtilde] = InexactLSR(A, b, s)
% This is the Inexact Least Square Regression algorithm
% ======================================
% Input: A, b: the LSR problem is min|Ax - b|_2
%        s: the parameter to generate CountSketch
% ======================================

d = size(A, 2);
sketch = ( CountSketch([A,b]', s) )';
Asketch = sketch(:, 1:d);
bsketch = sketch(:, end);
xtilde = Asketch\bsketch;

end

