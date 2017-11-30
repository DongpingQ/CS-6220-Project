function [xtilde] = InexactLSR(A, b, s, type)
% This is the Inexact Least Square Regression algorithm
% ======================================
% Input: A, b: the LSR problem is min|Ax - b|_2
%        s: the parameter to generate CountSketch
% ======================================

d = size(A, 2);
if type == 1
    sketch = ( GaussianProjection([A,b]', s) )';
elseif type == 2
    sketch = ( CountSketch([A,b]', s) )';
end
Asketch = sketch(:, 1:d);
bsketch = sketch(:, end);
xtilde = Asketch\bsketch;

end

