function F = CalcFFromParameters(parameters)
%F = CalcFFromParameters(parameters)
%
%Calculate the fundamental matrix F from the parameters vector <parameters>
%of the form [[K(1) K(2) K(3)] [R1, R2, R3]'  [t1, t2, t3]']

K = CalcKFromParameters(parameters(1), parameters(2), parameters(3));

R = RotationMatrix(parameters(4), parameters(5), parameters(6));

tx = SkewMatrix(parameters(7), parameters(8), parameters(9));

K1 = K;
K2 = K;

F = (K1^-1)'*tx*R*(K2^-1);
end

